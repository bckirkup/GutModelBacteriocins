#include "chemical_field.h"
#include "chem_environment_config.h"
#include "agent.h"
#include "domain.h"
#include "error.h"
#include "input_parser.h"
#include "species_names.h"
#include "vbf.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>
#include <string_view>
#include <vector>

using namespace gutibm;

namespace {

Domain make_domain(Int nx, Int ny, Int nz) {
  DomainConfig cfg;
  cfg.lo = {0.0, 0.0, 0.0};
  cfg.hi = {nx * 5.0e-6, ny * 5.0e-6, nz * 5.0e-6};
  cfg.grid_dx = 5.0e-6;
  cfg.hash_cell_size = 10.0e-6;
  Domain domain;
  domain.init(cfg);
  return domain;
}

ChemicalSpec diffusing_species(Real diffusion, Real initial, Real boundary) {
  ChemicalSpec spec;
  spec.name = species::OXYGEN;
  spec.diff_coeff = diffusion;
  spec.retardation = 1.0;
  spec.initial_conc = initial;
  spec.boundary_conc = boundary;
  spec.diffusion_enabled = true;
  return spec;
}

ChemicalSpec delivery_species(EpithelialBoundaryMode mode, Real initial,
                              Real boundary, Real transfer_coeff,
                              Real flux) {
  ChemicalSpec spec = diffusing_species(2.1e-9, initial, boundary);
  spec.epithelial_boundary_mode = mode;
  spec.epithelial_transfer_coeff = transfer_coeff;
  spec.epithelial_flux = flux;
  spec.z_gradient_enabled = false;
  return spec;
}

Real inventory(const ChemicalField& chem, const Domain& domain) {
  const Real cell_volume = domain.cell_volume();
  Real total = 0.0;
  for (Int cell = 0; cell < chem.ncells(); ++cell) {
    total += chem.conc(0, cell) * cell_volume;
  }
  return total;
}

Real total_inventory(const ChemicalField& chem, const Domain& domain,
                     Int species_index) {
  Real total = 0.0;
  for (Int cell = 0; cell < chem.global_ncells(); ++cell) {
    if (!chem.owns_global_cell(cell)) continue;
    total += chem.conc_global(species_index, cell)
        * domain.cell_volume();
  }
  return total;
}

struct DeliveryClosureResult {
  Real initial = 0.0;
  Real final = 0.0;
  Real boundary = 0.0;
  Real gradient_source = 0.0;
  Real vbf_source = 0.0;
  Real agent_uptake = 0.0;
  Real maintenance = 0.0;
  Real vbf_sink = 0.0;
  Real reaction_clip = 0.0;
  Real residual = 0.0;
  Real relative_residual = 0.0;
};

DeliveryClosureResult run_delivery_closure(
    std::string_view species_name, bool gradient_enabled) {
  constexpr Real dt = 60.0;
  constexpr Int steps = 8;
  Domain domain = make_domain(1, 1, 6);
  ChemicalSpec spec = diffusing_species(1.0e-20, 1.0e-3, 1.0e-3);
  spec.name = std::string(species_name);
  spec.delivery_enabled = true;
  spec.z_gradient_enabled = gradient_enabled;
  spec.z_gradient_lambda = 5.0e-6;

  ChemicalField chem;
  chem.init(domain, {spec});
  const Int species_index = chem.find(spec.name);
  const Int agent_cell = domain.cell_index(0, 0, 2);
  Agent agent = Agent::create_default(
      1, 1, {2.5e-6, 2.5e-6, 12.5e-6}, 0.0);
  agent.grid_cell = agent_cell;
  const Real cell_volume = domain.cell_volume();
  VBFConfig vbf_config;
  vbf_config.mucin_liberation = 5.0e-5;
  vbf_config.carbon_sink_vmax = 0.0;
  vbf_config.nutrient_sink = 0.0;
  VBF vbf;
  vbf.init(vbf_config, domain);
  OxygenConfig oxygen;
  AcetateConfig acetate;
  MucinConfig mucin;

  DeliveryClosureResult result;
  result.initial = total_inventory(chem, domain, species_index);
  for (Int step = 0; step < steps; ++step) {
    chem.zero_reactions();
    const Real concentration = chem.conc_global(species_index, agent.grid_cell);
    const Real demand = 0.25 * concentration * cell_volume;
    const Real sink_rate = demand / (concentration * cell_volume * dt);
    chem.flux_accounting().add_uptake_demand(species_index, demand);
    chem.add_sink_rate_global(species_index, agent.grid_cell, sink_rate);
    VbfFluxTotals vbf_totals;
    vbf.apply_nutrient_coupling(
        chem, domain, dt, oxygen, acetate, mucin, &vbf_totals);
    chem.flux_accounting().add_interval(
        species_index, 0.0, vbf_totals.carbon_source,
        vbf_totals.carbon_sink, 0.0);
    for (Int cell = 0; cell < chem.global_ncells(); ++cell) {
      if (!chem.owns_global_cell(cell)) continue;
      const Real updated = chem.conc_global(species_index, cell)
          + chem.reac_global(species_index, cell) * dt;
      if (updated < 0.0) {
        chem.flux_accounting().add_reaction_clip(
            species_index, -updated * cell_volume);
      }
      chem.conc_global(species_index, cell) = std::max(updated, 0.0);
    }

    chem.apply_diffusion(domain, dt);
    Real realized = 0.0;
    for (Int cell = 0; cell < chem.global_ncells(); ++cell) {
      if (chem.owns_global_cell(cell)) {
        realized += chem.sink_realized_global(species_index, cell);
      }
    }
    chem.flux_accounting().add_agent_uptake(species_index, realized);
    chem.flux_accounting().add_uptake_shortfall(
        species_index, std::max(demand - realized, 0.0));
    if (realized + 1.0e-30 < demand) {
      chem.flux_accounting().add_uptake_limited(species_index, 1.0);
    }
    chem.flux_accounting().commit_agent_uptake_step();
    chem.flux_accounting().commit_boundary_and_reaction_step();
    chem.flux_accounting().close_interval();
  }

  result.final = total_inventory(chem, domain, species_index);
  const auto& flux = chem.flux_accounting();
  result.boundary = flux.boundary_cumulative[
      static_cast<size_t>(species_index)];
  result.gradient_source = flux.gradient_source_cumulative[
      static_cast<size_t>(species_index)];
  result.vbf_source = flux.vbf_source_cumulative[
      static_cast<size_t>(species_index)];
  result.agent_uptake = flux.agent_uptake_cumulative[
      static_cast<size_t>(species_index)];
  result.maintenance = flux.maintenance_cumulative[
      static_cast<size_t>(species_index)];
  result.vbf_sink = flux.vbf_sink_cumulative[
      static_cast<size_t>(species_index)];
  result.reaction_clip = flux.reaction_clip_cumulative[
      static_cast<size_t>(species_index)];
  const Real lhs = result.initial + result.boundary + result.gradient_source
      + result.vbf_source
      - result.agent_uptake
      - result.maintenance - result.vbf_sink + result.reaction_clip;
  result.residual = lhs - result.final;
  const Real scale = std::max(
      {std::abs(lhs), std::abs(result.final), 1.0e-300});
  result.relative_residual = std::abs(result.residual) / scale;
  return result;
}

void test_delivery_mass_closure_gradient_parameterization() {
  const DeliveryClosureResult no_gradient =
      run_delivery_closure(species::CARBON, false);
  const DeliveryClosureResult carbon_gradient =
      run_delivery_closure(species::CARBON, true);
  const DeliveryClosureResult oxygen_no_gradient =
      run_delivery_closure(species::OXYGEN, false);
  const DeliveryClosureResult oxygen_gradient =
      run_delivery_closure(species::OXYGEN, true);

  const auto report = [](std::string_view label,
                         const DeliveryClosureResult& result) {
    std::cout << "  " << label << " initial=" << result.initial
              << " final=" << result.final
              << " boundary=" << result.boundary
              << " gradient_source=" << result.gradient_source
              << " vbf_source=" << result.vbf_source
              << " agent_uptake=" << result.agent_uptake
              << " maintenance=" << result.maintenance
              << " vbf_sink=" << result.vbf_sink
              << " reaction_clip=" << result.reaction_clip
              << " residual=" << result.residual
              << " relative_residual=" << result.relative_residual << "\n";
  };
  report("carbon_gradient_off", no_gradient);
  report("carbon_gradient_on", carbon_gradient);
  report("oxygen_gradient_off", oxygen_no_gradient);
  report("oxygen_gradient_on", oxygen_gradient);

  constexpr Real tolerance = 1.0e-6;
  assert(no_gradient.relative_residual <= tolerance);
  assert(carbon_gradient.relative_residual <= tolerance);
  assert(oxygen_no_gradient.relative_residual <= tolerance);
  assert(oxygen_gradient.relative_residual <= tolerance);
}

void test_delivery_step_mass_closure_dirichlet_refill() {
  constexpr Real dt = 60.0;
  constexpr Real diffusion = 4.166666666666667e-9;
  constexpr Real initial = 1.0e-3;
  constexpr Real boundary = 1.0e-3;
  constexpr Real gradient_lambda = 5.0e-6;
  constexpr Real sink_rate = 1.0e3;
  Domain domain = make_domain(1, 1, 6);
  ChemicalSpec spec = diffusing_species(
      diffusion, initial, boundary);
  spec.delivery_enabled = true;
  spec.z_gradient_enabled = true;
  spec.z_gradient_lambda = gradient_lambda;

  ChemicalField chem;
  chem.init(domain, {spec});
  const Int species_index = chem.find(species::OXYGEN);
  const Int sink_cell = domain.cell_index(0, 0, 1);
  chem.zero_reactions();
  const Real before = inventory(chem, domain);
  chem.add_sink_rate_global(species_index, sink_cell, sink_rate);
  chem.apply_diffusion(domain, dt);

  Real realized = 0.0;
  for (Int cell = 0; cell < chem.global_ncells(); ++cell) {
    if (chem.owns_global_cell(cell)) {
      realized += chem.sink_realized_global(species_index, cell);
    }
  }
  chem.flux_accounting().add_agent_uptake(species_index, realized);
  chem.flux_accounting().commit_agent_uptake_step();
  chem.flux_accounting().commit_boundary_and_reaction_step();

  const Real after = inventory(chem, domain);
  const auto& flux = chem.flux_accounting();
  const auto index = static_cast<size_t>(species_index);
  const Real residual = before + flux.boundary_last_step[index]
      + flux.gradient_source_last_step[index]
      + flux.vbf_source_last_step[index]
      - flux.agent_uptake_last_step[index]
      - flux.maintenance_last_step[index]
      - flux.vbf_sink_last_step[index]
      + flux.reaction_clip_last_step[index] - after;
  const Real scale = std::max(
      {std::abs(before), std::abs(after), std::abs(realized), 1.0e-300});
  const Real relative_residual = std::abs(residual) / scale;

  std::cout << "  test_delivery_step_mass_closure_dirichlet_refill:"
            << " alpha=" << diffusion * dt
                / (domain.dx_z() * domain.dx_z())
            << " boundary_step=" << flux.boundary_last_step[index]
            << " vbf_source_step=" << flux.vbf_source_last_step[index]
            << " realized=" << realized
            << " residual=" << residual
            << " relative_residual=" << relative_residual << "\n";
  assert(relative_residual <= 1.0e-12);
}

void test_delivery_step_mass_closure_vbf_source() {
  constexpr Real dt = 60.0;
  constexpr Real diffusion = 4.166666666666667e-9;
  constexpr Real initial = 1.0e-3;
  constexpr Real boundary = 1.0e-3;
  constexpr Real gradient_lambda = 5.0e-6;
  constexpr Real sink_rate = 1.0e3;
  Domain domain = make_domain(1, 1, 6);
  ChemicalSpec spec = diffusing_species(
      diffusion, initial, boundary);
  spec.name = species::CARBON;
  spec.delivery_enabled = true;
  spec.z_gradient_enabled = true;
  spec.z_gradient_lambda = gradient_lambda;

  ChemicalField chem;
  chem.init(domain, {spec});
  const Int species_index = chem.find(species::CARBON);
  const Int sink_cell = domain.cell_index(0, 0, 1);
  VBFConfig vbf_config;
  vbf_config.mucin_liberation = 5.0e-5;
  vbf_config.carbon_sink_vmax = 0.0;
  vbf_config.nutrient_sink = 0.0;
  VBF vbf;
  vbf.init(vbf_config, domain);
  OxygenConfig oxygen;
  AcetateConfig acetate;
  MucinConfig mucin;

  chem.zero_reactions();
  const Real before = inventory(chem, domain);
  chem.add_sink_rate_global(species_index, sink_cell, sink_rate);
  VbfFluxTotals vbf_totals;
  vbf.apply_nutrient_coupling(
      chem, domain, dt, oxygen, acetate, mucin, &vbf_totals);
  chem.flux_accounting().add_interval(
      species_index, 0.0, vbf_totals.carbon_source,
      vbf_totals.carbon_sink, 0.0);
  for (Int cell = 0; cell < chem.global_ncells(); ++cell) {
    const Real updated = chem.conc_global(species_index, cell)
        + chem.reac_global(species_index, cell) * dt;
    if (updated < 0.0) {
      chem.flux_accounting().add_reaction_clip(
          species_index, -updated * domain.cell_volume());
    }
    chem.conc_global(species_index, cell) = std::max(updated, 0.0);
  }
  chem.apply_diffusion(domain, dt);

  Real realized = 0.0;
  for (Int cell = 0; cell < chem.global_ncells(); ++cell) {
    if (chem.owns_global_cell(cell)) {
      realized += chem.sink_realized_global(species_index, cell);
    }
  }
  chem.flux_accounting().add_agent_uptake(species_index, realized);
  chem.flux_accounting().commit_agent_uptake_step();
  chem.flux_accounting().commit_boundary_and_reaction_step();

  const Real after = inventory(chem, domain);
  const auto& flux = chem.flux_accounting();
  const auto index = static_cast<size_t>(species_index);
  const Real vbf_source = flux.vbf_source_last_step[index];
  const Real residual = before + flux.boundary_last_step[index]
      + flux.gradient_source_last_step[index] + vbf_source
      - flux.agent_uptake_last_step[index]
      - flux.maintenance_last_step[index]
      - flux.vbf_sink_last_step[index]
      + flux.reaction_clip_last_step[index] - after;
  const Real residual_without_vbf_source = residual - vbf_source;
  const Real scale = std::max(
      {std::abs(before), std::abs(after), std::abs(realized),
       std::abs(vbf_source), 1.0e-300});
  const Real relative_residual = std::abs(residual) / scale;

  std::cout << "  test_delivery_step_mass_closure_vbf_source:"
            << " vbf_source_step=" << vbf_source
            << " pre_fix_residual=" << residual_without_vbf_source
            << " residual=" << residual
            << " relative_residual=" << relative_residual << "\n";
  assert(vbf_source > 0.0);
  assert(std::abs(residual_without_vbf_source) > 1.0e-12 * scale);
  assert(relative_residual <= 1.0e-12);
}

void test_anisotropic_diffusion_invariants() {
  DomainConfig cfg;
  cfg.lo = {0.0, 0.0, 0.0};
  cfg.hi = {50.0e-6, 20.0e-6, 10.0e-6};
  cfg.grid_dx = 5.0e-6;
  cfg.chemistry_stride = {2, 1, 1};
  Domain domain;
  domain.init(cfg);
  ChemicalField chem;
  chem.init(domain, {diffusing_species(1.0e-12, 0.0, 0.0)});
  const Int source = domain.cell_index(2, 1, 1);
  chem.conc(0, source) = 1.0;
  const Real before = inventory(chem, domain);
  const Real dt = 2.5;
  chem.apply_periodic_x_diffusion(domain, dt);
  const Real after = inventory(chem, domain);
  assert(std::abs(after - before) <= 1.0e-12 * before);
  for (Int ix = 0; ix < domain.nx(); ++ix) {
    const Int cell = domain.cell_index(ix, 1, 1);
    const Int left = domain.cell_index(
        (ix + domain.nx() - 1) % domain.nx(), 1, 1);
    const Int right = domain.cell_index(
        (ix + 1) % domain.nx(), 1, 1);
    const Real value = chem.conc(0, cell);
    const Real laplacian = (chem.conc(0, left) + chem.conc(0, right)
        - 2.0 * value) / (domain.dx_x() * domain.dx_x());
    const Real residual = value - dt * 1.0e-12 * laplacian;
    assert(std::abs(residual - (ix == 2 ? 1.0 : 0.0)) < 1.0e-12);
    assert(value >= 0.0);
  }
  assert(std::abs(chem.conc(0, domain.cell_index(1, 1, 1))
                  - chem.conc(0, domain.cell_index(3, 1, 1)))
         < 1.0e-12);
  std::cout << "  test_anisotropic_diffusion_invariants: PASSED\n";
}

void test_stride_sensitivity_is_ordered() {
  std::array<Real, 3> neighbor_values{};
  const std::array<Int, 3> strides = {1, 2, 4};
  for (size_t i = 0; i < strides.size(); ++i) {
    DomainConfig cfg;
    cfg.hi = {80.0e-6, 10.0e-6, 10.0e-6};
    cfg.grid_dx = 5.0e-6;
    cfg.chemistry_stride = {strides[i], 1, 1};
    Domain domain;
    domain.init(cfg);
    ChemicalField chem;
    chem.init(domain, {diffusing_species(1.0e-12, 0.0, 0.0)});
    const Int source = domain.cell_index(domain.nx() / 2, 0, 0);
    chem.conc(0, source) = 1.0;
    chem.apply_periodic_x_diffusion(domain, 2.5);
    neighbor_values[i] = chem.conc(
        0, domain.cell_index(domain.nx() / 2 + 1, 0, 0));
  }
  std::cout << "  stride neighbors: " << neighbor_values[0] << ", "
            << neighbor_values[1] << ", " << neighbor_values[2] << "\n";
  assert(neighbor_values[0] > neighbor_values[1]);
  assert(neighbor_values[1] > neighbor_values[2]);
  std::cout << "  test_stride_sensitivity_is_ordered: PASSED\n";
}

void test_anisotropic_yz_residuals() {
  DomainConfig cfg;
  cfg.hi = {20.0e-6, 40.0e-6, 60.0e-6};
  cfg.grid_dx = 5.0e-6;
  cfg.chemistry_stride = {1, 2, 4};
  Domain domain;
  domain.init(cfg);
  ChemicalSpec spec = diffusing_species(1.0e-12, 0.0, 0.0);

  ChemicalField y_field;
  y_field.init(domain, {spec});
  const Int y_source = domain.cell_index(0, domain.ny() / 2, 0);
  y_field.conc(0, y_source) = 1.0;
  const Real y_before = inventory(y_field, domain);
  y_field.apply_periodic_y_diffusion(domain, 2.5, 0);
  assert(std::abs(inventory(y_field, domain) - y_before)
         <= 1.0e-12 * y_before);
  const Real y_value = y_field.conc(
      0, domain.cell_index(0, domain.ny() / 2, 0));
  const Real y_left = y_field.conc(
      0, domain.cell_index(0, domain.ny() / 2 - 1, 0));
  const Real y_right = y_field.conc(
      0, domain.cell_index(0, domain.ny() / 2 + 1, 0));
  const Real y_residual = y_value - 2.5e-12
      * (y_left + y_right - 2.0 * y_value)
      / (domain.dx_y() * domain.dx_y());
  assert(std::abs(y_residual - 1.0) < 1.0e-12);
  assert(std::abs(y_left - y_right) < 1.0e-12);

  ChemicalField z_field;
  z_field.init(domain, {spec});
  const Int z_source = domain.cell_index(0, 0, 1);
  z_field.conc(0, z_source) = 1.0;
  z_field.apply_bounded_z_diffusion(domain, 2.5, 0);
  for (Int iz = 1; iz < domain.nz(); ++iz) {
    assert(z_field.conc(0, domain.cell_index(0, 0, iz)) >= 0.0);
  }
  const Real z_value = z_field.conc(0, z_source);
  const Real z_below = spec.boundary_conc;
  const Real z_above = z_field.conc(0, domain.cell_index(0, 0, 2));
  const Real z_residual = z_value - 2.5e-12
      * (z_below + z_above - 2.0 * z_value)
      / (domain.dx_z() * domain.dx_z());
  assert(std::abs(z_residual - 1.0) < 1.0e-12);
  std::cout << "  test_anisotropic_yz_residuals: PASSED\n";
}

void test_slab_point_source_invariants() {
  DomainConfig cfg;
  cfg.lo = {0.0, 0.0, 0.0};
  cfg.hi = {25.0e-6, 5.0e-6, 10.0e-6};
  cfg.grid_dx = 5.0e-6;
  cfg.grid_halo_width = 2;
  Domain domain;
  domain.init(cfg);
  ChemicalField chem;
  chem.init(domain, {diffusing_species(1.0e-12, 0.0, 0.0)}, "slab");
  const Int source = domain.cell_index(2, 0, 1);
  chem.conc_global(0, source) = 1.0;
  const Real before = inventory(chem, domain);
  chem.apply_diffusion(domain, 2.5);
  const Real after = inventory(chem, domain);
  const Real alpha = 2.5 * 1.0e-12 / (domain.dx_x() * domain.dx_x());
  assert(std::abs(after * (1.0 + alpha) - before)
         <= 1.0e-12 * std::abs(before));
  for (Int ix = 0; ix < domain.nx(); ++ix) {
    const Int cell = domain.cell_index(ix, 0, 1);
    const Int left = domain.cell_index(
        (ix + domain.nx() - 1) % domain.nx(), 0, 1);
    const Int right = domain.cell_index((ix + 1) % domain.nx(), 0, 1);
    const Real value = chem.conc_global(0, cell);
    const Real laplacian = (chem.conc_global(0, left)
        + chem.conc_global(0, right) - 2.0 * value)
        / (domain.dx_x() * domain.dx_x());
    const Real residual = (1.0 + alpha)
        * (value - 2.5e-12 * laplacian);
    assert(std::abs(residual - (ix == 2 ? 1.0 : 0.0)) < 1.0e-12);
    assert(value >= 0.0);
  }
  std::cout << "  test_slab_point_source_invariants: PASSED\n";
}

void test_uniform_field_is_fixed_point() {
  Domain domain = make_domain(5, 4, 6);
  ChemicalField chem;
  chem.init(domain, {diffusing_species(2.1e-9, 0.25, 0.25)});

  chem.apply_diffusion(domain, 300.0);

  for (Int cell = 0; cell < chem.ncells(); ++cell) {
    assert(std::abs(chem.conc(0, cell) - 0.25) < 1.0e-10);
  }
  std::cout << "  test_uniform_field_is_fixed_point: PASSED\n";
}

void test_point_source_invariants_and_residual() {
  Domain domain = make_domain(5, 1, 2);
  ChemicalField chem;
  chem.init(domain, {diffusing_species(1.0e-12, 0.0, 0.0)});
  chem.conc(0, domain.cell_index(2, 0, 1)) = 1.0;
  const Real before = inventory(chem, domain);

  chem.apply_diffusion(domain, 2.5);

  const Real after = inventory(chem, domain);
  const Real alpha = 2.5 * 1.0e-12 / (domain.dx_x() * domain.dx_x());
  const Real periodic_ring_inventory = after * (1.0 + alpha);
  assert(std::abs(periodic_ring_inventory - before) <=
         1.0e-12 * std::abs(before));
  const Int source = 2;
  for (Int ix = 0; ix < domain.nx(); ++ix) {
    const Real value = chem.conc(0, domain.cell_index(ix, 0, 1));
    const Real left = chem.conc(
        0, domain.cell_index((ix + domain.nx() - 1) % domain.nx(), 0, 1));
    const Real right = chem.conc(
        0, domain.cell_index((ix + 1) % domain.nx(), 0, 1));
    assert(value > 0.0);
    if (ix < source) {
      assert(value < chem.conc(0, domain.cell_index(ix + 1, 0, 1)));
    } else if (ix > source) {
      assert(value < chem.conc(0, domain.cell_index(ix - 1, 0, 1)));
    }
    const Real initial = ix == source ? 1.0 : 0.0;
    const Real laplacian = (left + right - 2.0 * value)
        / (domain.dx_x() * domain.dx_x());
    // The directional z solve follows the periodic x solve and scales this
    // ring by (1 + alpha); undo that factor before checking the x operator.
    const Real residual = (1.0 + alpha)
        * (value - 2.5 * 1.0e-12 * laplacian);
    assert(std::abs(residual - initial) < 1.0e-12);
  }
  const Real center = chem.conc(0, domain.cell_index(source, 0, 1));
  assert(center > chem.conc(0, domain.cell_index(1, 0, 1)));
  assert(std::abs(chem.conc(0, domain.cell_index(1, 0, 1))
                  - chem.conc(0, domain.cell_index(3, 0, 1)))
         < 1.0e-12);
  assert(std::abs(chem.conc(0, domain.cell_index(0, 0, 1))
                  - chem.conc(0, domain.cell_index(4, 0, 1)))
         < 1.0e-12);
  std::cout << "  test_point_source_invariants_and_residual: PASSED\n";
}

ChemicalField diffuse_point_source(const Domain& domain, Real diffusion,
                                   bool enabled) {
  ChemicalSpec spec = diffusing_species(diffusion, 0.0, 0.0);
  spec.diffusion_enabled = enabled;
  ChemicalField chem;
  chem.init(domain, {spec});
  chem.conc(0, domain.cell_index(2, 0, 1)) = 1.0;
  chem.apply_diffusion(domain, 2.5);
  return chem;
}

void test_diffusion_enable_and_coefficient_sensitivity() {
  Domain domain = make_domain(5, 1, 2);
  const ChemicalField disabled = diffuse_point_source(domain, 1.0e-12, false);
  const ChemicalField slow = diffuse_point_source(domain, 5.0e-13, true);
  const ChemicalField fast = diffuse_point_source(domain, 2.0e-12, true);
  const Int center = domain.cell_index(2, 0, 1);
  const Int neighbor = domain.cell_index(1, 0, 1);

  assert(disabled.conc(0, center) == 1.0);
  assert(disabled.conc(0, neighbor) == 0.0);
  assert(fast.conc(0, center) < slow.conc(0, center));
  assert(fast.conc(0, neighbor) > slow.conc(0, neighbor));
  std::cout << "  test_diffusion_enable_and_coefficient_sensitivity: PASSED\n";
}

void test_dirichlet_neumann_boundary_gradient() {
  Domain domain = make_domain(4, 4, 8);
  ChemicalField chem;
  chem.init(domain, {diffusing_species(2.1e-9, 0.0, 1.0)});

  chem.apply_diffusion(domain, 60.0);
  chem.flux_accounting().commit_boundary_and_reaction_step();

  Real previous = chem.conc(0, domain.cell_index(0, 0, 0));
  assert(std::abs(previous - 1.0) < 1.0e-15);
  for (Int iz = 1; iz < domain.nz(); ++iz) {
    const Real current = chem.conc(0, domain.cell_index(0, 0, iz));
    assert(std::isfinite(current));
    assert(current > 0.0);
    assert(current <= previous);
    previous = current;
  }
  assert(previous < 1.0);
  std::cout << "  test_dirichlet_neumann_boundary_gradient: PASSED\n";
}

void test_configured_z_gradient_is_background_fixed_point() {
  Domain domain = make_domain(4, 3, 8);
  ChemicalSpec spec = diffusing_species(5.0e-10, 5.0e-3, 5.0e-3);
  spec.z_gradient_enabled = true;
  spec.z_gradient_lambda = 25.0e-6;
  ChemicalField chem;
  chem.init(domain, {spec});

  chem.apply_diffusion(domain, 60.0);

  assert(std::abs(chem.conc(0, domain.cell_index(0, 0, 0))
                  - spec.boundary_conc) < 1.0e-15);
  for (Int iz = 1; iz < domain.nz() - 1; ++iz) {
    const Real z_rel = (iz + 0.5) * domain.dx_z();
    const Real expected = spec.initial_conc
        * std::exp(-z_rel / spec.z_gradient_lambda);
    for (Int iy = 0; iy < domain.ny(); ++iy) {
      for (Int ix = 0; ix < domain.nx(); ++ix) {
        const Real actual = chem.conc(0, domain.cell_index(ix, iy, iz));
        assert(std::abs(actual - expected) < 1.0e-12);
      }
    }
  }
  const Real top = chem.conc(0, domain.cell_index(0, 0, domain.nz() - 1));
  const Real below = chem.conc(0, domain.cell_index(0, 0, domain.nz() - 2));
  assert(std::abs(top - below) < 1.0e-15);

  chem.conc(0, domain.cell_index(1, 1, 1)) = 0.0;
  chem.apply_diffusion(domain, 300.0);
  for (Int cell = 0; cell < chem.ncells(); ++cell) {
    assert(chem.conc(0, cell) >= 0.0);
  }
  std::cout << "  test_configured_z_gradient_is_background_fixed_point: PASSED\n";
}

void test_large_timestep_is_positive_and_finite() {
  Domain domain = make_domain(6, 5, 7);
  ChemicalField chem;
  chem.init(domain, {diffusing_species(2.1e-9, 0.0, 0.0)});
  for (Int iz = 1; iz < domain.nz(); ++iz) {
    for (Int iy = 0; iy < domain.ny(); ++iy) {
      for (Int ix = 0; ix < domain.nx(); ++ix) {
        chem.conc(0, domain.cell_index(ix, iy, iz)) =
            ((ix + iy + iz) % 2 == 0) ? 1.0 : 0.0;
      }
    }
  }

  chem.apply_diffusion(domain, 300.0);

  for (Int cell = 0; cell < chem.ncells(); ++cell) {
    const Real value = chem.conc(0, cell);
    assert(std::isfinite(value));
    assert(value >= 0.0);
    assert(value <= 1.0 + 1.0e-12);
  }
  std::cout << "  test_large_timestep_is_positive_and_finite: PASSED\n";
}

void test_boundary_accounting_closes_diffusion_only_inventory() {
  Domain domain = make_domain(5, 4, 6);
  ChemicalField chem;
  chem.init(domain, {diffusing_species(2.1e-9, 0.0, 1.0)});
  const Real before = inventory(chem, domain);

  chem.apply_diffusion(domain, 60.0);
  chem.flux_accounting().commit_boundary_and_reaction_step();

  const Real after = inventory(chem, domain);
  const Real recorded = chem.flux_accounting().boundary_interval[0];
  const Real tolerance =
      1.0e-12 * (std::abs(before) + std::abs(after));
  assert(std::abs((after - before) - recorded) < tolerance);
  std::cout << "  test_boundary_accounting_closes_diffusion_only_inventory: PASSED\n";
}

void test_boundary_accounting_is_boundary_concentration_sensitive() {
  Domain domain = make_domain(5, 4, 6);
  std::array<Real, 3> boundary_values = {0.25, 1.0, 4.0};
  std::array<Real, 3> fluxes{};

  for (size_t i = 0; i < boundary_values.size(); ++i) {
    ChemicalField chem;
    chem.init(domain, {diffusing_species(2.1e-9, 0.0, boundary_values[i])});
    chem.apply_diffusion(domain, 60.0);
    chem.flux_accounting().commit_boundary_and_reaction_step();
    fluxes[i] = chem.flux_accounting().boundary_interval[0];
  }

  assert(fluxes[0] < fluxes[1]);
  assert(fluxes[1] < fluxes[2]);
  assert(fluxes[1] - fluxes[0] > 1.0e-16);
  assert(fluxes[2] - fluxes[1] > 1.0e-16);
  std::cout << "  test_boundary_accounting_is_boundary_concentration_sensitive: PASSED\n";
}

void test_gradient_boundary_accounting_is_nonzero() {
  Domain domain = make_domain(5, 4, 6);
  ChemicalSpec spec = diffusing_species(5.0e-10, 1.0e-3, 1.0e-3);
  spec.z_gradient_enabled = true;
  spec.z_gradient_lambda = 10.0e-6;
  ChemicalField chem;
  chem.init(domain, {spec});
  for (Int iz = 1; iz < domain.nz(); ++iz) {
    for (Int iy = 0; iy < domain.ny(); ++iy) {
      for (Int ix = 0; ix < domain.nx(); ++ix) {
        chem.conc(0, domain.cell_index(ix, iy, iz)) = 0.0;
      }
    }
  }

  chem.apply_diffusion(domain, 60.0);
  chem.flux_accounting().commit_boundary_and_reaction_step();

  const Real recorded = chem.flux_accounting().boundary_interval[0];
  assert(std::abs(recorded) > 1.0e-20);
  std::cout << "  test_gradient_boundary_accounting_is_nonzero: PASSED\n";
}

void advance_delivery(ChemicalField& chem, const Domain& domain, Int steps,
                      Real dt) {
  for (Int step = 0; step < steps; ++step) {
    chem.apply_diffusion(domain, dt);
    chem.flux_accounting().commit_boundary_and_reaction_step();
    chem.flux_accounting().close_interval();
  }
}

Real bottom_average(const ChemicalField& chem, const Domain& domain) {
  Real total = 0.0;
  for (Int iy = 0; iy < domain.ny(); ++iy) {
    for (Int ix = 0; ix < domain.nx(); ++ix) {
      total += chem.conc(
          0, domain.cell_index(ix, iy, 0));
    }
  }
  return total / static_cast<Real>(domain.nx() * domain.ny());
}

void test_delivery_boundary_conservation_and_sensitivity() {
  Domain domain = make_domain(4, 3, 6);
  for (const auto mode : {EpithelialBoundaryMode::Robin,
                          EpithelialBoundaryMode::Flux}) {
    ChemicalField chem;
    chem.init(domain, {delivery_species(mode, 0.0, 1.0, 0.0, 0.0)});
    for (Int cell = 0; cell < chem.ncells(); ++cell) {
      chem.conc(0, cell) = 0.1 + 0.01 * static_cast<Real>(cell);
    }
    const Real before = inventory(chem, domain);
    chem.apply_diffusion(domain, 30.0);
    chem.flux_accounting().commit_boundary_and_reaction_step();
    const Real after = inventory(chem, domain);
    assert(std::abs(after - before) < 1.0e-12 * before);
    assert(chem.flux_accounting().boundary_cumulative[0] == 0.0);
  }

  std::array<Real, 3> robin_bottom{};
  std::array<Real, 3> robin_supply{};
  const std::array<Real, 3> transfer_coeffs = {0.0, 2.0e-5, 1.0e-4};
  for (size_t i = 0; i < transfer_coeffs.size(); ++i) {
    ChemicalField chem;
    chem.init(domain, {delivery_species(
        EpithelialBoundaryMode::Robin, 0.0, 1.0,
        transfer_coeffs[i], 0.0)});
    advance_delivery(chem, domain, 50, 1.0);
    robin_bottom[i] = bottom_average(chem, domain);
    robin_supply[i] = chem.flux_accounting().boundary_cumulative[0];
  }
  assert(robin_bottom[0] < robin_bottom[1]);
  assert(robin_bottom[1] < robin_bottom[2]);
  assert(robin_supply[0] < robin_supply[1]);
  assert(robin_supply[1] < robin_supply[2]);

  std::array<Real, 3> flux_bottom{};
  std::array<Real, 3> flux_supply{};
  const std::array<Real, 3> fixed_fluxes = {0.0, 1.0e-10, 2.0e-10};
  for (size_t i = 0; i < fixed_fluxes.size(); ++i) {
    ChemicalField chem;
    chem.init(domain, {delivery_species(
        EpithelialBoundaryMode::Flux, 0.0, 1.0, 0.0, fixed_fluxes[i])});
    advance_delivery(chem, domain, 50, 1.0);
    flux_bottom[i] = bottom_average(chem, domain);
    flux_supply[i] = chem.flux_accounting().boundary_cumulative[0];
  }
  assert(flux_bottom[0] < flux_bottom[1]);
  assert(flux_bottom[1] < flux_bottom[2]);
  assert(flux_supply[0] < flux_supply[1]);
  assert(flux_supply[1] < flux_supply[2]);
  std::cout << "  test_delivery_boundary_conservation_and_sensitivity: PASSED\n";
}

void test_delivery_boundary_limits_depletion_and_accounting() {
  Domain domain = make_domain(4, 3, 6);
  ChemicalField robin;
  robin.init(domain, {delivery_species(
      EpithelialBoundaryMode::Robin, 0.0, 1.0, 1.0e-1, 0.0)});
  robin.apply_diffusion(domain, 1.0);
  assert(std::abs(bottom_average(robin, domain) - 1.0) < 1.0e-3);

  const Real flux = 1.0e-10;
  ChemicalField flux_field;
  flux_field.init(domain, {delivery_species(
      EpithelialBoundaryMode::Flux, 0.0, 1.0, 0.0, flux)});
  advance_delivery(flux_field, domain, 7, 1.0);
  const Real expected = 7.0 * flux * domain.nx() * domain.ny()
      * domain.dx_x() * domain.dx_y();
  assert(std::abs(flux_field.flux_accounting().boundary_cumulative[0]
                  - expected) < 1.0e-24);

  ChemicalField control;
  control.init(domain, {delivery_species(
      EpithelialBoundaryMode::Flux, 0.1, 1.0, 0.0, 1.0e-14)});
  const Real control_before = inventory(control, domain);
  for (Int step = 0; step < 200; ++step) {
    control.apply_diffusion(domain, 1.0);
  }
  const Real control_after = inventory(control, domain);
  assert(control_after >= control_before);
  const Real control_bottom = bottom_average(control, domain);

  ChemicalField depleted;
  depleted.init(domain, {delivery_species(
      EpithelialBoundaryMode::Flux, 0.1, 1.0, 0.0, 1.0e-14)});
  const Int sink_cell = domain.cell_index(1, 1, 2);
  for (Int step = 0; step < 200; ++step) {
    depleted.conc(0, sink_cell) =
        std::max(0.0, depleted.conc(0, sink_cell) - 2.0e-2);
    depleted.apply_diffusion(domain, 1.0);
  }
  assert(control_bottom - bottom_average(depleted, domain) > 0.05);
  std::cout << "  test_delivery_boundary_limits_depletion_and_accounting: PASSED\n";
}

void test_delivery_boundary_slab_matches_replicated() {
  DomainConfig cfg;
  cfg.hi = {40.0e-6, 15.0e-6, 30.0e-6};
  cfg.grid_dx = 5.0e-6;
  cfg.grid_halo_width = 2;
  Domain domain;
  domain.init(cfg);
  for (const auto mode : {EpithelialBoundaryMode::Dirichlet,
                          EpithelialBoundaryMode::Robin,
                          EpithelialBoundaryMode::Flux}) {
    const ChemicalSpec spec = delivery_species(
        mode, 0.0, 1.0, mode == EpithelialBoundaryMode::Robin ? 2.0e-5 : 0.0,
        mode == EpithelialBoundaryMode::Flux ? 1.0e-10 : 0.0);
    ChemicalField replicated;
    ChemicalField slab;
    replicated.init(domain, {spec});
    slab.init(domain, {spec}, "slab");
    for (Int cell = 0; cell < domain.ncells(); ++cell) {
      const Real value = 0.01 * static_cast<Real>(cell + 1);
      replicated.conc_global(0, cell) = value;
      slab.conc_global(0, cell) = value;
    }
    replicated.apply_diffusion(domain, 2.0);
    slab.apply_diffusion(domain, 2.0);
    for (Int cell = 0; cell < domain.ncells(); ++cell) {
      assert(std::abs(replicated.conc_global(0, cell)
                      - slab.conc_global(0, cell)) < 1.0e-12);
    }
  }
  std::cout << "  test_delivery_boundary_slab_matches_replicated: PASSED\n";
}

void test_dirichlet_default_is_unchanged() {
  Domain domain = make_domain(4, 3, 6);
  ChemicalSpec implicit_default = diffusing_species(2.1e-9, 0.0, 1.0);
  ChemicalSpec explicit_dirichlet = implicit_default;
  explicit_dirichlet.epithelial_boundary_mode =
      EpithelialBoundaryMode::Dirichlet;
  ChemicalField default_field;
  ChemicalField explicit_field;
  default_field.init(domain, {implicit_default});
  explicit_field.init(domain, {explicit_dirichlet});
  for (Int cell = 0; cell < domain.ncells(); ++cell) {
    const Real value = 0.01 * static_cast<Real>(cell + 1);
    default_field.conc(0, cell) = value;
    explicit_field.conc(0, cell) = value;
  }
  default_field.apply_diffusion(domain, 2.0);
  explicit_field.apply_diffusion(domain, 2.0);
  for (Int cell = 0; cell < domain.ncells(); ++cell) {
    assert(default_field.conc(0, cell) == explicit_field.conc(0, cell));
  }
  std::cout << "  test_dirichlet_default_is_unchanged: PASSED\n";
}

void test_delivery_boundary_rejects_gradient() {
  Domain domain = make_domain(3, 2, 4);
  ChemicalSpec spec = delivery_species(
      EpithelialBoundaryMode::Robin, 0.0, 1.0, 1.0e-5, 0.0);
  spec.z_gradient_enabled = true;
  bool rejected = false;
  try {
    ChemicalField chem;
    chem.init(domain, {spec});
  } catch (const ConfigError&) {
    rejected = true;
  }
  assert(rejected);
  std::cout << "  test_delivery_boundary_rejects_gradient: PASSED\n";
}

void test_default_species_configuration() {
  const SimulationConfig cfg = InputParser::default_config();
  const auto diffusion_enabled = [&cfg](std::string_view name) {
    const auto it = std::ranges::find_if(
        cfg.chemicals, [name](const ChemicalSpec& spec) { return spec.name == name; });
    assert(it != cfg.chemicals.end());
    return it->diffusion_enabled;
  };

  assert(diffusion_enabled(species::CARBON));
  assert(diffusion_enabled(species::IRON));
  assert(diffusion_enabled(species::B12));
  assert(diffusion_enabled(species::ACETATE));
  assert(diffusion_enabled(species::ETHANOLAMINE));
  assert(!diffusion_enabled(species::BACTERIOCIN_BTUB));
  assert(!diffusion_enabled(species::BACTERIOCIN_FEPA));
  assert(!diffusion_enabled(species::BACTERIOCIN_CIRA));
  assert(!diffusion_enabled(species::BACTERIOCIN_FHUA));
  std::cout << "  test_default_species_configuration: PASSED\n";
}

}  // namespace

int main() {
  std::cout << "=== Nutrient Diffusion Tests ===\n";
  test_uniform_field_is_fixed_point();
  test_anisotropic_diffusion_invariants();
  test_stride_sensitivity_is_ordered();
  test_anisotropic_yz_residuals();
  test_slab_point_source_invariants();
  test_point_source_invariants_and_residual();
  test_diffusion_enable_and_coefficient_sensitivity();
  test_dirichlet_neumann_boundary_gradient();
  test_configured_z_gradient_is_background_fixed_point();
  test_large_timestep_is_positive_and_finite();
  test_boundary_accounting_closes_diffusion_only_inventory();
  test_boundary_accounting_is_boundary_concentration_sensitive();
  test_gradient_boundary_accounting_is_nonzero();
  test_delivery_boundary_conservation_and_sensitivity();
  test_delivery_boundary_limits_depletion_and_accounting();
  test_delivery_boundary_slab_matches_replicated();
  test_dirichlet_default_is_unchanged();
  test_delivery_boundary_rejects_gradient();
  test_delivery_mass_closure_gradient_parameterization();
  test_delivery_step_mass_closure_dirichlet_refill();
  test_delivery_step_mass_closure_vbf_source();
  test_default_species_configuration();
  std::cout << "All nutrient diffusion tests passed.\n";
  return 0;
}
