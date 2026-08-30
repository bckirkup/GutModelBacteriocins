#include "chemical_field_gpu.h"
#include "chemical_field.h"
#include "cuda_aware_mpi.h"
#include "diffusion_gpu.h"
#include "domain.h"
#include "dispatch.h"
#include "gpu_kernels.h"
#include "device_memory.h"
#include "species_names.h"
#include "carbon_maintenance.h"

#include <cstdlib>
#include <iostream>

#ifdef GUTIBM_CUDA
#include <cuda_runtime.h>
#endif

#ifdef GUTIBM_MPI
#include <mpi.h>
#endif

namespace gutibm {

void ChemicalFieldGpu::init(ChemicalField& field) {
  active_ = gpu_runtime_enabled();
  nspec_ = field.num_species();
  ncells_ = field.ncells();
  global_nx_ = field.global_nx();
  global_ny_ = field.global_ny();
  global_nz_ = field.global_nz();
  owned_x_begin_ = field.owned_global_x_begin();
  owned_x_end_ = field.owned_global_x_end();
  halo_width_ = field.grid_halo_width();
  storage_nx_ = field.storage_nx();
  slab_mode_ = field.slab_mode();
  diffusion_fallback_warning_emitted_ = false;
  if (!active_ || nspec_ <= 0 || ncells_ <= 0) return;

  d_conc_.resize(static_cast<size_t>(nspec_));
  d_reac_.resize(static_cast<size_t>(nspec_));
  for (Int s = 0; s < nspec_; ++s) {
    d_conc_[static_cast<size_t>(s)].allocate(static_cast<size_t>(ncells_));
    d_reac_[static_cast<size_t>(s)].allocate(static_cast<size_t>(ncells_));
  }

  std::vector<double> bc(static_cast<size_t>(nspec_));
  for (Int s = 0; s < nspec_; ++s) {
    bc[static_cast<size_t>(s)] = field.spec(s).boundary_conc;
  }
  d_boundary_conc_.upload(bc);
  d_boundary_injected_.allocate(static_cast<size_t>(nspec_));
  d_agent_uptake_.allocate(4);
  d_maintenance_available_.allocate(static_cast<size_t>(ncells_));
  d_uptake_limit_totals_.allocate(5);
  d_vbf_totals_.allocate(4);
  d_reaction_clip_.allocate(static_cast<size_t>(nspec_));
  d_agent_counts_.allocate(static_cast<size_t>(ncells_));
  sync_to_device(field);
}

void ChemicalFieldGpu::reset_agent_uptake() {
  if (!active_) return;
  d_agent_uptake_.upload(std::vector(4, 0.0));
}

void ChemicalFieldGpu::download_agent_uptake(ChemicalField& field) const {
  if (!active_) return;
  gpu_sync_compute();
  std::vector values(4, 0.0);
  d_agent_uptake_.download(values);
  const Int carbon = field.find(species::CARBON);
  const Int iron = field.find(species::IRON);
  if (carbon >= 0) field.flux_accounting().add_agent_uptake(carbon, values[0]);
  if (iron >= 0) field.flux_accounting().add_agent_uptake(iron, values[1]);
  if (carbon >= 0) {
    field.flux_accounting().add_maintenance(carbon, values[2]);
    field.flux_accounting().add_maintenance_shortfall(carbon, values[3]);
  }
}

void ChemicalFieldGpu::prepare_maintenance(
    const ChemicalField& field, Int carbon, Real cell_volume) {
  if (!active_ || carbon < 0 || cell_volume <= 0.0) return;
  std::vector available(static_cast<size_t>(ncells_), 0.0);
  for (Int cell = 0; cell < ncells_; ++cell) {
    available[static_cast<size_t>(cell)] = carbon_maintenance::available(
        field.conc(carbon, cell), cell_volume);
  }
  d_maintenance_available_.upload(available);
}

void ChemicalFieldGpu::reset_uptake_limit_totals() {
  if (!active_) return;
  d_uptake_limit_totals_.upload(std::vector(5, 0.0));
}

void ChemicalFieldGpu::download_uptake_limit_totals(ChemicalField& field) const {
  if (!active_) return;
  gpu_sync_compute();
  std::vector values(5, 0.0);
  d_uptake_limit_totals_.download(values);
  const Int carbon = field.find(species::CARBON);
  const Int iron = field.find(species::IRON);
  if (carbon >= 0) {
    field.flux_accounting().add_uptake_demand(carbon, values[0]);
    field.flux_accounting().add_uptake_limited(carbon, values[2]);
    field.flux_accounting().add_maintenance_limited_agents(carbon, values[4]);
  }
  if (iron >= 0) {
    field.flux_accounting().add_uptake_demand(iron, values[1]);
    field.flux_accounting().add_uptake_limited(iron, values[3]);
  }
}

void ChemicalFieldGpu::reset_vbf_totals() {
  if (!active_) return;
  d_vbf_totals_.upload(std::vector(4, 0.0));
}

double* ChemicalFieldGpu::vbf_totals_device() {
  return active_ ? d_vbf_totals_.data() : nullptr;
}

void ChemicalFieldGpu::download_vbf_totals(
    std::vector<double>& values) const {
  if (!active_) return;
  d_vbf_totals_.download(values);
}

void ChemicalFieldGpu::reset_agent_counts() {
#ifndef GUTIBM_CUDA
  return;
#else
  if (!active_) return;
  cudaMemsetAsync(d_agent_counts_.data(), 0,
                  static_cast<size_t>(ncells_) * sizeof(int),
                  gpu_compute_stream());
#endif
}

int* ChemicalFieldGpu::agent_counts_device() {
  return active_ ? d_agent_counts_.data() : nullptr;
}

void ChemicalFieldGpu::reset_reaction_clip() {
#ifndef GUTIBM_CUDA
  return;
#else
  if (!active_) return;
  cudaMemset(d_reaction_clip_.data(), 0,
             static_cast<size_t>(nspec_) * sizeof(double));
#endif
}

void ChemicalFieldGpu::download_reaction_clip(ChemicalField& field) const {
#ifndef GUTIBM_CUDA
  (void)field;
  return;
#else
  if (!active_) return;
  gpu_sync_compute();
  std::vector values(static_cast<size_t>(nspec_), 0.0);
  d_reaction_clip_.download(values);
  for (Int s = 0; s < nspec_; ++s) {
    field.flux_accounting().add_reaction_clip(
        s, values[static_cast<size_t>(s)]);
  }
#endif
}

void ChemicalFieldGpu::sync_to_device(const ChemicalField& field) {
  sync_concentrations_to_device(field);
  sync_reactions_to_device(field);
}

void ChemicalFieldGpu::sync_to_host(ChemicalField& field) {
  sync_concentrations_to_host(field);
  sync_reactions_to_host(field);
}

void ChemicalFieldGpu::sync_concentrations_to_device(const ChemicalField& field) {
  if (!active_) return;
  for (Int s = 0; s < nspec_; ++s) {
    d_conc_[static_cast<size_t>(s)].upload(
        field.conc_data()[static_cast<size_t>(s)]);
  }
}

void ChemicalFieldGpu::sync_reactions_to_device(const ChemicalField& field) {
  if (!active_) return;
  for (Int s = 0; s < nspec_; ++s) {
    std::vector<double> host(static_cast<size_t>(ncells_));
    for (Int c = 0; c < ncells_; ++c) {
      host[static_cast<size_t>(c)] = field.reac(s, c);
    }
    d_reac_[static_cast<size_t>(s)].upload(host);
  }
}

void ChemicalFieldGpu::sync_concentrations_to_host(ChemicalField& field) {
  if (!active_) return;
  for (Int s = 0; s < nspec_; ++s) {
    std::vector<double> host(static_cast<size_t>(ncells_));
    d_conc_[static_cast<size_t>(s)].download(host);
    for (Int c = 0; c < ncells_; ++c) {
      field.conc(s, c) = host[static_cast<size_t>(c)];
    }
  }
}

void ChemicalFieldGpu::sync_reactions_to_host(ChemicalField& field) {
  if (!active_) return;
  for (Int s = 0; s < nspec_; ++s) {
    std::vector<double> host(static_cast<size_t>(ncells_));
    d_reac_[static_cast<size_t>(s)].download(host);
    for (Int c = 0; c < ncells_; ++c) {
      field.reac(s, c) = host[static_cast<size_t>(c)];
    }
  }
}

void ChemicalFieldGpu::accumulate_reactions_to_host(ChemicalField& field) {
  if (!active_) return;
  for (Int s = 0; s < nspec_; ++s) {
    std::vector<double> device_reactions(static_cast<size_t>(ncells_));
    d_reac_[static_cast<size_t>(s)].download(device_reactions);
    for (Int c = 0; c < ncells_; ++c) {
      field.reac(s, c) += device_reactions[static_cast<size_t>(c)];
    }
  }
  zero_reactions_on_device();
}

void ChemicalFieldGpu::sync_species_concentrations_to_host(ChemicalField& field,
                                                           Int spec) {
  if (!active_ || spec < 0 || spec >= nspec_) return;
  std::vector<double> host(static_cast<size_t>(ncells_));
  d_conc_[static_cast<size_t>(spec)].download(host);
  for (Int c = 0; c < ncells_; ++c) {
    field.conc(spec, c) = host[static_cast<size_t>(c)];
  }
}

void ChemicalFieldGpu::sync_species_concentrations_to_device(
    const ChemicalField& field, Int spec) {
  if (!active_ || spec < 0 || spec >= nspec_) return;
  d_conc_[static_cast<size_t>(spec)].upload(
      field.conc_data()[static_cast<size_t>(spec)]);
}

void ChemicalFieldGpu::zero_reactions_on_device() {
#ifndef GUTIBM_CUDA
  return;
#else
  if (!active_) return;
  for (Int s = 0; s < nspec_; ++s) {
    cudaMemset(d_reac_[static_cast<size_t>(s)].data(), 0,
               static_cast<size_t>(ncells_) * sizeof(double));
  }
#endif
}

bool ChemicalFieldGpu::apply_reactions(double dt, const Domain& domain) {
  (void)domain;
#ifndef GUTIBM_CUDA
  (void)dt;
  return false;
#else
  if (!active_) return false;
  reset_reaction_clip();
  for (Int s = 0; s < nspec_; ++s) {
    gpu::launch_field_update_kernel(
        d_conc_[static_cast<size_t>(s)].data(),
        d_reac_[static_cast<size_t>(s)].data(),
        ncells_, 1, dt,
        d_reaction_clip_.data() + static_cast<size_t>(s),
        domain.cell_volume(), storage_nx_,
        global_ny_, global_nz_,
        slab_mode_ ? halo_width_ : 0,
        slab_mode_ ? halo_width_ + (owned_x_end_ - owned_x_begin_) : global_nx_,
        gpu_compute_stream());
  }
  gpu_sync_compute();
  gpu_check_error("field_update_kernel");
  return true;
#endif
}

bool ChemicalFieldGpu::apply_diffusion(const Domain& domain,
                                       ChemicalField& field,
                                       Real dt) {
#ifndef GUTIBM_CUDA
  (void)domain;
  (void)field;
  (void)dt;
  return false;
#else
  if (!active_) return false;

  const int max_line = gpu::diffusion_max_line_length();
  if (!diffusion_all_species_within(domain, field, max_line)) {
    if (!diffusion_fallback_warning_emitted_) {
      diffusion_fallback_warning_emitted_ = true;
      if (domain.rank() == 0) {
        for (Int s = 0; s < nspec_; ++s) {
          const ChemicalSpec& spec = field.spec(s);
          if (!spec.diffuses() || diffusion_line_lengths_within(
                                    domain, spec.epithelial_boundary_mode,
                                    max_line)) {
            continue;
          }
          const int z_line = spec.epithelial_boundary_mode
                  == EpithelialBoundaryMode::Dirichlet
              ? domain.nz() - 1 : domain.nz();
          std::cerr << "Warning: GPU diffusion fallback for species '"
                    << spec.name << "': line length exceeds cap"
                    << " (nx=" << domain.nx()
                    << ", ny=" << domain.ny()
                    << ", z_line=" << z_line
                    << ", nz=" << domain.nz()
                    << ", cap=" << max_line
                    << ", boundary_mode=";
          switch (spec.epithelial_boundary_mode) {
            case EpithelialBoundaryMode::Dirichlet:
              std::cerr << "Dirichlet";
              break;
            case EpithelialBoundaryMode::Robin:
              std::cerr << "Robin";
              break;
            case EpithelialBoundaryMode::Flux:
              std::cerr << "Flux";
              break;
          }
          std::cerr << ")\n";
          break;
        }
      }
    }
    return false;
  }

  bool applied = false;
  std::vector zero(static_cast<size_t>(nspec_), 0.0);
  d_boundary_injected_.upload(zero);
  const Int owned_storage_begin = owned_storage_x_begin();
  const Int owned_storage_end = owned_storage_x_end();
  for (Int s = 0; s < nspec_; ++s) {
    bool species_applied = false;
    if (slab_mode_) {
      SlabDiffusionContext context{
          field, s, storage_nx_, owned_storage_begin, owned_storage_end};
      species_applied = gpu_apply_species_diffusion_slab_device(
          domain, field.spec(s), d_conc_[static_cast<size_t>(s)].data(),
          d_boundary_injected_.data() + s, dt, context);
    } else {
      species_applied = gpu_apply_species_diffusion_device(
          domain, field.spec(s), d_conc_[static_cast<size_t>(s)].data(),
          d_boundary_injected_.data() + s, dt);
    }
    if (species_applied) {
      applied = true;
    }
  }

  if (applied) {
    gpu_sync_compute();
    gpu_check_error("ChemicalFieldGpu::apply_diffusion");
    std::vector injected(static_cast<size_t>(nspec_), 0.0);
    d_boundary_injected_.download(injected);
    for (Int s = 0; s < nspec_; ++s) {
      field.flux_accounting().add_boundary(s, injected[static_cast<size_t>(s)]);
    }
  }
  return applied;
#endif
}

bool ChemicalFieldGpu::apply_boundaries(const Domain& domain,
                                        ChemicalField& field) {
#ifndef GUTIBM_CUDA
  (void)domain;
  (void)field;
  return false;
#else
  if (!active_) return false;

  const int nx = storage_nx_;
  const int ny = domain.ny();
  const int nz = domain.nz();
  const Int owned_storage_begin = owned_storage_x_begin();
  const Int owned_storage_end = owned_storage_x_end();

  std::vector zero(static_cast<size_t>(nspec_), 0.0);
  d_boundary_injected_.upload(zero);
  for (Int s = 0; s < nspec_; ++s) {
    const ChemicalSpec& spec = field.spec(s);
    double* d_conc = d_conc_[static_cast<size_t>(s)].data();
    if (spec.epithelial_boundary_mode
        == EpithelialBoundaryMode::Dirichlet) {
      gpu::launch_set_epithelial_boundary(
          d_conc, nx, ny, owned_storage_begin, owned_storage_end,
          spec.boundary_conc, domain.cell_volume(),
          d_boundary_injected_.data() + s, gpu_compute_stream());
    }
    if (!spec.diffusion_enabled && nz >= 2) {
      gpu::launch_set_luminal_neumann(
          d_conc, nx, ny, nz,
          owned_storage_begin, owned_storage_end,
          gpu_compute_stream());
    }
  }

  gpu_sync_compute();
  gpu_check_error("ChemicalFieldGpu::apply_boundaries");
  std::vector injected(static_cast<size_t>(nspec_), 0.0);
  d_boundary_injected_.download(injected);
  for (Int s = 0; s < nspec_; ++s) {
    field.flux_accounting().add_boundary(s, injected[static_cast<size_t>(s)]);
  }
  return true;
#endif
}

double* ChemicalFieldGpu::conc_device(Int spec) {
  if (!active_ || spec < 0 || spec >= nspec_) return nullptr;
  return d_conc_[static_cast<size_t>(spec)].data();
}

const double* ChemicalFieldGpu::conc_device(Int spec) const {
  if (!active_ || spec < 0 || spec >= nspec_) return nullptr;
  return d_conc_[static_cast<size_t>(spec)].data();
}

double* ChemicalFieldGpu::reac_device(Int spec) {
  if (!active_ || spec < 0 || spec >= nspec_) return nullptr;
  return d_reac_[static_cast<size_t>(spec)].data();
}

const double* ChemicalFieldGpu::reac_device(Int spec) const {
  if (!active_ || spec < 0 || spec >= nspec_) return nullptr;
  return d_reac_[static_cast<size_t>(spec)].data();
}

bool ChemicalFieldGpu::try_sum_reactions_on_device(ChemicalField& field) {
#ifndef GUTIBM_CUDA
  (void)field;
  return false;
#else
  (void)field;
  if (!active_) return false;
  if (slab_mode_) return false;

#ifdef GUTIBM_MPI
  int initialized = 0;
  int finalized = 0;
  MPI_Initialized(&initialized);
  MPI_Finalized(&finalized);
  if (!initialized || finalized) return false;

  int ranks = 1;
  MPI_Comm_size(MPI_COMM_WORLD, &ranks);
  if (ranks <= 1) return false;
#endif

  if (const char* env = std::getenv("GUTIBM_CUDA_AWARE_MPI");
      env == nullptr || (env[0] != '1' && env[0] != 't' && env[0] != 'T')) {
    return false;
  }
  if (!cuda_aware_mpi_runtime_available()) {
    return false;
  }

  for (Int s = 0; s < nspec_; ++s) {
    double* d_reac = reac_device(s);
    if (d_reac == nullptr) return false;
#ifdef GUTIBM_MPI
    MPI_Allreduce(MPI_IN_PLACE, d_reac, ncells_, MPI_DOUBLE, MPI_SUM,
                  MPI_COMM_WORLD);
#endif
  }

  return true;
#endif
}

}  // namespace gutibm
