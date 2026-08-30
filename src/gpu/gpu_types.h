#ifndef GUTIBM_GPU_TYPES_H
#define GUTIBM_GPU_TYPES_H

#include <array>

namespace gutibm::gpu {

struct DomainParams {
  int nx;
  int ny;
  int nz;
  double dx_x;
  double dx_y;
  double dx_z;
  std::array<double, 3> lo;
  std::array<double, 3> extent;
  std::array<bool, 3> periodic;
  double z_lo;
  double z_hi;
};

struct AdvectionParams {
  double v_radial_max;
  double v_distal_max;
  double h;
  double lo_z;
  double profile_alpha;
  bool crypts_enabled;
  double crypt_depth;
  bool peristaltic_enabled;
  double peristaltic_period;
  double peristaltic_amplitude;
  double peristaltic_wavelength;
  double current_time;
};

struct GfSourceParams {
  double diff_coeff;
  double source_rate;
  double retardation;
  double decay_rate;
  double lumen_transfer_length;
  double robin_cutoff;
  int lumen_transfer_basis_free;
  int robin_table_index;
  double image_series_relative_tolerance;
  int image_series_max_shells;
  int image_series_legacy_reflections;
};

struct VbfLaunchParams {
  int storage_nx = 0;
  int owned_x_begin = 0;
  int owned_x_end = 0;
  int global_nx = 0;
  int ny = 0;
  int nz = 0;
  int owned_global_x_begin = 0;
  int owned_global_x_end = 0;
  double dx_x = 0.0;
  double dx_y = 0.0;
  double dx_z = 0.0;
  double nutrient_sink = 0.0;
  double carbon_sink_vmax = 0.0;
  double carbon_sink_km = 0.0;
  double agent_carbon_coupling = 0.0;
  int* agent_counts = nullptr;
  int use_dynamic_mucin = 0;
  int mucin_z_gradient_enabled = 0;
  double mucin_z_gradient_lambda = 0.0;
  double mucin_liberation = 0.0;
  double vbf_density = 0.0;
  int oxygen_enabled = 0;
  int oxygen_delivery_enabled = 0;
  double oxygen_vbf_sink = 0.0;
  int acetate_enabled = 0;
  double acetate_vbf_production = 0.0;
  double acetate_vbf_consumption = 0.0;
  double acetate_epithelial_uptake = 0.0;
  int mucin_enabled = 0;
  double mucin_secretion_rate = 0.0;
  double mucin_Km_degradation = 0.0;
  double mucin_k_liberation = 0.0;
};

struct AgentUptakeDevice {
  double* totals = nullptr;
};

struct MechanicsLaunchParams {
  double hertz_k = 0.0;
  int hertzian_enabled = 0;
  int adhesion_enabled = 0;
  double adhesion_strength = 0.0;
  double adhesion_range = 0.0;
  double viscosity = 0.0;
  double max_displacement_fraction = 0.0;
  double dt = 0.0;
  double sim_time = 0.0;
  double corpse_persistence = 0.0;
  int cdi_enabled = 0;
  double lo0 = 0.0;
  double lo1 = 0.0;
  double lo2 = 0.0;
  double hi0 = 0.0;
  double hi1 = 0.0;
  double hi2 = 0.0;
  int periodic_x = 0;
  int periodic_y = 0;
  int periodic_z = 0;
  double cell_size = 0.0;
  int nx_cells = 0;
  int ny_cells = 0;
  int nz_cells = 0;
};

}  // namespace gutibm::gpu

#endif  // GUTIBM_GPU_TYPES_H
