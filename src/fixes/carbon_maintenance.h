#ifndef GUTIBM_CARBON_MAINTENANCE_H
#define GUTIBM_CARBON_MAINTENANCE_H

#ifdef __CUDACC__
#define GUTIBM_CARBON_MAINTENANCE_HOST_DEVICE __host__ __device__
#else
#define GUTIBM_CARBON_MAINTENANCE_HOST_DEVICE
#endif

namespace gutibm::carbon_maintenance {

GUTIBM_CARBON_MAINTENANCE_HOST_DEVICE inline double requested(
    double rate, double biomass, double dt) {
  if (rate <= 0.0 || biomass <= 0.0 || dt <= 0.0) return 0.0;
  return rate * biomass * dt;
}

GUTIBM_CARBON_MAINTENANCE_HOST_DEVICE inline double available(
    double concentration, double cell_volume) {
  if (concentration <= 0.0 || cell_volume <= 0.0) return 0.0;
  return concentration * cell_volume;
}

GUTIBM_CARBON_MAINTENANCE_HOST_DEVICE inline double realized(
    double requested_amount, double concentration, double cell_volume) {
  const double available_amount = available(concentration, cell_volume);
  return requested_amount < available_amount
      ? requested_amount : available_amount;
}

}  // namespace gutibm::carbon_maintenance

#undef GUTIBM_CARBON_MAINTENANCE_HOST_DEVICE

#endif  // GUTIBM_CARBON_MAINTENANCE_H
