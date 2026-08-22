#ifndef GUTIBM_METABOLIC_MODE_H
#define GUTIBM_METABOLIC_MODE_H

#ifdef __CUDACC__
#define GUTIBM_METABOLIC_MODE_HOST_DEVICE __host__ __device__
#else
#define GUTIBM_METABOLIC_MODE_HOST_DEVICE
#endif

#include <cmath>

namespace gutibm::metabolic_mode {

GUTIBM_METABOLIC_MODE_HOST_DEVICE inline double clamp01(double value) {
  return value < 0.0 ? 0.0 : (value > 1.0 ? 1.0 : value);
}

GUTIBM_METABOLIC_MODE_HOST_DEVICE inline double monod(
    double concentration, double half_saturation) {
  if (concentration <= 0.0) return 0.0;
  if (half_saturation <= 0.0) return 1.0;
  return concentration / (half_saturation + concentration);
}

GUTIBM_METABOLIC_MODE_HOST_DEVICE inline double oxygen_availability(
    double concentration, double half_saturation) {
  return monod(concentration, half_saturation);
}

GUTIBM_METABOLIC_MODE_HOST_DEVICE inline double respiratory_capacity(
    double availability, double mu, double mu_crit) {
  const double denominator = mu > mu_crit ? mu : mu_crit;
  if (denominator <= 0.0) return 0.0;
  return clamp01(availability * mu_crit / denominator);
}

GUTIBM_METABOLIC_MODE_HOST_DEVICE inline double fermentation_fraction(
    double oxygen_availability, double mu, double mu_crit) {
  return clamp01(1.0 - respiratory_capacity(
      oxygen_availability, mu, mu_crit));
}

GUTIBM_METABOLIC_MODE_HOST_DEVICE inline double interpolate(
    double aerobic, double anaerobic, double fraction) {
  const double f = clamp01(fraction);
  return aerobic + f * (anaerobic - aerobic);
}

GUTIBM_METABOLIC_MODE_HOST_DEVICE inline double relax(
    double realized, double instantaneous, double dt, double tau) {
  if (dt <= 0.0 || tau <= 0.0) return clamp01(realized);
  const double alpha = 1.0 - exp(-dt / tau);
  return clamp01(realized + alpha * (instantaneous - realized));
}

GUTIBM_METABOLIC_MODE_HOST_DEVICE inline double undissociated_fraction(
    double ph, double pka) {
  return clamp01(1.0 / (1.0 + exp(2.302585092994046 * (ph - pka))));
}

GUTIBM_METABOLIC_MODE_HOST_DEVICE inline double acid_inhibition(
    double total_acetate, double ph, double pka, double ki, double maximum) {
  const double undissociated = total_acetate
      * undissociated_fraction(ph, pka);
  if (undissociated <= 0.0 || ki <= 0.0) return 0.0;
  return clamp01(maximum * undissociated / (ki + undissociated));
}

}  // namespace gutibm::metabolic_mode

#undef GUTIBM_METABOLIC_MODE_HOST_DEVICE

#endif  // GUTIBM_METABOLIC_MODE_H
