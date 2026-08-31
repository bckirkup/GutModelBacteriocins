#ifndef GUTIBM_VBF_FIRST_ORDER_SINK_H
#define GUTIBM_VBF_FIRST_ORDER_SINK_H

#ifdef __CUDACC__
#define GUTIBM_VBF_HOST_DEVICE __host__ __device__
#else
#define GUTIBM_VBF_HOST_DEVICE
#endif

namespace gutibm::vbf {

GUTIBM_VBF_HOST_DEVICE inline double implicit_first_order_sink(
    double concentration, double rate, double dt) {
  if (concentration <= 0.0 || rate <= 0.0 || dt <= 0.0) return 0.0;
  return concentration * rate / (1.0 + rate * dt);
}

}  // namespace gutibm::vbf

#undef GUTIBM_VBF_HOST_DEVICE

#endif  // GUTIBM_VBF_FIRST_ORDER_SINK_H
