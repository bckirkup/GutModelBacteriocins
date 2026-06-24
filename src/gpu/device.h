#ifndef GUTIBM_DEVICE_H
#define GUTIBM_DEVICE_H

#include <string>

namespace gutibm {

class DeviceContext {
 public:
  DeviceContext() = default;

  // Select and initialize a CUDA device. Returns false if unavailable.
  bool init(int device_id);

  bool active() const { return active_; }
  int  device_id() const { return device_id_; }

  static int device_count();
  static std::string last_error();

 private:
  bool active_    = false;
  int  device_id_ = 0;
};

}  // namespace gutibm

#endif  // GUTIBM_DEVICE_H
