/* -----------------------------------------------------------------------
   GutIBM – shared scalar/vector HDF5 checkpoint helpers

   Small read/write utilities shared by the Layer-2 segment checkpoint
   (src/layer2/segment_io.cpp) and the Layer-3 chain checkpoint
   (src/layer3/chain_io.cpp). Header-only so both layers persist the
   same layouts without duplicated boilerplate. Only compiled under
   GUTIBM_HDF5.
   ----------------------------------------------------------------------- */

#ifndef GUTIBM_H5_UTIL_H
#define GUTIBM_H5_UTIL_H

#include "error.h"
#include "random.h"

#ifdef GUTIBM_HDF5
extern "C" {
#include <hdf5.h>
}

#include <algorithm>
#include <sstream>
#include <string>
#include <vector>

namespace gutibm::h5 {

constexpr hsize_t kRngStrLen = 16384;

inline std::string rng_to_string(const RNG& rng) {
  std::ostringstream oss;
  oss << rng.engine();
  return oss.str();
}

inline void rng_from_string(const char* buf, RNG& rng) {
  std::istringstream iss{std::string(buf)};
  iss >> rng.engine();
}

struct GroupGuard {
  hid_t g;
  ~GroupGuard() { H5Gclose(g); }
};

struct FileGuard {
  hid_t f;
  ~FileGuard() { H5Fclose(f); }
};

// Creates `name` plus any missing parent groups so a driver can nest
// state under paths like /regions/<i>/segment without pre-building
// the tree.
inline void make_group(hid_t fid, const std::string& name) {
  std::string cur;
  std::istringstream parts(name);
  std::string part;
  while (std::getline(parts, part, '/')) {
    if (part.empty()) continue;
    cur += "/" + part;
    if (H5Lexists(fid, cur.c_str(), H5P_DEFAULT) > 0) continue;
    GroupGuard g{H5Gcreate2(fid, cur.c_str(), H5P_DEFAULT, H5P_DEFAULT,
                            H5P_DEFAULT)};
    if (g.g < 0) {
      throw HDF5Error("checkpoint cannot create group: " + cur);
    }
  }
}

inline void write_scalar_i64(hid_t fid, const char* name, int64_t v) {
  hid_t space = H5Screate(H5S_SCALAR);
  hid_t ds = H5Dcreate2(fid, name, H5T_NATIVE_INT64, space, H5P_DEFAULT,
                        H5P_DEFAULT, H5P_DEFAULT);
  if (ds < 0 || H5Dwrite(ds, H5T_NATIVE_INT64, H5S_ALL, H5S_ALL,
                         H5P_DEFAULT, &v) < 0) {
    throw HDF5Error(std::string("checkpoint write failed: ") + name);
  }
  H5Dclose(ds);
  H5Sclose(space);
}

inline void write_scalar_f64(hid_t fid, const char* name, double v) {
  hid_t space = H5Screate(H5S_SCALAR);
  hid_t ds = H5Dcreate2(fid, name, H5T_NATIVE_DOUBLE, space, H5P_DEFAULT,
                        H5P_DEFAULT, H5P_DEFAULT);
  if (ds < 0 || H5Dwrite(ds, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
                         H5P_DEFAULT, &v) < 0) {
    throw HDF5Error(std::string("checkpoint write failed: ") + name);
  }
  H5Dclose(ds);
  H5Sclose(space);
}

template <typename T>
void write_vec(hid_t fid, const char* name, hid_t type,
               const std::vector<T>& v) {
  hsize_t dims[1] = {v.empty() ? 0 : v.size()};
  hid_t space = H5Screate_simple(1, dims, nullptr);
  hid_t ds = H5Dcreate2(fid, name, type, space, H5P_DEFAULT, H5P_DEFAULT,
                        H5P_DEFAULT);
  if (ds < 0) {
    H5Sclose(space);
    throw HDF5Error(std::string("checkpoint write failed: ") + name);
  }
  if (!v.empty() &&
      H5Dwrite(ds, type, H5S_ALL, H5S_ALL, H5P_DEFAULT, v.data()) < 0) {
    H5Dclose(ds);
    H5Sclose(space);
    throw HDF5Error(std::string("checkpoint write failed: ") + name);
  }
  H5Dclose(ds);
  H5Sclose(space);
}

inline void write_str_vec(hid_t fid, const char* name,
                          const std::vector<std::string>& v) {
  hsize_t dims[1] = {v.empty() ? 0 : v.size()};
  hid_t type = H5Tcopy(H5T_C_S1);
  H5Tset_size(type, kRngStrLen);
  hid_t space = H5Screate_simple(1, dims, nullptr);
  hid_t ds = H5Dcreate2(fid, name, type, space, H5P_DEFAULT, H5P_DEFAULT,
                        H5P_DEFAULT);
  if (ds < 0) {
    H5Sclose(space);
    H5Tclose(type);
    throw HDF5Error(std::string("checkpoint write failed: ") + name);
  }
  if (!v.empty()) {
    std::vector<char> flat(v.size() * kRngStrLen, '\0');
    for (size_t i = 0; i < v.size(); ++i) {
      std::copy_n(v[i].data(),
                  std::min(v[i].size(), static_cast<size_t>(kRngStrLen)),
                  flat.data() + i * kRngStrLen);
    }
    H5Dwrite(ds, type, H5S_ALL, H5S_ALL, H5P_DEFAULT, flat.data());
  }
  H5Dclose(ds);
  H5Sclose(space);
  H5Tclose(type);
}

inline int64_t read_scalar_i64(hid_t fid, const char* name) {
  hid_t ds = H5Dopen2(fid, name, H5P_DEFAULT);
  if (ds < 0) {
    throw HDF5Error(std::string("checkpoint missing dataset: ") + name);
  }
  int64_t v = 0;
  H5Dread(ds, H5T_NATIVE_INT64, H5S_ALL, H5S_ALL, H5P_DEFAULT, &v);
  H5Dclose(ds);
  return v;
}

// Optional dataset: absent in checkpoints written before the field
// existed (the field then takes `fallback`).
inline int64_t read_scalar_i64_opt(hid_t fid, const std::string& name,
                                 int64_t fallback) {
  if (H5Lexists(fid, name.c_str(), H5P_DEFAULT) <= 0) {
    return fallback;
  }
  return read_scalar_i64(fid, name.c_str());
}

inline double read_scalar_f64(hid_t fid, const char* name) {
  hid_t ds = H5Dopen2(fid, name, H5P_DEFAULT);
  if (ds < 0) {
    throw HDF5Error(std::string("checkpoint missing dataset: ") + name);
  }
  double v = 0;
  H5Dread(ds, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &v);
  H5Dclose(ds);
  return v;
}

template <typename T>
std::vector<T> read_vec(hid_t fid, const char* name, hid_t type) {
  hid_t ds = H5Dopen2(fid, name, H5P_DEFAULT);
  if (ds < 0) {
    throw HDF5Error(std::string("checkpoint missing dataset: ") + name);
  }
  hid_t space = H5Dget_space(ds);
  hsize_t dims[1] = {0};
  H5Sget_simple_extent_dims(space, dims, nullptr);
  std::vector<T> v(dims[0]);
  if (dims[0] > 0) {
    H5Dread(ds, type, H5S_ALL, H5S_ALL, H5P_DEFAULT, v.data());
  }
  H5Sclose(space);
  H5Dclose(ds);
  return v;
}

inline std::vector<std::string> read_str_vec(hid_t fid,
                                             const char* name) {
  hid_t ds = H5Dopen2(fid, name, H5P_DEFAULT);
  if (ds < 0) {
    throw HDF5Error(std::string("checkpoint missing dataset: ") + name);
  }
  hid_t space = H5Dget_space(ds);
  hsize_t dims[1] = {0};
  H5Sget_simple_extent_dims(space, dims, nullptr);
  hid_t type = H5Tcopy(H5T_C_S1);
  H5Tset_size(type, kRngStrLen);
  std::vector<std::string> out(dims[0]);
  if (dims[0] > 0) {
    std::vector<char> flat(dims[0] * kRngStrLen);
    H5Dread(ds, type, H5S_ALL, H5S_ALL, H5P_DEFAULT, flat.data());
    for (size_t i = 0; i < dims[0]; ++i) {
      out[i] = std::string(flat.data() + i * kRngStrLen);
      const auto nul = out[i].find('\0');
      if (nul != std::string::npos) {
        out[i].resize(nul);
      }
    }
  }
  H5Sclose(space);
  H5Dclose(ds);
  H5Tclose(type);
  return out;
}

}  // namespace gutibm::h5
#endif  // GUTIBM_HDF5

#endif  // GUTIBM_H5_UTIL_H
