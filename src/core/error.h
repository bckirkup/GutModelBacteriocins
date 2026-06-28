/* -----------------------------------------------------------------------
   GutIBM – Domain-specific exception hierarchy
   ----------------------------------------------------------------------- */

#ifndef GUTIBM_ERROR_H
#define GUTIBM_ERROR_H

#include <stdexcept>
#include <string>

namespace gutibm {

class Error : public std::runtime_error {
 public:
  using std::runtime_error::runtime_error;
};

class ConfigError : public Error {
 public:
  using Error::Error;
};

class IOError : public Error {
 public:
  using Error::Error;
};

class HDF5Error : public IOError {
 public:
  using IOError::IOError;
};

class PathError : public IOError {
 public:
  using IOError::IOError;
};

class CheckpointError : public Error {
 public:
  using Error::Error;
};

class SimulationError : public Error {
 public:
  using Error::Error;
};

}  // namespace gutibm

#endif  // GUTIBM_ERROR_H
