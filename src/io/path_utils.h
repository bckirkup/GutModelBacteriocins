/* -----------------------------------------------------------------------
   GutIBM – Filesystem path validation helpers
   ----------------------------------------------------------------------- */

#ifndef GUTIBM_PATH_UTILS_H
#define GUTIBM_PATH_UTILS_H

#include <string>

namespace gutibm {

// Reject null bytes and parent-directory traversal components.
void validate_path_syntax(const std::string& path);

// Input files must exist as regular files; symlinks in world-writable
// directories are rejected. Returns a canonical path for I/O.
std::string validate_input_file_path(const std::string& path);

// Output paths must have a usable parent directory; symlinks in
// world-writable directories are rejected.
void validate_output_file_path(const std::string& path);

// Create a private temporary file with an unpredictable name (0600).
std::string secure_temp_file(const std::string& prefix);

// Resolve an HDF5 test output path from an env override or secure temp.
std::string resolve_test_h5_path(const char* env_var, const std::string& tag);

}  // namespace gutibm

#endif  // GUTIBM_PATH_UTILS_H
