/* -----------------------------------------------------------------------
   GutIBM – Filesystem path validation helpers
   ----------------------------------------------------------------------- */

#include "path_utils.h"

#include <algorithm>
#include <cerrno>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include "error.h"
#include <string>
#include <string_view>
#include <system_error>
#include <unistd.h>
#include <vector>

#include <sys/stat.h>

namespace gutibm {

namespace {

namespace fs = std::filesystem;

bool path_has_null_byte(const std::string_view path) {
  return path.find('\0') != std::string_view::npos;
}

bool path_has_parent_traversal(const std::string& path) {
  const fs::path p(path);
  return std::ranges::find(p, fs::path("..")) != p.end();
}

bool is_world_writable_directory(const fs::path& dir) {
  struct stat st {};
  if (stat(dir.c_str(), &st) != 0) return false;
  if (!S_ISDIR(st.st_mode)) return false;
  return (st.st_mode & S_IWOTH) != 0;
}

bool has_sticky_bit(const fs::path& dir) {
  struct stat st {};
  if (stat(dir.c_str(), &st) != 0) return false;
  return (st.st_mode & S_ISVTX) != 0;
}

std::string validate_temp_directory(const fs::path& dir) {
  validate_path_syntax(dir.string());
  if (!fs::exists(dir) || !fs::is_directory(dir)) {
    throw PathError("temporary directory is not a directory: " + dir.string());
  }
  if (fs::is_symlink(fs::symlink_status(dir))) {
    throw PathError("refusing to use symlinked temporary directory: " +
                    dir.string());
  }

  const fs::path canon = fs::weakly_canonical(dir);
  if (is_world_writable_directory(canon) && !has_sticky_bit(canon)) {
    throw PathError(
        "world-writable temporary directory lacks sticky bit: " +
        canon.string());
  }
  return canon.string();
}

fs::path parent_directory_for(const fs::path& path) {
  if (path.has_parent_path()) return path.parent_path();
  return fs::path(".");
}

void reject_symlink_in_world_writable_parent(const fs::path& path,
                                             const char* operation) {
  const fs::path parent = parent_directory_for(path);
  if (!fs::exists(parent)) return;

  if (const fs::path parent_canon = fs::weakly_canonical(parent);
      !is_world_writable_directory(parent_canon)) {
    return;
  }

  if (auto status = fs::symlink_status(parent); fs::is_symlink(status)) {
    throw PathError(
        std::string("refusing to ") + operation +
        " via symlinked parent in world-writable directory: " + parent.string());
  }

  if (fs::exists(path) && fs::is_symlink(fs::symlink_status(path))) {
    throw PathError(
        std::string("refusing to ") + operation +
        " through symlink in world-writable directory: " + path.string());
  }
}

std::string temp_directory() {
  fs::path dir;
  if (const char* tmpdir = std::getenv("TMPDIR"); tmpdir && tmpdir[0] != '\0') {  // NOSONAR cpp:S5443 — validated below
    dir = fs::path(tmpdir);
  } else {
    dir = fs::temp_directory_path();  // NOSONAR cpp:S5443 — validated below; mkstemp creates private file
  }
  return validate_temp_directory(dir);
}

}  // namespace

void validate_path_syntax(const std::string& path) {
  if (path.empty()) {
    throw PathError("empty path");
  }
  if (path_has_null_byte(path)) {
    throw PathError("path contains null byte");
  }
  if (path_has_parent_traversal(path)) {
    throw PathError("path contains parent-directory traversal ('..')");
  }
}

std::string validate_input_file_path(const std::string& path) {
  validate_path_syntax(path);

  const fs::path p(path);
  if (!fs::exists(p)) {
    throw PathError("input file not found: " + path);
  }

  reject_symlink_in_world_writable_parent(p, "read");

  if (!fs::is_regular_file(p)) {
    throw PathError("input path is not a regular file: " + path);
  }

  return fs::weakly_canonical(p).string();
}

void validate_output_file_path(const std::string& path) {
  validate_path_syntax(path);

  const fs::path p(path);
  const fs::path parent = parent_directory_for(p);

  if (!fs::exists(parent)) {
    throw PathError("output directory does not exist: " + parent.string());
  }
  if (!fs::is_directory(parent)) {
    throw PathError("output parent is not a directory: " + parent.string());
  }

  reject_symlink_in_world_writable_parent(p, "write");

  if (fs::exists(p)) {
    const fs::file_status status = fs::symlink_status(p);
    if (fs::is_symlink(status)) {
      throw PathError("refusing to overwrite symlink: " + path);
    }
    if (!fs::is_regular_file(status) && !fs::is_regular_file(p)) {
      throw PathError("output path exists and is not a regular file: " + path);
    }
  }
}

std::string secure_temp_file(const std::string& prefix) {
  validate_path_syntax(prefix);

  const std::string dir = temp_directory();
  std::string pattern = dir;
  if (!pattern.empty() && pattern.back() != '/') pattern += '/';
  pattern += prefix;
  pattern += "XXXXXX";

  std::vector<char> buf(pattern.begin(), pattern.end());
  buf.push_back('\0');

  const int fd = mkstemp(buf.data());
  if (fd < 0) {
    throw PathError(
        std::string("failed to create secure temporary file in ") + dir +
        ": " + std::strerror(errno));
  }
  close(fd);

  return std::string(buf.data());
}

std::string resolve_test_h5_path(const char* env_var, const std::string& tag) {
  validate_path_syntax(tag);

  if (env_var) {
    if (const char* env_path = std::getenv(env_var); env_path && env_path[0] != '\0') {
      validate_output_file_path(env_path);
      return env_path;
    }
  }

  return secure_temp_file("gutibm_" + tag + "_");
}

}  // namespace gutibm
