/* -----------------------------------------------------------------------
   GutIBM – Graceful shutdown signal handlers (Spec 4)
   ----------------------------------------------------------------------- */

#include "stop_signal.h"

#include <csignal>

namespace gutibm {

namespace {

volatile std::sig_atomic_t g_stop_requested = 0;

void on_stop_signal(int /*sig*/) { g_stop_requested = 1; }

}  // namespace

void install_stop_signal_handlers() {
  std::signal(SIGTERM, on_stop_signal);
  std::signal(SIGINT, on_stop_signal);
}

bool gutibm_stop_requested() { return g_stop_requested != 0; }

void gutibm_reset_stop_request() { g_stop_requested = 0; }

void gutibm_request_stop() { g_stop_requested = 1; }

}  // namespace gutibm
