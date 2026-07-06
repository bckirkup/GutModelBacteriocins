/* -----------------------------------------------------------------------
   GutIBM – Graceful shutdown on SIGTERM/SIGINT (Spec 4)
   ----------------------------------------------------------------------- */

#ifndef GUTIBM_STOP_SIGNAL_H
#define GUTIBM_STOP_SIGNAL_H

namespace gutibm {

void install_stop_signal_handlers();
bool gutibm_stop_requested();
void gutibm_reset_stop_request();
void gutibm_request_stop();  // test hook / cooperative shutdown

}  // namespace gutibm

#endif  // GUTIBM_STOP_SIGNAL_H
