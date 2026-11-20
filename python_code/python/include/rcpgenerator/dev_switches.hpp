// Developer switches (BUILD.md section 15, Y6; Ken, 2026-11-19: "hidden by construct, only a developer who is willing
// to recompile from source can turn them back on"). The environment variables of the study tools are read only in a
// build configured with -DRCP_DEV_SWITCHES=ON; in a normal build rcp_dev_getenv() returns null, whatever is set.
#pragma once
#include <cstdlib>

namespace rcpgenerator {
inline const char* rcp_dev_getenv(const char* name) {
#ifdef RCP_DEV_SWITCHES
    return std::getenv(name);
#else
    (void)name;
    return nullptr;
#endif
}
inline bool rcp_dev_build() {
#ifdef RCP_DEV_SWITCHES
    return true;
#else
    return false;
#endif
}
}  // namespace rcpgenerator
