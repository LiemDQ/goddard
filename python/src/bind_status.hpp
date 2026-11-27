#pragma once
#include <string>

#include "goddard/status.hpp"

/** Text of a status for a `__repr__`: "ok", or "FAILED (<reason>)". */
inline std::string status_tag(const Goddard::SolveStatus& status) {
    if (status.ok()) {
        return "ok";
    }
    return "FAILED (" + std::string(Goddard::to_string(status.reason)) + ")";
}
