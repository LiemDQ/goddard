#pragma once
#include <format>
#include <iostream>
#include <string>
#include <string_view>
#include <vector>
#include "goddard/moc.hpp"
#include "goddard/moc_thermo.hpp"
#include "goddard/profile.hpp"

namespace Goddard {

/** Messages collected during one solve(); surfaces as MocResult::messages. */
class MocLog {
public:
    explicit MocLog(MocLogLevel level = MocLogLevel::NORMAL) : m_level(level) {}

    template <typename... Args>
    void warning(std::string_view fmt, Args&&... args) {
        messages.push_back("Warning: " + std::vformat(fmt, std::make_format_args(args...)));
    }
    void warning(const std::string& msg) {
        messages.push_back("Warning: " + msg);
    }

    template <typename... Args>
    void info(std::string_view fmt, Args&&... args) {
        messages.push_back("Info: " + std::vformat(fmt, std::make_format_args(args...)));
    }
    void info(const std::string& msg) {
        messages.push_back("Info: " + msg);
    }

    // Verbose kernel/initialization trace, only active when level() == MocLogLevel::DEBUG.
    // Collected in `messages` like warning()/info() (so it surfaces via MocResult::messages
    // and the Python binding with no extra plumbing), and additionally echoed live to stderr
    // immediately as each call happens -- useful for a hang or a solve that never returns
    // (maxiter reached), where messages collected only in the returned MocResult would never
    // be seen.
    template <typename... Args>
    void debug(std::string_view fmt, Args&&... args) {
        if (m_level == MocLogLevel::DEBUG) {
            debug(std::vformat(fmt, std::make_format_args(args...)));
        }
    }
    void debug(const std::string& msg) {
        if (m_level != MocLogLevel::DEBUG) return;
        std::string full = "Debug: " + msg;
        messages.push_back(full);
        std::cerr << full << std::endl;
    }

    MocLogLevel level() const { return m_level; }

    std::vector<std::string> messages;

private:
    MocLogLevel m_level;
};

/**
 * Everything a kernel or unit process needs for one solve(), besides the points it works on.
 * Immutable for the solve: built once solve() has the thermo object and the resolved contour,
 * and passed by const reference everywhere.
 */
struct MocSolveContext {
    const MocOptions& options;
    const NozzleProfile& wall;    // the contour the solve marches against (options.nozzle_profile, or the generated Rao contour)
    const MocThermo& thermo;
    MocLog& log;
};

} // namespace Goddard
