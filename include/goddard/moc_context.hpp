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
    /** @param level MocLogLevel::DEBUG additionally enables debug() and echoes it to stderr. */
    explicit MocLog(MocLogLevel level = MocLogLevel::NORMAL) : m_level(level) {}

    /** Record a formatted warning; always collected. */
    template <typename... Args>
    void warning(std::string_view fmt, Args&&... args) {
        messages.push_back("Warning: " + std::vformat(fmt, std::make_format_args(args...)));
    }
    /** Record a warning; always collected. */
    void warning(const std::string& msg) {
        messages.push_back("Warning: " + msg);
    }

    /** Record a formatted informational message; always collected. */
    template <typename... Args>
    void info(std::string_view fmt, Args&&... args) {
        messages.push_back("Info: " + std::vformat(fmt, std::make_format_args(args...)));
    }
    /** Record an informational message; always collected. */
    void info(const std::string& msg) {
        messages.push_back("Info: " + msg);
    }

    /**
     * Record a formatted verbose kernel/initialization trace entry -- a no-op unless
     * level() == MocLogLevel::DEBUG.
     *
     * Collected in `messages` like warning()/info() (so it surfaces via MocResult::messages
     * and the Python binding with no extra plumbing), and additionally echoed live to stderr
     * immediately as each call happens -- useful for a hang or a solve that never returns
     * (maxiter reached), where messages collected only in the returned MocResult would never
     * be seen.
     */
    template <typename... Args>
    void debug(std::string_view fmt, Args&&... args) {
        if (m_level == MocLogLevel::DEBUG) {
            debug(std::vformat(fmt, std::make_format_args(args...)));
        }
    }
    /** Record a verbose trace entry -- a no-op unless level() == MocLogLevel::DEBUG. */
    void debug(const std::string& msg) {
        if (m_level != MocLogLevel::DEBUG) return;
        std::string full = "Debug: " + msg;
        messages.push_back(full);
        std::cerr << full << std::endl;
    }

    /** The verbosity this log was constructed with. */
    MocLogLevel level() const { return m_level; }

    /** Every message recorded so far, in call order, each prefixed by its kind. */
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
    const MocOptions& options; ///< The options the solve was invoked with.
    /** The contour the solve marches against: `options.nozzle_profile` in analysis mode, or
     *  the freshly generated contour in a design mode. */
    const NozzleProfile& wall;
    const MocThermo& thermo;   ///< Chemistry dispatch shared by every unit process this solve.
    MocLog& log;                ///< Sink for warnings, info, and (at DEBUG level) trace entries.
};

} // namespace Goddard
