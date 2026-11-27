#pragma once
/**
 * @file status_test_helpers.hpp
 * @brief Assertions on `SolveStatus` values and `SolveError` exceptions.
 *
 * Use them as `EXPECT_TRUE(failed_with(result.status, FailureReason::NO_SOLUTION, "below 1"))`, so
 * that a failing assertion is reported at the line of the test and prints the status it found.
 */

#include "goddard/error.hpp"
#include "goddard/status.hpp"

#include "gtest/gtest.h"

#include <exception>
#include <string>
#include <string_view>
#include <utility>

namespace Goddard {
namespace test_helpers {

/** Success if `status` is ok; otherwise a failure that prints what the status describes. */
inline ::testing::AssertionResult is_ok(const SolveStatus& status) {
    if (status.ok()) {
        return ::testing::AssertionSuccess();
    }
    return ::testing::AssertionFailure() << "status is not ok: " << status.describe();
}

/**
 * Success if `status` is a failure with the given reason whose message contains `message_part`
 * (any message if empty).
 */
inline ::testing::AssertionResult failed_with(const SolveStatus& status, FailureReason reason,
                                              std::string_view message_part = "")
{
    if (status.reason != reason) {
        return ::testing::AssertionFailure()
            << "expected reason " << to_string(reason) << ", got " << to_string(status.reason)
            << " (" << status.describe() << ")";
    }
    if (status.ok()) {
        return ::testing::AssertionFailure() << "expected a failure, got an ok status";
    }
    if (status.message.empty()) {
        return ::testing::AssertionFailure() << "a failed status must carry a message";
    }
    if (status.message.find(message_part) == std::string::npos) {
        return ::testing::AssertionFailure()
            << "expected a message containing \"" << message_part << "\", got \"" << status.message << "\"";
    }
    return ::testing::AssertionSuccess();
}

/**
 * Success if `call()` throws a `SolveError` (or a class derived from it) with the given reason
 * whose message contains `message_part`. Any other exception, and no exception, is a failure.
 */
template <typename Call>
::testing::AssertionResult throws_solve_error(Call&& call, FailureReason reason,
                                              std::string_view message_part = "")
{
    try {
        std::forward<Call>(call)();
    } catch (const SolveError& error) {
        return failed_with(error.status(), reason, message_part);
    } catch (const std::exception& error) {
        return ::testing::AssertionFailure()
            << "expected a SolveError, got another exception: " << error.what();
    }
    return ::testing::AssertionFailure() << "expected a SolveError, nothing was thrown";
}

} // namespace test_helpers
} // namespace Goddard
