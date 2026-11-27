#include "goddard/error.hpp"
#include "goddard/solve_guard.hpp"
#include "goddard/status.hpp"

#include "cantera/base/ctexceptions.h"
#include "cantera/core.h"
#include "gtest/gtest.h"
#include "status_test_helpers.hpp"

#include <cmath>
#include <stdexcept>
#include <string>
#include <type_traits>

using namespace Goddard;
using Goddard::test_helpers::failed_with;
using Goddard::test_helpers::is_ok;

// ============================================================
//  SolveStatus
// ============================================================

TEST(SolveStatus, DefaultIsOkWithNoDetails) {
    const SolveStatus status;
    EXPECT_TRUE(status.ok());
    EXPECT_EQ(status.reason, FailureReason::NONE);
    EXPECT_TRUE(status.message.empty());
    EXPECT_TRUE(status.stage.empty());
    EXPECT_EQ(status.iterations, -1);
    EXPECT_TRUE(std::isnan(status.tolerance));
    EXPECT_TRUE(std::isnan(status.residual));
    EXPECT_TRUE(status.warnings.empty());
    EXPECT_EQ(status.describe(), "ok");
}

TEST(SolveStatus, WarningsDoNotMakeAStatusFail) {
    SolveStatus status;
    status.warnings.push_back("resolution");
    EXPECT_TRUE(status.ok());
    EXPECT_NO_THROW(status.raise_if_failed());
}

TEST(SolveStatus, DescribeListsTheWarnings) {
    SolveStatus status;
    status.warnings = {"extrapolated", "resolution"};
    EXPECT_EQ(status.describe(), "ok; warning: extrapolated; warning: resolution");

    SolveStatus failed = SolveStatus::failure(FailureReason::NO_SOLUTION, "detached");
    failed.warnings = {"extrapolated"};
    EXPECT_EQ(failed.describe(), "detached (NO_SOLUTION); warning: extrapolated");
}

TEST(SolveStatus, FailureCarriesReasonAndMessage) {
    const SolveStatus status = SolveStatus::failure(FailureReason::NO_SOLUTION, "detached");
    EXPECT_TRUE(failed_with(status, FailureReason::NO_SOLUTION, "detached"));
    EXPECT_FALSE(status.ok());
}

TEST(SolveStatus, DescribeIncludesStageMessageReasonAndIterationData) {
    SolveStatus status = SolveStatus::failure(FailureReason::NOT_CONVERGED, "no root");
    status.stage = "throat";
    status.iterations = 12;
    status.tolerance = 1e-6;
    status.residual = 2.5e-3;
    const std::string text = status.describe();
    EXPECT_NE(text.find("throat"), std::string::npos) << text;
    EXPECT_NE(text.find("no root"), std::string::npos) << text;
    EXPECT_NE(text.find("NOT_CONVERGED"), std::string::npos) << text;
    EXPECT_NE(text.find("12 iterations"), std::string::npos) << text;
    EXPECT_NE(text.find("residual 0.0025"), std::string::npos) << text;
    EXPECT_NE(text.find("tolerance 1e-06"), std::string::npos) << text;
    EXPECT_LT(text.find("throat"), text.find("no root"));

    // Values that were not reported are left out
    const std::string bare = SolveStatus::failure(FailureReason::NO_SOLUTION, "none").describe();
    EXPECT_EQ(bare.find("iterations"), std::string::npos) << bare;
    EXPECT_EQ(bare.find("residual"), std::string::npos) << bare;
    EXPECT_EQ(bare.find("tolerance"), std::string::npos) << bare;
}

TEST(SolveStatus, ToStringNamesEveryReason) {
    EXPECT_EQ(to_string(FailureReason::NONE), "NONE");
    EXPECT_EQ(to_string(FailureReason::NOT_CONVERGED), "NOT_CONVERGED");
    EXPECT_EQ(to_string(FailureReason::NO_SOLUTION), "NO_SOLUTION");
    EXPECT_EQ(to_string(FailureReason::OUT_OF_DATA_RANGE), "OUT_OF_DATA_RANGE");
    EXPECT_EQ(to_string(FailureReason::THERMO_FAILURE), "THERMO_FAILURE");
    EXPECT_EQ(to_string(FailureReason::INCOMPLETE), "INCOMPLETE");
    EXPECT_EQ(to_string(FailureReason::NOT_COMPUTED), "NOT_COMPUTED");
}

TEST(SolveStatus, NotComputedNamesTheUpstreamStageAndMessage) {
    SolveStatus upstream = SolveStatus::failure(FailureReason::NOT_CONVERGED, "no root");
    upstream.stage = "throat";
    const SolveStatus downstream = SolveStatus::not_computed(upstream);
    EXPECT_TRUE(failed_with(downstream, FailureReason::NOT_COMPUTED, "throat"));
    EXPECT_NE(downstream.message.find("no root"), std::string::npos) << downstream.message;

    // Without a stage the message still repeats the upstream message
    const SolveStatus anonymous = SolveStatus::not_computed(SolveStatus::failure(FailureReason::NO_SOLUTION, "choked"));
    EXPECT_TRUE(failed_with(anonymous, FailureReason::NOT_COMPUTED, "choked"));
}

TEST(SolveStatus, NotComputedOfANotComputedStatusIsPassedThroughUnchanged) {
    SolveStatus root = SolveStatus::failure(FailureReason::NO_SOLUTION, "choked");
    root.stage = "finite-area chamber";
    const SolveStatus once = SolveStatus::not_computed(root);
    const SolveStatus twice = SolveStatus::not_computed(once);
    EXPECT_EQ(twice.reason, FailureReason::NOT_COMPUTED);
    EXPECT_EQ(twice.message, once.message);
    EXPECT_EQ(twice.stage, once.stage);
    // The reason is named once, not as "not computed because ... not computed because ...".
    const std::string needle = "not computed because";
    const size_t first = twice.message.find(needle);
    ASSERT_NE(first, std::string::npos) << twice.message;
    EXPECT_EQ(twice.message.find(needle, first + 1), std::string::npos) << twice.message;
    EXPECT_NE(twice.message.find("choked"), std::string::npos) << twice.message;
}

TEST(SolveStatus, RaiseIfFailedDoesNothingForAnOkStatus) {
    EXPECT_NO_THROW(SolveStatus{}.raise_if_failed());
}

TEST(SolveStatus, RaiseIfFailedThrowsConvergenceErrorOnlyForNotConverged) {
    SolveStatus status = SolveStatus::failure(FailureReason::NOT_CONVERGED, "stalled");
    status.iterations = 30;
    status.tolerance = 1e-8;
    status.residual = 1e-3;
    try {
        status.raise_if_failed();
        FAIL() << "expected a ConvergenceError";
    } catch (const ConvergenceError& error) {
        EXPECT_EQ(error.reason(), FailureReason::NOT_CONVERGED);
        EXPECT_EQ(error.iterations(), 30);
        EXPECT_EQ(error.tolerance(), 1e-8);
        EXPECT_EQ(error.residual(), 1e-3);
        EXPECT_EQ(error.status().message, "stalled");
    }

    for (FailureReason reason : {FailureReason::NO_SOLUTION, FailureReason::OUT_OF_DATA_RANGE,
                                 FailureReason::THERMO_FAILURE, FailureReason::INCOMPLETE,
                                 FailureReason::NOT_COMPUTED}) {
        try {
            SolveStatus::failure(reason, "x").raise_if_failed();
            FAIL() << "expected a SolveError for " << to_string(reason);
        } catch (const ConvergenceError&) {
            FAIL() << "only NOT_CONVERGED gives a ConvergenceError, got one for " << to_string(reason);
        } catch (const SolveError& error) {
            EXPECT_EQ(error.reason(), reason);
        }
    }
}

// ============================================================
//  Exceptions
// ============================================================

TEST(SolveErrorTypes, HierarchyAndWhat) {
    static_assert(std::is_base_of_v<std::runtime_error, SolveError>);
    static_assert(std::is_base_of_v<SolveError, ConvergenceError>);

    const SolveError error(FailureReason::OUT_OF_DATA_RANGE, "T outside the data");
    EXPECT_EQ(error.reason(), FailureReason::OUT_OF_DATA_RANGE);
    EXPECT_EQ(error.status().message, "T outside the data");
    EXPECT_EQ(std::string(error.what()), error.status().describe());
}

TEST(SolveErrorTypes, ConvergenceErrorIsASolveErrorAndReportsWhatItWasGiven) {
    try {
        throw ConvergenceError("Newton", 25, 1e-9, 4e-4);
    } catch (const SolveError& error) {
        EXPECT_EQ(error.reason(), FailureReason::NOT_CONVERGED);
        const auto* convergence = dynamic_cast<const ConvergenceError*>(&error);
        ASSERT_NE(convergence, nullptr);
        EXPECT_EQ(convergence->iterations(), 25);
        EXPECT_EQ(convergence->tolerance(), 1e-9);
        EXPECT_EQ(convergence->residual(), 4e-4);
    }
}

TEST(SolveErrorTypes, UnreportedConvergenceDataIsMinusOneIterationsAndNaN) {
    const ConvergenceError error("Newton");
    EXPECT_EQ(error.iterations(), -1);
    EXPECT_TRUE(std::isnan(error.tolerance()));
    EXPECT_TRUE(std::isnan(error.residual()));
}

// ============================================================
//  Capturing failures at a solver boundary
// ============================================================

TEST(CaptureFailure, ReturnsAnOkStatusWhenTheBodyReturns) {
    int calls = 0;
    const SolveStatus status = capture_failure("stage", [&] { calls++; });
    EXPECT_TRUE(is_ok(status));
    EXPECT_EQ(calls, 1);
}

TEST(CaptureFailure, SolveErrorBecomesItsStatusWithTheStage) {
    const SolveStatus status = capture_failure("exit station 2", [] {
        throw SolveError(FailureReason::NO_SOLUTION, "choked");
    });
    EXPECT_TRUE(failed_with(status, FailureReason::NO_SOLUTION, "choked"));
    EXPECT_EQ(status.stage, "exit station 2");
}

TEST(CaptureFailure, ConvergenceDataSurvivesTheCapture) {
    const SolveStatus status = capture_failure("throat", [] { throw ConvergenceError("Newton", 7, 1e-6, 0.5); });
    EXPECT_TRUE(failed_with(status, FailureReason::NOT_CONVERGED, "Newton"));
    EXPECT_EQ(status.iterations, 7);
    EXPECT_EQ(status.tolerance, 1e-6);
    EXPECT_EQ(status.residual, 0.5);
}

TEST(CaptureFailure, AStageAlreadyNamedIsKept) {
    SolveStatus named = SolveStatus::failure(FailureReason::NO_SOLUTION, "x");
    named.stage = "inner";
    const SolveStatus status = capture_failure("outer", [&] { throw SolveError(named); });
    EXPECT_EQ(status.stage, "inner");
}

TEST(CaptureFailure, CanteraErrorBecomesAThermoFailure) {
    const SolveStatus status = capture_failure("equilibrium", [] {
        throw Cantera::CanteraError("Test::procedure", "no convergence for T");
    });
    EXPECT_TRUE(failed_with(status, FailureReason::THERMO_FAILURE, "no convergence for T"));
    EXPECT_EQ(status.stage, "equilibrium");
}

TEST(CaptureFailure, InvalidArgumentAndOtherExceptionsPropagate) {
    EXPECT_THROW(capture_failure("stage", [] { throw std::invalid_argument("bad"); }), std::invalid_argument);
    EXPECT_THROW(capture_failure("stage", [] { throw NotImplementedError("later"); }), NotImplementedError);
    EXPECT_THROW(capture_failure("stage", [] { throw std::out_of_range("index"); }), std::out_of_range);
}

TEST(SetStageIfUnset, SetsOnlyAFailedStatusThatHasNoStage) {
    SolveStatus ok;
    set_stage_if_unset(ok, "stage");
    EXPECT_TRUE(ok.stage.empty());

    SolveStatus failed = SolveStatus::failure(FailureReason::NO_SOLUTION, "x");
    set_stage_if_unset(failed, "first");
    set_stage_if_unset(failed, "second");
    EXPECT_EQ(failed.stage, "first");
}

// ============================================================
//  Translating backend errors
// ============================================================

TEST(TranslateBackendErrors, ReturnsTheValueOfASuccessfulCall) {
    EXPECT_EQ(translate_backend_errors("op", [] { return 42; }), 42);
    EXPECT_EQ(translate_backend_usage_errors("op", [] { return std::string("ok"); }), "ok");
}

TEST(TranslateBackendErrors, CanteraErrorBecomesASolveErrorThatNamesTheOperation) {
    try {
        translate_backend_errors("Gas::set_state_TP(T = 0 K)", []() -> int {
            throw Cantera::CanteraError("ThermoPhase::setState_TP", "temperature must be positive");
        });
        FAIL() << "expected a SolveError";
    } catch (const SolveError& error) {
        EXPECT_EQ(error.reason(), FailureReason::THERMO_FAILURE);
        const std::string message = error.status().message;
        EXPECT_EQ(message.rfind("Gas::set_state_TP(T = 0 K)", 0), 0u) << message;
        EXPECT_NE(message.find("temperature must be positive"), std::string::npos) << message;
    }
}

TEST(TranslateBackendErrors, OtherExceptionsPassThroughUnchanged) {
    EXPECT_THROW(translate_backend_errors("op", [] { throw std::invalid_argument("bad"); }), std::invalid_argument);
    EXPECT_THROW(translate_backend_errors("op", [] { throw SolveError(FailureReason::NO_SOLUTION, "x"); }),
                 SolveError);
}

TEST(TranslateBackendUsageErrors, CanteraErrorBecomesAnInvalidArgument) {
    try {
        translate_backend_usage_errors("Gas::create", []() -> int {
            throw Cantera::CanteraError("newSolution", "no such phase");
        });
        FAIL() << "expected std::invalid_argument";
    } catch (const std::invalid_argument& error) {
        const std::string message = error.what();
        EXPECT_NE(message.find("Gas::create"), std::string::npos) << message;
        EXPECT_NE(message.find("no such phase"), std::string::npos) << message;
    }
}

TEST(TranslateBackendUsageErrors, OtherExceptionsPassThroughUnchanged) {
    EXPECT_THROW(translate_backend_usage_errors("op", [] { throw SolveError(FailureReason::NO_SOLUTION, "x"); }),
                 SolveError);
}

TEST(TranslateBackendErrors, MessageCallableRunsOnlyWhenTheCallFails) {
    int message_calls = 0;
    const auto message = [&] {
        message_calls++;
        return std::string("Gas::set_state_TP(T = 0 K)");
    };

    EXPECT_EQ(translate_backend_errors(message, [] { return 7; }), 7);
    EXPECT_EQ(translate_backend_usage_errors(message, [] { return 8; }), 8);
    EXPECT_EQ(message_calls, 0) << "the message was built although nothing failed";

    try {
        translate_backend_errors(message, []() -> int {
            throw Cantera::CanteraError("ThermoPhase::setState_TP", "temperature must be positive");
        });
        FAIL() << "expected a SolveError";
    } catch (const SolveError& error) {
        EXPECT_EQ(error.status().message.rfind("Gas::set_state_TP(T = 0 K)", 0), 0u);
    }
    EXPECT_EQ(message_calls, 1);

    EXPECT_THROW(translate_backend_usage_errors(message, []() -> int {
        throw Cantera::CanteraError("newSolution", "no such phase");
    }), std::invalid_argument);
    EXPECT_EQ(message_calls, 2);
}

TEST(ValidateEquilibriumRequest, AcceptsExactlyThePairsCanteraKnows) {
    for (const char* pair : {"TP", "TV", "HP", "SP", "SV", "UV"}) {
        EXPECT_NO_THROW(validate_equilibrium_request(pair, "auto")) << pair;
    }
    for (const char* pair : {"PT", "VT", "PH", "TT", "TH", "HS", "UP", "tp", "T", "TPX", ""}) {
        EXPECT_THROW(validate_equilibrium_request(pair, "auto"), std::invalid_argument) << pair;
    }
    for (const char* solver : {"auto", "element_potential", "gibbs", "vcs"}) {
        EXPECT_NO_THROW(validate_equilibrium_request("TP", solver)) << solver;
    }
    EXPECT_THROW(validate_equilibrium_request("TP", "gibs"), std::invalid_argument);
}
