// Nozzle expansion with condensed products (work package D).
//
// The equilibrium and frozen nozzle stations that carry condensed phases are not implemented yet;
// this file exists so the suite is wired into CMake ahead of that work.

#include "goddard/gas.hpp"

#include "gtest/gtest.h"

using namespace Goddard;

TEST(NozzleCondensed, GasOnlyNozzleGasHasNoCondensedCandidates) {
    Gas gas("h2o2.yaml", "ohmech", GasChemistry::EQUILIBRIUM);
    EXPECT_FALSE(gas.has_condensed_candidates());
    EXPECT_FALSE(gas.has_condensed_phases());
}
