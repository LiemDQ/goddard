#pragma once

#include <variant>
#include <vector>

namespace Goddard {


enum class CombustorType {
    NONE,
    INFINITE_AREA,
    FINITE_MASS_FLUX,
    FINITE_CONTRACTION_RATIO
};

struct CombustorOptions {
    CombustorType type;
    std::vector<double> pressures;
    double mass_flux;
    double contraction_ratio;
};

enum class NozzleChemistryType {
    NONE,
    FROZEN,
    EQUILIBRIUM,
    KINETIC
};

enum class ExpansionType {
    SUPERSONIC_AREA_RATIO,
    SUBSONIC_AREA_RATIO,
    PRESSURE_RATIO
};


struct NozzleOptions {
    NozzleChemistryType chemistry;
    ExpansionType expansion_type;
    std::vector<double> expansion_ratios;
    unsigned int frozen_NFZ = 1;
};


} //namespace goddard