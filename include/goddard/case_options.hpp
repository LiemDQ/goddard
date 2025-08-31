#pragma once

#include <variant>

namespace Goddard {


enum class CombustorType {
    NONE,
    INFINITE_AREA,
    FINITE_MASS_FLUX,
    FINITE_CONTRACTION_RATIO
};

struct CombustorOptions {
    CombustorType type;    
    double mass_flux;
    double contraction_ratio;
};

struct FrozenRxOpts {
    unsigned int NFZ = 1;
};

struct EquilibriumRxOpts {
    //no options
};

struct NoRxOpts {
    //no options
};

struct KineticRxOpts {
    //no options
};

using ReactionOptions = std::variant<NoRxOpts, FrozenRxOpts, EquilibriumRxOpts, KineticRxOpts>;

struct RocketCaseOptions {
    CombustorOptions combustor_options;
    ReactionOptions rx_options;
    bool include_ionized = false;
};

} //namespace goddard