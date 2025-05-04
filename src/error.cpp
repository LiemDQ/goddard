#include "goddard/error.hpp"
#include <string>

namespace Goddard {

const char * NotImplementedError::what() const noexcept {
    std::string msg = "Functionality not yet implemented: ";
    msg += std::logic_error::what();

    return msg.c_str();
}



}