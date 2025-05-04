#pragma once
#include <stdexcept>
#include <string>

namespace Goddard {

class NotImplementedError: public std::logic_error {
    using std::logic_error::logic_error;

    public:
    virtual const char * what() const noexcept override;

};


}