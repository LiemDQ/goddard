#include "goddard/global.hpp"
#include "goddard/config.h"

namespace Goddard {

    
void setup_defaults(){
    static bool is_initialized;
    if (!is_initialized){
        add_directory(DATA_DIR);

        is_initialized = true;
    }
}

void add_directory(const std::string& dir){
    Cantera::addDataDirectory(dir);
}


}