#include "goddard/global.hpp"
#include "goddard/config.h"

namespace Goddard {

    
void setup_defaults(){
    
    add_directory(DATA_DIR);
}

void add_directory(const std::string& dir){
    Cantera::addDirectory(dir);
}


}