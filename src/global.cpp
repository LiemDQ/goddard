#include "goddard/global.hpp"
#include "goddard/config.h"

#include <filesystem>
#include <system_error>

namespace Goddard {

    
void setup_defaults(){
    static bool is_initialized;
    if (!is_initialized){
        // DATA_DIR is empty in wheels, and absent when the source tree or install prefix it was
        // configured with has moved. The Python package then registers its own data directory.
        const std::filesystem::path data_dir(DATA_DIR);
        std::error_code error;
        if (!data_dir.empty() && std::filesystem::is_directory(data_dir, error)){
            add_directory(data_dir.string());
        }

        is_initialized = true;
    }
}

void add_directory(const std::string& dir){
    Cantera::addDataDirectory(dir);
}


}
