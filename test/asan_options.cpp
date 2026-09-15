// Default AddressSanitizer options for the test binary, read by the ASan runtime at startup.
// ASan reports an ODR violation for libstdc++'s inline variable `std::piecewise_construct`, which
// is defined in both the instrumented test executable and the instrumented goddard_lib shared
// library. The report is a known false positive for inline variables and would abort every run.
// Setting the option here applies however the binary is launched (ctest, test discovery, or
// directly). The function is unused in builds without ASan.
extern "C" const char* __asan_default_options();
extern "C" const char* __asan_default_options() {
    return "detect_odr_violation=0";
}
