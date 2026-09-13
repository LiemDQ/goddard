# Diagnostic tools. Off by default: these are for investigation, not for the library's
# consumers, and building them adds link time to every ordinary build.
add_executable(moc_sweep moc_sweep.cpp)
target_link_libraries(moc_sweep PRIVATE goddard_lib)
