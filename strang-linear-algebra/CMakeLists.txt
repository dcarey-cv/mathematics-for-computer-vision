add_library(strang_deps INTERFACE)
target_link_libraries(strang_deps INTERFACE Eigen3::Eigen)
target_link_libraries(strang_deps INTERFACE fmt::fmt)

if(STRICT_WARNINGS)
  target_compile_options(strang_deps INTERFACE -Wall -Wextra -Werror -Wpedantic)
endif()

add_executable(GaussianElimination GaussianElimination.cpp)
target_link_libraries(GaussianElimination PRIVATE strang_deps)
target_link_options(GaussianElimination PRIVATE -Wl,--subsystem,console -Wl,-e,mainCRTStartup)
target_precompile_headers(GaussianElimination PRIVATE <Eigen/Dense> <fmt/format.h>)

add_executable(SystemDecomposition SystemDecomposition.cpp)
target_link_libraries(SystemDecomposition PRIVATE strang_deps)
target_link_options(SystemDecomposition PRIVATE -Wl,--subsystem,console -Wl,-e,mainCRTStartup)
target_precompile_headers(SystemDecomposition REUSE_FROM GaussianElimination)