add_library(strang_deps INTERFACE)
target_link_libraries(strang_deps INTERFACE Eigen3::Eigen)
target_link_libraries(strang_deps INTERFACE fmt::fmt)

add_library(strang_viz_deps INTERFACE)
target_link_libraries(strang_viz_deps INTERFACE Eigen3::Eigen pango_display pango_opengl pango_vars pango_windowing pango_plot pango_image pango_video)

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

add_executable(chapter_03 chapter_03.cpp)
target_link_libraries(chapter_03 PRIVATE strang_viz_deps)
target_link_options(chapter_03 PRIVATE -Wl,--subsystem,console -Wl,-e,mainCRTStartup)
target_precompile_headers(chapter_03 PRIVATE <Eigen/Dense> <pangolin/pangolin.h>)