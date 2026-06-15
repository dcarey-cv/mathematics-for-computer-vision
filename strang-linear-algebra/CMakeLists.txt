add_library(strang_deps INTERFACE)
target_link_libraries(strang_deps INTERFACE Eigen3::Eigen)
target_link_libraries(strang_deps INTERFACE fmt::fmt)

if(STRICT_WARNINGS)
  target_compile_options(strang_deps INTERFACE -Wall -Wextra -Werror -Wpedantic)
endif()

add_executable(chapter_01 chapter_01.cpp)
target_link_libraries(chapter_01 PRIVATE strang_deps)
target_link_options(chapter_01 PRIVATE -Wl,--subsystem,console -Wl,-e,mainCRTStartup)

add_executable(chapter_02 chapter_02.cpp)
target_link_libraries(chapter_02 PRIVATE strang_deps)
target_link_options(chapter_02 PRIVATE -Wl,--subsystem,console -Wl,-e,mainCRTStartup)