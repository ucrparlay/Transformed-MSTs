add_executable(tmst_bench bench.cpp)
target_include_directories(tmst_bench PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_compile_options(tmst_bench PRIVATE -Wall -Wextra -Wpedantic -march=native)
