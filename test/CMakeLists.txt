add_executable(tmst_test test.cpp)
target_include_directories(tmst_test PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_compile_options(tmst_test PRIVATE -Wall -Wextra -Wpedantic)
add_test(NAME tmst_test COMMAND tmst_test)
# A broken link can create a cycle and loop forever; count that as a failure.
set_tests_properties(tmst_test PROPERTIES TIMEOUT 600)
