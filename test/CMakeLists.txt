


add_executable(test_genData test_genData.f90)
target_link_libraries(test_genData FitLogNormal )
add_test(NAME test_genData_test COMMAND test_genData )
set_tests_properties(test_genData_test
    PROPERTIES
    LABELS "EXAMPLE;SHORT"
    #ctest will look for "PASS" in the stdout and consider it a pass if it finds it
    PASS_REGULAR_EXPRESSION "PASS"
) 
 