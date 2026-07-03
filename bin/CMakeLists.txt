 



add_executable(iceberg ${CMAKE_SOURCE_DIR}/src/iceberg.f90)
target_link_libraries(iceberg FitLogNormal)

#Install the exec to bin folder of install path
install(TARGETS iceberg
    RUNTIME
    DESTINATION bin)

