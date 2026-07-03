#Change this to whatever cmake version you have
cmake_minimum_required(VERSION 3.18) 

#Change the name and the version as you like
project(ICEBERG
    VERSION 0.0.2
    #Add C, CXX etc if interface moves to the same project
    LANGUAGES C Fortran)


    


include(CMakePrintHelpers)
cmake_print_variables(PROJECT_NAME)
cmake_print_variables(CMAKE_INSTALL_PREFIX)
cmake_print_variables(PROJECT_SOURCE_DIR)
cmake_print_variables(CMAKE_BUILD_TYPE)



#Add the option to turn tests off i.e. cmake -DENABLE_TESTS=OFF ....
option(ENABLE_TESTS "Build unit tests" ON)

set(CMAKE_POSITION_INDEPENDENT_CODE ON)

if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    #Set some compiler flags
    message(STATUS "setting compile/link options for GNU fortran compiler")    
    set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -ffree-line-length-256")
    add_compile_options(
        #Added regardless of build options
        "-I/usr/include;-cpp;-Wall;-Wno-uninitialized"
        #For CMAKE_BUILD_TYPE=Release
        "$<$<CONFIG:Release>:-O2>"
        #For general profiling, prevents inlining, and adds debug info but still optimised
        "$<$<CONFIG:Perf>:-O2;-fno-inline;-fno-omit-frame-pointer;-g;-p;-pg>"
        #For coverage testing
        "$<$<CONFIG:Coverage>:-Og;-g;-coverage>"
        #For gprof, prevents inlining, and adds debug info but still optimised
        "$<$<CONFIG:Gprof>:-Og;-fno-inline;-fno-omit-frame-pointer;-g>"
        #debug
        "$<$<CONFIG:Debug>:-Og;-g;-fno-inline;-fbacktrace;-fbounds-check;-ffpe-trap=invalid,zero,overflow>"
        )
    add_link_options( 
        "$<$<CONFIG:Debug>:-fbacktrace;-fbounds-check;-ffpe-trap=invalid,zero,overflow;-g;-O0;-Wsurprising>"
        "$<$<CONFIG:Gprof>:-p;-pg>"
        "$<$<CONFIG:Coverage>:-coverage>"
         )
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")    
    #Add flags for intel (am assuming the old and LLVM compilers have the same flags 
    #I think can use STREQUAL "Intel" and STREQUAL "IntelLLVM" to distinguish
    #if needed
endif()

if(CMAKE_C_COMPILER_ID STREQUAL "GNU")
    #Set some compiler flags
    message(STATUS "setting compile/link options for GNU C compiler")    
endif()



#Normally you want to add -DCMAKE_INSTALL_PREFIX=/path/to/install/location
#but this will default to the build directory for convenience
if (CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    set(CMAKE_INSTALL_PREFIX ${CMAKE_BINARY_DIR} CACHE PATH "..." FORCE)
endif()


include_directories(src)



add_subdirectory(src)
add_subdirectory(bin)


#   set an environment variable so that the tests can find the necessary data files

if (ENABLE_TESTS)
    include(CTest)
    enable_testing()
    add_subdirectory(test)
endif()


#This stuff is boiler plate to export the details of
#the project to the install directory
#once done another cmake project can
#easily link to the library
include(CMakePackageConfigHelpers)
write_basic_package_version_file(
    ${PROJECT_NAME}ConfigVersion.cmake
    VERSION ${PACKAGE_VERSION}
    COMPATIBILITY AnyNewerVersion)

install(FILES 
    ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake
    DESTINATION lib/cmake)

install(EXPORT ${PROJECT_NAME}_Targets
    FILE ${PROJECT_NAME}Config.cmake
    NAMESPACE
    ${PROJECT_NAME}::
    DESTINATION
    lib/cmake)
export(TARGETS FitLogNormal FILE ${PROJECT_NAME}Config.cmake)
