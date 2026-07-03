

#Configure the install
 
find_package(LAPACK REQUIRED)
 

#Split the command line parsing code into its own lib
#And do the same thing
add_library(FitLogNormal STATIC
    NBAX_StringTokenizers.f90
    Lib_CommandLineArguments.f90
    Lib_ColouredTerminal.f90
    Lib_LogNormal.f90
    Lib_Linint.f90
    Lib_ChiSquare.f90
    Lib_SafeExp.F90
    Lib_LogisticFunction.f90
    Lib_RandomSeed.f90
    Lib_Quicksort.f90
    Lib_DownhillSimplex.f90
    )
set_target_properties(FitLogNormal PROPERTIES 
    Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/include)
target_include_directories(FitLogNormal
    PUBLIC
    $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/include>
    $<INSTALL_INTERFACE:include>)

    


#Copy libraries to the lib directory of
#the install directory

install(TARGETS FitLogNormal
    EXPORT ${PROJECT_NAME}_Targets
    ARCHIVE
    DESTINATION lib
    INCLUDES 
    DESTINATION
    include)
     

#...and the fortran module to the include directory
install(DIRECTORY
    ${CMAKE_BINARY_DIR}/include
    DESTINATION
    ${CMAKE_INSTALL_PREFIX})


 