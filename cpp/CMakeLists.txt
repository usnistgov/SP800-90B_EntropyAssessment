include(GNUInstallDirs)

find_package(OpenMP COMPONENTS CXX REQUIRED)
find_package(PkgConfig REQUIRED)

pkg_check_modules(BASE_DEPS REQUIRED IMPORTED_TARGET bzip2 jsoncpp libcrypto libdivsufsort libdivsufsort64)
pkg_check_modules(EXTRA_COND_DEPS REQUIRED IMPORTED_TARGET gmp mpfr)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2")
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED On)

add_executable(ea_conditioning conditioning_main.cpp)
set_target_properties(ea_conditioning PROPERTIES LINKER_LANGUAGE CXX)
target_link_libraries(ea_conditioning PRIVATE PkgConfig::BASE_DEPS PkgConfig::EXTRA_COND_DEPS OpenMP::OpenMP_CXX)

add_executable(ea_iid iid_main.cpp)
set_target_properties(ea_iid PROPERTIES LINKER_LANGUAGE CXX)
target_link_libraries(ea_iid PRIVATE PkgConfig::BASE_DEPS OpenMP::OpenMP_CXX)

add_executable(ea_non_iid non_iid_main.cpp)
set_target_properties(ea_non_iid PROPERTIES LINKER_LANGUAGE CXX)
target_link_libraries(ea_non_iid PRIVATE PkgConfig::BASE_DEPS OpenMP::OpenMP_CXX)

add_executable(ea_restart restart_main.cpp)
set_target_properties(ea_restart PROPERTIES LINKER_LANGUAGE CXX)
target_link_libraries(ea_restart PRIVATE PkgConfig::BASE_DEPS OpenMP::OpenMP_CXX)

add_executable(ea_transpose transpose_main.cpp)
set_target_properties(ea_transpose PROPERTIES LINKER_LANGUAGE CXX)
target_link_libraries(ea_transpose PRIVATE PkgConfig::BASE_DEPS OpenMP::OpenMP_CXX)

install(
	TARGETS ea_conditioning ea_iid ea_non_iid ea_restart ea_transpose
	RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)
