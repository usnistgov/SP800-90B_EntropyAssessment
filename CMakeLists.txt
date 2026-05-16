cmake_minimum_required(VERSION 3.14)

find_package(PkgConfig)

project(SP800-90B_EntropyAssessment CXX)
set(PROJECT_VENDOR "NIST CMVP")
set(PROJECT_CONTACT "esv-prod@nist.gov")
set(PROJECT_URL "https://github.com/usnistgov/SP800-90B_EntropyAssessment")
set(PROJECT_DESCRIPTION "The SP800-90B_EntropyAssessment C++package implements the min-entropy assessment methods included in Special Publication 800-90B.")

add_subdirectory(bin)
add_subdirectory(cpp)
