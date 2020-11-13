# --------------------------------------------------------------------------- #
#    CMake file for unit tests                                                #
# --------------------------------------------------------------------------- #

# ----- Fetch GoogleTest ---------------------------------------------------- #
include(FetchContent)
set(INSTALL_GTEST OFF)
fetchcontent_declare(
        googletest
        GIT_REPOSITORY https://github.com/google/googletest.git
        GIT_TAG release-1.10.0
)

fetchcontent_getproperties(googletest)
if (NOT googletest_POPULATED)
    if (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
        set(CMAKE_MACOSX_RPATH 1)
    endif ()
    fetchcontent_populate(googletest)
    add_subdirectory(${googletest_SOURCE_DIR} ${googletest_BINARY_DIR})
endif ()

# ----- tests_ThermalUnitBlock ---------------------------------------------- #
add_executable(tests_ThermalUnitBlock tests_ThermalUnitBlock.cpp)
target_link_libraries(tests_ThermalUnitBlock PRIVATE SMS++::UCBlock gtest)

add_test(
        NAME tests_ThermalUnitBlock
        COMMAND tests_ThermalUnitBlock
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)
