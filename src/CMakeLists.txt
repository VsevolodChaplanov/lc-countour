set(PDE_EXE pde)

set(PDE_EXE_ROOT_DIR ${CMAKE_CURRENT_LIST_DIR})

file(GLOB_RECURSE SOURCES ${PDE_EXE_ROOT_DIR}/include/*)

add_executable(${PDE_EXE} ${SOURCES} main.cpp)
add_executable(test ${SOURCES} test.cpp)

target_include_directories(test PUBLIC ${PDE_EXE_ROOT_DIR}/include/)
target_include_directories(${PDE_EXE} PUBLIC ${PDE_EXE_ROOT_DIR}/include/)