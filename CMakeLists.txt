project(ThermoPoroElasticity)

cmake_minimum_required(VERSION 2.6)

# Add local modules
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH}
                      ${PROJECT_SOURCE_DIR}/../../cmake/Modules)

if(NOT IFEM_CONFIGURED)
  find_package(IFEM REQUIRED)
endif()

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${IFEM_CXX_FLAGS}")

# Generate regtest script with correct paths
configure_file(${PROJECT_SOURCE_DIR}/../Common/scripts/regtest.sh.in
               regtest.sh)

set(ThermoPoroElasticity_SRCS main_ThermoPoroElasticity.C
             ThermoPoroElasticity.C PoroMaterial.C)

include_directories(${IFEM_INCLUDES} ${PROJECT_SOURCE_DIR}
                    ../Common)

set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR}/bin)

if(NOT TARGET CommonIFEM)
  add_subdirectory(../Common CommonIFEM.dir)
endif()

add_executable(ThermoPoroElasticity ${ThermoPoroElasticity_SRCS})
target_link_libraries(ThermoPoroElasticity CommonIFEM ${IFEM_LIBRARIES})

# Installation
install(TARGETS ThermoPoroElasticity DESTINATION bin COMPONENT bin)

enable_testing()

# Regression tests
if(IFEM_USE_PARALLEL_PETSC)
  # Add parallel tests here
else()
  # Add serial tests here
endif()

list(APPEND TEST_APPS ThermoPoroElasticity)

# Unit tests
#IFEM_add_test_app(${PROJECT_SOURCE_DIR}/Test/*.C
#                  ${PROJECT_SOURCE_DIR}/Test
#                  Poisson
#                  CommonPoisson ${IFEM_LIBRARIES})

if(IFEM_COMMON_APP_BUILD)
  set(TEST_APPS ${TEST_APPS} PARENT_SCOPE)
  set(UNIT_TEST_NUMBER ${UNIT_TEST_NUMBER} PARENT_SCOPE)
else()
  add_check_target()
endif()
