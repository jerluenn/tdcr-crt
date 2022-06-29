cmake_minimum_required(VERSION 3.10)
project(tdcr-crt)

set (CMAKE_CXX_STANDARD 11)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()
# set(CMAKE_BUILD_TYPE Debug)
SET(CMAKE_CXX_FLAGS "-Wall -Wextra -fPIC -fopenmp  -DBOOST_LOG_DYN_LINK")
SET(CMAKE_CXX_FLAGS_DEBUG "-g")
SET(CMAKE_CXX_FLAGS_RELEASE "-O3")

SET( EXECUTABLE_OUTPUT_PATH ${CMAKE_CURRENT_SOURCE_DIR}/bin )
SET( LIB_PATH ${CMAKE_CURRENT_SOURCE_DIR}/lib/shared )
SET( SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src )
SET( TEST_DIR ${CMAKE_CURRENT_SOURCE_DIR}/tests )
SET( CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/lib/static)
SET( CMAKE_LIBRARY_OUTPUT_DIRECTORY ${LIB_PATH} ) 
SET( acados_include "~/acados/include" )
SET( acados_lib "~/acados/lib")
SET( c_functions_controller "${SOURCE_DIR}/tdcr_controller/c_generated_code_tdcr_lmpc")
SET( c_functions_integration "${SOURCE_DIR}/tdcr_model/c_generated_code_multistage_straight")

MACRO(SUBDIRLIST result curdir)
  FILE(GLOB children RELATIVE ${curdir} ${curdir}/*)
  SET(dirlist "")
  FOREACH(child ${children})
    IF(IS_DIRECTORY ${curdir}/${child})
      LIST(APPEND dirlist ${child})
    ENDIF()
  ENDFOREACH()
  SET(${result} ${dirlist})
ENDMACRO()

find_package(Eigen3 3.4 REQUIRED)

# find_package(Eigen3)
if(NOT EIGEN3_FOUND)
  # Fallback to cmake_modules
  find_package(cmake_modules REQUIRED)
  find_package(Eigen REQUIRED)
  set(EIGEN3_INCLUDE_DIRS ${EIGEN_INCLUDE_DIRS})
  set(EIGEN3_LIBRARIES ${EIGEN_LIBRARIES})
endif()
find_package(PythonLibs REQUIRED)
find_package(catkin REQUIRED)
## System dependencies are found with CMake's conventions
find_package(Boost REQUIRED COMPONENTS system log_setup log REQUIRED)

include_directories(

  ${catkin_INCLUDE_DIRS}
	${SOURCE_DIR}/tdcr_model
	${acados_include}
	${acados_include}/blasfeo/include/
	${acados_include}/hpipm/include/
	${acados_include}/acados/
	${acados_include}/qpOASES_e/
	${PROJECT_SOURCE_DIR}/utils
	${LIB_PATH}
	${PYTHON_INCLUDE_DIRS}

)


SUBDIRLIST(C_FUNC_INT_LIST ${c_functions_integration})
SUBDIRLIST(C_FUNC_CON_LIST ${c_functions_controller})
FILE( GLOB TEST_FILES ${TEST_DIR}/*.cpp)

FOREACH(SUBDIR ${C_FUNC_INT_LIST})

  SUBDIRLIST(SUBSUBDIRS ${c_functions_integration}/${SUBDIR})

  FOREACH(SUBSUBDIR ${SUBSUBDIRS})

    INCLUDE_DIRECTORIES(${c_functions_integration}/${SUBDIR}/${SUBSUBDIR})
    INCLUDE_DIRECTORIES(${c_functions_integration}/${SUBDIR})

  ENDFOREACH()

ENDFOREACH()

FOREACH(SUBDIR ${C_FUNC_CON_LIST})

  SUBDIRLIST(SUBSUBDIRS ${c_functions_controller}/${SUBDIR})

  FOREACH(SUBSUBDIR ${SUBSUBDIRS})

    INCLUDE_DIRECTORIES(${c_functions_controller}/${SUBDIR}/${SUBSUBDIR})
    INCLUDE_DIRECTORIES(${c_functions_controller}/${SUBDIR})

  ENDFOREACH()

ENDFOREACH()

ADD_LIBRARY(MathUtils ${CMAKE_CURRENT_SOURCE_DIR}/utils/MathUtils.cpp)
TARGET_LINK_LIBRARIES(MathUtils 
PRIVATE Eigen3::Eigen 
${PYTHON_LIBRARIES} 
${acados_lib}/libacados.so)

ADD_LIBRARY(ControllerInterface ${CMAKE_CURRENT_SOURCE_DIR}/utils/ControllerInterface.cpp)
TARGET_LINK_LIBRARIES(ControllerInterface 
PUBLIC MathUtils
PRIVATE Eigen3::Eigen 
${PYTHON_LIBRARIES} 
${acados_lib}/libacados.so)

ADD_LIBRARY(IntegratorInterface ${CMAKE_CURRENT_SOURCE_DIR}/utils/IntegratorInterface.cpp)
TARGET_LINK_LIBRARIES(IntegratorInterface 
PUBLIC MathUtils
PRIVATE Eigen3::Eigen 
${PYTHON_LIBRARIES} 
${acados_lib}/libacados.so)

ADD_LIBRARY(MultistageTDCR_Solver ${CMAKE_CURRENT_SOURCE_DIR}/src/tdcr_model/MultistageTDCR_Solver.cpp)
TARGET_LINK_LIBRARIES(MultistageTDCR_Solver 
PUBLIC MathUtils
PUBLIC IntegratorInterface 
PUBLIC ControllerInterface 
PRIVATE Eigen3::Eigen 
${PYTHON_LIBRARIES} 
${acados_lib}/libacados.so)

ADD_LIBRARY(LevenbergMarquardtFunctor ${CMAKE_CURRENT_SOURCE_DIR}/utils/LevenbergMarquardtFunctor.cpp)
TARGET_LINK_LIBRARIES(LevenbergMarquardtFunctor
PUBLIC MultistageTDCR_Solver 
PRIVATE Eigen3::Eigen 
${PYTHON_LIBRARIES} 
${acados_lib}/libacados.so)

ADD_LIBRARY(TDCR_Interface ${CMAKE_CURRENT_SOURCE_DIR}/src/tdcr_model/TDCR_Interface.cpp)
TARGET_LINK_LIBRARIES(TDCR_Interface
PUBLIC LevenbergMarquardtFunctor
PUBLIC MultistageTDCR_Solver 
PRIVATE Eigen3::Eigen 
${PYTHON_LIBRARIES} 
${acados_lib}/libacados.so)


FOREACH(TEST_FILE ${TEST_FILES})

  FILE( GLOB LIBFILES ${LIB_PATH}/*.so)
	GET_FILENAME_COMPONENT( EXEC_NAME ${TEST_FILE} NAME_WE)
	ADD_EXECUTABLE(${EXEC_NAME} ${TEST_FILE})
	TARGET_LINK_LIBRARIES(${EXEC_NAME} 
  PUBLIC TDCR_Interface
  PUBLIC LevenbergMarquardtFunctor
  PUBLIC IntegratorInterface 
  PUBLIC ControllerInterface 
  PUBLIC MultistageTDCR_Solver 
  PRIVATE Eigen3::Eigen 
  ${PYTHON_LIBRARIES} 
  ${acados_lib}/libacados.so
  ${LIBFILES}
)

ENDFOREACH(TEST_FILE ${TEST_FILES})

## Mark executables and/or libraries for installation
INSTALL(TARGETS tdcr_test MathUtils TDCR_Interface LevenbergMarquardtFunctor IntegratorInterface ControllerInterface MultistageTDCR_Solver
  ARCHIVE DESTINATION ${CATKIN_PACKAGE_LIB_DESTINATION}
  LIBRARY DESTINATION ${CATKIN_PACKAGE_LIB_DESTINATION}
  RUNTIME DESTINATION ${CATKIN_PACKAGE_BIN_DESTINATION}
)