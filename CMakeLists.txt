cmake_minimum_required(VERSION 3.10)
project(tdcr-crt)

set (CMAKE_CXX_STANDARD 11)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

SET(CMAKE_CXX_FLAGS "-Wall -Wextra -fPIC -fopenmp  -DBOOST_LOG_DYN_LINK")
SET(CMAKE_CXX_FLAGS_DEBUG "-g")
# SET(CMAKE_CXX_FLAGS_RELEASE "-O3")

SET( EXECUTABLE_OUTPUT_PATH ${CMAKE_CURRENT_SOURCE_DIR}/bin )
SET( LIB_PATH ${CMAKE_CURRENT_SOURCE_DIR}/lib )
SET( SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src )
SET( TEST_DIR ${CMAKE_CURRENT_SOURCE_DIR}/tests )
SET( CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/lib/static)
SET( CMAKE_LIBRARY_OUTPUT_DIRECTORY ${LIB_PATH} ) 
SET( acados_include "${CMAKE_CURRENT_SOURCE_DIR}/../acados/include" )
SET( acados_lib "${CMAKE_CURRENT_SOURCE_DIR}/../acados/lib")
SET( c_functions "${SOURCE_DIR}/tendon_robots" )


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
find_package(ifopt)

## System dependencies are found with CMake's conventions
find_package(Boost REQUIRED COMPONENTS system log_setup log REQUIRED)

include_directories(

	${SOURCE_DIR}/tendon_robots
	${acados_include}
	${acados_include}/blasfeo/include/
	${acados_include}/hpipm/include/
	${acados_include}/acados/
	${acados_include}/qpOASES_e/
	${PROJECT_SOURCE_DIR}/utils
	${LIB_PATH}
	${PYTHON_INCLUDE_DIRS}

)


FILE( GLOB SOURCES ${SOURCE_DIR}/tendon_robots/*.cpp )
FILE( GLOB UTIL_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/utils/*.cpp )
FILE( GLOB TEST_FILES ${TEST_DIR}/*.cpp)
SUBDIRLIST(SUBDIRS ${c_functions})

FOREACH(SUBDIR ${SUBDIRS})

	SUBDIRLIST(SUBSUBDIRS ${c_functions}/${SUBDIR})
	
	FOREACH(SUBSUBDIR ${SUBSUBDIRS})

		SUBDIRLIST(SUBSUBSUBDIRS ${c_functions}/${SUBDIR}/${SUBSUBDIR})
	
		include_directories(${c_functions}/${SUBDIR}/${SUBSUBDIR}/${SUBSUBSUBDIR})		

		FOREACH(SUBSUBSUBDIR ${SUBSUBSUBDIRS})

			include_directories(${c_functions}/${SUBDIR}/${SUBSUBDIR}/${SUBSUBSUBDIR})
			# message(STATUS ${c_functions}/${SUBDIR}/${SUBSUBDIR}/${SUBSUBSUBDIR})

		ENDFOREACH()

	ENDFOREACH()

ENDFOREACH(SUBDIR ${c_functions})

# Compiling executables

FOREACH(UTIL_SRC ${UTIL_SOURCES})

	FILE( GLOB LIBFILES ${LIB_PATH}/shared/*.so)
	FILE( GLOB LIBS ${CMAKE_ARCHIVE_OUTPUT_DIRECTORY}/*.a )
	GET_FILENAME_COMPONENT( EXEC_NAME ${UTIL_SRC} NAME_WE )
	ADD_LIBRARY(${EXEC_NAME} ${UTIL_SRC})	
	TARGET_LINK_LIBRARIES(	${EXEC_NAME} 
	PRIVATE Eigen3::Eigen 
	${PYTHON_LIBRARIES} 
	${acados_lib}/libacados.so
	${LIBFILES}
	${LIBS}

)

ENDFOREACH(UTIL_SRC ${UTIL_SOURCES})

FOREACH(SRC ${SOURCES})

	FILE( GLOB LIBFILES ${LIB_PATH}/shared/*.so)
	GET_FILENAME_COMPONENT( EXEC_NAME ${SRC} NAME_WE )
	FILE( GLOB LIBS ${CMAKE_ARCHIVE_OUTPUT_DIRECTORY}/*.a )
	ADD_LIBRARY(${EXEC_NAME} ${SRC})	
	TARGET_LINK_LIBRARIES(${EXEC_NAME} 
	PRIVATE Eigen3::Eigen 
	${PYTHON_LIBRARIES} 
	${acados_lib}/libacados.so
	${LIBFILES}
	${LIBS}

)

ENDFOREACH(SRC ${SOURCES})

FOREACH(TEST_FILE ${TEST_FILES})

	FILE( GLOB LIBFILES ${LIB_PATH}/shared/*.so)
	message(STATUS ${LIBFILES})
	FILE( GLOB LIBS ${CMAKE_ARCHIVE_OUTPUT_DIRECTORY}/*.a )
	message(STATUS ${LIBS})
	GET_FILENAME_COMPONENT( EXEC_NAME ${TEST_FILE} NAME_WE)
	ADD_EXECUTABLE(${EXEC_NAME} ${TEST_FILE})
	TARGET_LINK_LIBRARIES(${EXEC_NAME} 
	PRIVATE Eigen3::Eigen 
	${PYTHON_LIBRARIES} 
	MultistageTDCR_Solver
	${acados_lib}/libacados.so
	${LIBFILES}
	${LIBS}
	Boost::log 
	Boost::log_setup
 
)

ENDFOREACH(TEST_FILE ${TEST_FILES})


