CMAKE_MINIMUM_REQUIRED(VERSION 3.17)

SET(CMAKE_C_COMPILER "gcc")
SET(CMAKE_CXX_COMPILER "clang++")
PROJECT(CMakeATAT
	VERSION 1.0
	LANGUAGES C CXX
	)
ADD_COMPILE_OPTIONS(-Wall -Wpedantic)

SET(CMAKE_CXX_STANDARD 17)
SET(CMAKE_EXPORT_COMPILE_COMMANDS ON )

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

INCLUDE(FetchContent)
SET(FETCHCONTENT_QUIET FALSE)
FIND_PACKAGE(Git REQUIRED)
FIND_PACKAGE(GSL REQUIRED)

FIND_PACKAGE(Eigen3 NO_MODULE)
IF(${Eigen3_FOUND})
	MESSAGE(STATUS "Found Eigen")
ELSE()
	MESSAGE(STATUS "Eigen3 Not Found. Downloading from Source")
	FETCHCONTENT_DECLARE(Eigen3
		GIT_REPOSITORY    https://gitlab.com/libeigen/eigen.git
		GIT_PROGRESS		  TRUE
		)
	IF(NOT eigen3_POPULATED)
		FETCHCONTENT_POPULATE(Eigen3)
		FETCHCONTENT_MAKEAVAILABLE(Eigen3)
		ADD_SUBDIRECTORY(${eigen3_SOURCE_DIR} ${eigen3_BUILD_DIR})
	ENDIF()
ENDIF()

MESSAGE(STATUS "Downloading IFOPT (Eigen wrapper with IPOPT)")
FETCHCONTENT_DECLARE(ifopt
	GIT_REPOSITORY    https://github.com/ethz-adrl/ifopt
	GIT_PROGRESS		  TRUE
	)
IF(NOT ifopt_POPULATED)
	FETCHCONTENT_POPULATE(ifopt)
	FETCHCONTENT_MAKEAVAILABLE(ifopt)
	ADD_SUBDIRECTORY(${ifopt_SOURCE_DIR} ${ifopt_BUILD_DIR})
ENDIF()

MESSAGE(STATUS "This is ${ifopt_SOURCE_DIR}")
IF( EXISTS "${CMAKE_CURRENT_BINARY_DIR}/compile_commands.json" )
  EXECUTE_PROCESS( COMMAND ${CMAKE_COMMAND} -E copy_if_different
    ${CMAKE_CURRENT_BINARY_DIR}/compile_commands.json
    ${CMAKE_CURRENT_SOURCE_DIR}/compile_commands.json
  )
ENDIF()

INCLUDE_DIRECTORIES(${PROJECT_SOURCE_DIR}/include ${ifopt_SOURCE_DIR}
	${ifopt_BUILD_DIR})

ADD_LIBRARY(cvmclushelp SHARED ${PROJECT_SOURCE_DIR}/src/cvmclushelp.c++)
ADD_LIBRARY(corrdumphelp SHARED ${PROJECT_SOURCE_DIR}/src/corrdumphelp.c++)
ADD_LIBRARY(gensqshelp SHARED ${PROJECT_SOURCE_DIR}/src/gensqshelp.c++)
ADD_LIBRARY(mcsqshelp SHARED ${PROJECT_SOURCE_DIR}/src/mcsqshelp.c++)
ADD_LIBRARY(mapshelp SHARED ${PROJECT_SOURCE_DIR}/src/mapshelp.c++)
ADD_LIBRARY(phbhelp SHARED ${PROJECT_SOURCE_DIR}/src/phbhelp.c++)
ADD_LIBRARY(apbhelp SHARED ${PROJECT_SOURCE_DIR}/src/apbhelp.c++)

ADD_LIBRARY(linearops SHARED 
	${PROJECT_SOURCE_DIR}/src/linalg.c++
	${PROJECT_SOURCE_DIR}/src/linsolve.c++
	${PROJECT_SOURCE_DIR}/src/lstsqr.c++
	)

ADD_LIBRARY(curvefit SHARED 
	${PROJECT_SOURCE_DIR}/src/curvefit.c++
	${PROJECT_SOURCE_DIR}/src/thermofunctions.c++
	)
TARGET_LINK_LIBRARIES(curvefit PUBLIC GSL::gsl	GSL::gslcblas)

ADD_LIBRARY(cvm SHARED
	${PROJECT_SOURCE_DIR}/src/CVMDataHolder.c++
	${PROJECT_SOURCE_DIR}/src/CVMModel.c++
	${PROJECT_SOURCE_DIR}/src/CVMLogger.c++
	)
TARGET_LINK_LIBRARIES(cvm PUBLIC Eigen3::Eigen ifopt_core ifopt_ipopt findsym calccorr linearops
	clus_str curvefit)

ADD_LIBRARY(parseops SHARED
	${PROJECT_SOURCE_DIR}/src/stringo.c++
	${PROJECT_SOURCE_DIR}/src/integer.c++
	${PROJECT_SOURCE_DIR}/src/parse.c++
	${PROJECT_SOURCE_DIR}/src/getvalue.c++
	)

ADD_LIBRARY(xtalutil SHARED ${PROJECT_SOURCE_DIR}/src/xtalutil.c++)

ADD_LIBRARY(findsym SHARED ${PROJECT_SOURCE_DIR}/src/findsym.c++)

ADD_LIBRARY(calccorr SHARED ${PROJECT_SOURCE_DIR}/src/calccorr.c++)
TARGET_LINK_LIBRARIES(calccorr PUBLIC linearops clus_str)

ADD_LIBRARY(gstate SHARED ${PROJECT_SOURCE_DIR}/src/gstate.c++)

ADD_LIBRARY(refine SHARED ${PROJECT_SOURCE_DIR}/src/refine.c++)
TARGET_LINK_LIBRARIES(refine PUBLIC calccorr gstate)

ADD_LIBRARY(clus_str SHARED ${PROJECT_SOURCE_DIR}/src/clus_str.c++)

ADD_LIBRARY(lattype SHARED ${PROJECT_SOURCE_DIR}/src/lattype.c++)

ADD_EXECUTABLE(cvmclus ${PROJECT_SOURCE_DIR}/src/cvmclus.c++)
TARGET_LINK_LIBRARIES(cvmclus PUBLIC 
	cvmclushelp
	ifopt_core
	ifopt_ipopt
	parseops
	xtalutil
	cvm
	)

ADD_EXECUTABLE(corrdump ${PROJECT_SOURCE_DIR}/src/corrdump.c++)
TARGET_LINK_LIBRARIES(corrdump PUBLIC
	corrdumphelp
	parseops
	xtalutil
	findsym
	calccorr
	)
	
ADD_EXECUTABLE(gensqs ${PROJECT_SOURCE_DIR}/src/gensqs.c++)
TARGET_LINK_LIBRARIES(gensqs PUBLIC 
	gensqshelp
	calccorr
	parseops
	xtalutil
	findsym
	lattype
	)

ADD_EXECUTABLE(mcsqs ${PROJECT_SOURCE_DIR}/src/mcsqs.c++)
TARGET_LINK_LIBRARIES(mcsqs PUBLIC 
	mcsqshelp
	parseops
	xtalutil
	findsym
	calccorr
	lattype
	)

ADD_EXECUTABLE(maps ${PROJECT_SOURCE_DIR}/src/maps.c++)
TARGET_LINK_LIBRARIES(maps PUBLIC
	mapshelp
	parseops
	xtalutil
	findsym
	lattype
	calccorr
	refine
	)
	
ADD_EXECUTABLE(checkcell ${PROJECT_SOURCE_DIR}/src/checkcell.c++)
TARGET_LINK_LIBRARIES(checkcell PUBLIC parseops xtalutil)

ADD_EXECUTABLE(fixcell ${PROJECT_SOURCE_DIR}/src/fixcell.c++)
TARGET_LINK_LIBRARIES(fixcell PUBLIC parseops xtalutil lattype)

ADD_EXECUTABLE(nntouch ${PROJECT_SOURCE_DIR}/src/nntouch.c++)
TARGET_LINK_LIBRARIES(nntouch PUBLIC parseops xtalutil)

ADD_EXECUTABLE(cellcvrt ${PROJECT_SOURCE_DIR}/src/nntouch.c++)
TARGET_LINK_LIBRARIES(cellcvrt PUBLIC parseops xtalutil lattype findsym)

ADD_EXECUTABLE(calces ${PROJECT_SOURCE_DIR}/src/calces.c++)
TARGET_LINK_LIBRARIES(calces PUBLIC parseops xtalutil)

ADD_EXECUTABLE(apb ${PROJECT_SOURCE_DIR}/src/apb.c++)
TARGET_LINK_LIBRARIES(apb PUBLIC apbhelp parseops xtalutil)

ADD_EXECUTABLE(plotcurv ${PROJECT_SOURCE_DIR}/src/plotcurv.c++)
TARGET_LINK_LIBRARIES(plotcurv PUBLIC linearops parseops xtalutil)

ENABLE_TESTING()

