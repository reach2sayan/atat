CMAKE_MINIMUM_REQUIRED(VERSION 3.17)

SET(CMAKE_C_COMPILER "gcc")
SET(CMAKE_CXX_COMPILER "clang++")
PROJECT(ATAT
	VERSION 1.0
	LANGUAGES C CXX
	)
ADD_COMPILE_OPTIONS(-Wall -Wpedantic)

SET(CMAKE_CXX_STANDARD 17)
SET(CMAKE_EXPORT_COMPILE_COMMANDS ON )

SET(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
SET(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
SET(CMAKE_PYTHON_DIRECTORY ${PROJECT_SOURCE_DIR}/python)


INCLUDE(FetchContent)
SET(FETCHCONTENT_QUIET FALSE)

FIND_PACKAGE(Git REQUIRED)

OPTION(USEGSL "Build with Python" OFF)
IF(USEGSL)
	ADD_DEFINITIONS(-DUSE_GSL)
	FIND_PACKAGE(GSL REQUIRED)
ENDIF()

OPTION(USEPYTHON "Build with Python" OFF)
IF(USEPYTHON)
	ADD_DEFINITIONS(-DUSE_PYTHON)
	FIND_PACKAGE(PythonLibs REQUIRED)
	FIND_PACKAGE(Python COMPONENTS Interpreter Development)
	INCLUDE_DIRECTORIES(${PYTHON_INCLUDE_DIRS})
	FIND_PROGRAM(PYTHON "python")

	IF(PYTHON)
		SET(SETUP_PY_IN "${PROJECT_SOURCE_DIR}/python/setup.py")
		SET(SETUP_PY    "${CMAKE_CURRENT_BINARY_DIR}/setup.py")
		SET(DEPS        "${PROJECT_SOURCE_DIR}/python/__init__.py")
		SET(OUTPUT      "${CMAKE_BINARY_DIR}")
		
		CONFIGURE_FILE(${SETUP_PY_IN} ${SETUP_PY})
		ADD_CUSTOM_COMMAND(OUTPUT ${OUTPUT}
			COMMAND ${PYTHON}
			ARGS setup.py
			build DEPENDS ${DEPS})

		EXECUTE_PROCESS(COMMAND ${PYTHON} ${SETUP_PY} install)
	ENDIF()

	FIND_PACKAGE(pybind11 CONFIG REQUIRED)
ENDIF()

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
ADD_LIBRARY(mcsqshelp SHARED ${PROJECT_SOURCE_DIR}/src/mcsqshelp.c++)
ADD_LIBRARY(gensqshelp SHARED ${PROJECT_SOURCE_DIR}/src/gensqshelp.c++)
ADD_LIBRARY(mapshelp SHARED ${PROJECT_SOURCE_DIR}/src/mapshelp.c++)
ADD_LIBRARY(phbhelp SHARED	${PROJECT_SOURCE_DIR}/src/phbhelp.c++)
ADD_LIBRARY(apbhelp SHARED	${PROJECT_SOURCE_DIR}/src/apbhelp.c++)
ADD_LIBRARY(mmapshelp SHARED ${PROJECT_SOURCE_DIR}/src/mmapshelp.c++)
ADD_LIBRARY(emc2help SHARED ${PROJECT_SOURCE_DIR}/src/emc2help.c++)
ADD_LIBRARY(gcehelp SHARED ${PROJECT_SOURCE_DIR}/src/gcehelp.c++)

ADD_LIBRARY(linearops SHARED
	${PROJECT_SOURCE_DIR}/src/linalg.c++
	${PROJECT_SOURCE_DIR}/src/linsolve.c++
	${PROJECT_SOURCE_DIR}/src/lstsqr.c++
	)

IF(USEGSL)
	ADD_LIBRARY(gslcurvefit SHARED
		${PROJECT_SOURCE_DIR}/src/gslcurvefit.c++
		${PROJECT_SOURCE_DIR}/src/thermofunctions.c++
		)
	TARGET_LINK_LIBRARIES(gslcurvefit PUBLIC GSL::gsl GSL::gslcblas)

ADD_LIBRARY(cvm SHARED
	${PROJECT_SOURCE_DIR}/src/CVMDataHolder.c++
	${PROJECT_SOURCE_DIR}/src/CVMModel.c++
	${PROJECT_SOURCE_DIR}/src/CVMLogger.c++
	)
TARGET_LINK_LIBRARIES(cvm PUBLIC 
	Eigen3::Eigen
	ifopt_core
	ifopt_ipopt
	findsym
	calccorr
	linearops
	clus_str
	gslcurvefit)
ENDIF()

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

ADD_LIBRARY(mrefine SHARED ${PROJECT_SOURCE_DIR}/src/mrefine.c++)

ADD_LIBRARY(clus_str SHARED ${PROJECT_SOURCE_DIR}/src/clus_str.c++)

ADD_LIBRARY(lattype SHARED ${PROJECT_SOURCE_DIR}/src/lattype.c++)

ADD_LIBRARY(chull SHARED ${PROJECT_SOURCE_DIR}/src/chull.c++)

ADD_LIBRARY(mpiinterf SHARED ${PROJECT_SOURCE_DIR}/src/mpiinterf.c++)

ADD_LIBRARY(eci SHARED ${PROJECT_SOURCE_DIR}/src/teci.c++
	${PROJECT_SOURCE_DIR}/src/keci.c++)

ADD_LIBRARY(drawpd SHARED ${PROJECT_SOURCE_DIR}/src/drawpd.c++)
ADD_LIBRARY(mclib SHARED ${PROJECT_SOURCE_DIR}/src/mclib.c++)
ADD_LIBRARY(fft SHARED ${PROJECT_SOURCE_DIR}/src/anyfft.c++ ${PROJECT_SOURCE_DIR}/src/fftn.c++)

ADD_LIBRARY(kspacees SHARED ${PROJECT_SOURCE_DIR}/src/kspacees.c++)

ADD_LIBRARY(pred SHARED
	${PROJECT_SOURCE_DIR}/src/predcs.c++
	${PROJECT_SOURCE_DIR}/src/predes.c++
	${PROJECT_SOURCE_DIR}/src/predrs.c++
	)

ADD_EXECUTABLE(cvmclus ${PROJECT_SOURCE_DIR}/src/cvmclus.c++)
SET(CVMCLUSLIBS cvmclushelp ifopt_core ifopt_ipopt parseops xtalutil cvm)
IF(USEPYTHON)
	LIST(APPEND CVMCLUSLIBS pybind11::embed)
ENDIF()
TARGET_LINK_LIBRARIES(cvmclus PUBLIC ${CVMCLUSLIBS})

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

ADD_EXECUTABLE(genstr ${PROJECT_SOURCE_DIR}/src/genstr.c++)
TARGET_LINK_LIBRARIES(genstr PUBLIC 
	parseops
	xtalutil
	findsym
	clus_str
	lattype
	mpiinterf
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
	mpiinterf
	)

ADD_EXECUTABLE(mmaps ${PROJECT_SOURCE_DIR}/src/mmaps.c++)
TARGET_LINK_LIBRARIES(mmaps PUBLIC
	mmapshelp
	chull
	parseops
	xtalutil
	calccorr
	lattype
	mrefine
	clus_str
	linearops
	pred
	findsym
	mpiinterf
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


ADD_EXECUTABLE(emc2 ${PROJECT_SOURCE_DIR}/src/emc2.c++)
TARGET_LINK_LIBRARIES(emc2 PUBLIC
	emc2help
	parseops
	findsym
	linearops
	clus_str
	lattype
	xtalutil
	calccorr
	drawpd
	mclib
	eci
	fft
	kspacees
	)

ADD_LIBRARY(gceutil SHARED ${PROJECT_SOURCE_DIR}/src/gceutil.c++
	${PROJECT_SOURCE_DIR}/src/tensorsym.c++
	${PROJECT_SOURCE_DIR}/src/phonlib.c++
	${PROJECT_SOURCE_DIR}/src/multipoly.c++
	)

ADD_EXECUTABLE(viewgce ${PROJECT_SOURCE_DIR}/src/viewgce.c++)
TARGET_LINK_LIBRARIES(viewgce PUBLIC
	corrdumphelp
	parseops
	linearops
	gceutil
	xtalutil
	findsym
	calccorr
	clus_str
	)

ADD_EXECUTABLE(gce ${PROJECT_SOURCE_DIR}/src/gce.c++)
TARGET_LINK_LIBRARIES(gce PUBLIC
	gcehelp
	linearops
	gceutil
	parseops
	xtalutil
	calccorr
	findsym
	clus_str)

ADD_EXECUTABLE(gencs ${PROJECT_SOURCE_DIR}/src/gencs.c++)
TARGET_LINK_LIBRARIES(gencs PUBLIC
	linearops
	parseops
	xtalutil
	findsym
	gceutil
	)

ENABLE_TESTING()

