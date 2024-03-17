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

OPTION(MAKENDVIEWER "Build ndviewer? Requires OpenGL" ON)
IF(MAKENDVIEWER)
	ADD_DEFINITIONS(-DMAKE_NDVIEWER)
	FIND_PACKAGE(OpenGL REQUIRED)
	FIND_PACKAGE(GLUT REQUIRED)
	INCLUDE_DIRECTORIES(${OPENGL_INCLUDE_DIRS} ${GLUT_INCLUDE_DIRS})
ENDIF()

INCLUDE(FetchContent)
SET(FETCHCONTENT_QUIET FALSE)

FIND_PACKAGE(Git REQUIRED)

OPTION(USEGSL "Build with GSL" ON)
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


file(GLOB CPP_LIST "${CMAKE_CURRENT_SOURCE_DIR}/*.cpp")

ADD_LIBRARY(itertools INTERFACE)
TARGET_INCLUDE_DIRECTORIES(itertools INTERFACE
	$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/itertools>
	$<INSTALL_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/lib>)

ADD_LIBRARY(cvmclushelp SHARED ${PROJECT_SOURCE_DIR}/src/cvmclushelp.c++)
ADD_LIBRARY(cvhelp SHARED ${PROJECT_SOURCE_DIR}/src/cvhelp.c++)
ADD_LIBRARY(corrdumphelp SHARED ${PROJECT_SOURCE_DIR}/src/corrdumphelp.c++)
ADD_LIBRARY(mcsqshelp SHARED ${PROJECT_SOURCE_DIR}/src/mcsqshelp.c++)
ADD_LIBRARY(gensqshelp SHARED ${PROJECT_SOURCE_DIR}/src/gensqshelp.c++)
ADD_LIBRARY(mapshelp SHARED ${PROJECT_SOURCE_DIR}/src/mapshelp.c++)
ADD_LIBRARY(phbhelp SHARED	${PROJECT_SOURCE_DIR}/src/phbhelp.c++)
ADD_LIBRARY(apbhelp SHARED	${PROJECT_SOURCE_DIR}/src/apbhelp.c++)
ADD_LIBRARY(mmapshelp SHARED ${PROJECT_SOURCE_DIR}/src/mmapshelp.c++)
ADD_LIBRARY(emc2help SHARED ${PROJECT_SOURCE_DIR}/src/emc2help.c++)
ADD_LIBRARY(memc2help SHARED ${PROJECT_SOURCE_DIR}/src/memc2help.c++)
ADD_LIBRARY(gcehelp SHARED ${PROJECT_SOURCE_DIR}/src/gcehelp.c++)
ADD_LIBRARY(fitsvslhelp SHARED ${PROJECT_SOURCE_DIR}/src/fitsvslhelp.c++)
ADD_LIBRARY(svslhelp SHARED ${PROJECT_SOURCE_DIR}/src/svslhelp.c++)
ADD_LIBRARY(fitfchelp SHARED ${PROJECT_SOURCE_DIR}/src/fitfchelp.c++)
ADD_LIBRARY(felechelp SHARED ${PROJECT_SOURCE_DIR}/src/felechelp.c++)
ADD_LIBRARY(fempmaghelp SHARED ${PROJECT_SOURCE_DIR}/src/fempmaghelp.c++)
ADD_LIBRARY(predcshelp SHARED ${PROJECT_SOURCE_DIR}/src/predcshelp.c++)
ADD_LIBRARY(infdethelp SHARED ${PROJECT_SOURCE_DIR}/src/infdethelp.c++)
ADD_LIBRARY(mindisthelp SHARED ${PROJECT_SOURCE_DIR}/src/mindisthelp.c++)
ADD_LIBRARY(epichelp SHARED ${PROJECT_SOURCE_DIR}/src/epichelp.c++)

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
ENDIF()

ADD_LIBRARY(cvm SHARED
	${PROJECT_SOURCE_DIR}/src/CVMDataHolder.c++
	${PROJECT_SOURCE_DIR}/src/CVMModel.c++
	${PROJECT_SOURCE_DIR}/src/CVMLogger.c++
	)
SET(CVMLIBS Eigen3::Eigen
	ifopt_core
	ifopt_ipopt
	findsym
	calccorr
	linearops
	clus_str
	)
IF(USEGSL)
	LIST(APPEND CVMLIBS gslcurvefit)
ENDIF()
TARGET_LINK_LIBRARIES(cvm PUBLIC ${CVMLIBS})

ADD_EXECUTABLE(cvmclus ${PROJECT_SOURCE_DIR}/src/cvmclus.c++)
SET(CVMCLUSLIBS cvmclushelp parseops xtalutil cvm	Eigen3::Eigen)
IF(USEPYTHON)
	LIST(APPEND CVMCLUSLIBS pybind11::embed)
	PYBIND11_ADD_MODULE(pycvm ${PROJECT_SOURCE_DIR}/src/CVMDataHolder.c++)
	TARGET_LINK_LIBRARIES(pycvm PUBLIC cvm)
ENDIF()
TARGET_LINK_LIBRARIES(cvmclus PUBLIC ${CVMCLUSLIBS})

ADD_LIBRARY(parseops SHARED
	${PROJECT_SOURCE_DIR}/src/stringo.c++
	${PROJECT_SOURCE_DIR}/src/integer.c++
	${PROJECT_SOURCE_DIR}/src/parse.c++
	${PROJECT_SOURCE_DIR}/src/getvalue.c++
	)

ADD_LIBRARY(strinterf SHARED ${PROJECT_SOURCE_DIR}/src/strinterf.c++)

ADD_LIBRARY(xtalutil SHARED ${PROJECT_SOURCE_DIR}/src/xtalutil.c++)

ADD_LIBRARY(mesh SHARED ${PROJECT_SOURCE_DIR}/src/meshutil.c++
	${PROJECT_SOURCE_DIR}/src/meshcalc.c++)

ADD_LIBRARY(findsym SHARED ${PROJECT_SOURCE_DIR}/src/findsym.c++)

ADD_LIBRARY(opti SHARED ${PROJECT_SOURCE_DIR}/src/opti.c++)

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
ADD_LIBRARY(normal SHARED ${PROJECT_SOURCE_DIR}/src/normal.c++)
ADD_LIBRARY(mclib SHARED ${PROJECT_SOURCE_DIR}/src/mclib.c++)
ADD_LIBRARY(mmclib SHARED ${PROJECT_SOURCE_DIR}/src/mmclib.c++)
ADD_LIBRARY(calcmf SHARED ${PROJECT_SOURCE_DIR}/src/calcmf.c++)
ADD_LIBRARY(equil SHARED ${PROJECT_SOURCE_DIR}/src/equil.c++)
ADD_LIBRARY(fft SHARED ${PROJECT_SOURCE_DIR}/src/anyfft.c++ ${PROJECT_SOURCE_DIR}/src/fftn.c++)

ADD_LIBRARY(kspacees SHARED ${PROJECT_SOURCE_DIR}/src/kspacees.c++)
ADD_LIBRARY(morsepot SHARED ${PROJECT_SOURCE_DIR}/src/morsepot.c++
	${PROJECT_SOURCE_DIR}/src/morsepotc3.c++)
TARGET_LINK_LIBRARIES(morsepot PUBLIC linearops)
ADD_LIBRARY(kspacecs SHARED ${PROJECT_SOURCE_DIR}/src/kspacecs.c++)
ADD_LIBRARY(tlambda SHARED ${PROJECT_SOURCE_DIR}/src/tlambda.c++)
ADD_LIBRARY(ridge SHARED ${PROJECT_SOURCE_DIR}/src/ridge.c++)

ADD_LIBRARY(pred SHARED
	${PROJECT_SOURCE_DIR}/src/predcs.c++
	${PROJECT_SOURCE_DIR}/src/predes.c++
	${PROJECT_SOURCE_DIR}/src/predrs.c++
	)

ADD_EXECUTABLE(fitsvsl ${PROJECT_SOURCE_DIR}/src/fitsvsl.c++)
TARGET_LINK_LIBRARIES(fitsvsl PUBLIC
	fitsvslhelp
	xtalutil
	linearops
	findsym
	parseops
	clus_str
	lattype
	gceutil
	calccorr
	)

ADD_EXECUTABLE(svsl ${PROJECT_SOURCE_DIR}/src/svsl.c++)
TARGET_LINK_LIBRARIES(svsl PUBLIC
	morsepot
	svslhelp
	xtalutil
	linearops
	parseops
	calccorr
	findsym
	clus_str
	gceutil
	)


ADD_EXECUTABLE(fitfc ${PROJECT_SOURCE_DIR}/src/fitfc.c++)
TARGET_LINK_LIBRARIES(fitfc PUBLIC
	fitfchelp
	xtalutil
	linearops
	parseops
	gceutil
	findsym
	lattype
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
	pred
	tlambda
	ridge
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

ADD_EXECUTABLE(analrelax ${PROJECT_SOURCE_DIR}/src/analrelax.c++)
TARGET_LINK_LIBRARIES(analrelax PUBLIC parseops xtalutil linearops)

ADD_EXECUTABLE(fixcell ${PROJECT_SOURCE_DIR}/src/fixcell.c++)
TARGET_LINK_LIBRARIES(fixcell PUBLIC parseops xtalutil lattype)

ADD_EXECUTABLE(nntouch ${PROJECT_SOURCE_DIR}/src/nntouch.c++)
TARGET_LINK_LIBRARIES(nntouch PUBLIC parseops xtalutil)

ADD_EXECUTABLE(cellcvrt ${PROJECT_SOURCE_DIR}/src/nntouch.c++)
TARGET_LINK_LIBRARIES(cellcvrt PUBLIC parseops xtalutil lattype findsym)

ADD_EXECUTABLE(csfit ${PROJECT_SOURCE_DIR}/src/csfit.c++)
TARGET_LINK_LIBRARIES(csfit PUBLIC predcshelp parseops xtalutil linearops)

ADD_EXECUTABLE(calces ${PROJECT_SOURCE_DIR}/src/calces.c++)
TARGET_LINK_LIBRARIES(calces PUBLIC parseops xtalutil)

ADD_EXECUTABLE(apb ${PROJECT_SOURCE_DIR}/src/apb.c++)
TARGET_LINK_LIBRARIES(apb PUBLIC apbhelp parseops xtalutil)

ADD_EXECUTABLE(plotcurv ${PROJECT_SOURCE_DIR}/src/plotcurv.c++)
TARGET_LINK_LIBRARIES(plotcurv PUBLIC linearops parseops xtalutil)

ADD_EXECUTABLE(cv ${PROJECT_SOURCE_DIR}/src/cv.c++)
TARGET_LINK_LIBRARIES(cv PUBLIC cvhelp xtalutil linearops parseops)

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

ADD_EXECUTABLE(memc2 ${PROJECT_SOURCE_DIR}/src/memc2.c++)
TARGET_LINK_LIBRARIES(memc2 PUBLIC
	memc2help
	parseops
	linearops
	equil
	eci
	drawpd
	clus_str
	gceutil
	chull
	calcmf
	fft
	calccorr
	findsym
	mmclib
	xtalutil
	)


ADD_EXECUTABLE(phb ${PROJECT_SOURCE_DIR}/src/phb.c++)
TARGET_LINK_LIBRARIES(phb PUBLIC
	phbhelp
	parseops
	linearops
	kspacees
	calccorr
	xtalutil
	drawpd
	fft
	eci
	mclib
	findsym
	clus_str
	)

ADD_EXECUTABLE(kmesh ${PROJECT_SOURCE_DIR}/src/kmesh.c++)
TARGET_LINK_LIBRARIES(kmesh PUBLIC
	findsym
	parseops
	xtalutil
	)

ADD_LIBRARY(gceutil SHARED ${PROJECT_SOURCE_DIR}/src/gceutil.c++
	${PROJECT_SOURCE_DIR}/src/tensorsym.c++
	${PROJECT_SOURCE_DIR}/src/phonlib.c++
	${PROJECT_SOURCE_DIR}/src/multipoly.c++
	)

#
#ADD_EXECUTABLE(skelread ${PROJECT_SOURCE_DIR}/src/skelread.c++)
#TARGET_LINK_LIBRARIES(skelread PUBLIC
#	parseops
#	xtalutil
#	findsym
#	calccorr
#	clus_str
#	linearops
#	)
#

ADD_EXECUTABLE(skelreadstr ${PROJECT_SOURCE_DIR}/src/skelreadstr.c++)
TARGET_LINK_LIBRARIES(skelreadstr PUBLIC xtalutil parseops)

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

ADD_EXECUTABLE(calcelas ${PROJECT_SOURCE_DIR}/src/calcelas.c++)
TARGET_LINK_LIBRARIES(calcelas PUBLIC
	xtalutil
	linearops
	findsym
	parseops
	gceutil
	)

ADD_EXECUTABLE(strainder ${PROJECT_SOURCE_DIR}/src/strainder.c++)
TARGET_LINK_LIBRARIES(strainder PUBLIC
	xtalutil
	linearops
	findsym
	parseops
	gceutil
	)

ADD_EXECUTABLE(loctst ${PROJECT_SOURCE_DIR}/src/loctst.c++)
TARGET_LINK_LIBRARIES(loctst PUBLIC xtalutil parseops)

ADD_EXECUTABLE(iwulff ${PROJECT_SOURCE_DIR}/src/iwulff.c++)
TARGET_LINK_LIBRARIES(iwulff PUBLIC xtalutil parseops)

#ADD_EXECUTABLE(clusover ${PROJECT_SOURCE_DIR}/src/clusover.c++)
#TARGET_LINK_LIBRARIES(clusover PUBLIC parseops linearops)

ADD_EXECUTABLE(tryclus ${PROJECT_SOURCE_DIR}/src/tryclus.c++)
TARGET_LINK_LIBRARIES(tryclus PUBLIC xtalutil parseops)

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

ADD_EXECUTABLE(csl ${PROJECT_SOURCE_DIR}/src/csl.c++)
TARGET_LINK_LIBRARIES(csl PUBLIC xtalutil parseops findsym)

ADD_EXECUTABLE(testatat ${PROJECT_SOURCE_DIR}/src/test.c++)
TARGET_LINK_LIBRARIES(testatat PUBLIC linearops parseops chull xtalutil)

ADD_EXECUTABLE(nnshell ${PROJECT_SOURCE_DIR}/src/nnshell.c++)
TARGET_LINK_LIBRARIES(nnshell PUBLIC xtalutil parseops)

ADD_EXECUTABLE(pdef ${PROJECT_SOURCE_DIR}/src/pdef.c++)
TARGET_LINK_LIBRARIES(pdef PUBLIC xtalutil findsym linearops parseops lattype)

ADD_EXECUTABLE(felec ${PROJECT_SOURCE_DIR}/src/felec.c++)
TARGET_LINK_LIBRARIES(felec PUBLIC felechelp xtalutil parseops)

ADD_EXECUTABLE(fmag ${PROJECT_SOURCE_DIR}/src/fmag.c++)
TARGET_LINK_LIBRARIES(fmag PUBLIC xtalutil parseops)

ADD_EXECUTABLE(fempmag ${PROJECT_SOURCE_DIR}/src/fempmag.c++)
TARGET_LINK_LIBRARIES(fempmag PUBLIC fempmaghelp xtalutil parseops)

ADD_EXECUTABLE(cna ${PROJECT_SOURCE_DIR}/src/cna.c++)
TARGET_LINK_LIBRARIES(cna PUBLIC xtalutil parseops)

ADD_EXECUTABLE(cmpstr ${PROJECT_SOURCE_DIR}/src/cmpstr.c++)
TARGET_LINK_LIBRARIES(cmpstr PUBLIC xtalutil parseops)

ADD_EXECUTABLE(ehmass ${PROJECT_SOURCE_DIR}/src/ehmass.c++)
TARGET_LINK_LIBRARIES(ehmass PUBLIC xtalutil parseops linearops findsym)

ADD_EXECUTABLE(seebeck ${PROJECT_SOURCE_DIR}/src/seebeck.c++)
TARGET_LINK_LIBRARIES(seebeck PUBLIC xtalutil parseops linearops findsym)

ADD_EXECUTABLE(bgsdqmc ${PROJECT_SOURCE_DIR}/src/bgsdqmc.c++)
TARGET_LINK_LIBRARIES(bgsdqmc PUBLIC normal)

ADD_EXECUTABLE(vctriv ${PROJECT_SOURCE_DIR}/src/vctriv.c++)
TARGET_LINK_LIBRARIES(vctriv PUBLIC normal xtalutil parseops)

ADD_EXECUTABLE(hessian ${PROJECT_SOURCE_DIR}/src/hessian.c++)
TARGET_LINK_LIBRARIES(hessian PUBLIC xtalutil linearops parseops)

ADD_EXECUTABLE(triph ${PROJECT_SOURCE_DIR}/src/triph.c++)
TARGET_LINK_LIBRARIES(triph PUBLIC xtalutil parseops)

ADD_EXECUTABLE(simplexize ${PROJECT_SOURCE_DIR}/src/simplexize.c++)
TARGET_LINK_LIBRARIES(simplexize PUBLIC xtalutil linearops parseops mesh)

ADD_EXECUTABLE(mkaxes ${PROJECT_SOURCE_DIR}/src/mkaxes.c++)
TARGET_LINK_LIBRARIES(mkaxes PUBLIC xtalutil linearops parseops mesh)

ADD_EXECUTABLE(pathint ${PROJECT_SOURCE_DIR}/src/pathint.c++)
TARGET_LINK_LIBRARIES(pathint PUBLIC xtalutil linearops parseops)

ADD_EXECUTABLE(greadyneb ${PROJECT_SOURCE_DIR}/src/greadyneb.c++)
TARGET_LINK_LIBRARIES(greadyneb PUBLIC linearops parseops xtalutil)

ADD_EXECUTABLE(mpitest ${PROJECT_SOURCE_DIR}/src/mpitest.c++)
TARGET_LINK_LIBRARIES(mpitest PUBLIC parseops xtalutil mpiinterf)

ADD_EXECUTABLE(equivclus ${PROJECT_SOURCE_DIR}/src/equivclus.c++)
TARGET_LINK_LIBRARIES(equivclus PUBLIC parseops findsym xtalutil)

ADD_EXECUTABLE(3d2gnu ${PROJECT_SOURCE_DIR}/src/3d2gnu.c++)
TARGET_LINK_LIBRARIES(3d2gnu PUBLIC xtalutil linearops parseops mesh)

ADD_EXECUTABLE(infdet ${PROJECT_SOURCE_DIR}/src/infdet.c++)
TARGET_LINK_LIBRARIES(infdet PUBLIC 
	infdethelp
	xtalutil
	linearops
	parseops
	opti
	strinterf
	)

ADD_EXECUTABLE(mindist ${PROJECT_SOURCE_DIR}/src/mindist.c++)
TARGET_LINK_LIBRARIES(mindist PUBLIC 
	mindisthelp
	xtalutil
	parseops
	)

ADD_EXECUTABLE(epic ${PROJECT_SOURCE_DIR}/src/epic.c++)
TARGET_LINK_LIBRARIES(epic PUBLIC 
	epichelp
	strinterf
	opti
	parseops
	xtalutil
	linearops
	)

IF(MAKENDVIEWER)
	ADD_EXECUTABLE(ndviewer ${PROJECT_SOURCE_DIR}/src/ndviewer.c++)
	TARGET_LINK_LIBRARIES(ndviewer PUBLIC 
		parseops
		mesh
		linearops
		xtalutil
		${OPENGL_LIBRARIES}
		${GLUT_LIBRARY}
		)
ENDIF()

FETCHCONTENT_DECLARE(
  Catch2
  GIT_REPOSITORY https://github.com/catchorg/Catch2.git
  GIT_TAG        v3.4.0 # or a later release
)
FETCHCONTENT_MAKEAVAILABLE(Catch2)

ADD_EXECUTABLE(itertoolstests
	${PROJECT_SOURCE_DIR}/tests/itertools_test_enumerate.c++
	${PROJECT_SOURCE_DIR}/tests/itertools_test_zip.c++
	#	${PROJECT_SOURCE_DIR}/tests/itertools_test_filter.c++
  ${PROJECT_SOURCE_DIR}/tests/itertools_test_range.c++)
TARGET_LINK_LIBRARIES(itertoolstests PRIVATE itertools Catch2::Catch2WithMain)

ENABLE_TESTING()
LIST(APPEND CMAKE_MODULE_PATH ${catch2_SOURCE_DIR}/extras)

INCLUDE(CTest)
INCLUDE(Catch)
CATCH_DISCOVER_TESTS(itertoolstests)
