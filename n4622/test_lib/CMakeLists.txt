INCLUDE_DIRECTORIES(${SPIRAL_BINARY_DIR}/libs)
LINK_DIRECTORIES(${SPIRAL_BINARY_DIR}/libs)

set(CMAKE_MODULE_PATH "${SPIRAL_SOURCE_DIR}/cmake_modules/")

#SET(SLATEC "${SPIRAL_SOURCE_DIR}/libs/libslatec.a")

ADD_EXECUTABLE(smooth_func smooth_func.f)
TARGET_LINK_LIBRARIES(smooth_func SPIRAL_LIB pgplot)
#
#ADD_EXECUTABLE(bar_force bar_force.f)
#TARGET_LINK_LIBRARIES(bar_force SPIRAL_LIB pgplot)
#


FIND_PACKAGE(PGPLOT REQUIRED)
INCLUDE_DIRECTORIES(${PGPLOT_INCLUDES})
LINK_LIBRARIES(${PGPLOT_LIBRARY})


ADD_LIBRARY(slatec STATIC IMPORTED)
SET_PROPERTY(TARGET slatec PROPERTY IMPORTED_LOCATION "${SPIRAL_SOURCE_DIR}/libs/libslatec.a" )
ADD_EXECUTABLE(companion_force companion_force.f)
TARGET_LINK_LIBRARIES(companion_force SPIRAL_LIB slatec pgplot)
SET_PROPERTY(TARGET companion_force PROPERTY LINKER_LANGUAGE Fortran)
SET_PROPERTY(TARGET companion_force PROPERTY Fortran_MODULE_DIRECTORY "${SPIRAL_BINARY_DIR}/libs")


#
#ADD_EXECUTABLE(companion_potential_fourier companion_potential_fourier.f)
#TARGET_LINK_LIBRARIES(companion_potential_fourier SPIRAL_LIB pgplot ${SLATEC})
#
ADD_EXECUTABLE(plot_rotation plot_rotation.f)
TARGET_LINK_LIBRARIES(plot_rotation SPIRAL_LIB pgplot slatec)
#
#ADD_EXECUTABLE(tlaplace_coefficient tlaplace_coefficient.f)
#TARGET_LINK_LIBRARIES(tlaplace_coefficient SPIRAL_LIB pgplot ${SLATEC})


#set_property(TARGET 
#	smooth
#	bar_force
#	companion_force
#	companion_potential_fourier
#	plot_rotation
#	tlaplace_coefficient
#
#PROPERTY LINKER_LANGUAGE Fortran)
