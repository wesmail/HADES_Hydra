cmake_minimum_required(VERSION 3.10)

if(NOT WIN32)
  string(ASCII 27 Esc)
  set(ColourReset "${Esc}[m")
  set(ColourBold  "${Esc}[1m")
  set(Red         "${Esc}[31m")
  set(Green       "${Esc}[32m")
  set(Yellow      "${Esc}[33m")
  set(Blue        "${Esc}[34m")
  set(Magenta     "${Esc}[35m")
  set(Cyan        "${Esc}[36m")
  set(White       "${Esc}[37m")
  set(BoldRed     "${Esc}[1;31m")
  set(BoldGreen   "${Esc}[1;32m")
  set(BoldYellow  "${Esc}[1;33m")
  set(BoldBlue    "${Esc}[1;34m")
  set(BoldMagenta "${Esc}[1;35m")
  set(BoldCyan    "${Esc}[1;36m")
  set(BoldWhite   "${Esc}[1;37m")
endif()

project(Analysis)

set(TORCHLIB_PATH $ENV{TORCHLIB_PATH})
list(APPEND CMAKE_PREFIX_PATH ${TORCHLIB_PATH} $ENV{ROOTSYS} $ENV{ROOTSYS}/etc/cmake $ENV{ROOTSYS}/cmake/modules/ ${CMAKE_CURRENT_SOURCE_DIR})
list(APPEND CMAKE_MODULE_PATH $ENV{ROOTSYS} $ENV{ROOTSYS}/etc/cmake $ENV{ROOTSYS}/cmake/modules/ ${CMAKE_CURRENT_SOURCE_DIR})

message("--${BoldGreen} ROOT${ColourReset} at " $ENV{ROOTSYS})
message ("--${BoldGreen} LibTorch${ColourReset} at " ${TORCHLIB_PATH})
message ("--${BoldGreen} Hydra${ColourReset} at " $ENV{HADDIR})

find_package(ROOT
    REQUIRED COMPONENTS Core RIO Net Hist Graf Graf3d Gpad Tree Rint Postscript Matrix Physics MathCore Thread
    OPTIONAL_COMPONENTS Cint Cling)
include(RootNewMacros)

find_package(HYDRA REQUIRED)
find_package(Torch REQUIRED)

set(HYDRA_LIBRARIES $ENV{HADDIR}/lib/)
set(HYDRA_INCLUDE_DIR $ENV{HADDIR}/include)

include_directories(AFTER
    ${ROOT_INCLUDE_DIRS}
    ${HYDRA_INCLUDE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}
)

message("--${BoldGreen} Hydra${ColourReset} at " ${HYDRA_LIBRARIES})

ROOT_GENERATE_DICTIONARY(G__HRefitCand HRefitCand.h LINKDEF HRefitCandLinkDef.h)
ROOT_GENERATE_DICTIONARY(G__HFitter HFitter.h LINKDEF HFitterLinkDef.h)

add_library(HRefitCand SHARED HRefitCand.cxx G__HRefitCand.cxx)
add_library(HFitter SHARED HFitter.cxx G__HFitter.cxx)
target_link_libraries(HRefitCand PRIVATE ${ROOT_LIBRARIES} ${TORCH_LIBRARIES})
target_link_libraries(HFitter PRIVATE ${ROOT_LIBRARIES} ${TORCH_LIBRARIES})
set_target_properties(HRefitCand HFitter PROPERTIES CXX_STANDARD 11)
