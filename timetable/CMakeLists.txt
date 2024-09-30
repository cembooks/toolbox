#******************************************************************************
# Copyright (C) Siarhei Uzunbajakau, 2023.
#
# This program is free software. You can use, modify, and redistribute it under
# the terms of the GNU Lesser General Public License as published by the Free
# Software Foundation, either version 3 or (at your option) any later version.
# This program is distributed without any warranty.
#
# Refer to COPYING.LESSER for more details.
#*****************************************************************************/

set(PRJ_NAME "timetable")

cmake_minimum_required(VERSION 2.8.12)

include_directories(${DEAL_II_DIR} /usr/include/suitesparse ./include
../shared/include)

set(SOURCE_FILES "src/main.cpp")

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

find_package(deal.II 9.6.0 QUIET HINTS ~/dealii-9.6.0-timetable)

if(NOT ${deal.II_FOUND})
message(FATAL_ERROR "\n"
	"*** Could not locate a (sufficiently recent) version of deal.II. ***\n\n"
	"You may want to either pass a flag -DDEAL_II_DIR=/path/to/deal.II to cmake\n"
	"or set an environment variable \"DEAL_II_DIR\" that contains this path."
)
endif()
DEAL_II_INITIALIZE_CACHED_VARIABLES()

project(${PRJ_NAME})

set(TARGET "${PRJ_NAME}")

add_executable(${TARGET} ${SOURCE_FILES})
DEAL_II_SETUP_TARGET(${TARGET})

set_target_properties(${TARGET}
	PROPERTIES RUNTIME_OUTPUT_DIRECTORY
	"${PROJECT_SOURCE_DIR}/bin/$<CONFIG>")

message(STATUS "TARGET=${TARGET}")

