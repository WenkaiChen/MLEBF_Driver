cmake_minimum_required (VERSION 2.8)
project (MLEBF_Driver)
add_definitions(-std=c++11)

AUX_SOURCE_DIRECTORY(. SRC_LIST)
ADD_SUBDIRECTORY(src)

INCLUDE_DIRECTORIES(include)

ADD_EXECUTABLE(MLEBF_Driver ${SRC_LIST})

TARGET_LINK_LIBRARIES( MLEBF_Driver src )
