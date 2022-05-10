cmake_minimum_required(VERSION 3.10)

PROJECT(final_project)
SET(VTK_DIR /Users/luyingcai/Desktop/W410/build)

find_package(VTK REQUIRED)

add_executable(final_project MACOSX_BUNDLE final_project)

SET(CMAKE_CXX_LINK_FLAGS "-framework OpenGL")
if(VTK_LIBRARIES)
  target_link_libraries(final_project ${VTK_LIBRARIES})
else()
  target_link_libraries(final_project vtkHybrid)
endif()

