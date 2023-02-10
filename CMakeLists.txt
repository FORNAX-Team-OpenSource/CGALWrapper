cmake_minimum_required(VERSION 3.1...3.20)
project(CGALGeometryFunction)

add_definitions(-DCGAL_MESH_3_NO_DEPRECATED_SURFACE_INDEX
                -DCGAL_MESH_3_NO_DEPRECATED_C3T3_ITERATORS)

if(CGAL_MESH_3_VERBOSE)
  add_definitions(-DCGAL_MESH_3_VERBOSE)
endif()

find_package(CGAL REQUIRED COMPONENTS ImageIO)
find_package(Boost)

option(CGAL_ACTIVATE_CONCURRENT_MESH_3 "Activate parallelism in Mesh_3" OFF)

if(CGAL_ACTIVATE_CONCURRENT_MESH_3 OR ENV{CGAL_ACTIVATE_CONCURRENT_MESH_3})
  add_definitions(-DCGAL_CONCURRENT_MESH_3)
  find_package(TBB REQUIRED)
  include(CGAL_TBB_support)
endif()

# Use Eigen
find_package(Eigen3 3.1.0 QUIET) #(3.1.0 or greater)
include(CGAL_Eigen3_support)
if(NOT TARGET CGAL::Eigen3_support)
  message( STATUS "NOTICE: All examples need the Eigen3 library, and will not be compiled." )
  return()
endif() #CGAL::Eigen_3_support

create_single_source_cgal_program(
  "mesh_polyhedral_domain_with_surface_inside.cpp")
target_link_libraries(mesh_polyhedral_domain_with_surface_inside
                      PUBLIC CGAL::Eigen3_support)