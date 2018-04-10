cmake_minimum_required(VERSION 3.6)
project (mesh_imp_vorocrust)

file(GLOB SRC_FILES src/*.cpp src/*.h src/operator/*.cpp src/operator/*.h src/constraint/*.cpp src/constraint/*.h src/util/*.cpp src/util/*.h)
add_executable( mesh_imp ${SRC_FILES} )