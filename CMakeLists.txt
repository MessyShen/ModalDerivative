cmake_minimum_required(VERSION 3.16)
project(example)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

list(PREPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

# Libigl
include(libigl)

# Enable the target igl::glfw
igl_include(glfw)
# Other modules you could enable
#igl_include(embree)
igl_include(imgui)
#igl_include(opengl)
#igl_include(stb)
#igl_include(predicates)
#igl_include(xml)
#igl_include(copyleft cgal)
#igl_include(copyleft comiso)
#igl_include(copyleft core)
#igl_include(copyleft cork)
#igl_include(copyleft tetgen)
#igl_include(restricted matlab)
#igl_include(restricted mosek)
#igl_include(restricted triangle)

file(GLOB SRCFILES 
    "src/*.cpp"
    "src/ModalDerivative/*.cpp"
    "src/LoboDynamic/*.cpp"
    "src/Utils/*.cpp"
    "src/Functions/*.cpp"
    "src/LoboDynamic/LoboVolumtricMesh/*.cpp"
    "src/LoboDynamic/materials/*.cpp"
    "src/Utils/SparseMatrix/*.cpp"
    "src/Utils/BasicIO/*.cpp"
)

FOREACH(item ${SRCFILES})
    IF(${item} MATCHES "main.cpp")
        LIST(REMOVE_ITEM SRCFILES ${item})
    ENDIF(${item} MATCHES "main.cpp")
ENDFOREACH(item)

add_library(${PROJECT_NAME}_dev ${SRCFILES})

target_include_directories(${PROJECT_NAME}_dev PUBLIC
    "include"
    "src"
    "src/ModalDerivative"
    "src/Functions"
    "src/Utils"
    "src/Utils/SparseMatrix"
    "src/Utils/BasicIO"
    "src/LoboDynamic/LoboVolumtricMesh"
    "src/LoboDynamic/materials"
)


# Link igl (and the glfw module) to your project
target_link_libraries(${PROJECT_NAME}_dev PUBLIC 
  igl::glfw
  ## Other modules you could link to
  # igl::embree
  igl::imgui
  # igl::opengl
  # igl::stb
  # igl::predicates
  # igl::xml
  # igl_copyleft::cgal
  # igl_copyleft::comiso
  # igl_copyleft::core
  # igl_copyleft::cork
  # igl_copyleft::tetgen
  # igl_restricted::matlab
  # igl_restricted::mosek
  # igl_restricted::triangle
)

# Spectra
set(SPECTRA_ROOT ${CMAKE_CURRENT_SOURCE_DIR}/external/spectra-0.8.1/include)
target_include_directories(${PROJECT_NAME}_dev PUBLIC ${SPECTRA_ROOT})



add_executable(${PROJECT_NAME}_bin "src/main.cpp")

# Link against the library
target_link_libraries(${PROJECT_NAME}_bin PUBLIC ${PROJECT_NAME}_dev)