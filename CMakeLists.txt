cmake_minimum_required(VERSION 3.10)
project(mwe_particle_emitter)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

set(CMAKE_CXX_STANDARD 23)

file(GLOB_RECURSE SOURCES "src/*.cpp")
add_executable(${PROJECT_NAME} ${SOURCES})

find_package(glfw3)
find_package(glad)
find_package(stb)
find_package(glm)
find_package(spdlog)
target_link_libraries(${PROJECT_NAME} glfw glad::glad stb::stb glm::glm spdlog::spdlog)
