include_directories(${CMAKE_BINARY_DIR})


add_executable(experiment main.cpp)
target_compile_definitions(
        experiment
        PRIVATE
        DATA_PATH="${PROJECT_SOURCE_DIR}/data/"
)
target_link_libraries(experiment PUBLIC Fields2Cover)

