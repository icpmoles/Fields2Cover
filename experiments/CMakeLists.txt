include_directories(${CMAKE_BINARY_DIR})


add_executable(poli poli.cpp)
target_compile_definitions(
        poli
        PRIVATE
        DATA_PATH="${PROJECT_SOURCE_DIR}/data/"
)
target_link_libraries(poli PUBLIC Fields2Cover)


add_executable(poli_gml poli_gml.cpp)
target_compile_definitions(
        poli_gml
        PRIVATE
        DATA_PATH="${PROJECT_SOURCE_DIR}/data/"
)
target_link_libraries(poli_gml PUBLIC Fields2Cover)
