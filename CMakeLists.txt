cmake_minimum_required(VERSION 3.16)

# ESP32 enviroment
if(IDF_TARGET)
    idf_component_register(
        SRCS "mtxLib.c"
        INCLUDE_DIRS "."
    )
else()
    # Generic build Linux/Ubuntu/Windows
    project(mtxLib C)
    add_library(mtxLib STATIC mtxLib.c)
    target_include_directories(mtxLib PUBLIC .)
endif()
