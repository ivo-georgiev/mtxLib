cmake_minimum_required(VERSION 3.16)

# Проверка за ESP32 среда
if(IDF_TARGET)
    idf_component_register(
        SRCS "mtxLib.c"
        INCLUDE_DIRS "."
    )
else()
    # Общ случай за Linux/Ubuntu/Windows
    project(mtxLib C)
    add_library(mtxLib STATIC mtxLib.c)
    target_include_directories(mtxLib PUBLIC .)
endif()
