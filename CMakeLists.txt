cmake_minimum_required(VERSION 3.16)
project(ALNview LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTORCC ON)

find_package(Qt6 6.9 REQUIRED COMPONENTS Widgets)
find_package(ZLIB REQUIRED)
find_package(Threads REQUIRED)

add_executable(ALNview
    main.cpp main_window.cpp open_window.cpp
    sticks.c doter.c alncode.c align.c gene_core.c ONElib.c GDB.c hash.c select.c
    main_window.h open_window.h viewer.qrc
)
target_link_libraries(ALNview PRIVATE Qt6::Widgets ZLIB::ZLIB Threads::Threads)
if(UNIX)
    target_link_libraries(ALNview PRIVATE m)
endif()
if(CMAKE_C_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(ALNview PRIVATE -Wall -Wno-unused-result)
endif()
install(TARGETS ALNview RUNTIME DESTINATION bin)
