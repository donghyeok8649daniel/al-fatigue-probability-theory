cmake_minimum_required(VERSION 3.16)
project(fatigue_firmware_core C)

set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED ON)

add_library(fatigue_controller
    src/fatigue_controller.c
)
target_include_directories(fatigue_controller PUBLIC include)
target_link_libraries(fatigue_controller PUBLIC m)

add_executable(firmware_host_test src/host_test.c)
target_link_libraries(firmware_host_test PRIVATE fatigue_controller)
