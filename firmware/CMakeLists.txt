# === 한국어 빌드 파일 안내 시작 ===
# - 파일 역할: host-side firmware scaffold와 테스트 실행 파일을 빌드하는 CMake 설정이다.
# - 주요 기능: firmware library/host executable의 source, include path, link 구성을 지정한다.
# === 한국어 빌드 파일 안내 끝 ===

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
