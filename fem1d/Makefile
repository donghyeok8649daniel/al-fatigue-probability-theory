# === 한국어 파일 안내 시작 ===
# - 파일 역할: 1D bar FEM C 코드를 빌드하고 실행파일 bin/fem1d_solver를 생성한다.
# - 주요 타깃: all, clean, self-test, demo
# - 주의: 확률이론과는 독립적으로 빌드된다.
# === 한국어 파일 안내 끝 ===

CC ?= cc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra -Wpedantic
CPPFLAGS ?= -Iinclude
LDLIBS ?= -lm

BIN_DIR := bin
TARGET := $(BIN_DIR)/fem1d_solver
SOURCES := src/main.c src/fem1d.c

.PHONY: all clean self-test demo

all: $(TARGET)

$(TARGET): $(SOURCES) include/fem1d.h
	mkdir -p $(BIN_DIR)
	$(CC) $(CPPFLAGS) $(CFLAGS) $(SOURCES) -o $(TARGET) $(LDLIBS)

self-test: $(TARGET)
	./$(TARGET) --self-test

demo: $(TARGET)
	./$(TARGET) --outdir ../results/data/fem1d_demo

clean:
	rm -rf $(BIN_DIR)
