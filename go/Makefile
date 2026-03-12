# Makefile for scMicrobe Go components

.PHONY: all build clean test install

BINARY_NAME=scmicro-go
BUILD_DIR=build

# Detect OS
ifeq ($(OS),Windows_NT)
	BINARY_NAME=scmicro-go.exe
	RM=del /Q
	RMDIR=rmdir /S /Q
else
	BINARY_NAME=scmicro-go
	RM=rm -f
	RMDIR=rm -rf
endif

all: build

build:
	go build -o $(BINARY_NAME) .

build-linux:
	GOOS=linux GOARCH=amd64 go build -o $(BINARY_NAME)-linux-amd64 .

build-darwin:
	GOOS=darwin GOARCH=amd64 go build -o $(BINARY_NAME)-darwin-amd64 .
	GOOS=darwin GOARCH=arm64 go build -o $(BINARY_NAME)-darwin-arm64 .

build-windows:
	GOOS=windows GOARCH=amd64 go build -o $(BINARY_NAME)-windows-amd64.exe .

clean:
ifeq ($(OS),Windows_NT)
	-$(RM) $(BINARY_NAME) 2>nul
	-$(RM) *.exe 2>nul
else
	-$(RM) $(BINARY_NAME)
	-$(RM) *.exe
endif

test:
	go test -v ./...

deps:
	go mod download
	go mod tidy

install: build
ifeq ($(OS),Windows_NT)
	copy $(BINARY_NAME) $(GOPATH)\bin\
else
	cp $(BINARY_NAME) $(GOPATH)/bin/
endif

.DEFAULT_GOAL := build
