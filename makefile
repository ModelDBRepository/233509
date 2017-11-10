# created using the help provided at the below URL
# https://stackoverflow.com/questions/2908057/can-i-compile-all-cpp-files-in-src-to-os-in-obj-then-link-to-binary-in


all: Net.exe

clean:
	rm Net.exe *.o
	
#Net.exe: Net.o IN.o RE.o TC.o
#	g++ -o Net Net.o IN.o Re.o TC.o Net.cpp
	
#Tell make to make one .out file for each .cpp file found in the current directory
all: $(patsubst %.cpp, %.out, $(wildcard *.cpp))

#Rule how to create arbitary .out files. 
#First state what is needed for them e.g. additional headers, .cpp files in an include folder...
#Then the command to create the .out file, probably you want to add further options to the g++ call.
%.out: %.cpp Makefile
	g++ $< -o $@ -std=c++0x

SRC_DIR := ./src
OBJ_DIR := ./obj
SRC_FILES := $(wildcard $(SRC_DIR)/*.cpp)
OBJ_FILES := $(patsubst $(SRC_DIR)/%.cpp,$(OBJ_DIR)/%.o,$(SRC_FILES))

TOP_SRC_DIR := ./
TOP_OBJ_DIR := ./
TOP_SRC_FILES := $(wildcard $(TOP_SRC_DIR)/*.cpp)
TOP_OBJ_FILES := $(patsubst $(TOP_SRC_DIR)/%.cpp,$(TOP_OBJ_DIR)/%.o,$(TOP_SRC_FILES))
 
#LDFLAGS := ...
#CPPFLAGS := ...
#CXXFLAGS := ...

Net.exe: $(OBJ_FILES) $(TOP_OBJ_FILES) 
	g++ $(LDFLAGS) -o $@ $^

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp
	g++ $(CPPFLAGS) $(CXXFLAGS) -c -o $@ $<
$(TOP_OBJ_DIR)/%.o: $(TOP_SRC_DIR)/%.cpp
	g++ $(CPPFLAGS) $(CXXFLAGS) -c -o $@ $<

	