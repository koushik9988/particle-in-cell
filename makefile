# Name of the executable
TARGET = pic

# Location of header files relative of makefile location.
INCLUDE_DIR = include

# Include directories
INCLUDES = -I/usr/include/hdf5/serial/

# Library directories
LIB_DIRS = -L/usr/lib/x86_64-linux-gnu/hdf5/serial/

# Libraries to link
LIBS = -lhdf5 -lhdf5_cpp
FFLAGS = -lfftw3 -lm

# Location of source files relative to makefile location.
SRC_DIR = src

OBJS = $(SRC_DIR)/main.o $(SRC_DIR)/init.o $(SRC_DIR)/species.o $(SRC_DIR)/field.o $(SRC_DIR)/domain.o  $(SRC_DIR)/iniparser.o $(SRC_DIR)/output.o  # Add more as needed

# Compiler and C++ standard and optimization flags
CC = g++
#CC = clang
CFLAGS = -Wall -O2 -std=c++17 -I$(INCLUDE_DIR)

# Default target
all: $(TARGET)

# Rule to build the executable
$(TARGET): $(OBJS)
	@echo "linking object file to make $(TARGET)"
	$(CC) $(CFLAGS) $(LIB_DIRS) -o $(TARGET) $(OBJS) $(LIBS) $(FFLAGS)
#rm -f $(OBJS)  

# Rule to build object files from source files
$(SRC_DIR)/%.o: $(SRC_DIR)/%.cpp
	$(CC) $(CFLAGS) $(INCLUDES) -c -o $@ $<  

# Rule to clean up object files and executable
clean:
	@echo "Cleaning executables and object files"
	rm -f $(OBJS) $(TARGET)

