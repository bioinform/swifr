CC=g++
CPPFLAGS=-std=c++11 -O3 -w
CPPFLAGS2=-std=c++0x
BIN=bin
SRC=src/main.cpp

TESTSRC=tests/*.cpp
# include directories -- where to search for header files
# include for tclap, io-lib, etc...

INC=-I include/ -I /home/dannebar/software_downloads/io_lib-1.14.6/ -I ./external/io_lib_wrapper/ -I ./external/tclap-1.2.1/include/ 
# lib directories -- whereto search for compiled static/dynamic libraries
LIB=-L /home/dannebar/software_downloads/io_lib-1.14.6/lib
IOLIBSO=$(BIN)/
# io_tools for operating on the bams/sams
# LIBS=lstaden-read
EXE=swifr

all: $(EXE) test

$(EXE):
	if [ ! -e $(BIN) ]; then mkdir $(BIN); fi
	$(CC) $(CPPFLAGS) $(INC) $(LIB) -o $(BIN)/$@ $(SRC) -lstaden-read -lpthread -Wl,-rpath,"/home/dannebar/software_downloads/io_lib-1.14.6/lib/"


test:
	 $(CC) $(CPPFLAGS) $(INC) $(LIB) -o $(BIN)/$@ $(TESTSRC) -lstaden-read -Wl,-rpath,"/home/dannebar/software_downloads/io_lib-1.14.6/lib/" 

clean:
	rm $(BIN)/$(EXE)
	# remove generated test files

