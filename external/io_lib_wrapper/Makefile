# compile examples
CC=g++
# CPPFLAGS=-std=c++11 -O3 -w
# quick compilation
CPPFLAGS=-std=c++11 -w
SRC=examples
TESTSRC=tests/test_*.cpp
BIN=.
INC=-I . -I /usr/local/include/
TESTINCLUDE=-I io_lib_wrapper -I tests
LIB=-L /usr/local/lib/
# io_tools for operating on the bams/sams
# LIBS=lstaden-read

all: create reading_bam writing_bam reading_fasta reading_reference_fasta test

create:
	if [ ! -e $(BIN) ]; then mkdir $(BIN); fi

reading_bam:
	$(CC) $(CPPFLAGS) $(INC) $(LIB) -o $(BIN)/$@ $(SRC)/reading_bam.cpp -lstaden-read

writing_bam:
	$(CC) $(CPPFLAGS) $(INC) $(LIB) -o $(BIN)/$@ $(SRC)/writing_bam.cpp -lstaden-read

reading_fasta:
	$(CC) $(CPPFLAGS) $(INC) $(LIB) -o $(BIN)/$@ $(SRC)/reading_fasta.cpp -lstaden-read

reading_reference_fasta:
	$(CC) $(CPPFLAGS) $(INC) $(LIB) -o $(BIN)/$@ $(SRC)/reading_reference_fasta.cpp -lstaden-read

# compile tests
test:
	$(CC) $(CPPFLAGS) $(INC) $(TESTINCLUDE) $(LIB) -o $(BIN)/$@ $(TESTSRC) -lstaden-read

clean:
	# remove examples, testing binary
	rm -f $(BIN)/test $(BIN)/reading_bam $(BIN)/writing_bam $(BIN)/reading_fasta $(BIN)/reading_reference_fasta
	rm *.sam
	rm *.bam	
