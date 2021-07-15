#include <iostream>
#include <memory>
#include <ctime>
#include "io_lib_wrapper/bam_reader.hpp"
#include "io_lib_wrapper/alignment.hpp"

using namespace std;

int main() {
	clock_t start = std::clock();
	string fname = "examples/test.bam";
	BamReader bamReader(fname);
	if (!bamReader.is_open()) {
		exit(1);
	}

	int primary_alignments = 0;
	shared_ptr<Alignment> alignment = bamReader.getNextAlignment();
	while ( alignment != nullptr ) {
		if (!alignment->is_unaligned() && alignment->is_primary() ) {
			primary_alignments++;
		}
		alignment = bamReader.getNextAlignment();
	}
	cerr << "Read " << bamReader.getLineCount() << " alignments," 
		" of them " << primary_alignments << " are unique best" << endl;
	cerr << "Time to read: " << 
		(std::clock() - start) / (double)(CLOCKS_PER_SEC / 1000) << "ms" << endl;
}