#include <iostream>
#include <memory>

#include "io_lib_wrapper/bam_reader.hpp"
#include "io_lib_wrapper/bam_writer.hpp"
#include "io_lib_wrapper/optional_tags.hpp"

int main() {
	clock_t start = std::clock();
	string f_input = "A01_1.ccs.bam";
	string fname = "a01_1.primary.sam";

	BamReader bamReader(f_input);
	BamWriter bamWriter(fname);
	cerr << "Writing a  " << fname << " file" << endl;
	if (!bamWriter.is_open()) {
		exit(1);
	}
	// fill out header to be the same as in test.bam
	bamWriter.write_header(bamReader);
	const uint32_t R = 100;
	// read from one bam and write to another, but add a new tag
	shared_ptr<Alignment> alignment = bamReader.getNextAlignment();
	int i = 0;
	while ( alignment != nullptr ) {
		// only write the even alingment
		if (i % 2 == 0) {
			MutableAlignment new_alignment(*alignment);
			new_alignment.set_read_id("read" + to_string(i));
			// array<int, 2> data;
			// data[0] = 3;
			// data[1] = 1;
			// new_alignment.add_tag("XX", data);
			string name = "AC";
			string category = "category1";
			new_alignment.add_tag(name, category);
			bamWriter.writeAlignment(new_alignment);
		}
		i++;
		alignment = bamReader.getNextAlignment();
	}
	

	// for (int i = 0; i < 10000; i++) {
	// 	// create an alignment
	// 	MutableAlignment sketch("read" + to_string(i), i, i + R);
	// 	sketch.set_sequence("ACGT");
	// 	sketch.set_flags(0x4 /* unmapped */);
	// 	bamWriter.writeAlignment(sketch);
	// 	// bamWriter << sketch << endl;
	// }
	bamWriter.flush();
	cerr << "Time to write 10000 alignments: " << 
		(std::clock() - start) / (double)(CLOCKS_PER_SEC / 1000) << "ms" << endl;
}