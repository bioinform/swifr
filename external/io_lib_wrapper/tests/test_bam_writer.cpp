#include "catch.hpp"

#include <unordered_set>
#include <cmath>

#include "bam_reader.hpp"
#include "bam_writer.hpp"
#include "optional_tags.hpp"

TEST_CASE( "Reading a header and a single alignment and writing to a SAM", "[bam_writer]" ) {
	string f_input = "bam/unaligned_ccs_read.bam";
	string fname = "bam/unaligned_ccs_read_copy.sam";

	BamReader bamReader(f_input);
	REQUIRE(bamReader.is_open());
	BamWriter bamWriter(fname);
	REQUIRE(bamWriter.is_open());

	// fill out header to be the same as in test.bam
	bamWriter.write_header(bamReader);

	shared_ptr<AlignmentCore> alignment = bamReader.getNextAlignment();
	REQUIRE(alignment != nullptr);
	bamWriter.writeAlignment(alignment);
}

TEST_CASE( "Reading, writing, and reading to validate that written info matched", "[bam_writer]" ) {
	string f_input = "bam/unaligned_ccs_read.bam";
	string fname = "bam/unaligned_ccs_read_copy2.sam";

	BamReader bamReader(f_input);
	REQUIRE(bamReader.is_open());

	BamWriter bamWriter(fname);
	REQUIRE(bamWriter.is_open());


	// fill out header to be the same as in test.bam
	bamWriter.write_header(bamReader);

	shared_ptr<AlignmentCore> alignment = bamReader.getNextAlignment();
	REQUIRE(alignment != nullptr);
	bamWriter.writeAlignment(alignment);
	bamWriter.flush();
	bamWriter.close();
	bamReader.close();

	// now read from the SAM that was just created
	BamReader bamReader2(fname);
	REQUIRE(bamReader2.is_open());
	alignment = bamReader2.getNextAlignment();
	REQUIRE(alignment != nullptr);
	// compare all fields
	REQUIRE(alignment->get_read_id() == "m54013_160609_014924/31392526/ccs");
	REQUIRE(alignment->get_flags() == 4);
	REQUIRE(alignment->get_reference_name() == "*");
	REQUIRE(alignment->get_alignment_start() == 0);
	REQUIRE(alignment->get_mapq() == 255);
	REQUIRE(alignment->get_cigar() == "*");
	REQUIRE(alignment->get_rnext() 	== "*");
	REQUIRE(alignment->get_pnext() 	== 0);
	REQUIRE(alignment->get_tlen() 	== 0);
	REQUIRE(alignment->get_sequence().size() == 2580);
	REQUIRE(alignment->get_qualities().size() == 2580);
	// check optional flags -- all present and values are equal
	shared_ptr<vector<OptionalTag*>> tags = alignment->get_optional_tags();
	REQUIRE(tags != nullptr);
	REQUIRE(tags->size() == 8);
	// now validate that tag names and values are as expected
	unordered_set<string> tag_names = {"RG", "np", "rq", "rs", "sn", "za", "zm", "zs"};
	for (auto * tag : *tags) {
		REQUIRE(tag_names.find(tag->get_name() ) != tag_names.end() );
		if (tag->get_name() == "rs") {
			ArrayTag<int> * int_tag = dynamic_cast<ArrayTag<int>*>(tag);
			REQUIRE(int_tag != NULL);
			vector<int> expected_values = {7,0,0,0,0,0};
			vector<int> values = int_tag->get_values();
			REQUIRE(values.size() == 6);
			for (int i = 0; i < expected_values.size(); i++)
				REQUIRE(values[i] == expected_values[i]);
		}
		else if (tag->get_name() == "sn" ) {
			vector<float> expected_values = {5.57006,10.0712,7.08381,10.9168};
			ArrayTag<float> * float_tag = dynamic_cast<ArrayTag<float>*>(tag);
			REQUIRE(float_tag != NULL);
			vector<float> values = float_tag->get_values();
			REQUIRE(values.size() == 4);
			for (int i = 0; i < expected_values.size(); i++)
				REQUIRE(values[i] == expected_values[i]);
		}
		else if (tag->get_name() == "zs" ) {
			vector<float> expected_values = {NAN,NAN,NAN,NAN,NAN,NAN,NAN};
			ArrayTag<float> * float_tag = dynamic_cast<ArrayTag<float>*>(tag);
			REQUIRE(float_tag != NULL);
			vector<float> values = float_tag->get_values();
			REQUIRE(values.size() == 7);
			for (int i = 0; i < expected_values.size(); i++)
				REQUIRE(std::isnan(values[i]) );
		}
	}
}

TEST_CASE( "Writing de novo BAM w/ default header", "[bam_writer]" ) {
	BamWriter writer("bam/de_novo.sam");
	REQUIRE(writer.is_open());
	writer.write_header();
	shared_ptr<AlignmentCore> alignment(new MutableAlignment("read_one", "AACCGGTT", "ZZZZZZZZ"));
	writer.writeAlignment(alignment);
	writer.close();

	// read and verify that what we read matches output
	BamReader reader("bam/de_novo.sam");
	REQUIRE(reader.is_open());
	shared_ptr<Alignment> parsed_alignment = reader.getNextAlignment();
	// cerr << *parsed_alignment == *alignment << endl;
	REQUIRE(*parsed_alignment == *alignment);
}