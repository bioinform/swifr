#include "catch.hpp"

#include <vector>
#include <unordered_set>

#include "bam_reader.hpp"
#include "optional_tags.hpp"

TEST_CASE( "Testing non-existent bam", "[bam_reader]" ) {
	BamReader reader("doesnotexist.bam");
	REQUIRE( !reader.is_open() );
}

TEST_CASE( "Testing empty bam", "[bam_reader]" ) {
	BamReader reader("bam/empty.bam");
	REQUIRE( !reader.is_open() );
}

// header is optional
TEST_CASE( "Testing no header sam", "[bam_reader]" ) {
	BamReader reader("bam/no_header.sam");
	REQUIRE( reader.is_open() );
	for (int i = 0; i < 10; i++) {
		auto alignment = reader.getNextAlignment();
		REQUIRE(alignment != nullptr);
		REQUIRE( reader.is_open() );
	}
}

TEST_CASE( "Testing no alignments bam", "[bam_reader]" ) {
	BamReader reader("bam/no_alignments.sam");
	REQUIRE( reader.is_open() );
	auto alignment = reader.getNextAlignment();
	REQUIRE( alignment == nullptr);
	REQUIRE( reader.getLineCount() == 0);
}

TEST_CASE( "Testing 1 alignment bam", "[bam_reader]" ) {
	BamReader reader("bam/one_alignment.sam");
	REQUIRE( reader.is_open() );
	auto alignment = reader.getNextAlignment();
	REQUIRE( alignment != nullptr);
	REQUIRE( reader.getLineCount() == 1);
	alignment = reader.getNextAlignment();
	REQUIRE( alignment == nullptr);
	REQUIRE( reader.getLineCount() == 1);
}

TEST_CASE( "Testing single unaligned CCS bam", "[bam_reader]" ) {
	BamReader reader("bam/unaligned_ccs_read.bam");
	REQUIRE( reader.is_open() );
	auto alignment = reader.getNextAlignment();
	REQUIRE( alignment != nullptr);
	REQUIRE( reader.getLineCount() == 1);
	// validate fields in the alignment
	REQUIRE(alignment->get_cigar() == "*");

	shared_ptr<vector<OptionalTag*>> tags = alignment->get_optional_tags();
	unordered_set<string> expected_tags = {"RG", "za", "np", "rq", "rs", "sn", "zm", "zs"};
	for (OptionalTag* tag : *tags) {
		REQUIRE(expected_tags.find(tag->get_name()) != expected_tags.end() );
		// validate that the values were read correctly
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