// test alignment class
#include "catch.hpp"

#include "edit_pair.hpp"

TEST_CASE( "Testing edit_pair constructor", "[edit_pair]" ) {
	EditPair ep(0, 0);
	CHECK(ep.is_mismatch());
	CHECK(!ep.is_insertion());
	CHECK(!ep.is_deletion());
}

TEST_CASE( "Testing edit_pair getters", "[edit_pair]" ) {
	auto code = 0;
	int position = 0;
	EditPair ep(code, position);
	CHECK(ep.get_code() == code);
	CHECK(ep.get_position() == position);
	code = 5;
	position = 10;
	EditPair ep2(code, position);
	CHECK(ep2.get_code() == code);
	CHECK(ep2.get_position() == position);
}

TEST_CASE( "Testing edit_pair get insertion code function", "[edit_pair]" ) {
	auto code = EditPair::get_insertion_code('A');
	CHECK(code == INSERTION_A);
	code = EditPair::get_insertion_code('C');
	CHECK(code == INSERTION_C);
	code = EditPair::get_insertion_code('G');
	CHECK(code == INSERTION_G);
	code = EditPair::get_insertion_code('T');
	CHECK(code == INSERTION_T);
	code = EditPair::get_insertion_code('N');
	CHECK(code == INSERTION_N);
	code = EditPair::get_insertion_code('X');
	CHECK(code == UNKNOWN_BASE);
}

TEST_CASE( "Testing edit_pair get mismatch code function", "[edit_pair]" ) {
	auto code = EditPair::get_mismatch_code('A');
	CHECK(code == MISMATCH_A);
	code = EditPair::get_mismatch_code('C');
	CHECK(code == MISMATCH_C);
	code = EditPair::get_mismatch_code('G');
	CHECK(code == MISMATCH_G);
	code = EditPair::get_mismatch_code('T');
	CHECK(code == MISMATCH_T);
	code = EditPair::get_mismatch_code('N');
	CHECK(code == MISMATCH_N);
	code = EditPair::get_mismatch_code('X');
	CHECK(code == UNKNOWN_BASE);
}

TEST_CASE( "Testing edit_pair get deletion code function", "[edit_pair]" ) {
	auto code = EditPair::get_deletion_code();
	CHECK(code == DELETION_N);
}

// TEST_CASE( "Testing alignment w/ no edits", "[alignment]" ) {
// 	BamReader reader("bam/no_edits_alignment.bam");
// 	REQUIRE( !reader.is_open() );
// }

// TEST_CASE( "Testing alignment w/ only mismatches in MD tag", "[alignment]" ) {
// 	BamReader reader("bam/only_md_tag_edits.bam");
// 	REQUIRE( !reader.is_open() );
// }

// TEST_CASE( "Testing alignment w/ only indels in cigar string", "[alignment]" ) {
// 	BamReader reader("bam/only_cigar_string_edits.bam");
// 	REQUIRE( !reader.is_open() );
// }

// TEST_CASE( "Testing alignment w/ indels in cigar and mismatches in MD tag", "[alignment]" ) {
// 	BamReader reader("bam/cigar_and_md_edits.bam");
// 	REQUIRE( !reader.is_open() );
// }