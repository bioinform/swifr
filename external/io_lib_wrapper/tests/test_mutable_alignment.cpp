// test alignment class
#include "catch.hpp"

#include "mutable_alignment.hpp"

TEST_CASE( "Testing de novo mutable_alignment constructor", "[mutable_alignment]" ) {
	string name = "read_one", seq = "AACCGGTT", quals = "MMMMMMMM";
	MutableAlignment alignment( name/* read_name */,
		seq /* sequence */,
		quals /* qualities */);
	REQUIRE(alignment.get_read_id() == name);
	REQUIRE(alignment.get_sequence() == seq);
	REQUIRE(alignment.get_qualities() == quals);
	REQUIRE(alignment.get_cigar() == "*");
	REQUIRE( (alignment.get_flags() & 4) > 0 );
}

TEST_CASE( "Testing mutable_alignment constructor", "[mutable_alignment]" ) {
	vector<EditPair> mutations;
	// mutations.emplace_back(EditPair::get_mismatch_code('A'), 105);
	// mutations.emplace_back(EditPair::get_mismatch_code('C'), 115);
	auto start_coordinate = 100;
	auto length = 30;
	MutableAlignment as(mutations, start_coordinate, length);
	CHECK(as.get_alignment_start() == start_coordinate);
	CHECK(as.length() == length);
}

TEST_CASE( "Testing mutable_alignment constructor -- from Alignment obj", "[mutable_alignment]" ) {
	vector<EditPair> mutations;
	// mutations.emplace_back(EditPair::get_mismatch_code('A'), 105);
	// mutations.emplace_back(EditPair::get_mismatch_code('C'), 115);
	auto start_coordinate = 100;
	auto length = 30;
	MutableAlignment as(mutations, start_coordinate, length);
	CHECK(as.get_alignment_start() == start_coordinate);
	CHECK(as.length() == length);
	// TODO
}

TEST_CASE( "Testing mutable_alignment constructor -- from AlignmentCore obj", "[mutable_alignment]" ) {
	auto start_coordinate = 100;
	auto length = 30;
	vector<EditPair> mutations;
	// mutations.emplace_back(EditPair::get_mismatch_code('A'), start_coordinate + 5);
	MutableAlignment as(mutations, start_coordinate, length);
	CHECK(as.get_alignment_start() == start_coordinate);
	CHECK(as.length() == length);

	
	MutableAlignment copy(as);

	REQUIRE(copy.get_alignment_start() == start_coordinate);
	REQUIRE(copy.length() == length);
	
	REQUIRE(copy.getMutations().size() == 0);
	REQUIRE(copy.get_read_id() == as.get_read_id() );
	REQUIRE(copy.get_flags() == as.get_flags() );
	REQUIRE(copy.get_reference_name() == as.get_reference_name() );
	REQUIRE(copy.get_mapq() == as.get_mapq() );
	
	REQUIRE(copy.get_cigar() == as.get_cigar() );
	REQUIRE(copy.get_rnext() == as.get_rnext() );
	REQUIRE(copy.get_pnext() == as.get_pnext() );
	REQUIRE(copy.get_tlen() == as.get_tlen() );
	REQUIRE(copy.get_sequence() == as.get_sequence() );
	REQUIRE(copy.get_qualities() == as.get_qualities() );
	REQUIRE(copy.get_optional_tags_string() == as.get_optional_tags_string() );
}

TEST_CASE( "Testing mutable_alignment has_mutation_at", "[mutable_alignment]" ) {
	vector<EditPair> mutations;
	mutations.emplace_back(EditPair::get_mismatch_code('A'), 105);
	mutations.emplace_back(EditPair::get_mismatch_code('C'), 115);
	auto start_coordinate = 100;
	auto length = 30;
	MutableAlignment as(mutations, start_coordinate, length);
	CHECK(as.has_mutation_at(105) );
	CHECK(as.has_mutation_at(115) );
}

TEST_CASE( "Testing mutable_alignment getMutations", "[mutable_alignment]" ) {
	vector<EditPair> mutations;
	mutations.emplace_back(EditPair::get_mismatch_code('A'), 105);
	mutations.emplace_back(EditPair::get_mismatch_code('C'), 115);
	auto start_coordinate = 100;
	auto length = 30;
	MutableAlignment as(mutations, start_coordinate, length);
	auto as_mutations = as.getMutations();
	CHECK(as_mutations.size() == mutations.size() );
}

//==============================================================================
// CIGAR::DELETIONS
//==============================================================================

TEST_CASE( "Testing cigar -- no indels", "[mutable_alignment]" ) {
	vector<EditPair> mutations;
	MutableAlignment alignment(mutations, 6, 100);
	string cigar = alignment.get_cigar();
	REQUIRE(cigar == "100M");
}

TEST_CASE( "Testing cigar -- single deletion surrounded by matches", "[mutable_alignment]" ) {
	vector<EditPair> mutations;
	mutations.emplace_back(EditPair::get_deletion_code(), 10);
	MutableAlignment alignment(mutations, 6, 100);
	string cigar = alignment.get_cigar();
	REQUIRE(cigar == "4M1D96M");
}

TEST_CASE( "Testing cigar -- single deletion at the end", "[mutable_alignment]" ) {
	vector<EditPair> mutations;
	mutations.emplace_back(EditPair::get_deletion_code(), 6 + 100 - 1);
	MutableAlignment alignment(mutations, 6, 100);
	string cigar = alignment.get_cigar();
	REQUIRE(cigar == "99M1D1M");
}

TEST_CASE( "Testing cigar -- double deletion at the end", "[mutable_alignment]" ) {
	vector<EditPair> mutations;
	mutations.emplace_back(EditPair::get_deletion_code(), 6 + 100 - 2);
	mutations.emplace_back(EditPair::get_deletion_code(), 6 + 100 - 1);
	MutableAlignment alignment(mutations, 6, 100);
	string cigar = alignment.get_cigar();
	REQUIRE(cigar == "98M2D2M");
}

//==============================================================================
// CIGAR::INSERTIONS
//==============================================================================

TEST_CASE( "Testing cigar -- single insertion surrounded by matches", "[mutable_alignment]" ) {
	vector<EditPair> mutations;
	mutations.emplace_back(EditPair::get_insertion_code('A'), 50);
	MutableAlignment alignment(mutations, 10, 100);
	string cigar = alignment.get_cigar();
	REQUIRE(cigar == "40M1I60M");
}

TEST_CASE( "Testing cigar -- single insertion at the end", "[mutable_alignment]" ) {
	vector<EditPair> mutations;
	mutations.emplace_back(EditPair::get_insertion_code('A'), 110);
	MutableAlignment alignment(mutations, 10, 100);
	string cigar = alignment.get_cigar();
	REQUIRE(cigar == "100M1I");
}

TEST_CASE( "Testing cigar -- double insertion in the middle", "[mutable_alignment]" ) {
	vector<EditPair> mutations;
	mutations.emplace_back(EditPair::get_insertion_code('A'), 50);
	mutations.emplace_back(EditPair::get_insertion_code('C'), 50);
	MutableAlignment alignment(mutations, 10, 100);
	string cigar = alignment.get_cigar();
	REQUIRE(cigar == "40M2I60M");
}
