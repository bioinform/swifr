#include <unordered_set>

// test alignment class
#include "catch.hpp"

#include "bam_reader.hpp"
#include "alignment.hpp"

TEST_CASE( "Testing for primary alignments", "[alignment]" ) {
	BamReader reader("bam/test.10.sam");
	int primary = 0;
	REQUIRE( reader.is_open() );
	for (int i = 0; i < 10; i++) {
		auto alignment = reader.getNextAlignment();
		primary += alignment->is_primary();
	}
	REQUIRE(primary == 8);
}

TEST_CASE( "Testing for unaligned reads", "[alignment]" ) {
	BamReader reader("bam/unaligned.sam");
	REQUIRE( reader.is_open() );
	auto alignment = reader.getNextAlignment();
	REQUIRE(alignment != nullptr);
	REQUIRE(alignment->is_unaligned());
}


TEST_CASE( "Testing alignment w/ no edits", "[alignment]" ) {
	BamReader reader("bam/no_edits_alignment.bam");
	REQUIRE( reader.is_open() );

	shared_ptr<Alignment> alignment = reader.getNextAlignment();
	REQUIRE(alignment != nullptr);
	REQUIRE(alignment->getCigarMutations().size() == 0);

	shared_ptr<LongSequence> reference(new LongSequence("bam/hxb2.fa") );
	REQUIRE(reference->read_sequence() );
	alignment->set_reference(reference);
	vector<EditPair> mutations = alignment->getMutations();
	REQUIRE(mutations.size() == 0); 
}

TEST_CASE( "Testing alignment w/ only mismatches in MD tag", "[alignment]" ) {
	BamReader reader("bam/only_md_tag_edits.bam");
	REQUIRE( reader.is_open() );
	shared_ptr<Alignment> alignment = reader.getNextAlignment();
	REQUIRE(alignment != nullptr);
	REQUIRE(alignment->getCigarMutations().size() == 0);

	// get MD mutations now
	shared_ptr<LongSequence> reference(new LongSequence("bam/hxb2.fa") );
	REQUIRE(reference->read_sequence() );
	alignment->set_reference(reference);
	vector<EditPair> mutations = alignment->getMutations();
	// for (auto m : mutations) cerr << m << " ";
	REQUIRE(mutations.size() == 3);
	REQUIRE(mutations[0].get_position() == 15);
	REQUIRE(mutations[1].get_position() == 20);
	REQUIRE(mutations[2].get_position() == 27);
}

TEST_CASE( "Testing alignment w/ only indels in cigar string", "[alignment]" ) {
	BamReader reader("bam/only_cigar_string_edits.bam");
	REQUIRE( reader.is_open() );
	shared_ptr<Alignment> alignment = reader.getNextAlignment();
	REQUIRE(alignment != nullptr);
	vector<EditPair> cigar_muts = alignment->getCigarMutations();
	REQUIRE(cigar_muts.size() == 2);

	// for (auto e : cigar_muts) cerr << e << " ";
		// cerr << endl;

	// TODO validate positions of all insertions and deletions
	CHECK(cigar_muts[0].get_position() == 18);
	CHECK(cigar_muts[0].get_code() == DELETION_N);
	CHECK(cigar_muts[1].get_position() == 23);
	CHECK(cigar_muts[1].get_code() == INSERTION_G);
}

TEST_CASE( "Testing alignment w/ indels in cigar and mismatches in MD tag", "[alignment]" ) {
	// cerr << " ===^^^^====" << endl;
	BamReader reader("bam/cigar_and_md_edits.bam");
	REQUIRE( reader.is_open() );
	shared_ptr<Alignment> alignment = reader.getNextAlignment();
	REQUIRE(alignment != nullptr);
	vector<EditPair> cigar_muts = alignment->getMutations();
	REQUIRE(cigar_muts.size() == 2);
	CHECK(cigar_muts[0].get_position() == 16);
	CHECK(cigar_muts[0].get_code() == DELETION_N);
	CHECK(cigar_muts[1].get_position() == 20);
	CHECK(cigar_muts[1].get_code() == INSERTION_G);

	shared_ptr<LongSequence> reference(new LongSequence("bam/hxb2.fa") );
	REQUIRE(reference->read_sequence() );
	alignment->set_reference(reference);
	vector<EditPair> all_mutations = alignment->getMutations();
	// for (auto m : all_mutations) cerr << m << " ";
	// cerr << endl;
	REQUIRE(all_mutations.size() == 5);
	// validate all locations
	CHECK(all_mutations[0].get_position() == 13);
	CHECK(all_mutations[0].get_code() == MISMATCH_G);
	// also accounts for the deletion
	CHECK(all_mutations[1].get_position() == 16);
	CHECK(all_mutations[1].get_code() == DELETION_N);
	// 2nd mismatch
	CHECK(all_mutations[2].get_position() == 19);
	CHECK(all_mutations[2].get_code() == MISMATCH_A);
	// insertion
	CHECK(all_mutations[3].get_position() == 20);
	CHECK(all_mutations[3].get_code() == INSERTION_G);
	// 3rd mismatch
	CHECK(all_mutations[4].get_position() == 22);
	CHECK(all_mutations[4].get_code() == MISMATCH_T);

}


TEST_CASE( "Testing alignment w/ 2 insertions and a mismatch", "[alignment]" ) {
	// cerr << " ===^^^^====" << endl;
	BamReader reader("bam/cigar_and_md_edits2.bam");
	REQUIRE( reader.is_open() );
	shared_ptr<Alignment> alignment = reader.getNextAlignment();
	REQUIRE(alignment != nullptr);
	vector<EditPair> cigar_muts = alignment->getMutations();
	REQUIRE(cigar_muts.size() == 2);
	CHECK(cigar_muts[0].get_position() == 14);
	CHECK(cigar_muts[0].get_code() == INSERTION_G);
	CHECK(cigar_muts[1].get_position() == 14);
	CHECK(cigar_muts[1].get_code() == INSERTION_G);

	shared_ptr<LongSequence> reference(new LongSequence("bam/hxb2.fa") );
	REQUIRE(reference->read_sequence() );
	alignment->set_reference(reference);
	vector<EditPair> all_mutations = alignment->getMutations();
	// for (auto m : all_mutations) cerr << m << " ";
	// cerr << endl;
	REQUIRE(all_mutations.size() == 3);
	// validate all locations
	CHECK(all_mutations[0].get_position() == 14);
	CHECK(all_mutations[0].get_code() == INSERTION_G);
	// also accounts for the deletion
	CHECK(all_mutations[1].get_position() == 14);
	CHECK(all_mutations[1].get_code() == INSERTION_G);
	// 2nd mismatch
	CHECK(all_mutations[2].get_position() == 18);
	CHECK(all_mutations[2].get_code() == MISMATCH_G);
}


TEST_CASE( "Testing along PacBio read with multiple insertions and deletions", "[alignment]" ) {
	// cerr << " ==========" << endl;
	BamReader reader("bam/single_pacbio_read.sam");
	REQUIRE( reader.is_open() );
	shared_ptr<Alignment> alignment = reader.getNextAlignment();
	REQUIRE(alignment != nullptr);
	// get indels
	vector<EditPair> cigar_muts = alignment->getMutations();

	// cerr << "Indels: ";
	// for (auto m : cigar_muts) cerr << m << " ";
	// cerr << endl;


	REQUIRE(cigar_muts.size() > 10);
	vector<EditPair> true_indels;
	true_indels.emplace_back(EditPair::get_insertion_code('T'), 598);
	true_indels.emplace_back(EditPair::get_insertion_code('A'), 607);
	true_indels.emplace_back(EditPair::get_insertion_code('T'), 622);
	true_indels.emplace_back(EditPair::get_insertion_code('C'), 672);
	true_indels.emplace_back(EditPair::get_insertion_code('A'), 672);
	true_indels.emplace_back(EditPair::get_deletion_code(), 972);
	true_indels.emplace_back(EditPair::get_insertion_code('A'), 1945);

	
	for (auto const & m : true_indels) {
		// cerr << "true indel: " << m << endl;
		bool found = false;
		for (auto const & n : cigar_muts) {
			if (n.get_position() == m.get_position() )
				if (n.get_code() == m.get_code() ) {
					found = true;
				}
		}
		REQUIRE(found);
	}

	// get indels and mismatches
	shared_ptr<LongSequence> reference(new LongSequence("bam/hxb2.fa") );
	REQUIRE(reference->read_sequence() );
	alignment->set_reference(reference);
	vector<EditPair> all_mutations = alignment->getMutations();

	// for (auto m : all_mutations) cerr << m << " ";
	// cerr << endl;

	// validate that a subset of the mutaitons is present
	vector<EditPair> true_mutations;
	true_mutations.emplace_back(EditPair::get_mismatch_code('G'), 	653);
	true_mutations.emplace_back(EditPair::get_deletion_code(), 		659);
	true_mutations.emplace_back(EditPair::get_deletion_code(), 		660);
	true_mutations.emplace_back(EditPair::get_deletion_code(), 		661);
	true_mutations.emplace_back(EditPair::get_deletion_code(), 		662);
	true_mutations.emplace_back(EditPair::get_mismatch_code('A'), 	681 );
	true_mutations.emplace_back(EditPair::get_mismatch_code('A'), 	832 );
	true_mutations.emplace_back(EditPair::get_mismatch_code('A'), 	833 );
	true_mutations.emplace_back(EditPair::get_mismatch_code('T'), 	954 );
	

	for (auto const & m : true_mutations) {
		bool found = false;
		for (auto const & n : all_mutations) {
			if (n.get_position() == m.get_position() )
				if (n.get_code() == m.get_code() ) {
					found = true;
				}
		}
		REQUIRE(found);
	}
}