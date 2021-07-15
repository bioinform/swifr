#define CATCH_CONFIG_MAIN 

#include <vector>
#include <string>
#include <memory>
#include <iostream>
#include <fstream>
#include <unordered_map>
#include <map>
#include <algorithm>
#include <sys/stat.h>
#include <unistd.h>
#include <thread>
#include <mutex>
#include <chrono>

#include "io_lib_wrapper/fasta_reader.h"
#include "alignment_report.hpp"
#include "alignment_parameters.hpp"
#include "sw_aligner.hpp"
#include "io_lib_wrapper/mutable_alignment.hpp"
#include "fastq_reader_wrapper.hpp"
#include "universal_sequence.hpp"
#include "basename.hpp"
#include "catch.hpp"

using namespace std;

/*
These test cases currently assume the preset:

struct alignment_parameters {
  int match = 1;
  int mismatch = -2;
  int insertion_open = -2;
  int insertion_extend = -1;
  int deletion_open = -1;
  int deletion_extend = -1;
  int min_aln_score = 10;
}; 
*/


TEST_CASE( "Testing exact match", "[sw_aligner]" ) {
  alignment_parameters default_settings;
  default_settings.min_aln_score = 9;
  bool debug = false;
  bool global_alignment = false;
  SWAligner aligner(default_settings, debug);
  string barcode_str = "TACGACGTC";
  string read_str = "ATATACGACGTCATA";
  shared_ptr<MutableAlignment> read(new MutableAlignment("test_read", read_str));
  shared_ptr<MutableAlignment> barcode(new MutableAlignment("test_barcode", barcode_str));
  shared_ptr< vector<alignment_report> > alignments(new vector<alignment_report>());

  aligner.align(barcode, read, global_alignment, alignments);

  REQUIRE(alignments->size() == 1);
  alignment_report aln = alignments->at(0);
  REQUIRE(aln.reference_start == 4);
  REQUIRE(aln.reference_end == 12);
  REQUIRE(aln.cigar == "9M");
  REQUIRE(aln.query_start == 1);
  REQUIRE(aln.query_end == 9);
  REQUIRE(aln.aln_score == 9);
  REQUIRE(aln.cigar_values == "MMMMMMMMM");
  REQUIRE(aln.strand == "+");
  REQUIRE(aln.reference_name == "test_read");
  REQUIRE(aln.query_name == "test_barcode");   
}

TEST_CASE( "Testing deletion, insertion, and mismatch ", "[sw_aligner]" ) {
  alignment_parameters default_settings;
  default_settings.min_aln_score = 6;
  bool debug = false;
  bool global_alignment = false;
  SWAligner aligner(default_settings, debug);
  string barcode_str = "TACGTTGACACACGTC";
  string read_str = "AATATACGATGAACACCGTCATA";
  shared_ptr<MutableAlignment> read(new MutableAlignment("test_read", read_str));
  shared_ptr<MutableAlignment> barcode(new MutableAlignment("test_barcode", barcode_str));
  shared_ptr< vector<alignment_report> > alignments(new vector<alignment_report>());

  aligner.align(barcode, read, global_alignment, alignments);
  
  REQUIRE(alignments->size() == 1);
  alignment_report aln = alignments->at(0);
  REQUIRE(aln.reference_start == 5);
  REQUIRE(aln.reference_end == 20);
  REQUIRE(aln.cigar == "4M1X2M1D4M1I4M");
  REQUIRE(aln.query_start == 1);
  REQUIRE(aln.query_end == 16);
  REQUIRE(aln.aln_score == 6);
  REQUIRE(aln.cigar_values == "MMMMXMMDMMMMIMMMM");
  REQUIRE(aln.strand == "+");
  REQUIRE(aln.reference_name == "test_read");
  REQUIRE(aln.query_name == "test_barcode");   
}


TEST_CASE( "Testing affine deletion penalities", "[sw_aligner]" ) {
  alignment_parameters default_settings;
  default_settings.min_aln_score = 6;
  default_settings.deletion_open = -2;
  bool debug = false;
  bool global_alignment = false;
  SWAligner aligner(default_settings, debug);
  string barcode_str = "ATCGAACGTCAT";
  string read_str = "ATCGAAAACGTCAT";
  shared_ptr<MutableAlignment> read(new MutableAlignment("test_read", read_str));
  shared_ptr<MutableAlignment> barcode(new MutableAlignment("test_barcode", barcode_str));
  shared_ptr< vector<alignment_report> > alignments(new vector<alignment_report>());

  aligner.align(barcode, read, global_alignment, alignments);

  REQUIRE(alignments->size() == 1);
  alignment_report aln = alignments->at(0);
  REQUIRE(aln.reference_start == 1);
  REQUIRE(aln.reference_end == 14);
  REQUIRE(aln.cigar == "4M2D8M");
  REQUIRE(aln.query_start == 1);
  REQUIRE(aln.query_end == 12);
  REQUIRE(aln.aln_score == 9);
  REQUIRE(aln.cigar_values == "MMMMDDMMMMMMMM");
  REQUIRE(aln.strand == "+");
  REQUIRE(aln.reference_name == "test_read");
  REQUIRE(aln.query_name == "test_barcode");   
}

TEST_CASE( "Testing affine insertion penalities", "[sw_aligner]" ) {
  alignment_parameters default_settings;
  default_settings.min_aln_score = 6;
  default_settings.insertion_open = -2;
  bool debug = false;
  bool global_alignment = false;
  SWAligner aligner(default_settings, debug);
  string barcode_str = "ATCGAAACGTCAT";
  string read_str = "ATCGAACGTCAT";
  shared_ptr<MutableAlignment> read(new MutableAlignment("test_read", read_str));
  shared_ptr<MutableAlignment> barcode(new MutableAlignment("test_barcode", barcode_str));
  shared_ptr< vector<alignment_report> > alignments(new vector<alignment_report>());

  aligner.align(barcode, read, global_alignment, alignments);

  REQUIRE(alignments->size() == 1);
  alignment_report aln = alignments->at(0);
  REQUIRE(aln.reference_start == 1);
  REQUIRE(aln.reference_end == 12);
  REQUIRE(aln.cigar == "4M1I8M");
  REQUIRE(aln.query_start == 1);
  REQUIRE(aln.query_end == 13);
  REQUIRE(aln.aln_score == 10);
  REQUIRE(aln.cigar_values == "MMMMIMMMMMMMM");
  REQUIRE(aln.strand == "+");
  REQUIRE(aln.reference_name == "test_read");
  REQUIRE(aln.query_name == "test_barcode");   
}

TEST_CASE( "Testing - strand alignment", "[sw_aligner]" ) {
  alignment_parameters default_settings;
  default_settings.min_aln_score = 6;
  bool debug = false;
  bool global_alignment = false;
  SWAligner aligner(default_settings, debug);
  string barcode_str = "ACTGACA";
  string read_str = "TGTCAGT";
  shared_ptr<MutableAlignment> read(new MutableAlignment("test_read", read_str));
  shared_ptr<MutableAlignment> barcode(new MutableAlignment("test_barcode", barcode_str));
  shared_ptr< vector<alignment_report> > alignments(new vector<alignment_report>());

  aligner.align(barcode, read, global_alignment, alignments);

  REQUIRE(alignments->size() == 1);
  alignment_report aln = alignments->at(0);
  REQUIRE(aln.reference_start == 1);
  REQUIRE(aln.reference_end == 7);
  REQUIRE(aln.cigar == "7M");
  REQUIRE(aln.query_start == 1);
  REQUIRE(aln.query_end == 7);
  REQUIRE(aln.aln_score == 7);
  REQUIRE(aln.cigar_values == "MMMMMMM");
  REQUIRE(aln.strand == "-");
  REQUIRE(aln.reference_name == "test_read");
  REQUIRE(aln.query_name == "test_barcode");   
}


TEST_CASE( "Testing +/- strand, multi alignments", "[sw_aligner]" ) {
  alignment_parameters default_settings;
  default_settings.min_aln_score = 7;
  bool debug = false;
  bool global_alignment = false;
  SWAligner aligner(default_settings, debug);
  string barcode_str = "TGTCAGT";
  string read_str = "ACTGACACGGCTGTCAGT";
  shared_ptr<MutableAlignment> read(new MutableAlignment("test_read", read_str));
  shared_ptr<MutableAlignment> barcode(new MutableAlignment("test_barcode", barcode_str));
  shared_ptr< vector<alignment_report> > alignments(new vector<alignment_report>());

  aligner.align(barcode, read, global_alignment, alignments);

  REQUIRE(alignments->size() == 2);
  alignment_report aln = alignments->at(0);
  REQUIRE(aln.reference_start == 12);
  REQUIRE(aln.reference_end == 18);
  REQUIRE(aln.cigar == "7M");
  REQUIRE(aln.query_start == 1);
  REQUIRE(aln.query_end == 7);
  REQUIRE(aln.aln_score == 7);
  REQUIRE(aln.cigar_values == "MMMMMMM");
  REQUIRE(aln.strand == "+");
  REQUIRE(aln.reference_name == "test_read");   
  REQUIRE(aln.query_name == "test_barcode");
 
  aln = alignments->at(1);
  REQUIRE(aln.reference_start == 1);
  REQUIRE(aln.reference_end == 7);
  REQUIRE(aln.cigar == "7M");
  REQUIRE(aln.query_start == 1);
  REQUIRE(aln.query_end == 7);
  REQUIRE(aln.aln_score == 7);
  REQUIRE(aln.cigar_values == "MMMMMMM");
  REQUIRE(aln.strand == "-");
  REQUIRE(aln.reference_name == "test_read");
  REQUIRE(aln.query_name == "test_barcode");   

}

/*
TEST_CASE( "Testing global alignment", "[sw_aligner]" ) {
  alignment_parameters default_settings;
  default_settings.min_aln_score = 6;
  bool debug =false;
  bool global_alignment = true;
  SWAligner aligner(default_settings, debug);
  string barcode_str = "TGTCAGA";
  string read_str = "ACGAGCATGTCAGTACGAGCA";
  shared_ptr<MutableAlignment> read(new MutableAlignment("test_read", read_str));
  shared_ptr<MutableAlignment> barcode(new MutableAlignment("test_barcode", barcode_str));
  shared_ptr< vector<alignment_report> > alignments(new vector<alignment_report>());

  cerr << "testing" << endl;
  aligner.align(barcode, read, global_alignment, alignments);
  
  REQUIRE(alignments->size() == 1);
  alignment_report aln = alignments->at(0);
  cerr << aln.aln_score << endl;
  
  REQUIRE(aln.reference_start == 8);
  REQUIRE(aln.reference_end == 15);
  REQUIRE(aln.cigar == "6M1D1M");
  REQUIRE(aln.query_start == 1);
  REQUIRE(aln.query_end == 7);
  REQUIRE(aln.aln_score == 6);
  REQUIRE(aln.cigar_values == "MMMMMMDM");
  REQUIRE(aln.strand == "+");
  REQUIRE(aln.reference_name == "test_read");
  REQUIRE(aln.query_name == "test_barcode");   
}
*/

/*
TEST_CASE( "Testing soft-clipped alignment", "[sw_aligner]" ) {
  alignment_parameters default_settings;
  default_settings.min_aln_score = 10;
  bool debug = true;
  bool global_alignment = false;
  SWAligner aligner(default_settings, debug);
  string barcode_str = "TATCTCTCTACTGACTGTCCTC";
  string read_str = "CAAAAAACACGTCAGTAGAGAGATA";
  shared_ptr<MutableAlignment> read(new MutableAlignment("test_read", read_str));
  shared_ptr<MutableAlignment> barcode(new MutableAlignment("test_barcode", barcode_str));
  shared_ptr< vector<alignment_report> > alignments(new vector<alignment_report>());

  cerr << "testing" << endl;
  aligner.align(barcode, read, global_alignment, alignments);
  
  REQUIRE(alignments->size() == 1);
  alignment_report aln = alignments->at(0);
  cerr << aln.aln_score << endl;
    
  REQUIRE(aln.reference_start == 7);
  REQUIRE(aln.reference_end == 15);
  REQUIRE(aln.cigar == "6M1D1M");
  REQUIRE(aln.query_start == 1);
  REQUIRE(aln.query_end == 7);
  REQUIRE(aln.aln_score == 6);
  REQUIRE(aln.cigar_values == "MMMMMMDM");
  REQUIRE(aln.strand == "+");
  REQUIRE(aln.reference_name == "test_read");
  REQUIRE(aln.query_name == "test_barcode");   
}
*/



