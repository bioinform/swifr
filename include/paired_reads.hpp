#ifndef PAIRED_READS_HPP
#define PAIRED_READS_HPP

using namespace std;

//alignment object check against
struct paired_reads {
  string read_name;
  string ref_index;
  shared_ptr<MutableAlignment> read1_ptr;
  string read1_trim_seq;
  string read1_trim_qual;
  string read2_trim_seq;
  string read2_trim_qual;
};


#endif
