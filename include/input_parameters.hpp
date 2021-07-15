#ifndef INPUT_PARAMETERS_HPP
#define INPUT_PARAMETERS_HPP
#include "alignment_parameters.hpp"

using namespace std;

////////////////////////////////////////
//  Command-line parameters structure
////////////////////////////////////////
struct input_parameters {
  string read_path;
  string query_path;
  int max_report = 5;
  string output_basename; 
  string output_file_type;
  alignment_parameters align_params; 
  bool debug_mode = false;
  bool global_alignment = false;
  bool alignment_report = false;
  bool complete_search = false;
  bool verbose = false;
  int n_reads = -1;
  int n_threads = 1;
  int aln_len = -1;
  int kmer_size;
  int kmer_mismatches;
  float kmer_freq = 0.5;
  
};

#endif
