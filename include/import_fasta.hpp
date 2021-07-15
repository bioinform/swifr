#ifndef IMPORT_FASTA_HPP
#define IMPORT_FASTA_HPP

#include <vector>
#include <memory>
#include "io_lib_wrapper/fasta_reader.h"
#include "input_parser.hpp"
#include "universal_sequence.hpp"
#include "fasta_reader_wrapper.hpp"

using namespace std;

vector< shared_ptr< MutableAlignment > > import_fasta(string primer_path){
  shared_ptr<InputParser> seq_reader;
  shared_ptr< MutableAlignment > seq_ptr;
  vector< shared_ptr< MutableAlignment > > seqs;
  if(primer_path == ""){
    return seqs;
  }
  seq_reader = shared_ptr<InputParser>(new FastaReaderWrapper(primer_path) );
  while ( true ){
    seq_ptr = seq_reader->getNextSequence();
    if(seq_ptr == nullptr){
      break;
    }
    seqs.push_back(seq_ptr);
  }
  //confirm not empty
  if(seqs.size() == 0){
    cerr << "no sequences found in fasta file " << primer_path << endl;
    exit(1);
  }
  return seqs;
}

#endif
