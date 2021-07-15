#ifndef FASTQ_READER_WRAPPER_LIB
#define FASTQ_READER_WRAPPER_LIB

#include <memory>
#include "io_lib_wrapper/fasta_reader.h"
#include "input_parser.hpp"
#include "universal_sequence.hpp"
#include "basename.hpp"

class FastqReaderWrapper : public InputParser {
private:
  FastaReader reader_;

public:
  
  
  FastqReaderWrapper(const string & path):
    reader_(path.c_str()) {
  }
  
  shared_ptr<MutableAlignment> getNextSequence() {
    kseq_t * s = reader_.nextSequence();
    if (s != nullptr){

      string seqName = s->name.s;
      string seqComment = "*";
      string seqString = s->seq.s;
      string seqQual = s->qual.s;
      
      shared_ptr<MutableAlignment> mutAln(new MutableAlignment(seqName, seqString, seqQual, seqComment));
      return mutAln;
    }
    else{
      return nullptr;
    }
  }  
  
  
};

#endif
