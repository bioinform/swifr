#ifndef FASTA_READER_WRAPPER_LIB
#define FASTA_READER_WRAPPER_LIB

#include <memory>
#include "io_lib_wrapper/fasta_reader.h"
#include "input_parser.hpp"
#include "universal_sequence.hpp"

class FastaReaderWrapper : public InputParser {
private:
  FastaReader reader_;

public:


  FastaReaderWrapper(const string & path):
    reader_(path.c_str()) {
  }

  shared_ptr<MutableAlignment> getNextSequence() {
    kseq_t * s = reader_.nextSequence();
    if (s != nullptr){

      string seqName = s->name.s;
      string seqIndex = "*";
      string seqString = s->seq.s;
      string seqQual = "*";
      shared_ptr<MutableAlignment> mutAln(new MutableAlignment(seqName, seqString, seqQual, seqIndex));

      return mutAln;
    }
    else{
      return nullptr;
    }
  }


};

#endif
