#ifndef OUTPUT_PARSER_LIB
#define OUTPUT_PARSER_LIB

#include <iostream>
#include <memory>
#include "universal_sequence.hpp"
#include "input_parameters.hpp"
#include "read_segment.hpp"
#include "paired_reads.hpp"

class OutputParser {

public:
  OutputParser() {}

  virtual void writeSequence(input_parameters ip, shared_ptr< paired_reads > read_ptr) = 0;

  virtual void close_files() = 0;

};

#endif 

