#ifndef INPUT_PARSER_LIB
#define INPUT_PARSER_LIB

#include <iostream>
#include <memory>
#include "io_lib_wrapper/mutable_alignment.hpp"

class InputParser {
public:
  InputParser() {}

  virtual shared_ptr<MutableAlignment> getNextSequence() = 0;
  
};

#endif 
