#ifndef READ_SEGMENT_HPP
#define READ_SEGMENT_HPP

#include "io_lib_wrapper/mutable_alignment.hpp"
using namespace std;

struct read_segment {
  string sequence;
  string qualities;
  string read_name;
  string start_pos = "";
  string end_pos = "";
  shared_ptr<MutableAlignment> read_ptr;
};

#endif
