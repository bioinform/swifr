#ifndef GET_REPORT_FILENAME_HPP
#define GET_REPORT_FILENAME_HPP

#include "basename.hpp"
using namespace std;

string get_report_filename(string out_dir, string filename, string extension){
  filename = basename(filename);
  size_t lastindex = filename.find_last_of(".");
  int replace_len = filename.size() - lastindex;
  string out_file = filename.replace(lastindex, replace_len, extension);
  string out_path = out_dir + "/" + out_file;
  return out_path;
}

#endif
