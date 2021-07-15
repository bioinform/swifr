#ifndef FASTQ_WRITER_WRAPPER_LIB
#define FASTQ_WRITER_WRAPPER_LIB

#include <memory>
#include <iostream>
#include <fstream>
#include "output_parser.hpp"
#include "universal_sequence.hpp"
#include "basename.hpp"
#include "paired_reads.hpp"

using namespace std;

class FastqWriterWrapper : public OutputParser {

private:
  map< string, ofstream* > file_handler_;
  
public:
  string out_dir_;
  
  FastqWriterWrapper(string output_dir){
    out_dir_ = output_dir;
  }

  
  void writeSequence(input_parameters ip, shared_ptr< paired_reads > read_ptr){

    string outFile;
    string inputName1 = ip.read_path;
    size_t lastindex1 = inputName1.find_last_of(".");
    string baseName1 = basename(inputName1.substr(0, lastindex1));
    outFile1 = out_dir_ + "/" + baseName1 + "_swift.fastq";
    
    string outFile2;
    string inputName2 = ip.read2_path;
    size_t lastindex2 = inputName2.find_last_of(".");
    string baseName2 = basename(inputName2.substr(0, lastindex2));
    outFile2 = out_dir_ + "/" + baseName2 + "_swift.fastq";
    
    ////////////////////////////////////////////////////
    // write fastq lines, keep output files dynamically
    ////////////////////////////////////////////////////
    auto file_search1 = file_handler_.find(outFile1);
    if(file_search1 == file_handler_.end()){
      file_handler_[outFile1] = new ofstream(outFile1);     
    }
    auto myfile1 = file_handler_.find(outFile1);        
    *myfile1->second << "@" + read_ptr->read_name + "\n";
    *myfile1->second << read_ptr->read1_trim_seq + "\n";
    *myfile1->second << "+\n";
    *myfile1->second << read_ptr->read1_trim_qual + "\n";    

    //write read2...gonna imporve this later
    auto file_search2 = file_handler_.find(outFile2);
    if(file_search2 == file_handler_.end()){
      file_handler_[outFile2] = new ofstream(outFile2);     
    }
    auto myfile2 = file_handler_.find(outFile2);

    *myfile2->second << "@" + read_ptr->read_name + "\n";
    *myfile2->second << read_ptr->read2_trim_seq + "\n";
    *myfile2->second << "+\n";
    *myfile2->second << read_ptr->read2_trim_qual + "\n";    

  }

  void close_files(){
    for(auto& file_handle : file_handler_){
      delete file_handle.second;
      file_handle.second = 0;
    }
  }

    
};

#endif
