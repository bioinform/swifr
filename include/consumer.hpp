#ifndef CONSUMER_HPP
#define CONSUMER_HPP

#include <vector>
#include <string>
#include <memory>
#include <iostream>
#include <fstream>
#include <unordered_map>
#include <map>
#include <algorithm>
#include <sys/stat.h>
#include <unistd.h>
#include <thread>
#include <mutex>
#include <chrono>
#include <regex>
#undef max
#undef min

//arg parsing
#include "arg_parsing.hpp"
#include "input_parameters.hpp"
//input-output
#include "io_lib_wrapper/mutable_alignment.hpp"
#include "io_lib_wrapper/fasta_reader.h"
#include "fasta_reader_wrapper.hpp"
#include "fasta_writer_wrapper.hpp"
#include "fastq_reader_wrapper.hpp"
#include "fastq_writer_wrapper.hpp"
#include "paired_reads.hpp"
#include "universal_sequence.hpp"
//alignment
#include "sw_aligner.hpp"
#include "alignment_parameters.hpp"
#include "alignment_report.hpp"
#include "align_primers.hpp"
#include "filter_alignments.hpp"
#include "compare_aln_scores.hpp"
//reporting
#include "basename.hpp"
#include "write_aln_report.hpp"
#include "get_report_filename.hpp"

#include "do_work.hpp"

using namespace std;

void consumer(input_parameters ip,
	      shared_ptr< vector < paired_reads_ptr > > read_queue,
	      vector< shared_ptr< MutableAlignment > > vseqs,
	      vector< shared_ptr< MutableAlignment > > vprobes,
	      vector< shared_ptr< MutableAlignment > > jprobes,
	      shared_ptr< vector<alignment_report> > vhits,
	      shared_ptr< vector<alignment_report> > jhits,
	      shared_ptr< vector<alignment_report> > targetHits,
	      shared_ptr< bool> has_data,
	      shared_ptr< OutputParser > writer,
	      shared_ptr< int > counter){
  
  //////////////////////////////////////////////////////
  // wait for reads in queue, break when queue is empty 
  //////////////////////////////////////////////////////
  while(true){
    readBarrier.lock(); // prevent collisions
    if(*has_data){
      if(!read_queue->empty()){
	paired_reads_ptr read = read_queue->back();
	read_queue->pop_back();
	readBarrier.unlock(); // concurrently align reads
	do_work(ip, read, vseqs, vprobes, jprobes, vhits, jhits, targetHits, writer, counter);
      }
      else{
	readBarrier.unlock();
      }
    }
    else{
      readBarrier.unlock();
      break;
    }
  }
}

#endif
