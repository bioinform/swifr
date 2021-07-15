#ifndef DO_WORK_HPP
#define DO_WORK_HPP

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

using namespace std;

void do_work(input_parameters ip,
	     paired_reads_ptr read_ptr,
	     vector< shared_ptr< MutableAlignment > > vseqs,
	     vector< shared_ptr< MutableAlignment > > vprobes,
	     vector< shared_ptr< MutableAlignment > > jprobes,
	     shared_ptr< vector<alignment_report> > vhits,
	     shared_ptr< vector<alignment_report> > jhits,
	     shared_ptr< vector<alignment_report> > targetHits,
	     shared_ptr< OutputParser > writer,
	     shared_ptr< int > counter){
  
  SWAligner aligner(ip.align_params, ip.debug_mode);
  shared_ptr< vector<alignment_report> > vseq_alignments(new vector<alignment_report>());
  shared_ptr< vector<alignment_report> > jprobe_alignments(new vector<alignment_report>());
  shared_ptr< vector<alignment_report> > vprobe_alignments(new vector<alignment_report>());

  //look into making a decorator to search index with kmers from read pointer,
  // returning a smaller set of probes to align against the read
  // V primers
  align_primers(ip, aligner, read_ptr->read2_ptr,
		vprobes, vprobe_alignments, "+", 70);
  
  // J primers
  align_primers(ip, aligner, read_ptr->read1_ptr,
		jprobes, jprobe_alignments, "+", 40);

  /*
  // V gene region 
  align_primers(ip, aligner, read_ptr->read2_ptr,
		vseqs, vseq_alignments, "+",
		read_ptr->read2_ptr->get_sequence().size());
  */

  /////////////////////
  // filter and report
  /////////////////////  
  resultBarrier.lock();
  
  *counter = *counter + 1;
  if(*counter % 1000 == 0){
    if(*counter == 1000){
      cerr << "aligned " << *counter/1000 << "K"; //<< endl;
    }
    else{
      cerr << "\r" << "aligned " << *counter/1000 << "K";
    }
  }

  filter_alignments(vprobe_alignments, vhits, ip.vprobe_start);
  filter_alignments(jprobe_alignments, jhits, ip.jprobe_start);
    
  if(vseq_alignments->size() > 0){
    for(int i=0; i < vseq_alignments->size(); ++i){
      targetHits->push_back(vseq_alignments->at(i));
    }
  }
  
  resultBarrier.unlock();
}

#endif
