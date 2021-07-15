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
#include "fastq_reader_wrapper.hpp"
//#include "fastq_writer_wrapper.hpp"
#include "paired_reads.hpp"
#include "universal_sequence.hpp"
#include "import_fasta.hpp"
//alignment
#include "sw_aligner.hpp"
#include "alignment_parameters.hpp"
#include "alignment_report.hpp"
#include "align_primers.hpp"
#include "alignment_reporter.hpp"
#include "compare_aln_scores.hpp"
#include "kmer_index.hpp"
//reporting
#include "basename.hpp"
#include "write_aln_report.hpp"
#include "get_report_filename.hpp"

#include "paired_reads.hpp"
#include "reverse_complement.hpp"
#include <numeric>

using namespace std;

mutex readBarrier;
mutex resultBarrier;

typedef shared_ptr< MutableAlignment > shared_read_ptr;
typedef shared_ptr< paired_reads > paired_reads_ptr;


void do_work(input_parameters ip,
	     paired_reads_ptr read_pair_ptr,
	     AlignmentReporter reporter,
	     shared_ptr< int > counter,
	     shared_ptr< int > workCounter,
	     shared_ptr< int > alnCounter){
  
  SWAligner aligner(ip.align_params, ip.debug_mode);
    
  string read1 = read_pair_ptr->read1_ptr->get_sequence();
  string read2 = read_pair_ptr->read2_ptr->get_sequence();

  aligner.align_pairs(read_pair_ptr->read1_ptr,
		      read_pair_ptr->read2_ptr
		      ip.global_alignment, 2);


  strign read1_rc = reverse_complement(read1)

  cerr << "read1 size = " << read1.size() << endl;
  exit(0);
  
  // align probes to reads
  /*
  if(filtQuery.size() > 0){
    align_primers(ip, aligner, read_ptr->,
		  filtQuery, query_alignments, "+",
		  aln_len);
  }
  */
  
  //prevent collisions in reporting
  resultBarrier.lock();
    
  //keep track of count
  *counter = *counter + 1;
  if(*counter % 1000 == 0){
    if(*counter == 1000){
      cerr << "aligned " << *counter/1000 << "K"; 
    }
    else{
      cerr << "\r" << "aligned " << *counter/1000 << "K";
    }
  }
  resultBarrier.unlock();
}


void consumer(input_parameters ip,
	      shared_ptr< vector < paired_reads_ptr > > read_queue,
	      AlignmentReporter reporter,
	      shared_ptr< bool> has_data,
	      shared_ptr< int > counter,
	      shared_ptr< int > workCounter,
	      shared_ptr< int > alnCounter){
  
  //////////////////////////////////////////////////////
  // wait for reads in queue, break when queue is empty 
  //////////////////////////////////////////////////////
  while(true){
    readBarrier.lock(); // prevent collisions
    if(*has_data){
      if(!read_queue->empty()){
	paired_reads_ptr read_pair = read_queue->back();	
	read_queue->pop_back();
	readBarrier.unlock(); // concurrently align reads
	do_work(ip, read_pair, reporter,
		counter, workCounter, alnCounter);
	
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


int main(int argc, char * argv []) {
  
  //////////////////////////////////
  // Initialize classes & variables
  //////////////////////////////////  
  ArgParser arg_parser(argc, argv);
  input_parameters ip = arg_parser.parse_arguments();
  //shared_ptr< vector < shared_read_ptr > > read_queue(new vector< shared_read_ptr >());
  shared_ptr< bool > has_data(new bool(true));
  vector<thread> threads;
  shared_ptr< int > counter(new int(0));
  shared_ptr< int > workCounter(new int(0));
  shared_ptr< int > alnCounter(new int(0));
  shared_ptr<MutableAlignment> read_ptr;
  //shared_ptr<InputParser> reader;
  shared_ptr<InputParser> reader1;
  shared_ptr<InputParser> reader2;
  shared_ptr< MutableAlignment > read1_ptr;
  shared_ptr< MutableAlignment > read2_ptr;
  reader1 = shared_ptr<InputParser>(new FastqReaderWrapper(ip.read1_path));
  reader2 = shared_ptr<InputParser>(new FastqReaderWrapper(ip.read2_path));

  shared_ptr< vector < paired_reads_ptr > > read_queue(new vector< paired_reads_ptr >());
  
  auto time_start = chrono::system_clock::now();
  string out_file = get_report_filename("./", ip.read1_path, "_alignments.tsv");
  AlignmentReporter reporter(ip.max_report, out_file);
  
  //////////////////////
  // Initialize threads
  //////////////////////  
  for(int i = 0; i < ip.n_threads; ++i){
    threads.push_back( thread(consumer, ip, read_queue,
			      reporter, 
			      has_data, counter, workCounter, alnCounter) );
  }

  //////////////////////////////////////////////////////
  // Loop through the input file and perform alignments
  //////////////////////////////////////////////////////
  cerr << "aligning..." << endl;
  int count = 0;
  int total_aln = 0;
  while ( true ){
    shared_ptr< paired_reads > read_pair(new paired_reads());
    read1_ptr = reader1->getNextSequence();
    read2_ptr = reader2->getNextSequence();
    if(read1_ptr != nullptr and read2_ptr != nullptr){
      if(ip.n_reads > 0){
	if(count > ip.n_reads){
	  break;
	}
      }
    }
    else{
      break;
    }
    ++count;
    read_pair->read1_ptr = read1_ptr;
    read_pair->read2_ptr = read2_ptr;
    read_pair->read_name = read1_ptr->get_read_id();
    //add read pairs to queue
    readBarrier.lock();
    read_queue->emplace_back(read_pair); 
    readBarrier.unlock();
    // stall if queue back up
    if(read_queue->size() > 10000){
     usleep(10000); // 0.1 seconds
    }
  }
  
  ////////////////////////////////
  // wait for the queue to empty
  ////////////////////////////////
  while(true){
   if(read_queue->size() == 0){
     break;
   }
   else{
     usleep(1000);
   }
  }
  *has_data = false;

  //////////////////////////
  // join and close threads
  //////////////////////////
  for(int i; i<threads.size(); ++i){
    threads[i].join();
  }
    
  //////////////////////////////
  // write report, record time
  /////////////////////////////
  cerr << endl << "total reads with alignment = " << *alnCounter << endl;
  reporter.close_files();
  //string out_file = get_report_filename("./", ip.read_path, "_alignments.tsv");
  //write_aln_report(out_file, queryHits);
  cerr << "total query alignments performed = " << *workCounter/1000 << "K" << endl;
  auto time_end = chrono::system_clock::now();
  auto time_diff = time_end-time_start;
  cout << "total time (seconds)= " << time_diff.count()/1000000000.0f << endl;	
  
  return 0;
}

