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
#include <set>

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
#include "bloom_filter.hpp"

#include <numeric>

using namespace std;

mutex readBarrier;
mutex resultBarrier;

typedef shared_ptr< MutableAlignment > shared_read_ptr;
//typedef tuple <shared_ptr< MutableAlignment >,char> indexType;
typedef tuple <shared_ptr< MutableAlignment >,char, set<int>> indexType;
typedef shared_ptr< KmerIndex > index_ptr;

void do_work(input_parameters ip,
	     shared_read_ptr read_ptr,
	     index_ptr queryIndex,
	     AlignmentReporter reporter,
	     shared_ptr< vector<alignment_report> > queryHits,
	     shared_ptr< int > counter,
	     shared_ptr< int > workCounter,
	     shared_ptr< int > alnCounter){
  
  SWAligner aligner(ip.align_params, ip.debug_mode);
  shared_ptr< vector<alignment_report> > query_alignments(new vector< alignment_report >());

  vector< indexType > filtQuery;

  // first check for exact match
    
  // returning a smaller set of probes to align the read against
  if(read_ptr->get_sequence().size() > ip.kmer_size+5){

    if(ip.kmer_size > 0){
      filtQuery = queryIndex->filter_by_kmers(read_ptr->get_sequence(), false);
    }
    //pass all seqs in +/- orientation
    else{
      filtQuery = queryIndex->all_seqs();
    }
    
    if(filtQuery.size() == 0 && ip.complete_search){
      filtQuery = queryIndex->all_seqs();
    }
    
    int aln_len = read_ptr->get_sequence().size();
    
    if(ip.aln_len > 0){
      if(aln_len > ip.aln_len){
	aln_len = ip.aln_len;
      }
    }
    
    // align probes to reads
    if(filtQuery.size() > 0){
      align_primers(ip, aligner, read_ptr,
		    filtQuery, query_alignments, 
		    aln_len);
    }
  }  
  //prevent collisions in reporting
  resultBarrier.lock();
  
  //store alignments in table
  *workCounter = *workCounter + query_alignments->size();   
  if( query_alignments->size() > 0){
    *alnCounter = *alnCounter + 1;   
  }
  reporter.report_alignments(query_alignments, read_ptr);

  //keep track of count
  *counter = *counter + 1;
  if(ip.verbose == true){
    if(*counter % 1000 == 0){
      if(*counter == 1000){
	cerr << "aligned " << *counter/1000 << "K"; 
      }
      else{
	cerr << "\r" << "aligned " << *counter/1000 << "K";
      }
    }
  }
  resultBarrier.unlock();
}


void consumer(input_parameters ip,
	      shared_ptr< vector < shared_read_ptr > > read_queue,
	      index_ptr queryIndex,
	      AlignmentReporter reporter,
	      shared_ptr< vector<alignment_report> > queryHits,
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
	shared_ptr<MutableAlignment> read = read_queue->back();
	read_queue->pop_back();
	readBarrier.unlock(); // concurrently align reads
	do_work(ip, read, queryIndex, reporter,
		queryHits, counter, workCounter, alnCounter);

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
  shared_ptr< vector < shared_read_ptr > > read_queue(new vector< shared_read_ptr >());
  shared_ptr< bool > has_data(new bool(true));
  vector<thread> threads;
  shared_ptr< int > counter(new int(0));
  shared_ptr< int > workCounter(new int(0));
  shared_ptr< int > alnCounter(new int(0));
  vector< shared_read_ptr > query_seqs = import_fasta(ip.query_path);
  shared_ptr< vector<alignment_report> > queryHits(new vector<alignment_report>());
  shared_ptr<MutableAlignment> read_ptr;
  shared_ptr<InputParser> reader;
  reader = shared_ptr<InputParser>(new FastqReaderWrapper(ip.read_path));    
  //KmerIndex queryIndex(query_seqs, ip.kmer_size, ip.kmer_mismatches, ip.kmer_freq);  
  //shared_ptr<KmerIndex> index_ptr(new KmerIndex(query_seqs, ip.kmer_size, ip.kmer_mismatches, ip.kmer_freq));

  
  vector<shared_ptr<KmerIndex>> index_ptrs;
  cerr << "building index..." << endl;
  for(int i = 0; i < ip.n_threads; ++i){
    shared_ptr<KmerIndex> index_ptr(new KmerIndex(query_seqs, ip.kmer_size, ip.kmer_mismatches, ip.kmer_freq));
    index_ptrs.push_back(index_ptr);
  }

  
  auto time_start = chrono::system_clock::now();
  //string out_file = get_report_filename("./", ip.read_path, "_alignments.sam");
  string out_file = ip.output_basename + ".sam";
  AlignmentReporter reporter(ip.max_report, out_file, query_seqs);
  
  //////////////////////
  // Initialize threads
  //////////////////////  
  for(int i = 0; i < ip.n_threads; ++i){
    threads.push_back( thread(consumer, ip, read_queue,
			      index_ptrs[i], reporter, queryHits,
			      has_data, counter, workCounter, alnCounter) );
  }
  
  //////////////////////
  // Add Reads to Queue
  //////////////////////
  cerr << "aligning..." << endl;
  int count = 0;
  int total_aln = 0;
  while ( true ){
    shared_ptr< paired_reads > read_pair(new paired_reads());
    read_ptr = reader->getNextSequence();
    if(read_ptr != nullptr){
      if(ip.n_reads > 0){
	if(count >= ip.n_reads){
	  break;
	}
      }
    }
    else{
      break;
    }
    ++count;
    readBarrier.lock();
    read_queue->emplace_back(read_ptr); 
    readBarrier.unlock();
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

