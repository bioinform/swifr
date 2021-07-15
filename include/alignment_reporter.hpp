#ifndef ALIGNMENT_REPORTER_HPP
#define ALIGNMENT_REPORTER_HPP

#include "reverse_complement.hpp"

typedef shared_ptr< MutableAlignment > shared_read_ptr;

using namespace std;

class AlignmentReporter{

private:
  string out_file_;
  int max_report_;
  map< string, ofstream* > file_handler_;  
  
public:
  
  AlignmentReporter(int max_report,
		    string out_file,
		    vector< shared_ptr< MutableAlignment > > sequences){
    max_report_ = max_report;
    out_file_ = out_file;
    auto file_search = file_handler_.find(out_file_);
    if(file_search == file_handler_.end()){
      file_handler_[out_file_] = new ofstream(out_file_);
      auto output = file_handler_.find(out_file);        
      for(auto & refSeq : sequences){
	*output->second << "@SQ\tSN:" + refSeq->get_read_id()
	  + "\tLN:" + to_string(refSeq->get_sequence().size()) + "\n";      
      }
      //*output->second << "name\tref_len\taln_start\taln_end\tstrand\tquery_start\tquery_end\tscore\tcigar\tquery_name\n"; // \tseq\tqual\n";
      //*output->second << "name\tref_len\taln_start\taln_end\tstrand\tquery_start\tquery_end\tscore\tcigar\tquery_name\tseq\tqual\n";      
    }
  }

  void report_alignments(shared_ptr< vector<alignment_report> > alignments,
			 shared_read_ptr read_ptr){

    //filter alignments based on user input
    int count = 0;
    vector<alignment_report> hits;
    if(alignments->size() > 0){
      for(int i=0; i < alignments->size(); ++i){
	++ count;
	if(count > max_report_){
	  break;
	}
	hits.push_back(alignments->at(i));
      }
    }
    
    auto output = file_handler_.find(out_file_);        
    int aln_count = 0;
    int total_aln = hits.size();
    for (auto & aln : hits){
      string seq = read_ptr->get_sequence();
      string qual = read_ptr->get_qualities();
      string reportSeq;
      ++aln_count;
      int bitflag = 0;
      string strand = aln.strand;
      if(strand == "-"){
	bitflag = 16;
	// reverse quality string
	reverse(qual.begin(), qual.end());
	//reverse complement seq string
	reportSeq = reverse_complement(seq);	
      }
      else{
	reportSeq = seq;	
      }
      if(aln_count > 1 and strand == "+"){
	bitflag = 256;	  	
      }
      if(aln_count > 1 and strand == "-"){
	bitflag = 272;	  	
      }
      //*output->second << aln.reference_name + "\t" + to_string(aln.reference_length) + "\t" + to_string(aln.reference_start) + "\t" + to_string(aln.reference_end) + "\t" + aln.strand + "\t" + to_string(aln.query_start) + "\t" + to_string(aln.query_end) + "\t" + to_string(aln.aln_score) + "\t" + aln.cigar + "\t" + aln.query_name + "\t" + seq + "\t" + qual + "\n";

      *output->second << aln.query_name + "\t"	
	+ to_string(bitflag) + "\t"		
	+ aln.reference_name + "\t"		
	+ to_string(aln.reference_start) + "\t"				
	+ "0" + "\t"  // setting MapQ to zero for good alignment
	+ aln.cigar + "\t"				
	+ "*" + "\t"  
	+ "0" + "\t"  
	+ "0" + "\t"  
	+ reportSeq + "\t"  
	+ qual + "\t"			    
	+ "AS:i:" + to_string(aln.aln_score) + "\t"	
	+ "NM:i:" + to_string(aln.edit_distance) + "\n";           
      
	}
  } 
  
  void close_files(){
    for(auto& file_handle : file_handler_){
      delete file_handle.second;
      file_handle.second = 0;
    }
  }
};


#endif
