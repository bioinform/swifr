/**
This is an object designed to perform smith waterman alignment between two strings 
and produce a list of alignments, if any are found.

Primary Functionality:
1) allocate memory for score and traceback matrices
2) score matrices using the smith-waterman method (parameters: match, mismatch, gap open, gap extend)
3) find local maxima in the score matrix satisfying a minimum alignment score.
4) trace all alignments producing cigars
5) optionally identify end-to-end alignments
6) deallocate memory 
7) report all aspects of the alignment:

reference start= 291
reference end= 337
query start= 1
query end= 45
cigar= 2M1I8M2D7M1X17M1D9M1S
alignment score= 28
strand = +
GTGTAAAACGA--CGGCCAGTAAAAACGGAGGAGGAGG-ACAGTCAGTA
|| ||||||||  |||||||*||||||||||||||||| ||||||||| 
GT-TAAAACGACTCGGCCAGAAAAAACGGAGGAGGAGGAACAGTCAGT 
 
*/

/*
TODO
consider global alignments, right now I am searching for maxima at the end of the query sequence
however, the trace can just as easily terminate in the middle of the read
I could score globally and trace locally? (usefuly for demux)
Should I force the trace globally?
Should I permit 1-2 bases soft clip for global?
Considering barcode sequences, with a global score, local trace might still be ok.

define a function which also takes coordinates of a read
which will specify which part of the reference sequence to align against.

*/
#ifndef SWALIGNER_HPP
#define SWALIGNER_HPP

#include <iostream>
#include <string>
#include <algorithm>
#include <tuple>

#include "alignment_parameters.hpp"
#include "alignment_report.hpp"
#include "io_lib_wrapper/mutable_alignment.hpp"


using namespace std;

struct maxima_coords {
  int query_pos = 0;
  int ref_pos = 0;
  int aln_score = 0;
};
  
//////////////////////////////////////////////////////////////////////
// an object for performing Smith-Waterman alignments on two strings
//////////////////////////////////////////////////////////////////////
class SWAligner{
  
  alignment_parameters aln_settings_;
  string query_;
  string reference_;
  int * reference_maxima_;
  int * query_maxima_index_;
  int ** aln_array_;
  int ** traceback_matrix_;
  bool debug_;
  bool global_aln_;
  int search_dist_;
  string strand_;

  
  string reverse_complement_(string sequence){
    reverse(sequence.begin(), sequence.end());
    string reverseComp;
    for(int i=0; i < sequence.size(); ++i){
      char letter = sequence[i];
      if(letter == 'A'){
	reverseComp.append("T");
      }
      if(letter == 'T'){
	reverseComp.append("A");
      }
      if(letter == 'C'){
	reverseComp.append("G");
      }
      if(letter == 'G'){
	reverseComp.append("C");
      }
      if(letter == 'N'){
	reverseComp.append("N");
      }

    }
    return reverseComp;
  }
  
  
  ////////////////////////////////////////////////////
  //  allocate memory for alignment/traceback array
  ///////////////////////////////////////////////////
  int ** allocate_matrix_(){
    int que_pos = 0;
    int ref_pos = 0;
    // allocate new memory for array, with extra row & col for SW 
    int ** maloc_array = new int*[(query_.size()+1)];
    for(int i = 0; i < (query_.size()+1); ++i){
      maloc_array[i] = new int[(reference_.size()+1)](); //initialize to set to zero
    }
    return maloc_array;
  }
  
  ////////////////////////////////////////////////////////////
  //  allocate memory to index the maximum scoring alignments
  ////////////////////////////////////////////////////////////
  int * allocate_array_(){
    int * reference_array = new int[reference_.size()]();
    return reference_array;
  }
  
  ////////////////////////////////////////////////////
  //  allocate memory for alignment/traceback array
  ///////////////////////////////////////////////////  
  void allocate_memory_(){
    aln_array_ = allocate_matrix_();
    traceback_matrix_ = allocate_matrix_();
    query_maxima_index_ = allocate_array_();
    reference_maxima_ = allocate_array_();
  }

  ////////////////////////////////////////////////
  //  fill the positions in the alignment array
  ////////////////////////////////////////////////
  void score_matrices_(){
    string ref_base, que_base;
    int que_pos = 1;
    while(que_pos <= query_.size()){
      int ref_pos = 1;
      que_base = query_.substr((que_pos-1),1);
      while(ref_pos <= reference_.size()){
        ref_base = reference_.substr((ref_pos-1),1);
    	int score = calc_score_(ref_base, que_base, ref_pos, que_pos);
	aln_array_[que_pos][ref_pos] = score;
	// keeping track of maxima for tracing alignments
	if(score > reference_maxima_[ref_pos - 1]){
	  reference_maxima_[ref_pos - 1] = score;
	  query_maxima_index_[ref_pos - 1] = que_pos;
	}
	++ref_pos;
      }
      ++que_pos;
    }
  }
  
  ///////////////////////////////////////////////////////////////////////
  // define a function to calculate scores for positions in the matrices
  ///////////////////////////////////////////////////////////////////////
  int calc_score_(const string & ref_base,
		 const string & que_base,
		 int ref_pos,
		 int que_pos) {
    int score = 0;
    int all_scores[4] = {0,0,0,0};
    bool bases_match = false;
    //the 4th position remains zero, we reset the score to zero to support local alignments
    if(ref_base == que_base){
      bases_match = true;
      score = aln_settings_.match;
    }
    else{
      score = aln_settings_.mismatch;
    }
    // match or mismatch
    all_scores[0] = (aln_array_[que_pos-1][ref_pos-1] + score);
    // moving forward in the query = insertion
    all_scores[1] = (aln_array_[que_pos-1][ref_pos] + aln_settings_.insertion_open);
    // extended insertion
    if(que_pos > 2){
      // as long as previous base does not match, score extension
      if(traceback_matrix_[que_pos-1][ref_pos] != 0){ 
     	all_scores[1] = (aln_array_[que_pos-1][ref_pos] + aln_settings_.insertion_extend);
      }
    }
    // moving forward in the reference = deletion
    all_scores[2] = (aln_array_[que_pos][ref_pos-1] + aln_settings_.deletion_open);
    // extended deletion
    if(ref_pos > 2){
      // as long as previous base does not match, score extension
      if(traceback_matrix_[que_pos][ref_pos-1] != 0){ 
	all_scores[2] = (aln_array_[que_pos][ref_pos-1] + aln_settings_.deletion_extend);
      }
    }
    // trace the alignment
    trace_position_(que_pos, ref_pos, bases_match, all_scores);
    // if global, search for max in first 3 elements, since I don't want to restart to zero
    // it is necessary to support negative values to score the string end-to-end or globally
    int best_score = *max_element(all_scores, all_scores + search_dist_);
    return best_score;
  }

  ////////////////////////////////////////////////////////////////////
  //  Define the alignment path that should be taken given the score
  ////////////////////////////////////////////////////////////////////
  int trace_position_(int que_pos,
		      int ref_pos,
		      bool bases_match,
		      int all_scores[4]){
    int trace = 0;
    trace = distance(all_scores, max_element(all_scores, all_scores + 3));
    if(!bases_match and trace == 0){
      traceback_matrix_[que_pos][ref_pos] = -1;
    }
    else{
      traceback_matrix_[que_pos][ref_pos] = trace;
    }
  }

  /////////////////////////////////////////////////////////////////////////////////////////
  // find all alignments satisfying the minimum score, filter overlaping query alignments
  /////////////////////////////////////////////////////////////////////////////////////////
  vector< maxima_coords > find_local_maxima_(){
    vector< maxima_coords > local_maxima;
    if(debug_){
      cerr << "find local maxima" << endl;
    }
    // find all local maxima satisfying minimum score
    int score;
    int next_score;
    int last_score;
    bool is_maxima;
    for(int r = 2; r < reference_.size(); ++r){
      //starting from position 2 is fine, who would align sequences < 3nt anyway? (barcode or other)
      is_maxima = false;
      last_score = reference_maxima_[(r - 2)];
      score = reference_maxima_[(r - 1)];
      next_score = reference_maxima_[r];
      if(score >= aln_settings_.min_aln_score and score > last_score and score >= next_score){ 
	is_maxima = true;
      }
      else if(next_score > score and r == (reference_.size()-1) and next_score >= aln_settings_.min_aln_score){
	is_maxima = true;
	r = r + 1;
	score = reference_maxima_[(r - 1)];
      }
      //if point is maxima, keep it
      if(is_maxima){
	maxima_coords coords;
	int q = query_maxima_index_[(r - 1)];
	coords.query_pos = q;
	coords.ref_pos = r;
	coords.aln_score = score;
	local_maxima.push_back(coords);
      }
    }    
    if(debug_){
      cerr << "total maxima= " << local_maxima.size() << endl;
    }
    return filter_maxima_(local_maxima);
  }  
  
  /////////////////////////////////////////////////////////////////////////////////////////
  // find all alignments satisfying the minimum score, filter overlaping query alignments
  /////////////////////////////////////////////////////////////////////////////////////////
  vector< maxima_coords > find_global_maxima_(){
    vector< maxima_coords > global_maxima;
    if(debug_){
      cerr << "find global maxima" << endl;
    }
    // find all global maxima satisfying minimum score
    int score;
    int next_score;
    int last_score;
    bool is_maxima;
    for(int r = 3; r < reference_.size(); ++r){
      is_maxima = false;
      last_score = aln_array_[query_.size()][(r-2)];
      score = aln_array_[query_.size()][(r - 1)];
      next_score = aln_array_[query_.size()][r];
      if(score >= aln_settings_.min_aln_score and score > last_score and score >= next_score){ 
	is_maxima = true;
      }
      else if(next_score > score and r == (reference_.size()-1) and next_score >= aln_settings_.min_aln_score){
	is_maxima = true;
	r = r + 1;
	score = aln_array_[query_.size()][(r - 1)];
      }
      //if point is maxima, keep it
      if(is_maxima){
	maxima_coords coords;
	coords.query_pos = query_.size();
	coords.ref_pos = r - 1; // shift relative to string, aln_array_ has an extra row
	coords.aln_score = score;
	global_maxima.push_back(coords);
      }
    }  
    if(debug_){
      cerr << "total maxima= " << global_maxima.size() << endl;
    }
    return filter_maxima_(global_maxima);
  }  
  
  //////////////////////////////////////////////////////////////////
  // filter maxima to find best alignment of the query at each loci
  //////////////////////////////////////////////////////////////////
  vector< maxima_coords > filter_maxima_(vector< maxima_coords > maxima){
    vector< maxima_coords > filtered_maxima;
    //TODO: look at CLRS book, stock price algo.
    //maximum subarray problem
    //I should also consider to cost of a simple filter and then follow up with
    //tracing alignment and comparing start-end points together
    for(int i = 0; i < maxima.size(); ++i){
      bool keep = true;
      for(int j = 0; j < maxima.size(); ++j){
	if( i != j){
	  int dist = abs( (maxima[i].ref_pos - maxima[j].ref_pos) );
	  //alignments are filtered if they are contained in the same query space.
	  if( dist < ( query_.size() * 1.25 ) ){ 
	    if( maxima[j].aln_score > maxima[i].aln_score ){
	      keep = false;
	    }
	    else if( maxima[j].aln_score == maxima[i].aln_score ){
	      // if two equal scoring alignments exist, choose the one farther along the query sequence
	      if( maxima[j].query_pos > maxima[i].query_pos ){
		keep = false;
	      }
	      //otherwise choose the more concise alignments
	      else if( maxima[j].query_pos == maxima[i].query_pos ){
		if(maxima[j].ref_pos < maxima[i].ref_pos){
		  keep = false;
		}
	      }	      
	    }
	  }
	}
      }
      if(keep){
	filtered_maxima.push_back(maxima[i]);
      }
    }
    return filtered_maxima;
 }

  ////////////////////////////////////////
  // make a cigar from the traceback map
  ////////////////////////////////////////
  string cigar_compress_(string map){
    string cigar;
    int counter = 1;
    for(int cig = 0; cig < map.size(); ++ cig){
      if(map[cig]==map[cig+1]){
        counter = counter + 1;
      }
      else{
        string count, letter;
        count = to_string(counter);
        letter = map[cig];
        cigar = cigar.append(count);
        cigar = cigar.append(letter);
        counter = 1;
      }
    }
    return cigar;
  }
  
  /////////////////////////////////////////////////////
  // traceback the alignment from some given position
  /////////////////////////////////////////////////////
  // if global is true, add arguments for 
  void trace_alignment_(maxima_coords position, alignment_report * aln_report){
    if(debug_){
      cerr << "trace alignment" << endl;
    }
    int ref_pos = position.ref_pos;
    int que_pos = position.query_pos;
    aln_report->reference_end = ref_pos;
    aln_report->query_end = que_pos;
    string cigar_values;
    int trace_end = 0;
    int edit_dist = 0;
    // soft clip end if necessary
    if(que_pos != query_.size()){
      for(int i = 0; i < (query_.size() - que_pos); ++i){
	cigar_values = cigar_values.insert(0, "S");
      }
    }
    if(debug_){
      cerr << "tracing alignment" << endl;
    }
    while(trace_end < 1) {
      // 0 diag, 1 que, 2 ref
      int direction = traceback_matrix_[que_pos][ref_pos];
      // TODO: use a switch/case/default construct here
      /*
	if(debug_){
	cerr << "q= " << que_pos << " r=" << ref_pos << endl;
      }
      */
      if(direction <= 0) {
	que_pos = que_pos - 1;
	ref_pos = ref_pos - 1;
	if(direction == 0){
	  cigar_values = cigar_values.insert(0, "M");
	}
	else{
	  cigar_values = cigar_values.insert(0, "X");
	  ++edit_dist;
	}
      }
      else if(direction == 1){
	cigar_values = cigar_values.insert(0, "I");
	que_pos = que_pos - 1;
	++edit_dist;
      }
      else if(direction == 2){
	cigar_values = cigar_values.insert(0, "D");
	ref_pos = ref_pos - 1;
	++edit_dist;
      }
      //if we reach the end of a sequence we finished the alignment
      if( que_pos == 0 || ref_pos == 0){
	trace_end = 1;
      }
      if(!global_aln_){
	if(aln_array_[que_pos][ref_pos] == 0){
	  trace_end = 1;
	}
      }
      if(trace_end == 1){
	aln_report->reference_start = ref_pos+1;
	aln_report->query_start = que_pos+1;
      }
    }
    //soft clip the beginning if necessary
    if(que_pos != 0){
      for(int i = 0; i < que_pos; ++i){
	cigar_values = cigar_values.insert(0, "S");
      }
    }
    aln_report->cigar_values = cigar_values;
    aln_report->cigar = cigar_compress_(cigar_values);
    aln_report->aln_score = position.aln_score;
    aln_report->edit_distance = edit_dist;
    if(debug_){
      cerr << "printing alignment" << endl;
      show_alignment_(aln_report);
    }
  }
  
  //////////////////////////////////////////////////////
  // de-allocate memory for alignment/traceback arrays
  //////////////////////////////////////////////////////
  void delete_allocations_(){
    for(int i = 0; i < (query_.size()+1); ++i){
      delete aln_array_[i];
      delete traceback_matrix_[i];
    }
    delete aln_array_;
    delete traceback_matrix_;
    delete [] query_maxima_index_;
    delete [] reference_maxima_;
  }

  ////////////////////////////////////////////////////////////////////
  //  a function to print either the alignment or traceback matrix
  ////////////////////////////////////////////////////////////////////
  void show_matrix_(int ** allocated_matrix_){
    //print the query sequence
    cerr << endl;
    cerr << "    -  ";
    for(int i = 0; i < query_.size(); ++i)
      cerr << query_.substr(i, 1) << "  ";
    cerr << endl;
    //print the reference and scores
    string refScore, queScore;
    int val;
    int ref_pos = 0;
    while(ref_pos <= reference_.size()){
      int que_pos = 0;
      if (ref_pos == 0){
	refScore = "-";
      }
      else{
	refScore = reference_.substr((ref_pos-1),1);
      }
      cerr << refScore << "  ";
      while(que_pos <= query_.size()){
	val = allocated_matrix_[que_pos][ref_pos];
	if(val < 10 ){
	  if(val < 0){
	    cerr << val << " ";
	  }
	  else{
	    cerr << " " << val << " ";
	  }
	}
	else{
	  cerr << " " << val;
	}
	++que_pos;
      }
      cerr << endl;
      ++ref_pos;
    }
  }

  //function to print contents of arrays
    void show_array_(int * allocated_array_){
      //print the reference sequence
      cerr << endl;
      for(int i = 0; i < reference_.size(); ++i)
	cerr << reference_.substr(i, 1) << "  ";
      cerr << endl;
      //print the reference and scores
      string refScore, queScore;
      int val;
      for(int i = 0; i < reference_.size(); ++i){
	cerr << allocated_array_[i] << "  ";
      }
      cerr << endl;
    }
  
  void show_alignment_(alignment_report * aln_report){
    if(debug_){
      cerr << "show alignment" << endl;
    }
    string cig = aln_report->cigar_values;
    cerr << "reference name= " << aln_report->reference_name << endl;
    cerr << "reference start= " << aln_report->reference_start << endl;
    cerr << "reference end= " << aln_report->reference_end << endl;
    cerr << "query start= " << aln_report->query_start << endl;
    cerr << "query end= " << aln_report->query_end << endl;
    cerr << "strand= " << aln_report->strand << endl;
    cerr << "cigar= " << aln_report->cigar << endl;
    cerr << "alignment score= " << aln_report->aln_score << endl;

    //look at the barcode
    int que_pos = 0;
    for(int i=0; i < cig.size(); ++i){
      if(cig[i]=='D'){
	cerr << '-';
      }
      else{
	cerr << query_.substr(que_pos, 1);
	++que_pos;
      }
    }
    cerr << endl;

    //look at the mapping
    for(int i=0; i < cig.size(); ++i){
      string map_str;
      if(cig[i]=='M'){map_str = '|';}
      if(cig[i]=='X'){map_str = '*';}
      if(cig[i]=='I'){map_str = ' ';}
      if(cig[i]=='D'){map_str = ' ';}
      if(cig[i]=='S'){map_str = ' ';}
      cerr << map_str;
    }
    cerr << endl;

    //look at reference map
    int ref_pos = aln_report->reference_start - 1;

    if(aln_report->query_start > 1){
      for (int i = (aln_report->query_start - 1); i >= 1; --i) {
	if((ref_pos - i) >= 0){
	  cerr << reference_.substr( (ref_pos - i) ,1);
	}
	else{
	  cerr << " ";
	}
      }
    }
    
    for (int i=0; i < cig.size(); ++i) {
      if(cig[i]=='I'){
        cerr << '-';
      }
      else if(cig[i]=='S'){
	cerr << "";
      }
      else{
        cerr << reference_.substr(ref_pos,1);
        ++ref_pos;
      }
    }

    /*
    if(aln_report->query_end < query_.size()){
      for (int i = aln_report->query_end; i < query_.size(); ++i) {
	if((ref_pos + i) < reference_.size()){
	  cerr << reference_.substr( (ref_pos + i) ,1);
	}
      }
    }
    */
    
    cerr << endl << endl;
  }

  ////////////////////////////////////
  // find the maximum alignment score
  ////////////////////////////////////
  int max_alignment_score_(){
    int maxVal = 0;
    for(int j = 1; j <= reference_.size(); ++j){
      int value = aln_array_[query_.size()][j];
      if(value > maxVal){
	maxVal = value;
      }
    }
    return maxVal;
  }
  
  void sw_alignment_(shared_ptr< MutableAlignment > query,
		     shared_ptr<MutableAlignment> read,
		     string strand, bool global_aln,
		     shared_ptr< vector<alignment_report> > & alignments){
    global_aln_ = global_aln;
    strand_ = strand;
    if(global_aln_){
      search_dist_ = 3;
    }
    else{
      search_dist_ = 4;
    }
    // set some shared variables
    if(debug_){
      cerr << "smith waterman alignment on strand:" << strand << endl;
    }
    string ref_seq = read->get_sequence();
    int ref_len = ref_seq.size();
    reference_ = ref_seq;
    if(strand == "+"){
      query_ = query->get_sequence();
    }
    else if(strand == "-"){
      string rc = reverse_complement_(query->get_sequence()); 
      query_ = rc;
    }
    run_sw_(query, read, alignments, 0, ref_len);
  }

  void sw_alignment_begin_(shared_ptr< MutableAlignment > query,
			    shared_ptr<MutableAlignment> read,
			    string strand, bool global_aln,
			    shared_ptr< vector<alignment_report> > & alignments,
			    int aln_len){
    global_aln_ = global_aln;
    strand_ = strand;
    if(global_aln_){
      search_dist_ = 3;
    }
    else{
      search_dist_ = 4;
    }
    // set some shared variables
    if(debug_){
      cerr << "smith waterman alignment on strand:" << strand << endl;
    }
    string ref_seq = read->get_sequence().substr(0, aln_len);
    int ref_len = ref_seq.size();
    reference_ = ref_seq;
    if(strand == "+"){
      query_ = query->get_sequence();
    }
    else if(strand == "-"){
      string rc = reverse_complement_(query->get_sequence()); 
      query_ = rc;
    }
    run_sw_(query, read, alignments, 0, ref_len);
  }

  void sw_trim_(shared_ptr< MutableAlignment > query,
  		shared_ptr<MutableAlignment> read,
  		string strand, bool global_aln, int trim_len,
  		shared_ptr< vector<alignment_report> > & alignments){
    search_dist_ = 4;
    strand_ = strand;
    if(debug_){
      cerr << "smith waterman alignment on strand:" << strand << endl;
    }
    if(strand == "+"){
      query_ = query->get_sequence();
    }
    else if(strand == "-"){
      string rc = reverse_complement_(query->get_sequence()); 
      query_ = rc;
    }
    string ref_seq = read->get_sequence();
    int ref_len = ref_seq.size();
    if(2*trim_len < ref_seq.size()){
      reference_ = ref_seq.substr(0, trim_len);
      run_sw_(query, read, alignments, 0, ref_len);
      reference_ = ref_seq.substr(ref_len-trim_len, trim_len);
      run_sw_(query, read, alignments, ref_len-trim_len, ref_len);
    }
    else{
      reference_ = ref_seq;
      run_sw_(query, read, alignments, 0, ref_len);
    }
  }

  void run_sw_(shared_ptr< MutableAlignment > query,
  	       shared_ptr<MutableAlignment> read,
  	       shared_ptr< vector<alignment_report> > & alignments,
  	       int adj_pos, int ref_len){
    // perform key alignment steps 
    allocate_memory_();
    score_matrices_();    
    /*
    if(debug_){
      show_scores_();
      show_trace_();
      show_maxima_();
      show_maxima_index_();
    }
    */
    vector< maxima_coords > maxima;
    if(global_aln_){
      maxima = find_global_maxima_();
    }
    else{
      maxima = find_local_maxima_();
    }
    for (auto & max_pos : maxima){
      alignment_report alignment;
      alignment.strand = strand_;
      alignment.aln_score = max_pos.aln_score;
      alignment.query_end = max_pos.query_pos;
      alignment.reference_name = read->get_read_id();
      alignment.query_name = query->get_read_id();
      trace_alignment_(max_pos, &alignment);
      alignment.reference_start = alignment.reference_start + adj_pos;
      alignment.reference_end = alignment.reference_end + adj_pos;
      alignment.query_end = alignment.query_end + adj_pos;
      alignment.reference_length = ref_len;
      alignments->push_back(alignment);      
    }
    delete_allocations_();
  }
  
  void show_scores_(){
    show_matrix_( aln_array_ );
  }
  
  void show_trace_(){
    show_matrix_( traceback_matrix_ );
  }
  void show_maxima_(){
    show_array_( reference_maxima_ );
  }

  void show_maxima_index_(){
    show_array_( query_maxima_index_ );
  }

  ////////////////////////////
  // Public Facing Attributes
  ////////////////////////////
  
  public:

  ///////////////////////
  // Object constructor
  ///////////////////////
  SWAligner(alignment_parameters command_line_input, bool debug){
    aln_settings_ = command_line_input;
    debug_ = debug;
  }

  ////////////////////////////////////////////////////
  // step by step alignment performed by this object
  ////////////////////////////////////////////////////
  void align(shared_ptr< MutableAlignment > query,
	     shared_ptr< MutableAlignment > read,
	     bool global_aln,
	     shared_ptr< vector<alignment_report> >  alignments){
    sw_alignment_(query, read, "+", global_aln, alignments);
    sw_alignment_(query, read, "-", global_aln, alignments);
    //order the alignments by alignment start
    //I should also filter any overlaping +/- alignments (rare cases)
  }
  
  void align_strand(shared_ptr< MutableAlignment > query,
		    shared_ptr< MutableAlignment > read,
		    bool global_aln,
		    shared_ptr< vector<alignment_report> >  alignments,
		    string strand,
		    int aln_len){
    sw_alignment_begin_(query, read, strand, global_aln, alignments, aln_len);
  }

  void trim(shared_ptr< MutableAlignment > query,
	    shared_ptr<MutableAlignment> read,
	    bool global_aln, int trim_len,
	    shared_ptr< vector<alignment_report> > & alignments){
    sw_trim_(query, read, "+", global_aln, trim_len, alignments);
    sw_trim_(query, read, "-", global_aln, trim_len, alignments);
    //order the alignments by alignment start
    //I should also filter any overlaping +/- alignments (rare cases)
  }
  
  /*
  void align_set(vector< shared_ptr< MutableAlignment > > query_set,
		 shared_ptr< MutableAlignment > read,
		 bool global_aln,
		 shared_ptr< vector<alignment_report> >  alignments){
    for(int q=0; q < query_set.size(); ++q){
      align(query_set[q], read, global_aln, alignments);
    }
  }
  */
  
  ///////////////////////////////////////////////////
  // compare two string by computing max align score
  ///////////////////////////////////////////////////
  // int max_aln_score(string query, string reference, string strand){
  //   // this is faster than completing the alignment
  //   if(strand == "+"){
  //     query_ = query;
  //   }
  //   else if(strand == "-"){
  //     query_ = reverse_complement_(query);
  //   }
  //   int max_aln_score;
  //   allocate_memory_();
  //   score_matrices_();
  //   max_aln_score = max_alignment_score_();
  //   delete_allocations_();
  //   return max_aln_score;
  // }

};


#endif
