#ifndef KMER_INDEX_HPP
#define KMER_INDEX_HPP

#include "io_lib_wrapper/mutable_alignment.hpp"
#include <numeric>
#include <math.h>
#include <tuple>
#include "unordered_map.hpp"
#include "bloom_filter.hpp"
#include "reverse_complement.hpp"
//#include <hopscotch_map.h>
#include "bloom_filter.hpp"

typedef tuple <string, char> countType;
typedef tuple <shared_ptr< MutableAlignment >,char, set<int>> indexType;

class KmerIndex{

private:
  int kmer_size_;
  int mismatches_;
  float kmer_freq_;
  string dna_string_ = "ACGT";
  vector< shared_ptr< MutableAlignment > > seq_univ_;
  
  ska::flat_hash_map<int, vector< indexType> >  seq_index_;

  map<countType, int> query_sizes_;
  map<countType, shared_ptr< MutableAlignment >> query_seqs_;
  map<char, int> dna_vals_ {{'A', 0}, {'C', 1}, {'G', 2}, {'T', 3}};
    
  void build_index_(){    
    for(int i = 0; i < seq_univ_.size(); ++i){      
      string seq = seq_univ_[i]->get_sequence();
      string seqName = seq_univ_[i]->get_read_id();      
      vector < string > kmers = kmerize_(seq);
      // create keys to count kmers by reference
      tuple<string, char> countForward = make_tuple(seqName, '+');
      tuple<string, char> countReverse = make_tuple(seqName, '-');
      /// count number of kmers in each reference, on each strand
      query_sizes_[countForward] = kmers.size();
      query_seqs_[countForward] = seq_univ_[i];
      query_sizes_[countReverse] = kmers.size();
      query_seqs_[countReverse] = seq_univ_[i];      
      //vector < string > kmers2 = mismatch_kmers_(kmers);
      int count = 0;     
      for(auto &kmer : kmers){
	int kmerInt = kmer_2_int_(kmer);
	string seqName = seq_univ_[i]->get_read_id();
	string kmerRc = reverse_complement(kmer);
	int kmerRcInt = kmer_2_int_(kmerRc);
	// store forward and reverse kmer in index
	set<int> kmerRange;	
	for(int Kidx = count; Kidx < (count + kmer_size_); ++Kidx){ 
	  kmerRange.insert(Kidx);
	}
	tuple <shared_ptr<MutableAlignment>,char, set<int>> forward = make_tuple(seq_univ_[i], '+', kmerRange);
	tuple <shared_ptr<MutableAlignment>,char, set<int>> reverse = make_tuple(seq_univ_[i], '-', kmerRange);
	seq_index_[kmerInt].push_back(forward);
	seq_index_[kmerRcInt].push_back(reverse);
	count++;
      }
    }    
  }

  void add_kmer_to_index(shared_ptr<MutableAlignment> refSeq, string kmer, char strand){
    tuple <shared_ptr< MutableAlignment >,char> mapKey = make_tuple(refSeq, strand);    
  }
    
  int kmer_2_int_(string kmer){
    //check out bit-wise operations
    int val = 0;
    for(int i = 0; i < kmer.size(); ++i){
      int power = kmer.size() - i - 1;
      char base = kmer[i];
      if(power >0){
	val = val + pow(4, power)*dna_vals_[base];
      }
      else{
	val = val +dna_vals_[base];
      }
    }
    return val;
  }

  vector< string > kmerize_(string sequence){
    vector < string > kmers;
    for(int k = 0; k < ( (sequence.size() - kmer_size_) +1); ++k){
      string kmer = sequence.substr(k, kmer_size_);
      // do not duplicate kmer information
      //if(find(kmers.begin(), kmers.end(), kmer) == kmers.end()){
	kmers.push_back(kmer);
      //}
    }
    return kmers;
  }

  vector< string > mismatch_kmers_(vector<string> kmers){
    //add kmers of some edit distance away
    vector < string > kmer_universe = kmers;
    vector < string > mismatch_kmers;
    // each iteration adds a new mismatch level
    int iters = 0;
    while(true){
      ++iters;
      if(iters > mismatches_){
	break;
      }
      // for every kmer in our set, identify all possible 1nt snps
      for(auto &kmer : kmer_universe){
	add_mismatches_(kmer, mismatch_kmers);
      }
      // add mismatched kmers to universe
      for(auto &mismatch_kmer : mismatch_kmers){
	if(find(kmer_universe.begin(), kmer_universe.end(), mismatch_kmer) == kmer_universe.end()){
	  kmer_universe.push_back(mismatch_kmer);
	}
      }
      // if we loop again, 1nt errors will become 2nt errors, and so on
      // this is definitely not an optimal solution
    }
    return kmer_universe;
  }
			   
  
  void add_mismatches_(string ref_kmer, vector<string> & mismatch_kmers){
    //for each kmer position
    string mismatch_kmer = ref_kmer;
    for(int k = 0; k < mismatch_kmer.size(); ++k){
      //refresh kmer
      mismatch_kmer = ref_kmer;
      //add error for each position in the kmer
      for(int n = 0; n < dna_string_.size(); ++n){
	// can't match reference base
	if(mismatch_kmer[k] == dna_string_[n]){
	  continue;
	}
	mismatch_kmer[k] = dna_string_[n];
	// add kmer if it does not already exist
	if(find(mismatch_kmers.begin(), mismatch_kmers.end(), mismatch_kmer) == mismatch_kmers.end()){
	  mismatch_kmers.push_back(mismatch_kmer);
	}
      }
    }
  }
        
public:

  KmerIndex(vector< shared_ptr< MutableAlignment > > sequences,
	    int kmer_size,
	    int mismatches,
	    float kmer_freq){
    kmer_size_ = kmer_size;
    seq_univ_ = sequences;
    mismatches_ = mismatches;
    kmer_freq_ = kmer_freq;
    build_index_();
  }  
  
  vector< indexType > filter_by_kmers(string sequence,
				      bool search_hard){
    vector < string > subject_kmers = kmerize_(sequence);
    
    //map<countType, int> univ;
    map<countType, set<int>> univ;    
    vector< indexType > filt_query;
    
    map<countType, shared_ptr< MutableAlignment >>::iterator seqIter = query_seqs_.begin();

    while(seqIter != query_seqs_.end()){
      univ[seqIter->first] = set<int>(); 
      seqIter++;
    }
    //check index for kmers, returning all query seqs that match the kmers
    //I need to keep track of the kmer that matches, in case of repeating kmers
    int tot = 0;
    for(auto &kmer : subject_kmers){
      //if(!filter_.contains(kmer)){
      //continue;
      //}
      int kmerInt = kmer_2_int_(kmer);
      tot++;
      if(seq_index_.find(kmerInt) != seq_index_.end()){
	for(auto &x : seq_index_[kmerInt]){
	  string id = get<0>(x)->get_read_id();
	  tuple<string, char> toCount = make_tuple(id, get<1>(x));
	  //univ[toCount] = set_union(univ[toCount], get<2>(x));
	  for(auto &idx : get<2>(x)){
	    univ[toCount].insert(idx);
	  }
	  //univ[toCount]++;
	}
      }
    }
    
    // calculate average score by query seq
    vector<int> matches;
    map<countType, set<int>>::iterator iter = univ.begin();
    float total = 0.0;
    while(iter != univ.end()){
      matches.push_back(iter->second.size());
      iter++;
    }

    // find outliers, most likely hits
    total = accumulate( matches.begin(), matches.end(), 0.0);
    float avg = total/matches.size();
    // standard deviation
    double sqSum = std::inner_product(matches.begin(), matches.end(), matches.begin(), 0.0);
    double stdev = sqrt(sqSum / matches.size() - avg * avg);
    float cutoff1 = avg+(2*stdev);
    int maxMatch = *max_element(matches.begin(), matches.end());    
    //maxMatch - 1 std dev unit
    float cutoff2 = maxMatch - (1*stdev);    
    float cutoff = cutoff1;
    if(cutoff2 > cutoff1){
      cutoff = cutoff2;
    }   
    // choose query seqs to align, based off kmer identity
    iter = univ.begin();
    while(iter != univ.end()){
      float total_kmers = query_sizes_[iter->first];      
      if(iter->second.size() > cutoff){
	//cerr << "cutoff " << cutoff << " size " << iter->second.size() << endl;
	tuple <shared_ptr< MutableAlignment >,char, set<int>> seqStrand = make_tuple(query_seqs_[iter->first], get<1>(iter->first), set<int>());
	filt_query.push_back(seqStrand);
      }
      else if((iter->second.size() / total_kmers)  > kmer_freq_){
	//cerr << "fraction " << (iter->second.size() / total_kmers) << endl;
	tuple <shared_ptr< MutableAlignment >,char, set<int>> seqStrand = make_tuple(query_seqs_[iter->first], get<1>(iter->first), set<int>());
	filt_query.push_back(seqStrand);
      }
      iter++;
    }
    
    if(filt_query.size() == 0){
      vector<indexType> allSeqs;
      if(search_hard){	
	for(int i = 0; i < seq_univ_.size(); ++i){      
	  tuple<shared_ptr<MutableAlignment>, char, set<int>> forward = make_tuple(seq_univ_[i], '+', set<int>());
	  tuple<shared_ptr<MutableAlignment>, char, set<int>> reverse = make_tuple(seq_univ_[i], '-', set<int>());
	  allSeqs.push_back(forward);
	  allSeqs.push_back(reverse);
	}	      
	return allSeqs;
      }
    }    
    return filt_query;
  }

  vector< indexType > all_seqs(){
    vector<indexType> allSeqs;
    for(int i = 0; i < seq_univ_.size(); ++i){      
      tuple<shared_ptr<MutableAlignment>, char, set<int>> forward = make_tuple(seq_univ_[i], '+', set<int>());
      tuple<shared_ptr<MutableAlignment>, char, set<int>> reverse = make_tuple(seq_univ_[i], '-', set<int>());
      allSeqs.push_back(forward);
      allSeqs.push_back(reverse);
    }	      
    return allSeqs;
  }
  
};


#endif
