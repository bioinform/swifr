#ifndef ALIGN_PRIMERS_HPP
#define ALIGN_PRIMERS_HPP

#include "compare_aln_scores.hpp"


using namespace std;
//typedef tuple <shared_ptr< MutableAlignment >,char> indexType;
typedef tuple <shared_ptr< MutableAlignment >,char, set<int>> indexType;

void align_primers(input_parameters ip,
		   SWAligner aligner,
		   shared_ptr<MutableAlignment> read,
		   vector< indexType > probes,
		   shared_ptr< vector<alignment_report> > alignments,
		   int aln_len){
  for(int i=0; i < probes.size(); ++i){
    shared_ptr<MutableAlignment> probe = get<0>(probes[i]);
    //convert char to string
    char strandC = get<1>(probes[i]);    
    string strand(1, strandC);
    //aligner.align_strand(probe, read, ip.global_alignment,
    //alignments, strand, aln_len);
    aln_len = probe->get_sequence().size();
    aligner.align_strand(read, probe, ip.global_alignment,
			 alignments, strand, aln_len);

  }
  if(alignments->size() > 0){
    sort(alignments->begin(), alignments->end(), compare_aln_scores);
  }
}

#endif
