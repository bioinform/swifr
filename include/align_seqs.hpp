#ifndef ALIGN_SEQS_HPP
#define ALIGN_SEQS_HPP

#include "compare_aln_scores.hpp"

using namespace std;

void align_seqs(input_parameters ip,
		SWAligner aligner,
		shared_ptr<MutableAlignment> read,
		vector< shared_ptr< MutableAlignment > > probes,
		shared_ptr< vector<alignment_report> > alignments,
		string strand,
		int aln_len,
		int aln_side){
  for(int i=0; i < probes.size(); ++i){
    aligner.align_strand(probes[i], read, ip.global_alignment,
			 alignments, strand, aln_len);
  }
  if(alignments->size() > 0){
    sort(alignments->begin(), alignments->end(), compare_aln_scores);
  }
}

#endif
