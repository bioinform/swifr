#ifndef COMPARE_ALN_SCORES_HPP
#define COMPARE_ALN_SCORES_HPP

using namespace std;

bool compare_aln_scores(alignment_report & aln1, alignment_report & aln2){
  return aln1.aln_score > aln2.aln_score;
}

#endif
