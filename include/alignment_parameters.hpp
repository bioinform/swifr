#ifndef ALIGNMENT_PARAMETERS_HPP
#define ALIGNMENT_PARAMETERS_HPP

struct alignment_parameters {
  int match = 1;
  int mismatch = -2;
  int insertion_open = -4;
  int insertion_extend = -2;
  int deletion_open = -2;
  int deletion_extend = -1;
  int min_aln_score = 15;
}; 

#endif
