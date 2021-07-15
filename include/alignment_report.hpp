#ifndef ALIGNMENT_REPORT_HPP
#define ALIGNMENT_REPORT_HPP

using namespace std;

//alignment object check against
struct alignment_report {
  string query_name;
  string reference_name;
  int reference_length;
  int reference_start;
  int reference_end;
  int query_start;
  int query_end;
  string strand;
  string cigar_values;
  string cigar;
  int aln_score;
  int edit_distance;
};


#endif
