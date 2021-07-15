#ifndef WRITE_ALN_REPORT_HPP
#define WRITE_ALN_REPORT_HPP


void write_aln_report(string out_file, shared_ptr< vector<alignment_report> > hits){
  ofstream output;
  output.open(out_file);
  output << "name\tref_len\taln_start\taln_end\tstrand\tquery_start\tquery_end\tscore\tcigar\tquery_name\n";   
  for (auto & aln : *hits){
    output << aln.reference_name + "\t" + to_string(aln.reference_length) + "\t" + to_string(aln.reference_start) + "\t" + to_string(aln.reference_end) + "\t" + aln.strand + "\t" + to_string(aln.query_start) + "\t" + to_string(aln.query_end) + "\t" + to_string(aln.aln_score) + "\t" + aln.cigar + "\t" + aln.query_name + "\n"; 
  }
  output.close();
}

#endif
