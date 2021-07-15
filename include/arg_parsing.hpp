#ifndef ARG_PARSING_HPP
#define ARG_PARSING_HPP

#include "input_parameters.hpp"
#include <tclap/CmdLine.h>
#include <regex>

using namespace std;

/////////////////////////////////////////////
// object for parsing command line arguments
/////////////////////////////////////////////
class ArgParser {
public:
  //command line input
  int argc;
  char ** argv; // I am referencing a reference
  struct input_parameters ip;
  ArgParser(int argc_in, char * argv_in []){
    argc = argc_in;
    argv = argv_in;
  }

  input_parameters parse_arguments() {
    try {
      //cmd with a help message for the bottom of the screen
      TCLAP::CmdLine cmd("", ' ', "0.1");

      //add an argument for bam file
      TCLAP::SwitchArg debugArg("d", "debug", "debug mode: increase verbosity of alignments, including alignment matrices and traceback matrices, only recommended for a small set of reads.", cmd, false);

      //add an argument for output basename
      TCLAP::ValueArg<string> outputArg("o", "output", "specify an output file basename", false, "alignments", "alignments");
      cmd.add( outputArg );
      
      //add an argument for bam file
      TCLAP::SwitchArg verboseArg("v", "verbose", "verbose, print alignment status during alignment", cmd, false);
      
      //run alignment on the ends of the reads
      TCLAP::ValueArg<int> alnPosArg("l", "length", "Subject Alignment Length: Uses the beginning of the read for the alignment. Decreases alignment time.", false, -1, "int");
      cmd.add( alnPosArg );
      
      //add an argument for bam file
      TCLAP::SwitchArg globalArg("g", "global", "Score the adapter sequence end-to-end, allows negative integers in score matrix.", cmd, false);
      
      //run alignment on the ends of the reads
      TCLAP::ValueArg<int> nThreads("p", "processors", "The number of processors to use (default = 1)", false, 1, "int");
      cmd.add( nThreads );
      
      //run alignment on a subset of reads
      TCLAP::ValueArg<int> nreadArg("n", "nreads", "run alignment on 'n' number of reads", false, -1, "int");
      cmd.add( nreadArg );
      
      //add an argument for minScore
      TCLAP::ValueArg<int> scoreArg("s", "score", "minimum alignment score threshold (default 15)", false, 15, "int");
      cmd.add( scoreArg );

      //add an argument for max difference of expected start
      TCLAP::ValueArg<int> reportArg("m", "max_report", "maximum number of alignments that will be reported for a given read", false, 5, "int");
      cmd.add( reportArg );

      //add an argument for max difference of expected start
      TCLAP::ValueArg<float> freqArg("F", "kmer_frequency", "minumim kmer overlap required to seed alignment, if not outliers exist", false, 0.5, "int");
      cmd.add( freqArg );

      //argument for handeling query index searches
      TCLAP::SwitchArg indexArg("c", "complete", "If no query seqs are returned by index, align read against all seqs", cmd, false);

      /*
      //add an argument for kmer, mismatch
      TCLAP::ValueArg<string> kmerArg("k", "kmer_args", "Filter query by matching kmers with Read. By default, all query sequences are aligned to each read. Specifying the two integers < kmer size >, < mismatches> will set the size of kmers and total mismatches allowed (recommend 1)", false, "0,0", "int,int");
      cmd.add( kmerArg );
      */
      
      //add an argument for kmer, mismatch
      TCLAP::ValueArg<int> kmerArg("k", "kmer_args", "Filter query by matching kmers with Read. By default, all query sequences are aligned to each read. Specifying the kmer argument will filter the reference sequences by matching kmers ", false, 0, "int" );
      cmd.add( kmerArg );

      
      //add an argument for fasta file
      TCLAP::ValueArg<string> QueryArg("q", "query", "fasta file with query sequence(s) to be aligned against the reads", true, "missing", "query.fasta");
      cmd.add( QueryArg );
      
      //add an argument for bam file
      TCLAP::ValueArg<string> readArg("f", "fastq", "fastq input file of read sequences. Query sequences are aligned in the forward orientation.", true, "missing", "reads.fastq");
      cmd.add( readArg );
      
      //parse command line
      cmd.parse( argc, argv );

      //kmer index parameter parsing
      int kmer_size;
      //int kmer_mismatches;

      // parsing digits with comma separation
      // regex pattern2("(\\d+),(\\d+)");
      // string kmers = kmerArg.getValue();
      // if(regex_match(kmers, pattern2)){
      // 	size_t pos = kmers.find(",");
      // 	kmer_size = stoi(kmers.substr(0, pos));
      // 	kmer_mismatches = stoi(kmers.substr(pos + 1, kmers.size()));
      // } else{
      // 	cerr << "--kmer_args (-k) must be of the format '(\d+),(\d+)' for example: '12,1'" << endl;
      // 	// incorrect pattern
      // 	exit(1);
      // }
      
      //fill the input parameters
      ip.verbose = verboseArg.getValue();
      ip.kmer_freq = freqArg.getValue();
      ip.query_path = QueryArg.getValue();
      ip.complete_search = indexArg.getValue();
      ip.kmer_size = kmerArg.getValue();
      //ip.kmer_mismatches = kmer_mismatches;
      ip.max_report = reportArg.getValue();
      ip.read_path = readArg.getValue();
      ip.output_basename = outputArg.getValue();
      ip.align_params.min_aln_score = scoreArg.getValue();
      ip.n_reads = nreadArg.getValue();
      ip.n_threads =  nThreads.getValue();
      ip.alignment_report = true;
      ip.debug_mode = debugArg.getValue();
      ip.global_alignment = globalArg.getValue();
      ip.aln_len = alnPosArg.getValue();
      //ip.output_file_type = typeArg.getValue();
      
      return ip;
    }
    catch (TCLAP::ArgException &e)  // catch any exceptions
      { cerr << "error: " << e.error() << " for arg " << e.argId() << endl; }

  }


};

#endif


