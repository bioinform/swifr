## swifr
A Smith-Waterman implementation for fast read alignment 

### command-line interface

```
USAGE: 

   ./bin/swifr  -f <reads.fastq> -q <query.fasta> [-k <int>] [-c] [-F
                <int>] [-m <int>] [-s <int>] [-n <int>] [-p <int>] [-g] [-l
                <int>] [-v] [-o <alignments>] [-d] [--] [--version] [-h]


Where: 

   -f <reads.fastq>,  --fastq <reads.fastq>
     (required)  fastq input file of read sequences. Query sequences are
     aligned in the forward orientation.

   -q <query.fasta>,  --query <query.fasta>
     (required)  fasta file with query sequence(s) to be aligned against
     the reads

   -k <int>,  --kmer_args <int>
     Filter query by matching kmers with Read. By default, all query
     sequences are aligned to each read. Specifying the kmer argument will
     filter the reference sequences by matching kmers 

   -c,  --complete
     If no query seqs are returned by index, align read against all seqs

   -F <int>,  --kmer_frequency <int>
     minumim kmer overlap required to seed alignment, if not outliers exist

   -m <int>,  --max_report <int>
     maximum number of alignments that will be reported for a given read

   -s <int>,  --score <int>
     minimum alignment score threshold (default 15)

   -n <int>,  --nreads <int>
     run alignment on 'n' number of reads

   -p <int>,  --processors <int>
     The number of processors to use (default = 1)

   -g,  --global
     Score the adapter sequence end-to-end, allows negative integers in
     score matrix.

   -l <int>,  --length <int>
     Subject Alignment Length: Uses the beginning of the read for the
     alignment. Decreases alignment time.

   -v,  --verbose
     verbose, print alignment status during alignment

   -o <alignments>,  --output <alignments>
     specify an output file basename

   -d,  --debug
     debug mode: increase verbosity of alignments, including alignment
     matrices and traceback matrices, only recommended for a small set of
     reads.

   --,  --ignore_rest
     Ignores the rest of the labeled arguments following this flag.

   --version
     Displays version information and exits.

   -h,  --help
     Displays usage information and exits.




```

### command-line parameters
#### -f, --fastq
Fastq formatted file of read sequences. All query sequences are aligned to all reads, unless kmer indexing is specified.

#### -q, --query
Fasta file contain sequences to be searched for from the reads. Alignments are returned in order by alignment score. All alignments meeting the minimum *--score* threshold will be returned until the *--max_report* parameters is met or until there are no more alignments. 

#### -k, --kmer_args
The kmer size to use for the query index. Swifr produces a Map of kmers of size k across all query sequences. For each read, kmers of size k are produced and looked up in the index. For each read, the kmer coverage against each query sequence in the index is calculated. Query sequences with the highest coverage (Z-score >2) \are kept for the final alignment.Or, in the event that there are multiple similar query sequences, argument -F is used to align query sequences with some minimum coverage 

#### -c, --complete
In the event that the kmer index does not find a matching sequeunce, Align against all query sequences. If reads are noisy (pacbio, nanopore) this option could be useful to increase alginment sensitivity. This option will slow down alignments though, if there are many off target reads relative to the expected query sequences. 

#### -F, --kmer_frequency
After calculating coverage for each query sequence in the kmer-index, any query sequence with at least some fraction of coverage (i.e. 0.3 = 30%) is returned and aligned against using the smith-waterman algorithm. 


#### -s, --score
The minimum score, as calculated by the Smith-Waterman Algorithm
(see [Alignment Strategy](AlignmentStrategy.md),
necessary to qualify a sequence alignment. All alignments satisfying
this score will be used to split/trim the reads.

#### -m, --max_report
Maximum number of alignments that will be reported for a given read. Alignments are returned in order of highest score.

#### -n, --nreads
The total number of reads to run the alignment against. Useful for
debugging and testing parameters.

#### -p, --processors
The number of processors used to perform the alignments.

#### -g, --global
Optimize the alignment for global (end-to-end) alignments. Using this option will allow
negative values to the stored in the scoring matrix. Likewise, alignment maxima are only traced
from the end of the query sequence.  

#### -l, --length
Specify the amount of sequence at the *ends of the reads* to align barcodes against.
By default, adapters are aligned against the entire length of the read. Specifying --length,
for example 100, will align the adapter sequence against the first and last 100 bases of the reads.

#### -d, --debug
Increases verbosity during alignment. Prints highly detailed information about each alignment. Only useful for learning and exploring aspects of the alignments, not for use on large datasets.

