# SWIFR: Smith-Waterman Implementation for Full Read-identifaction.

SWIFR utilizes the Smith Waterman algorithm to align a set of query sequences against ngs reads

### Manual

Read the [Manual](manual/Manual.md) to learn about the command line parameters.
The [Alignment Strategy](manual/AlignmentStrategy.md) is also useful for understanding the alignment method.

### Dependencies

* `io_lib_wrapper` - (provided in external/io_lib_wrapper )
* `io_lib` - available at [SourceForge](http://sourceforge.net/projects/staden/files/io_lib/1.14.6/io_lib-1.14.6.tar.gz/download)
* `kseq.h` - For reading Fasta files
* `tclap` - Command Line Parsing available at [SourceForge](https://sourceforge.net/projects/tclap/files/tclap-1.2.1.tar.gz/download)


##### Installing io_lib

```
wget [link above] && tar -xvzf [io_lib-1.14.6.tar.gz]
./configure
make
make install
```

You may need to provide CCPFLAGS (additional directories to search for header) and CC (which compiler) variables to the configure:

```
./configure CPPFLAGS="-I /usr/local/include/" CC="gcc"
```

If LZMA is missing, you can install XZ utils with `homebrew install xz`.

##### Installing TCLAP

```
wget [link above] && tar -xvzf [tclap-1.2.1.tar.gz]
./configure
make
make install
```
 
#### compile swifr:

```
make
```

### Example of how to run swifr

```
/path/to/bin/swifr -f reads.fastq -q query_seqs.fasta -p 20 -s 20 -m 1 -k 10,1
```

