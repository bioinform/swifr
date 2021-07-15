////////////////////////////////////////////////////////////////////////////////
//
// LongSeqeunce --class to parse and whole reference sequences (long strings of 
// ACGT). If you reference is not a long string or is not ACGT, consider using
// other classes.
//
//
// Author: Darya Filippova, darya.filippova@gmail.com
//
// March 21, 2016
//
// Roche Sequencing Solutions
//
////////////////////////////////////////////////////////////////////////////////

#ifndef __REFERENCE_LIB
#define __REFERENCE_LIB

#include <unordered_map>
#include <string>
#include <iostream>
#include <sstream>

#include "fasta_reader.h"
#include "fai_entry.hpp"

using namespace std;

class LongSequence {

	unordered_map<string,FaiEntry> fai_index;

	unordered_map<string,size_t> ref_offsets;

	string reference_;

	string path_;

public:

	LongSequence() {
	}

	LongSequence(const string & path) : path_(path) {
		// TODO: parse input
		// create fai index file if needed (or bamtools index format -- faster indexing?)
		// doesn't matter for viruses since they are small anyway
	}

	bool read_sequence() {
		if (path_.size() == 0) return false;
		// cerr << "TODO: read the fasta" << endl;
		FastaReader reader(path_.c_str());
		if (!reader.is_open()) return false;
		stringstream ss;
		kseq_t * record;
		while ( (record = reader.nextSequence() ) != NULL ) {
			ss << string(record->seq.s, record->seq.l);
		}
		reference_ = ss.str();
		return true;
	}

	/* returns the reference sequence starting at position i (0-based coordinate system)
	and of length L*/
	string get_sequence(uint i, int L) const {
		// TODO: return a string
		return string('A', L);
	}
};

#endif