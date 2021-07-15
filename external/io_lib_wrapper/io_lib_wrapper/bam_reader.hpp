////////////////////////////////////////////////////////////////////////////////
//
// BAM Reader
//
// Author: Darya Filippova, darya.filippova@gmail.com
//
// March 21, 2016
//
// Roche Sequencing Solutions
//
////////////////////////////////////////////////////////////////////////////////

#ifndef __BAM_READER_LIB
#define __BAM_READER_LIB

#include <memory>
#include <io_lib/scram.h>
#include <io_lib/os.h>
#include <algorithm>
// #undef max
// #undef min

#include "alignment.hpp"
#include "io_tools.hpp"

class BamReader {

	scram_fd * _file_handler = NULL;

	SAM_hdr * _header = NULL;

	// maintain a counter of the number of lines read so far
	int _lines;

	void printReferences() {
		if (_header == NULL) return;
		cerr << _header->nref << " references in header" << endl;
		for (int i = 0; i < _header->nref; i++)
			cerr << _header->ref[i].name << endl;
	}

public:

	BamReader(const string & path) : _lines(0) {
		// check whether file exists
		// check if file is empty
		if ( !exists(path) ) {
			cerr << "[ERROR] File " << path << " does not eixsts" << endl; 
		}
		else if ( isEmpty(path) ) {
			cerr << "[ERROR] File " << path << " appears to be empty." << endl;
		}
		else if ( isSamFile(path) ) {
			_file_handler = scram_open(path.c_str(), "r");
		}
		else if ( isBamFile(path) ) {
			_file_handler = scram_open(path.c_str(), "rb");
		}
		else {
			cerr << "[ERROR] Unknown file type: " << path << endl;
			_file_handler = NULL;
		}

		if (_file_handler == NULL) {
			cerr << "[ERROR] Could not open the file for reading: " << path << endl;
		}
		else {
			// opened file successfully
			scram_set_option(_file_handler, CRAM_OPT_NTHREADS, 1);
			// if this is a headless BAM, then scram_get_header will read the 
			// first read instead of header
			_header = scram_get_header(_file_handler);
		}
	}

	~BamReader() {
		close();
	}

	/*
	 *
	 */
	SAM_hdr * get_header() const {
		if (_file_handler != NULL)
			return _header;
		else
			return NULL;
	}

	/*
	 * Returns true if the file is open and ready to be read, false otherwise.
	 */
	bool is_open() const {
		return _file_handler != NULL;
	}

	/*
	 * Returns the number of lines read so far
	 */
	int getLineCount() const {return _lines; }

	/*
	 * Read the next alignment from the file. Returns a NULL if read is not possible (for example, in the
	 * case that reached the end of file)
	 */
	shared_ptr<Alignment> getNextAlignment() {
		if (_file_handler == NULL) {
			cerr << "[INFO] File handler is NULL" << endl;
			return nullptr;
		}
		shared_ptr<Alignment> alignment;
		bam_seq_t * read = reinterpret_cast<bam_seq_t*>( calloc(1, sizeof(bam_seq_t) ) );
		int result = scram_get_seq(_file_handler, &read);
		if (result < 0) {
			cerr << "[INFO] Reached the end of file" << endl;
			return nullptr;
		}
		_lines++;
		alignment = shared_ptr<Alignment>(new Alignment(read, _header) );
		return alignment;
	}

	void close() {
		// if file was opened and associated w/ a file handler -- close it
		if (_file_handler != NULL) {
			scram_close(_file_handler);
			_file_handler = NULL;
		}
	}

};

#endif
