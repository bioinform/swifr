////////////////////////////////////////////////////////////////////////////////
//
// Write out alignments to the BAM/SAM file
//
// Author: Darya Filippova, darya.filippova@gmail.com
//
// March 21, 2016
//
// Roche Sequencing Solutions
//
////////////////////////////////////////////////////////////////////////////////

#ifndef __BAM_WRITER_LIB
#define __BAM_WRITER_LIB

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>

#include "io_lib/sam_header.h"

#include "optional_tags.hpp"
#include "alignment.hpp"
#include "mutable_alignment.hpp"

using namespace std;

/* a flag indicating that quality values are not written out */
const int SKIP_QUALITY = std::pow(2,10);

class BamWriter {

	string version_ = "0.1";

	bool stream_out_;

	ofstream out_;

	/*
	 *
	 */
	void write_header_version(ostream & stream, const string & version, const string & sorted) {
		stream << "@HD\tVN:" << version << "\tSO:" << sorted << "\tGO:none" << endl;
	}

public:

	/*
	 * // TODO: support BAM, output to stdout
	 */
	BamWriter(const string & name) : stream_out_(false) {
		cerr << "BamWriter only supports a SAM format for now" << endl;
		out_.open(name);
	}

	/*
	 *
	 */
	BamWriter() : stream_out_(true) {
		// TODO: write to STDOUT
		cerr << "Writing to STDOUT is unsupported right now" << endl;
		stream_out_ = true;
	}

	/*
	 * BamWriter destructor -- flush the streams
	 */
	~BamWriter() {
		flush();
		close();
	}

	/*
	 * Returns true if the object is ready to write: either into the standard 
	 * out or has successfully opened a file stream to write.
	 */
	bool is_open() {
		return stream_out_ ? true : out_.is_open();
	}

	/*
	 * Write SAM header HD tag
	 * TODO: write out sort order field if sort order is known
	 */
	void write_header() {
		if (stream_out_)
			write_header_version(cout, version_, "unknown");
		else
			write_header_version(out_, version_, "unknown");
	}

	/*
	 * Write SAM header out by copying the header of the provided bam file
	 */
	void write_header(const BamReader & bam_file_reader) {
		SAM_hdr * header = bam_file_reader.get_header();
		// iterate over entries in SAM
		auto version_type = sam_hdr_find(header, "HD", NULL, NULL);
		if (stream_out_)
			write_header_version(cout, version_type->tag->str, "unknown");
		else
			write_header_version(out_, version_type->tag->str, "unknown");
		
		// TODO(filippod): write out references
		// this could be hefty if there are a lot of references (as in 16S 
		// databases or metagenomic)
		auto num_ref = header->nref;
		for (auto i = 0; i < num_ref; i++) {
			// ref_seq_handler.setMapping(i, header->ref[i].name);
			if (stream_out_)
				cout << "@SQ\tSN:" << header->ref[i].name << "\tLN:" << header->ref[i].len << endl;
			else
				out_ << "@SQ\tSN:" << header->ref[i].name << "\tLN:" << header->ref[i].len << endl;
		}
		// TODO(filippod): write out command
		// TODO(filippod): write out comment
	}

	/*
	 *
	 */
	void writeAlignment(shared_ptr<AlignmentCore> & alignment) {
		if (stream_out_) {
			cout << alignment->get_read_id() << "\t" <<
					alignment->get_flags() << "\t" << 
					alignment->get_reference_name() << "\t" <<
					alignment->get_alignment_start() << "\t" << endl;
		}
		else {
			// out_ << alignment->get_read_id() << "\t" <<
			// 		alignment->get_flags() << "\t" << 
			// 		alignment->get_reference_name() << "\t" <<
			// 		alignment->get_alignment_start() << "\t" << endl;
			writeAlignment(*alignment);
		}
	}

	/*
	 *
	 */
	void writeAlignment(AlignmentCore & alignment, const uint8_t skip_fields = 0) {
		if (stream_out_) {
			cout << alignment.get_read_id() << "\t" <<
					alignment.get_flags() << "\t" << 
					alignment.get_reference_name() << "\t" <<
					alignment.get_alignment_start() << "\t" << endl;
		}
		else {
			// using namespace optional_tags;
			out_ << alignment.get_read_id() << "\t" <<
					alignment.get_flags() << "\t" <<
					alignment.get_reference_name() << "\t" <<
					alignment.get_alignment_start() << "\t" << 
					(int)alignment.get_mapq() << "\t" <<
					alignment.get_cigar() << "\t" << 
					alignment.get_rnext() << "\t" << 
					alignment.get_pnext() << "\t" << 
					alignment.get_tlen() << "\t" << 
					alignment.get_sequence() << "\t";
			if (skip_fields & SKIP_QUALITY) {
				out_ << "*\t";
			}
			else {
				out_ << alignment.get_qualities() << "\t";
			}
			out_ << alignment.get_optional_tags_string() << endl;
		}
	}

	/*
	 *
	 */
	void flush() {
		// cerr << "force flush" << endl;
		out_.flush();
	}

	/*
	 *
	 */
	void close() {
		// flush all the buffered bytes
		// release file handler
		if (!stream_out_ && out_)
				out_.close();
	}
};

#endif
