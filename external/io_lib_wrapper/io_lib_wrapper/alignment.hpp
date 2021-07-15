////////////////////////////////////////////////////////////////////////////////
//
// Alignment class -- holds everything we get from the io_tools. Wraps around bam_seq_t
// io_lib structure for convenience.
//
// Author: Darya Filippova, darya.filippova@gmail.com
//
// March 21, 2016
//
// Roche Sequencing Solutions
//
////////////////////////////////////////////////////////////////////////////////

#ifndef __ALIGNMENT_LIB
#define __ALIGNMENT_LIB

#include <vector>
#include <memory>
#include <iostream>
#include <io_lib/scram.h>
#include <io_lib/os.h>

#include "alignment_core.hpp"
#include "edit_pair.hpp"
#include "long_sequence.hpp"
#include "optional_tags.hpp"


/*
 * Alignment class wraps around bam_seq_t structure and provides logic to handle
 * mutations. The underlying alignment structure is immutable.
 */
class Alignment : public AlignmentCore {

	// reference to the ref. sequence
	shared_ptr<LongSequence> reference_ = nullptr;

	// holds the io_lib alignment struct associated w/ this alignment
	bam_seq_t * _read;

	SAM_hdr * header_;

	string cigar_;

	vector<EditPair> md_tag_mutations;

	vector<EditPair> cigar_string_mutations;

	vector<EditPair> merged_mutations;

	shared_ptr<vector<OptionalTag*>> tags_;

	// length for the left soft clipped region
	uint32_t left_soft_clip = 0;

	// length for the left hard clipped region
	unsigned char left_hard_clip = 0;

	// length for the right clipped region
	unsigned char right_soft_clip = 0;

	// length for the right hard clip
	unsigned char right_hard_clip = 0;


	/*
	 */
	void parseOptionalTags() {
		shared_ptr<vector<OptionalTag*>> tags(new vector<OptionalTag*> ());
		// a nice piece of C code to deal w/ io_lib
		char * iter_handle = NULL;
		char key[2];
		char type;
		bam_aux_t value;
		// while iterator is still pointing to a valid bam structure
		while (0 == bam_aux_iter(_read, &iter_handle, key, &type, &value) ) {
			string s_key(key, 2);
			switch (type) {
				case 'A': // printable char
					cerr << "Unsupported tag type: " << type << endl;
					break;
				case 'B': // integer or numeric array
					{
						if (value.B.t == 'i') {
							const int size = value.B.n;
							ArrayTag<int> * ints = new ArrayTag<int>(s_key, value.B.t, size);
							for(int i = 0; i < 4 * size; i += 4) {
								union {
								  uint32_t i;
								    unsigned char c[4];
								} u;
							    u.c[0] = value.B.s[i];
							    u.c[1] = value.B.s[i + 1];
							    u.c[2] = value.B.s[i + 2];
							    u.c[3] = value.B.s[i + 3];
							    ints->add_value(u.i);
							}
							tags->push_back(ints);
						}
						else if (value.B.t == 'C'){
						  const int size = value.B.n;
						  ArrayTag<int> * ints = new ArrayTag<int>(s_key, value.B.t, size);
						  for(int i = 0; i < size; i += 1) {
						    union {
						      uint8_t i; //see samtools tag specs
						      unsigned char c[1];
						    } u;
						    u.c[0] = value.B.s[i];
						    ints->add_value(u.i);
						  }
						  tags->push_back(ints);
						}
						else if (value.B.t == 'f') {
							const int size = value.B.n;
							ArrayTag<float> * floats = new ArrayTag<float>(s_key, value.B.t, size);

							for(int i = 0; i < 4 * size; i += 4) {
								union {
								    float f;
								    unsigned char c[4];
								} u;
							    u.c[0] = value.B.s[i];
							    u.c[1] = value.B.s[i + 1];
							    u.c[2] = value.B.s[i + 2];
							    u.c[3] = value.B.s[i + 3];
							    floats->add_value(u.f);
							}
							tags->push_back(floats);
						}
						else {
							cerr << "Unsupported array type tag: " << value.B.t << endl;
						}
					}
					break;
				case 'f': // float
					{
						tags->push_back( new FloatTag(s_key, value.f) );
					}
					break;
				case 'H': // byte array in HEX format
					cerr << "Unsupported tag type: " << type << endl;
					break;
				case 'i': // signed int
					{
						// IntTag tag(s_key, value.i);
						// tags.push_back(tag);
						tags->push_back( new IntTag(s_key, value.i) );
					}
					break;
				case 'Z': // printable string
					{
						// StringTag tag(s_key, value.s);
						// tags.push_back(tag);
						tags->push_back( new StringTag(s_key, value.s) );
					}
					break;
				default:
					cerr << "Unknown tag type: " << type << endl;
			}
		}
		tags_ = tags;
	}

	/*
	Parse the cigar string from the read structure; add all mutations found in the string
	as edit pairs to cigar_edits vector
	*/
	void parseCigarString() {
		stringstream cigar_stream;
		cigar_string_mutations.clear();
		// vector<EditPair> cigar_mutations;
		// get cigar string length
		auto cigar_len = bam_cigar_len(_read);
		// get binary cigar string encoding
	    auto cigar = bam_cigar(_read);
		if (cigar_len > 0) {
	    	size_t readIdx = 0;
	    	auto transcriptIdx = bam_pos(_read);
	    	// if the read starts before the beginning of the transcript,
		    // only consider the part overlapping the transcript
		    if (transcriptIdx < 0) {
		        readIdx = -transcriptIdx;
		        transcriptIdx = 0;
		    }
	    	size_t uTranscriptIdx = static_cast<size_t>(transcriptIdx);

	    	uint32_t opOffset = 0;
	    	// uint32_t left_clipped_bases = 0;
	    	uint32_t seen_insertions = 0;
	    	uint32_t seen_deletions = 0;
	    	for (uint32_t cigarIdx = 0; cigarIdx < cigar_len; ++cigarIdx) {
	            uint32_t operation_len = cigar[cigarIdx] >> BAM_CIGAR_SHIFT;
	            auto operation = cigar[cigarIdx] & BAM_CIGAR_MASK;
	            vector<EditPair> local_edits = parseSingleCigarOperation(
	            	operation, operation_len, opOffset, 
	            	left_soft_clip, seen_insertions, seen_deletions, cigar_stream);
	            for (const auto & edit : local_edits)
	            	cigar_string_mutations.push_back(edit);
	            // cerr << endl;
	        }
	    }
	    cigar_ = cigar_stream.str();
	    // return cigar_mutations;
	}

	/*
	 *
	 */
	vector<EditPair> parseSingleCigarOperation(	int32_t op, 
									uint32_t opLen, 
									uint32_t & offset_into_read,
									uint32_t & left_clipped_bases,
									uint32_t & seen_insertions,
									uint32_t & seen_deletions,
									stringstream & human_readable_cigar) const {
									// vector<EditPair> & cigar_edits) const {
		vector<EditPair> local_edits;
		// EditPair edit;
	    switch (op) {
	        case BAM_UNKNOWN:
	        	// skip the symbol
	            cerr << "[INFO] Unknown symbol in the cigar string." << endl;
	            break;

	        case BAM_CMATCH:
	        case BAM_CBASE_MATCH:
	           // do nothing
	        	// cerr << "MATCHES " << opLen << " ";
	        	human_readable_cigar << opLen << "M";
	        	offset_into_read += opLen;
	            break;

	        case BAM_CBASE_MISMATCH:
	            // TODO: account for mismatches (Blasr uses this code)
	        	offset_into_read += opLen;
	            break;

	        case BAM_CINS: {
	        		// cerr << "INS ";
	        		human_readable_cigar << opLen << "I";
	        		int reference_coordinate = get_alignment_start() + 
	        			offset_into_read - left_clipped_bases - seen_insertions + seen_deletions - 1;
	        		// cerr << reference_coordinate << " ";
		        	for (auto i = 0; i < opLen; i++) {
						seen_insertions++;
						auto seq = bam_seq(_read);
						// cerr << "offset into read: " << offset_into_read << " ";
						char base_insert = bit2char(bam_seqi(seq, offset_into_read) );

						local_edits.push_back( 
							EditPair(EditPair::get_insertion_code(base_insert), 
									 reference_coordinate)
							);
						offset_into_read++;	
					}
				}
	            break;

	        case BAM_CDEL: {
		        	// cerr << "DEL ";
	        		human_readable_cigar << opLen << "D";
					for (auto i = 0; i < opLen; i++) {
						seen_deletions++;
						int reference_coordinate = get_alignment_start() + 
		        			offset_into_read - left_clipped_bases - seen_insertions + seen_deletions - 1;
						// auto base_insert = bit2char(bam_seqi(seq, offset_into_read + i) );
						// cerr << "(" << offset_into_read << "," << i << "," << left_clipped_bases << ") ";
		        		// cerr << reference_coordinate << " ";

						local_edits.push_back( 
							EditPair(EditPair::get_deletion_code(), 
									 reference_coordinate ) );
						// return EditPair(EditPair::get_deletion_code(), 
									 // offset_into_read);
					}
				}
	            break;

	        case BAM_CREF_SKIP: {
		        	cerr << "[INFO] Did not expect a splicing event" << endl;
	        	}
	        	break;

	        case BAM_CSOFT_CLIP:
	        	// cerr << "SOFT CLIP ";
	        	human_readable_cigar << opLen << "S";
	        	if (offset_into_read > 0) { // S is the last character in the cigar string
	        		// TODO: is not captured
					// right_soft_clip = opLen;
					// do not include soft clips in the list of mutations
					// cigar_edits.push_back( EditPair('R', opLen) );
				}
				else {
					// TODO: not captured
					// left_soft_clip = opLen;
					offset_into_read += opLen;
					// cerr << offset_into_read << " ";
					left_clipped_bases = opLen;
				}
	        	break;

	        case BAM_CHARD_CLIP:
	        	human_readable_cigar << opLen << "H";
	        	if (offset_into_read > 0) {// not on the right side
	        		// TODO: is not captured
					// right_hard_clip = opLen;
					// do not add hard clips to the list of mutations
					// cigar_edits.push_back( EditPair('r', opLen) );
				}
				else {

					// TODO: value not captured
					// left_hard_clip = opLen;
					// cigar_edits.push_back( EditPair('l', opLen) );
				}
	            // offset_into_read += opLen;
	            break;

	        case BAM_CPAD: // ???
	        	cerr << "[INFO] Unexpected symbol in cigar string" << endl;
	            break;
	    }
	    return local_edits;
	}

	
    /*
     * 
     */
	void parseMdTag(const LongSequence & ref) {
		md_tag_mutations.clear();
		auto auxillary_fields = bam_aux(_read);
		auto md = bam_aux_find(_read, "MD");
		int seen_deletions = 0;
		if (md != NULL) {
			// skip 'Z' -- indicator of the printable string
			md++;
			auto mdLen = strlen(md);
			int read_pos = left_soft_clip;
			int offset = 0;
			for (auto i = 0; i < mdLen; i++) {
		        if ( isdigit(md[i]) ) {
		        	offset = offset * 10 + (md[i] - '0');
		        }
		        else {
		        	read_pos += offset;
		        	vector<EditPair> edits = parseSingleMdOperation(md, i, read_pos, seen_deletions, left_soft_clip);
		        	for (auto const & edit : edits)
		        		md_tag_mutations.push_back(edit);
					offset = 0;
		        }
		    }
		}
		// else {
		// 	// check cigar string for mismatches
		// 	md_edits = getMismatches(cigar_edits, ref);
		// }
	}

	
	vector<EditPair> parseSingleMdOperation(const char * md, int & i, int & read_pos, int & seen_deletions, const uint read_start_base) const {
		vector<EditPair> v;
		auto* seq = bam_seq(_read);
		switch (md[i]) {
			case '^':
	    		// this the same as D in CIGAR: this a duplicate record of the same mutation/error
				// "deletion from the reference" -- need to delete letters from the ref to get the read
				// e.g. gap in the read relative to the reference
				// scroll past the caret symbol
				i++;
				// skip all the bases -- they represent deleted sequence
				while ( isalpha(md[i]) ) {
					seen_deletions++;
					i++;
				}
				// HACK: scrol back one -- the for loop counter in the main 
				// function would otherwise eat up a digit and jump to an incorrect
				// position within the read
				i--;
			break;
			case 'A': case 'C': case 'G': case 'T': case 'N': {
				// mismatch between the reference and the read: get the character in the read
				// that mismatched with the reference
				// bam_seqi uses 0-coordinates
				char base_in_read = bit2char(bam_seqi(seq, read_pos) );
				read_pos++;
				// store the edit and the position relative to the start of aligned
				// bases (indexing from 1)
				int reference_coordinate = get_alignment_start() - read_start_base + read_pos + seen_deletions - 1;
				// cerr << "mismatch " << reference_coordinate << endl;
				// v.emplace_back(EditPair::get_mismatch_code(base_in_read), read_pos - read_start_base);
				v.emplace_back(EditPair::get_mismatch_code(base_in_read), reference_coordinate);
			}
			break;
			default:
				cerr << "[ERROR] Unexpected symbol (" << md[i] << ") in the MD string" << endl;
				// skip
    	}
    	return v;
	}

	/*
	merge inserions/deletions with mismatches accounting for offsets that 
	indels create for mismatches
	If relative_to_reference = true, the positions of the mismatched bases would
	be set while ignoring offsets imposed by the deletions and insertions.
	If relative_to_reference = false, the positions of the mismatch bases would
	be set relative to the first aligned base, shifting the mismatched bases
	downstream of the insertions by the number of insertions (i.e. mismatched
	bases would appear relative to the read, not the reference).
	*/
    void mergeCigarAndMdMutations(vector<EditPair> & cigar_mutations, 
    	vector<EditPair> & md_tag_mutations) {
    	merged_mutations.clear();
    	// vector<EditPair> merged_edits;
    	if (cigar_mutations.size() > 0 || md_tag_mutations.size() > 0) {
			if (cigar_mutations.size() > 0) {
				if (md_tag_mutations.size() > 0) {
					std::merge(	cigar_mutations.begin(), cigar_mutations.end(), 
							md_tag_mutations.begin(), md_tag_mutations.end(), 
							std::back_inserter(merged_mutations),
							[](EditPair const & a, EditPair const & b) {
								return a.get_position() < b.get_position();
							});
				}
				else {
					merged_mutations = cigar_mutations;
				}
			}
			else {
				merged_mutations = md_tag_mutations;
			}

			// adjust coordinates when we observe insertions/deletions
			int seen_insertions = 0, seen_deletions = 0;
			auto seq = bam_seq(_read);
			for (auto& p : merged_mutations) {
				// insertion
				if (p.is_insertion() ) seen_insertions++;
				if (p.is_deletion() ) seen_deletions++;
				bool seen_indels =  (seen_insertions > 0) || (seen_deletions > 0);
				if (p.is_mismatch() && seen_indels ) {
					// cerr << "adjusting the base at " << p.get_position() << 
					// 	" ins=" << seen_insertions << 
					// 	" del=" << seen_deletions << " ";
					int position_in_read = p.get_position() - get_alignment_start() + left_soft_clip + seen_insertions - seen_deletions;
					// cerr << "position in read " << position_in_read << endl;
					auto io_lib_base = bam_seqi(seq, position_in_read);
					char base = bit2char( io_lib_base );
					p.set_edit_code( EditPair::get_mismatch_code(base) );
					// if (!relative_to_reference) {
						// p.set_edit_position( p.get_position() + ins );
					// }
				}
			}
		}
	}

public:

	Alignment(bam_seq_t * read, SAM_hdr * header, const shared_ptr<LongSequence> reference_seq) :
		_read(read), 
		header_(header) {
		// TODO(filippod): set start, end for the interval
		if (is_unaligned() )
			cigar_ = "*";
		else
			parseCigarString();
		if (reference_seq != nullptr) {
			parseMdTag(*reference_seq);
			mergeCigarAndMdMutations(cigar_string_mutations, md_tag_mutations);	
		}
		parseOptionalTags();
	}

	Alignment(bam_seq_t * read, SAM_hdr * header) :
		_read(read), 
		header_(header) {
		if (is_unaligned() )
			cigar_ = "*";
		else
			parseCigarString();
		merged_mutations = cigar_string_mutations;
		parseOptionalTags();
	}

	bool has_reference() const {
		return reference_ != nullptr;
	}

	void set_reference(const shared_ptr<LongSequence> & reference_seq) {
		reference_ = reference_seq;
		// update MD mutations and merged
		if (reference_ != nullptr) {
			parseMdTag(*reference_seq);
			mergeCigarAndMdMutations(cigar_string_mutations, md_tag_mutations);	
		}
	}

	/*
	returns true if this is a primary alignment (not secondary and not a supplemental,
	in other words, if bits XXX are set to 0)
	*/
	bool is_primary() const {
		auto flag_value = bam_flag(_read);
		return ( (flag_value & 0x100) == 0) && ( (flag_value & 0x900) == 0);
	}

	/*
	Return true if the read is unaligned (bit XXX is set to 0 in the FLAGS value)
	*/
	bool is_unaligned() const {
		return (bam_flag(_read) & BAM_FUNMAP) > 0;
	}

	bool is_supplementary() const {
		return (bam_flag(_read) & 0x900) > 0;
	}

	/*
	get the reference sequence name to which this read aligned
	*/
	string get_reference_name() const {
		auto reference_index = bam_ref(_read);
		// find this reference index in the header
		if (header_ == NULL) {
			cerr << "[INFO] Headless BAM -- can not read reference name for index " <<
				reference_index << endl;
			return to_string(reference_index);
		}
		auto reference_count = header_->nref;
		for (int i = 0; i < reference_count; i++) {
			if (i == reference_index)
				return string(header_->ref[i].name);
		}
		// should never reach this part, unless BAM is headless
		// cerr << "[INFO] Reference name is not found in the header for index: headless BAM? " <<
			// reference_index << endl;
		return "*";
	}

	/*
	get the offset into the reference (1-based coordinate system) where the 
	first base of the read has aligned
	*/
	uint32_t get_alignment_start() const { 
		// bam_pos returns alignment coordinate in 0-based system, so we need to
		// convert to a 1-based
		return bam_pos(_read) + 1; 
	}

	/*
	 * return a string
	 */
	string get_read_id () const { return string(bam_name(_read)); }

	/*
	 * return original read length (includes unaligned and clipped bases)
	 */
	uint32_t length() const { return bam_seq_len(_read); }

	/*
	 * return the value for the flags (2nd) column of this alignment 
	 */
	uint16_t get_flags() const { return bam_flag(_read); }

	// int lsc() const { return left_soft_clip;}

	// int rsc() const { return right_soft_clip;}

	// int lhc() const { return left_hard_clip;}

	// int rhc() const { return right_hard_clip;}

	uint8_t get_mapq() const {
		return bam_map_qual(_read);
	}

	/*
	 *
	 */
	string get_cigar() const {
		return cigar_;
	}

	string get_rnext() const {return "*"; }// TODO}

	int get_pnext() const {return 0; }// TODO}

	int get_tlen() const {return 0; }// TODO}

	Sequence get_sequence() const {
		// return string(get_sequence());
		auto seq = bam_seq(_read);
		auto seq_len = bam_seq_len(_read);
		vector<uint8_t> char_seq(seq_len, 0);
		for (int i = 0; i < seq_len; i++) {
			char_seq[i] = bit2char(bam_seqi(seq, i));
		}
		return string(char_seq.begin(), char_seq.end());
	}

	string get_qualities() const {
		const char * qualities = bam_qual(_read);
		// map into visible characters space -- need to add 33 to io_lib's values
		int bam_qual_len = _read->len;
		char quals[bam_qual_len + 1];
		quals[bam_qual_len] = '\0';
		for (int i = 0; i < bam_qual_len; i++) {
			if (qualities[i] < 0) {
				// malformed quals or missing
				return "*";
			}
			quals[i] = qualities[i] + 33;
		}
		return string(quals);
	}

	string get_optional_tags_string() const {
		// using namespace optional_tags;
		stringstream s;
		for (int i = 0; i < tags_->size(); i++) {
			const auto tag_pair = (*tags_)[i];
		// for (const auto tag_pair : *tags_) {
			s << *tag_pair;
			if (i < tags_->size() - 1) s << "\t";
		}
		return s.str();
	}

	shared_ptr<vector<OptionalTag*>> get_optional_tags() const {
		return tags_;
		
	}

	/*
	 * TODO: return a shared_ptr
	 */
	vector<EditPair> getCigarMutations() const {
		return cigar_string_mutations;
		// uint left_clipped_bases = 0;
		// return parseCigarString(left_clipped_bases);
	};

	/*
	 * Not intended to be exposed w/o interleaved cigar mutations
	 */
	// vector<EditPair> getMdTagMutations(const LongSequence & reference_seq) const {
		// return parseMdTag(reference_seq);
	// }


	bool has_mutation_at(const genomic_coordinate_t x) const {
		cerr << "implement Alignment::has_mutation_at" << endl;
		// exit(1);
		// throw RuntimeError();
		return false;
	}

	shared_ptr<EditPair> getMutationAt(const genomic_coordinate_t x) const {
		cerr << "implement alignment::getMutationAt" << endl;
		// exit(1);
		return nullptr;
	}

	/*
	 * Return a vector w/ mutations present in the CIGAR string, with mutation 
	 * positions relative to the beinging of the read. Mutation locations are 
	 * computed in 1-based coordinate system. CLipped regions are not counted toward the 
	 * total aligned bases
	 */
	vector<EditPair> getMutations() const {
		return merged_mutations;
	};

	/*
	 * alphabet for remapping io_lib's two-bit encoding for sequence
	 */
	const uint8_t twoBitToChar[4] = {'A','C','G','T'};

	/*
	 * io_lib mapping for A/C/G/T nucleotides 
	 */
	const uint8_t samToTwoBit[16] = {
		0, 
		0/*A*/, 
		1/*C*/, 
		0, 
		2/*G*/, 
		0, 
		0, 
		0, 
		3,/*T*/
	    0, 
	    0, 
	    0, 
	    0, 
	    0, 
	    0, 
	    0};

	/*
	 * Returns a character A/C/G/T that corresponds to the right most two bits in the 
	 * insigned int c
	 */
	uint8_t bit2char(uint8_t c) const {
		return Alignment::twoBitToChar[Alignment::samToTwoBit[c]];
	}

	/*
	 * Operator overloading to enable pretty print for this object
	 */
	friend ostream & operator<<(ostream & os, const Alignment & a);
};

/*
 * Overloading of the stream output operator
 */
inline ostream & operator<<(ostream & os, const Alignment & a) {
	os << "Read: " << (int)a.get_alignment_start() << " [" << a.length() << "bp]" << endl;
	return os;
}

#endif // __ALIGNMENT_LIB
