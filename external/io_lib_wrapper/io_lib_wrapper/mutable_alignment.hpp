/*
 * contains only the start-end information about the read and locations of all mutations in it
 */

#ifndef __MUTATION_VECTOR_LIB
#define __MUTATION_VECTOR_LIB

#include <vector>
#include <string>
#include <sstream>
 #include <memory>
// #include <unordered_set>

#include "io_lib/bam.h"

#include "alignment_core.hpp"
#include "alignment.hpp"
#include "interval.hpp"
#include "edit_pair.hpp"
#include "sequence.hpp"
#include "optional_tags.hpp"
#include "definitions.hpp"

// using namespace optional_tags;

class MutableAlignment : public AlignmentCore {

	string read_id_ = "*";

	uint16_t flags_ = 0;

	string reference_name_ = "*";

	// start/end position -- in the Interval

	// TODO: should this be a set?
	vector<EditPair> mutations_;
	// unordered_set<EditPair> mutations_;

	string cigar_;

	uint8_t mapq_ = 0;

	// TODO(filippod): two-bit encoding for the string (array/vector of ints to 
	// allow variable length)
	Sequence sequence_;

        string qualities_ = "*";

        string read_comment_ = "*";

	shared_ptr<vector<OptionalTag*>> tags_;

public:

	// to suport multiple outputs from BAM
	/*
	 * Does not have to have an alignment coordinate
	 */
	MutableAlignment(const string & read_name,
			 const Sequence & sequence, 
			 const string & qualities = "",
			 const string & read_comment = "") :
		sequence_(sequence),
		qualities_(qualities),
		read_id_(read_name),
		read_comment_(read_comment),
		// read_length(sequence.size() ),
		tags_ (new vector<OptionalTag*>() ),
		cigar_("*") {
		// mark as an unaligned read
		flags_ = 4;
	}

	MutableAlignment(const vector<EditPair> & mutations, 
					const int start_coordinate, 
					const int read_length) : 
					AlignmentCore(start_coordinate, start_coordinate + read_length),
					mutations_(mutations),
					mapq_(0),
					tags_(new vector<OptionalTag*>()) {
		this->left_coordinate_ = start_coordinate;
		this->right_coordinate_ = start_coordinate + read_length - 1;	
	}

	// MutableAlignment(const AlignmentCore & alignment) {
	MutableAlignment(const AlignmentCore & alignment) {
		// TODO(filippod): instantiate all fields based on the alignment object
		read_id_ = alignment.get_read_id();
		flags_ = alignment.get_flags();
		reference_name_ = alignment.get_reference_name();
		left_coordinate_ = alignment.get_alignment_start();
		mapq_ = alignment.get_mapq();

		if (is_unaligned() ) {
			cigar_ = "*";
		}

		// cigar_ = alignment.get_cigar();
		sequence_ = alignment.get_sequence();
		qualities_ = alignment.get_qualities();
		// TODO: copy the tags

		// TODO: alignmentCore does not have mutations
		mutations_ = alignment.getMutations();
		tags_ = alignment.get_optional_tags();
	}

	MutableAlignment(const Alignment & alignment, const shared_ptr<LongSequence> reference) {
		// TODO(filippod): instantiate all fields based on the alignment object
		read_id_ = alignment.get_read_id();
		flags_ = alignment.get_flags();
		reference_name_ = alignment.get_reference_name();
		// cerr << "alignment start: " << alignment.get_alignment_start() << endl;
		// TODO: check if unaligned and set to 0
		if (alignment.is_unaligned() ) {
			left_coordinate_ = 0;
		}
		else
			left_coordinate_ = alignment.get_alignment_start();
		// right_coordinate_ = alignment.get
		mapq_ = alignment.get_mapq();

		if (is_unaligned() ) {
			cigar_ = "*";
		}
		else
			cigar_ = alignment.get_cigar();

		sequence_ = alignment.get_sequence();
		right_coordinate_ = left_coordinate_ + sequence_.size();
		qualities_ = alignment.get_qualities();
		tags_ = alignment.get_optional_tags();
		// if (reference != nullptr) {
			// if (!alignment.has_reference() )
				// alignment.set_reference(reference)
			mutations_ = alignment.getMutations();
		// }
		// else
			// mutations_ = alignment.getMutations();
	}

	/* */
	MutableAlignment(const string & read_name, 
					const uint32_t start_coordinate,
					const uint32_t read_length) : 
					AlignmentCore(start_coordinate, start_coordinate + read_length),
					read_id_(read_name),
					mapq_(0), 
					cigar_("*"), 
					tags_(new vector<OptionalTag*>()) {
	}

	vector<EditPair> getMutations() const {
		return mutations_;
	}

	/*
	 *
	 */
	bool has_mutation_at(const genomic_coordinate_t x) const {
		for (auto & m : mutations_) {
			if (m.get_position() == x) return true;
		}
		return false;
	}

	shared_ptr<EditPair> getMutationAt(const genomic_coordinate_t x) const {
		for (auto & m : mutations_) {
			if (m.get_position() == x) return make_shared<EditPair>(m);
		}
		return nullptr;
	}

	string get_read_id() const {
		return read_id_;
	}

  	void set_read_id(const string & id) {
		read_id_ = id;
	}

        string get_read_comment() const {
		return read_comment_;
	}

        void set_read_comment(const string & comment){
		 read_comment_ = comment;
	}

  	uint16_t get_flags() const { return flags_;}

	void set_flags(const uint16_t f) {
		// cerr << "settign to " << f << " ";
		flags_ = f;
	}

	bool is_primary() const {
		return flags_ & 0x100;
	}

	bool is_unaligned() const {
		return flags_ & BAM_FUNMAP;
	}

    bool is_supplementary() const {
    	return flags_ & 0x900;
    }

	string get_reference_name() const {
		return reference_name_;
	}

	/*
	 *
	 */
	uint32_t get_alignment_start() const {
		return this->left_coordinate_; 
	}

	void set_alignment_start(const uint32_t start) {
		this->left_coordinate_ = start;
	}

	uint8_t get_mapq() const {
		return mapq_;
	}

	string get_cigar() const {
		if (is_unaligned() ) {
			return cigar_;
		}

		// summarize cigar from the mutations
		if (cigar_ == "" || cigar_ == "*") {
			int last_indel_position = get_alignment_start();
			// TODO: account for soft clipping
			stringstream cigar_stream;
			if (mutations_.size() == 0) {
				cigar_stream << length() << "M";
			}
			else {
				// have at least 1 mutation
				EditPair previous = mutations_[0];
				int i = 0, count = 0;
				while (i < mutations_.size() ) {
					EditPair current = mutations_[i];
					if (current.is_deletion() ) {
						/*******************************************************
						 * collapse a run of deletions into a single cigar
						 ******************************************************/
						int matches = current.get_position() - last_indel_position;
						cigar_stream << matches << "M";
						int deletion_run = 1;
						while (current.is_deletion() && i < mutations_.size() ) {
							if (i < mutations_.size() - 1 ) {
								// look ahead: is a deletion && one base over?
								if (mutations_[i+1].is_deletion() && mutations_[i+1].get_position() == current.get_position() + 1) {
									deletion_run++;
									previous = mutations_[i];
									current = mutations_[i];
								}
								else
									break;
							}
							// else -- will quit the loop
							i++;
						}
						cigar_stream << deletion_run << "D";
						last_indel_position = current.get_position();
					}
					else if (current.is_insertion()) {
						/*******************************************************
						 * collapse a run of insertions into a single cigar
						 ******************************************************/
						int matches = current.get_position() - last_indel_position;
						cigar_stream << matches << "M";
						int insertion_run = 1;
						while (current.is_insertion() && i < mutations_.size() ) {
							if (i < mutations_.size() - 1 ) {
								// look ahead: is a deletion && one base over?
								if (mutations_[i+1].is_insertion() && mutations_[i+1].get_position() == current.get_position() ) {
									insertion_run++;
									previous = mutations_[i];
									current = mutations_[i];
								}
								else
									break;
							}
							// else -- will quit the loop
							i++;
						}
						cigar_stream << insertion_run << "I";
						last_indel_position = current.get_position();
					}
					else {
						// not an indel -- a mismatch or soft clip
					}
					previous = current;
					i++;
				}
				int trailing_matches = get_alignment_start() + length() - last_indel_position;
				if ( trailing_matches > 0 ) {
					cigar_stream << trailing_matches << "M";
				}
			}
			return cigar_stream.str();
		}
		return cigar_;
	}

	void set_cigar(const string & cigar) {
		cigar_ = cigar;
	}

	string get_rnext() const { return "*"; }// TODO}

	int get_pnext() const { return 0; }// TODO}

	int get_tlen() const { return 0; }// TODO}

	Sequence get_sequence() const {
		return sequence_;
	}

	void set_sequence(const string & seq) {
		sequence_ = seq;
	}

	virtual uint32_t length() const {
		if ( (flags_ & 4) != 0) {
			return sequence_.size();
		}
		else
			return AlignmentCore::length();
	}

	string get_qualities() const {
		return qualities_;
	}

	void set_qualities(const string & qualities) {
		qualities_ = qualities;
	}

	string get_optional_tags_string() const {
		if (tags_ == nullptr) return "";
		// using namespace optional_tags;
		stringstream s;
		for (int i = 0; i < tags_->size(); i++) {
			const auto tag_pair = (*tags_)[i];
		// for (const auto tag_pair : *tags_) {
			s << *tag_pair;
			if (i < tags_->size() - 1) s << "\t";
		}

		// string tags = s.str();
		// cerr << "All tags: " << tags << endl;
		// return tags;
		return s.str();
	}

	shared_ptr<vector<OptionalTag*>> get_optional_tags() const {
		return tags_;
	}

	void add_tag(const string & tag_name, const string & tag_value) {
		tags_->push_back(new StringTag(tag_name, tag_value) );
	}

	void add_tag(const string & tag_name, const float tag_value) {
		tags_->push_back(new FloatTag(tag_name, tag_value));
	}

	void add_tag(const string & tag_name, const int tag_value) {
		tags_->push_back(new IntTag(tag_name, tag_value));
	}

	template<typename T>
	void add_tag(const string & tag_name, const unsigned char * tag_value, const size_t size) {
		tags_->push_back(new ArrayTag<T>(tag_name, tag_value, size));
	}

	/*
	 * Friend operator overloading to be able to print the sketch in a pretty
	 * ofrmat
	 */
	friend ostream & operator<<(ostream & os, const MutableAlignment & p);
};

/*
 * Pretty print for the alignment sketch object. Prints:
 * "Sketch: <alignment_start> [<read_length>bp]"
 * where read_length is the total sequence length (including clipped regions) 
 * and alignment_start is the first base in (1-based? 0-based?) reference 
 * coordinates where the read sequence aligns to the reference.
 */
inline ostream & operator<<(ostream & os, const MutableAlignment & s) {
	os << "MutAlign [" << s.left_coordinate_ << "," << s.right_coordinate_ << "] (" << s.length() << "bp)";
	return os;
};

typedef shared_ptr<MutableAlignment> MutableAlignmentPtr;



#endif // __MUTATION_VECTOR_LIB
