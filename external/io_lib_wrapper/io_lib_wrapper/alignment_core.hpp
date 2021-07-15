////////////////////////////////////////////////////////////////////////////////
//
// Abstract base class for the different types of alignments
//
// Author: Darya Filippova, darya.filippova@gmail.com
//
// March 21, 2016
//
// Copyright (2016): Roche Sequencing Solutions
//
////////////////////////////////////////////////////////////////////////////////

#ifndef __ALIGNMENT_CORE_LIB
#define __ALIGNMENT_CORE_LIB

// TODO: list all includes
#include <vector>
#include <memory>

#include "sequence.hpp"
#include "long_sequence.hpp"
#include "interval.hpp"
#include "mutation_collection.hpp"
#include "edit_pair.hpp"
#include "optional_tags.hpp"

using namespace std;

class AlignmentCore : public MutationCollection {

public:

	AlignmentCore() : MutationCollection(0,0) {}

	AlignmentCore(const int start, const int end) : MutationCollection(start, end) {}

	/*
	 * Overload equality operator
	 */
	bool operator ==(const AlignmentCore &b) const {
		return this->get_read_id() == b.get_read_id() &&
			this->get_sequence() == b.get_sequence() &&
			this->get_qualities() == b.get_qualities() &&
			// this->alignment_start == b.alignment_start &&
			this->get_flags() == b.get_flags();
	}

	virtual string get_read_id() const =0;

	virtual uint16_t get_flags() const =0;

	virtual string get_reference_name() const =0;

	virtual uint32_t get_alignment_start() const =0;

	virtual uint8_t get_mapq() const =0;

	virtual string get_cigar() const =0;

	virtual string get_rnext() const =0;

	virtual int get_pnext() const =0;

	virtual int get_tlen() const =0;

	// call Interval's length
	// uint32_t length() const;

	virtual bool is_primary() const =0;

	virtual bool is_unaligned() const =0;

	virtual bool is_supplementary() const =0;

	// int get_rnext() =0;

	// int get_pnext() =0;

	// int get_tlen() =0;

	/*
	 * Return the whole read sequence, including the clipped regions
	 */
	virtual Sequence get_sequence() const =0;

	virtual string get_qualities() const =0;

	// virtual shared_ptr<string> get_qualities() const =0;

	virtual shared_ptr<vector<OptionalTag*>> get_optional_tags() const =0;

	virtual string get_optional_tags_string() const =0;

	// vector<EditPair> getCigarMutations() =0;

	// vector<EditPair> getMdTagMutations(const LongSequence & reference_seq) =0;

	// vector<EditPair> getAllMutations(const LongSequence & reference_seq) =0;

	

};

#endif // __ALIGNMENT_CORE_LIB