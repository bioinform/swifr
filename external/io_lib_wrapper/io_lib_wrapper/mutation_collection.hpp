////////////////////////////////////////////////////////////////////////////////
//
// Abstract base class for objects w/ mutations: tuples of (position,code)
//
// Author: Darya Filippova, darya.filippova@gmail.com
//
// March 21, 2016
//
// Copyright (2016): Roche Sequencing Solutions
//
////////////////////////////////////////////////////////////////////////////////

#ifndef __MUTATION_COLLECTION_LIB
#define __MUTATION_COLLECTION_LIB

#include <vector>
#include <memory>

#include "definitions.hpp"
#include "long_sequence.hpp"
#include "edit_pair.hpp"

class MutationCollection : public Interval {

public:

	MutationCollection() : Interval(0,0) {}

	MutationCollection(const genomic_coordinate_t start, const genomic_coordinate_t end) : Interval(start, end) {}

	virtual bool has_mutation_at(const genomic_coordinate_t x) const =0;

	virtual shared_ptr<EditPair> getMutationAt(const genomic_coordinate_t x) const =0;

	/*
	 * Return a vector of EditPairs w/ positions relative to the reference.
	 */
	virtual vector<EditPair> getMutations() const =0;

	/*
	 * Return a vector of EditPairs w/ positions relative to the reference
	 */
	// virtual vector<EditPair> getMutations(const LongSequence & reference) const =0;
};

#endif // __MUTATION_COLLECTION_LIB