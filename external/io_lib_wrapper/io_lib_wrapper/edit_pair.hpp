////////////////////////////////////////////////////////////////////////////////
//
// Edit_pair -- a convenience structure to keep track of individual mutations
//
// Author: Darya Filippova, darya.filippova@gmail.com
//
// March 22, 2016
//
// Roche Sequencing Solutions
//
////////////////////////////////////////////////////////////////////////////////

#ifndef __EDIT_PAIR_LIB
#define __EDIT_PAIR_LIB

#include <cassert>
#include <iostream>

// encoding of the edits
#define MISMATCH_A 0 // observed A in the read
#define MISMATCH_C 1
#define MISMATCH_G 2
#define MISMATCH_T 3
#define MISMATCH_N 4
#define INSERTION_A 5 // (was not present in the reference seq)
#define INSERTION_C 6
#define INSERTION_G 7
#define INSERTION_T 8
#define INSERTION_N 9
#define DELETION_A 10 // (was present in the reference, but not found in the read)
#define DELETION_C 11
#define DELETION_G 12
#define DELETION_T 13
#define DELETION_N 14

#define UNKNOWN_BASE 100

using namespace std;

class EditPair {
	
	// see above for the codes
	uint8_t code = 0;

	// mutation position relative to the reference sequence
	// TODO: how does this play with the deletions?
	int edit_pos = -1;

public:
	/*
	Default constructor: sets edit_op and edit_pos
	*/
	EditPair(unsigned char const & edit_op, int const & edit_pos): 
		code(edit_op),
		edit_pos(edit_pos) {};
	
	bool is_mismatch() const 	{ return code >= 0 && code <= 4; }

	bool is_insertion() const { return code >= 5 && code <= 9; }

	bool is_deletion() const	{ return code >= 10 && code <= 14; }

	/*
	 * Returns a mutation code (defined in edit_pair.hpp file).
	 */
	uint8_t get_code() const {
		return code;
	}

	void set_edit_code(uint8_t new_code) {
		assert(new_code >= MISMATCH_A && new_code <= DELETION_N);
		code = new_code;
	}

	int get_position() const {
		return edit_pos;
	}

	void set_edit_position(const int pos) {
		edit_pos = pos;
	}
	
	friend ostream & operator<<(ostream & os, const EditPair & p);

	static uint8_t get_insertion_code(char inserted_base) {
		switch (inserted_base) {
			case 'N': return INSERTION_N; break;
			case 'A': return INSERTION_A; break;
			case 'C': return INSERTION_C; break;
			case 'G': return INSERTION_G; break;
			case 'T': return INSERTION_T; break;
		}
		return UNKNOWN_BASE;
	}

	static uint8_t get_mismatch_code(char inserted_base) {
		switch (inserted_base) {
			case 'N': return MISMATCH_N; break;
			case 'A': return MISMATCH_A; break;
			case 'C': return MISMATCH_C; break;
			case 'G': return MISMATCH_G; break;
			case 'T': return MISMATCH_T; break;
		}
		return UNKNOWN_BASE;
	}

	static uint8_t get_deletion_code() {
		return DELETION_N;
	}
};

/*
 */
inline ostream & operator<<(ostream & os, const EditPair & p) {
    // write obj to stream
    os << (int)p.edit_pos << ": " << (int)p.code;
    return os;
}

#endif // __EDIT_PAIR_LIB
