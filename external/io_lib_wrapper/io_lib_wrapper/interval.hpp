#ifndef __INTERVAL_LIB
#define __INTERVAL_LIB

#undef max
#undef min

#include <iostream>
// #include <fstream>

#include "definitions.hpp"

using namespace std;

class Interval {

protected:

	genomic_coordinate_t left_coordinate_;

	genomic_coordinate_t right_coordinate_;

public:

	/*
	 *
	 */
	Interval(genomic_coordinate_t l, genomic_coordinate_t r) : 
		left_coordinate_(l), 
		right_coordinate_(r) {}

	/*
	 * Return the length (in bases) of this genomic interval.
	 * If interval starts at 10, end at 11, then its length is 2bp.
	 * If either of the interval ends at at UINT_MAX, the interval is considered
	 * empty and the function returns 0.
	 */
	virtual uint32_t length() const {
		// check for bounds -- may be an empty interval
		if (left_coordinate_ == UINT_MAX || right_coordinate_ == UINT_MAX) {
			return 0;
		}
		return right_coordinate_ - left_coordinate_ + 1;
	}

	/*
	 * Return true is intervals overlap and false otherwise. Implementation is
	 * answers a "not non-overlapping" question.
     * Two intervals do not overlap if only one of the following occurs:
     *    A       B    - A before B
     * --------  ---
     *
     *  B     A        - B before A
     * ---  --------
     * Therefore, we can return :
     *     not ( (self.end < B.start) or (B.end < self.start) )
	 */
	bool is_overlapping(const Interval & interval) const {
		return !no_overlap(interval);
	}

	/*
	 * Return true if intervals do not overlap and false otherwise.
	 */
	bool no_overlap(const Interval & interval) const {
		return 	(this->right_coordinate_ < interval.left_coordinate_) ||
				(interval.right_coordinate_ < this->left_coordinate_);
	}

	/*
	 *
	 */
	bool is_overlapping(const genomic_coordinate_t position) const {
		return position >= left_coordinate_ && position <= right_coordinate_;
	}

	/*
	 *
	 */
	uint overlap_length(const Interval & interval) const {
		if (no_overlap(interval)) 
			return 0;
		else {
			auto end_coordinate 	= std::min(right_coordinate_, interval.right_coordinate_);
			auto start_coordinate	= std::max(left_coordinate_, interval.left_coordinate_);
    		return end_coordinate - start_coordinate + 1;
		}
	}

	/*
	 *
	 */
	genomic_coordinate_t get_start() const {
		return this->left_coordinate_;
	}

	genomic_coordinate_t get_end() const {
		return this->right_coordinate_;
	}

	friend std::ostream & operator<<(std::ostream & os, const Interval & p);
};

inline std::ostream & operator<<(std::ostream & os, const Interval & s) {
	os << "gInt [" << s.left_coordinate_ << "," << s.right_coordinate_ << "] (" << s.length() << "bp)";
	return os;
}

#endif
