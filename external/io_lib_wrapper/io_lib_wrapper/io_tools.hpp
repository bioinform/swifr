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

#ifndef __IO_TOOLS_LIB
#define __IO_TOOLS_LIB

#include <string>
#include <fstream>
#include <vector>
#include <array>

using namespace std;
// TODO: define a namespace here for io_tools

// forward declarations
inline bool isOfType(const string & path, const string & expected_type);

inline bool isSamFile(const string & path) {
	return isOfType(path, ".sam");
}

inline bool isBamFile(const string & path) {
	return isOfType(path, ".bam");
}

inline bool isOfType(const string & path, const string & expected_type) {
	auto type_len = expected_type.size();
	auto index = path.rfind(expected_type);
	auto length = path.size();
	return (index != string::npos) && (index == (length - type_len) );
}

/*
 * Returns true if file has size 0
 */
inline bool isEmpty(const string & path) {
	ifstream in(path);
	// TODO: does this leave a dangling file handler?
	return (in.peek() == ifstream::traits_type::eof() );
}

inline bool exists(const string & path) {
	ifstream in(path.c_str());
    if (in.good()) {
        in.close();
        return true;
    } else {
        in.close();
        return false;
    }   
}

/*
 * searches for the maximum value in the vector and returns a pair object where
 * the first value is the index in the vector for the max value, and the second
 * value is the max value itself.
 */
template<typename T>
inline pair<int,T> max_element(const vector<T> & data) {
	assert(data.size() > 0);
	T max_value = data[0];
	int index_of_max = 0;
	for (int i = 0; i < data.size(); i++) {
		if (data[i] > max_value) {
			max_value = data[i];
			index_of_max = i;
		}
	}
	return make_pair(index_of_max, max_value);
}

// TODO: roll it into a single templatized function
template<typename T, unsigned long S>
inline pair<int,T> max_element(const array<T, S> & data) {
	assert(data.size() > 0);
	T max_value = data[0];
	int index_of_max = 0;
	for (int i = 0; i < data.size(); i++) {
		if (data[i] > max_value) {
			max_value = data[i];
			index_of_max = i;
		}
	}
	return make_pair(index_of_max, max_value);
}

#endif