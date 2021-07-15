#ifndef __FAI_ENTRY_LIB
#define __FAI_ENTRY_LIB

#include <unordered_map>
#include <string>

using namespace std;

class FaiEntry {

	unordered_map<string,FaiEntry> readFAI(string const & fname) {

	}

	void indexFile(string const & fname, unordered_map<string,FaiEntry> & map) {
	}

public:
	string ref_name;
	int64_t num_bases;
	int64_t byte_offset;
	int bases_per_line;
	int bytes_per_line;

	FaiEntry() {}

	FaiEntry(string s, int64_t n, int64_t b, int bp, int bytes):
		ref_name(s),
		num_bases(n),
		byte_offset(b),
		bases_per_line(bp),
		bytes_per_line(bytes) {}
};
#endif