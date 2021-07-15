#ifndef BASENAME_HPP
#define BASENAME_HPP

string basename( string const& file_path ){
  struct match_path{
    bool operator()(char ch) const{
      return ch == '/';
    }
  };
  return string(
		find_if(
			file_path.rbegin(),
			file_path.rend(),
			match_path()
			).base(),
		file_path.end()
		);
}

#endif
