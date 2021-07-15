#ifndef __OPTIONAL_TAGS_LIB
#define __OPTIONAL_TAGS_LIB

// define a namespace to keep these guys scoped
// namespace optional_tags {

	class OptionalTag {
	protected:
		string name;
	public:
		OptionalTag(const string & s): name(s) {}

		virtual void print(ostream & os) const { os << this->name; }

		string get_name() const {
			return name;
		}

		friend ostream & operator<<(ostream & os, const OptionalTag & p);
	};

	inline ostream & operator<<(ostream & os, const OptionalTag & p) {
		p.print(os);
		return os;
	}

	class FloatTag : public OptionalTag {
		float value;
	public:
		FloatTag(const string & n, const float f) : OptionalTag(n), value(f) {}
		// FloatTag(const char * n, const float f) : OptionalTag(n), value(f) {}

		virtual void print(ostream & os) const {
			os << name << ":f:" << value;
		}
	};
	
	class IntTag : public OptionalTag {
		int value;
	public:
		IntTag(const string & n, const int i) : OptionalTag(n), value(i) {}

		virtual void print(ostream & os) const {
			os << name << ":i:" << value; 
		}
	};

	class StringTag : public OptionalTag {
		string value;
	public:
		StringTag(const string & n, const string & s) : OptionalTag(n), value(s) {}

		virtual void print(ostream & os) const {
			os << name << ":Z:" << value; 
		}
	};

	template<typename T>
	class ArrayTag : public OptionalTag {
		vector<T> values;
	public:
		ArrayTag(const string & name, const char type, const size_t size) : 
			OptionalTag(name) {
		}

		void add_value(const T v) {
			values.push_back(v);
		}

		vector<T> get_values() const {
			return values;
		}

		virtual void print(ostream & os) const {
		        //name should include typedef.
		        os << name << ":B:"; //<< T << ":";
			if (std::is_same<T, int>::value)
			        os << "i,";
			else if (std::is_same<T, float>::value)
				os << "f,";
			else
				os << "?,";
			if (values.size() == 0) return;
			for (int i = 0; i < values.size() - 1; i++)
				os << values[i] << ",";
			os << values.back();
		}
	};

// }

#endif // __OPTIONAL_TAGS_LIB
