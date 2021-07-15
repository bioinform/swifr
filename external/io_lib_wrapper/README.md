# io_lib_wrapper

This is a header-only wrapper around `io_lib` library that will help you interact with `io_lib` without having to understand their complex API.

## Dependencies

#### Staden io_lib

`staden_io_lib` -- download the library from [SourceForge](http://sourceforge.net/projects/staden/files/io_lib/1.14.6/io_lib-1.14.6.tar.gz), unpack and follow the installation instructions:

```
wget [link above] && tar -xvzf [io_lib-1.14.6.tar.gz]
./configure CPPFLAGS="-I /usr/local/include/" CC="gcc"
make
make install
```

You may need to provide CCPFLAGS (additional directories to search for header) and CC (which compiler) variables to the configure (as shown above).

#### LZMA

Install LZMA through XZ utils:

```
homebrew install xz
```

#### TCLAP

C++ templatized header-only command line [parser](http://tclap.sourceforge.net/)

## How to use

#### Including in your projects

This library is header-only, so you can simply place it next to your include folder
and you are ready to go:

  `src/` -- your main and other *.cpp files

  `include/` -- your headers

  `io_lib_wrapper/` -- io_lib_wrapper
  
#### Reading a BAM file

```
  shared_ptr<Alignment> alignment = bamReader.getNextAlignment();
	while ( alignment != nullptr ) {
		alignment = bamReader.getNextAlignment();
	}
```
  
You will find more runnable examples in `examples/`.
