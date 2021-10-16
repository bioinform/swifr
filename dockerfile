## get ubuntu image
FROM gcc:5.4

RUN apt-get update

RUN apt-get install -y wget

RUN apt-get install ca-certificates

################
##install io_lib
################
RUN wget --no-check-certificate http://sourceforge.net/projects/staden/files/io_lib/1.14.6/io_lib-1.14.6.tar.gz

RUN tar xvfz io_lib-1.14.6.tar.gz

## RUN cd io_lib-1.14.6 && ./configure CPPFLAGS="-I /usr/local/include/" CC="gcc" && make && make install && cd ../

RUN cd io_lib-1.14.6 && ./configure  && make && make install && cd ../

################
##install tclap
################
RUN wget --no-check-certificate https://sourceforge.net/projects/tclap/files/tclap-1.2.1.tar.gz

RUN tar xvfz tclap-1.2.1.tar.gz

##RUN cd tclap-1.2.1 && ./configure CPPFLAGS="-I /usr/local/include/" CC="gcc" && make && make install && cd ../

RUN cd tclap-1.2.1 && ./configure && make && make install && cd ../

#################
## compile swifr
#################
ADD . /
RUN make swifr
RUN make test
ENV LD_LIBRARY_PATH=/usr/local/lib 

##CMD find / -name libstaden-read.so.11
##CMD swifr -h
CMD ./test













