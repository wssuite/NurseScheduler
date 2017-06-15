FROM nercury/cmake-cpp:gcc-5.2

# Create src directory
RUN mkdir -p build: ./containers/workers
WORKDIR /usr/src/inrc

# Install boost
RUN apt-get update
RUN apt-get install -y build-essential
RUN apt-get install -y libboost-all-dev
RUN apt-get install -y subversion
RUN apt-get install -y libbz2-dev libz-dev libblas-dev liblapack-dev

# Install BCP
RUN svn co https://projects.coin-or.org/svn/Bcp/stable/1.4 Bcp-1.4
RUN mkdir -p Bcp-1.4/build
RUN cd Bcp-1.4/build && ../configure --enable-shared=no --enable-static=yes --with-cpx=no && make && make test && make install

# Install CBC
RUN svn checkout https://projects.coin-or.org/svn/Cbc/releases/2.9.8 Cbc-2.9.8
RUN mkdir -p Cbc-2.9.8/build
RUN cd Cbc-2.9.8/build && ../configure --enable-shared=no --enable-static=yes --with-cpx=no && make && make test && make install

# Export path
ENV BCPDIROPT=/usr/src/inrc/Bcp-1.4/build
ENV CBCDIROPT=/usr/src/inrc/Cbc-2.9.8/build
ENV BOOST_DIR=/usr/include/boost

# Install java for the validator
RUN apt-get install -y openjdk-8-jre

# Copy src
COPY ./src/ /usr/src/inrc/src/
COPY ./Makefile /usr/src/inrc/
COPY ./make.bcp /usr/src/inrc/
COPY ./make.coin /usr/src/inrc/
COPY ./validator.jar /usr/src/inrc/
COPY ./validator.sh /usr/src/inrc/

# Compile inrc
RUN make

# Run static inrc
CMD ./bin/staticscheduler --dir datasets/ --instance $INST --weeks $WEEKS --his $HIS --param paramfiles/$PARAM --sol outfiles/$SOL --timeout $TIMEOUT && ./validator.sh ${INST} ${WEEKS} ${HIS} ${SOL}