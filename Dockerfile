FROM gcc:5.4

# Create src directory
RUN mkdir -p /ns/
WORKDIR /ns/

# Add depot for java
RUN echo "deb http://http.debian.net/debian jessie-backports main" > /etc/apt/sources.list.d/jessie-backports.list

# Install basics
RUN apt-get update
RUN apt-get install -y build-essential
RUN apt-get install -y subversion time
RUN apt-get install -y libbz2-dev libz-dev libblas-dev liblapack-dev

# Install java for the validator
RUN apt-get install -t jessie-backports -y openjdk-8-jre-headless

# Copy INSTALL.sh
COPY ./INSTALL.sh /ns/

# Install Boost and BCP
RUN ./INSTALL.sh

# Export path
ENV BCPDIROPT /ns/Bcp-1.4/build
ENV BOOST_DIR /ns/boost_1_64_0

# Copy src
COPY ./src/ /ns/src/
COPY ./Makefile /ns/
COPY ./make.* /ns/


# Compile nurse schedule
RUN make

# Copy src
COPY . /ns/

# Entrypoint for dynamic ns
ENTRYPOINT [ "./docker-entrypoint.sh" ]
