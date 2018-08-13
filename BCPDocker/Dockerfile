FROM gcc:5.4

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
COPY ./INSTALL.sh /

# Install Boost and BCP
ARG BUILD_TYPE=Debug
RUN mkdir -p /usr/local
RUN ./INSTALL.sh /usr/local $BUILD_TYPE
