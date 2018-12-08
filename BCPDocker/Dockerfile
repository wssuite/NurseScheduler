FROM legraina/cmake

# Add depot for java
RUN echo "deb http://http.debian.net/debian jessie-backports main" > /etc/apt/sources.list.d/jessie-backports.list && \
    # Install basics
    apt-get update && \
    apt-get install -y build-essential subversion time libbz2-dev libz-dev libblas-dev liblapack-dev && \
    # Install java for the validator
    apt-get install -t jessie-backports -y openjdk-8-jre-headless

# Copy INSTALL.sh
COPY ./INSTALL.sh /

# Install Boost and BCP
ARG BUILD_TYPE=Release
RUN mkdir -p /usr/local && \
    ./INSTALL.sh /usr/local $BUILD_TYPE
