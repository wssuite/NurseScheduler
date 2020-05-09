FROM legraina/cmake

# Add depot for java
RUN apt-get update && \
    # Install basics
    apt-get install -y \
        build-essential \
        libbz2-dev \
        libblas-dev \
        liblapack-dev \
        libz-dev \
        openjdk-11-jre-headless \
        subversion \
        time && \
    # Install valgrind
    apt-get install -y --force-yes --fix-missing valgrind

# Copy INSTALL.sh
COPY ./INSTALL.sh /

# Install Boost and BCP
ARG BUILD_TYPE=Release
RUN mkdir -p /usr/local && \
    ./INSTALL.sh /usr/local $BUILD_TYPE
