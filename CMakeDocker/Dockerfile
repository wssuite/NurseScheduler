FROM gcc:5.4

ENV cmake_version=3.5

# Install wget
RUN apt-get update && \
    apt-get install -y wget && \
    # Install cmake
    mkdir ~/temp && \
    cd ~/temp && \
    wget https://cmake.org/files/v$cmake_version/cmake-$cmake_version.0.tar.gz && \
    tar -xzvf cmake-$cmake_version.0.tar.gz && \
    cd cmake-$cmake_version.0/ && \
    ./bootstrap && \
    make && \
    make install && \
    rm -rf ~/temp