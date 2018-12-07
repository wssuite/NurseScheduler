FROM legraina/bcp

# Export path
ENV BCPDIROPT /usr/local/Bcp-1.4/build
ENV BCPDIRDBG /usr/local/Bcp-1.4/build
ENV BOOST_DIR /usr/local/include

# Create src directory
RUN mkdir -p /ns/
WORKDIR /ns/

# Copy src
COPY ./src /ns/src/
COPY ./CMakeLists.txt /ns/

ENTRYPOINT [ "sleep", "10000000000" ]


# Compile nurse schedule
RUN mkdir build && \\
    cd build && \\
    cmake ..

# Copy everything
COPY . /ns/

# Entrypoint for the nurse scheduler
ENTRYPOINT [ "./docker-entrypoint.sh" ]
CMD [ "-i" , "n005w4_0_1-2-3-3" ]