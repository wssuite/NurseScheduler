FROM legraina/bcp

# Export path
ENV BCPDIROPT /usr/local/Bcp-1.4/build
ENV BOOST_DIR /usr/local/include

# Create src directory
RUN mkdir -p /ns/
WORKDIR /ns/

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
