FROM legraina/bcp

# install valgrind

# create a user
RUN useradd -ms /bin/bash poly

# Change user
USER poly

# Copy sources
COPY --chown=poly ./CMakeLists.txt /home/poly/ns/CMakeLists.txt
COPY --chown=poly ./src/ /home/poly/ns/src/

# Set the working directory
WORKDIR /home/poly/ns/

# Compile the nurse scheduler
RUN echo "set(BCPDIROPT /usr/local/Bcp-1.4/build)" > CMakeDefinitionsLists.txt && \
    echo "set(BOOST_DIR /usr/local/include)" >> CMakeDefinitionsLists.txt && \
    mkdir build && \
    cd build && \
    cmake .. && \
    make

# Copy everything
COPY --chown=poly . /home/poly/ns/

# Entrypoint for the nurse scheduler
ENTRYPOINT [ "./docker-entrypoint.sh" ]
CMD [ "-h" ]
