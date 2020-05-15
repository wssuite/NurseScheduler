FROM legraina/bcp

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
ARG CMAKE_BUILD_TYPE=Release
RUN echo "set(BCPDIROPT /usr/local/Bcp-1.4/build)" > CMakeDefinitionsLists.txt && \
    echo "set(BCPDIRDBG /usr/local/Bcp-1.4/build)" >> CMakeDefinitionsLists.txt && \
    echo "set(BOOST_DIR /usr/local/include)" >> CMakeDefinitionsLists.txt && \
    mkdir build && \
    cd build && \
    cmake -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE} .. && \
    make

# Copy everything
COPY --chown=poly . /home/poly/ns/

# Entrypoint for the nurse scheduler
ENTRYPOINT [ "./docker-entrypoint.sh" ]
CMD [ "-h" ]
