FROM legraina/bcp

# create a user
RUN useradd -ms /bin/bash poly

# Change user
USER poly

# Copy everything
COPY --chown=poly . /home/poly/ns/

# Set the working directory
WORKDIR /home/poly/ns/

# Compile the nurse scheduler
RUN echo "set(BCPDIROPT /usr/local/Bcp-1.4/build)" > CMakeDefinitionsLists.txt && \
    echo "set(BOOST_DIR /usr/local/include)" >> CMakeDefinitionsLists.txt && \
    mkdir build && \
    cd build && \
    cmake .. && \
    make

# Entrypoint for the nurse scheduler
ENTRYPOINT [ "./docker-entrypoint.sh" ]
CMD [ "-h" ]
