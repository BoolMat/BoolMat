FROM gapsystem/gap-docker

MAINTAINER J. D. Mitchell <jdm3@st-andrews.ac.uk>

COPY --chown=1000:1000 ./src $HOME/src
COPY --chown=1000:1000 ./sh $HOME/sh
COPY --chown=1000:1000 Makefile $HOME
COPY --chown=1000:1000 README $HOME

RUN mkdir -p build/output

# Get libsemigroups
RUN cd $HOME/build \
    && git clone --branch enum --depth=1 \
       https://github.com/james-d-mitchell/libsemigroups 

# Get fmt for libsemigroups
RUN cd $HOME/build/libsemigroups/extern \
    && curl -L -O https://github.com/fmtlib/fmt/archive/5.3.0.tar.gz \
    && tar -xzf 5.3.0.tar.gz  \
    && rm -f 5.3.0.tar.gz 

# Get HPCombi for libsemigroups
RUN cd $HOME/build/libsemigroups/extern \
    && git clone --branch bmat8 https://github.com/james-d-mitchell/HPCombi.git \
    && cd $HOME/build/libsemigroups/extern/HPCombi \
    && git checkout 4ba8fef 

# Get bliss for libsemigroups
RUN cd $HOME/build/libsemigroups/extern/ \
    && git clone --branch master https://github.com/james-d-mitchell/bliss \
    && cd $HOME/build/libsemigroups/extern/bliss \
    && git checkout df8f056 \
    && cd $HOME/build/libsemigroups/extern/ \
    && mv bliss bliss-0.73 

# Copy new files into libsemigroups
RUN cp $HOME/src/test-bmat8-enum.cpp $HOME/build/libsemigroups/tests \
    && cp $HOME/src/Makefile.am $HOME/build/libsemigroups 
    
# Build libsemigroups
RUN cd $HOME/build/libsemigroups \
    && ./autogen.sh \
    && ./configure \
    && make test_bmat8_enum 

# Delete, get, and build GAP Digraphs package
RUN rm -rf $HOME/inst/gap-4.10.2/pkg/digraphs-* \
    && cd $HOME/inst/gap-4.10.2/pkg/ \
    && git clone --branch no-malloc-bliss --depth=1 \
       https://github.com/james-d-mitchell/Digraphs.git digraphs \
    && cd $HOME/inst/gap-4.10.2/pkg/digraphs \
    && ./autogen.sh \
    && ./configure \
    && make  

USER gap

WORKDIR $HOME
