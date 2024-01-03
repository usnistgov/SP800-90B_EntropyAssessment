FROM --platform=linux/amd64 ubuntu:latest

# Install dependencies
RUN apt-get update && DEBIAN_FRONTEND=noninteractive && apt-get install -y \
    g++ \
    libbz2-dev \
    libdivsufsort-dev \
    libjsoncpp-dev \
    libgmp-dev \
    libmpfr-dev \
    make \
    git \
    libssl-dev \
    gdb \
    && rm -rf /var/lib/apt/lists/*

# # Clone the repository
RUN git clone https://github.com/usnistgov/SP800-90B_EntropyAssessment.git

CMD ["/bin/bash"]