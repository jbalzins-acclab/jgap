FROM quay.io/pypa/manylinux_2_28_x86_64 AS manylinux

FROM ubuntu:24.04

ENV DEBIAN_FRONTEND=noninteractive

# Install base dependencies and GCC 15 toolchain
RUN apt-get update && apt-get install -y --no-install-recommends \
    software-properties-common \
    ca-certificates \
    curl \
    git \
    patchelf \
    && add-apt-repository -y ppa:ubuntu-toolchain-r/test \
    && apt-get update \
    && apt-get install -y --no-install-recommends \
    gcc-15 \
    g++-15 \
    cmake \
    ninja-build \
    libhdf5-dev \
    libopenblas-dev \
    && update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-15 100 \
    && update-alternatives --install /usr/bin/g++ g++ /usr/bin/g++-15 100 \
    && rm -rf /var/lib/apt/lists/*

# Copy pre-built Python environments from PyPA manylinux image
COPY --from=manylinux /opt/python /opt/python

# Install auditwheel into one of the Python environments and make it globally accessible
RUN /opt/python/cp312-cp312/bin/pip install --no-cache-dir auditwheel \
    && ln -s /opt/python/cp312-cp312/bin/auditwheel /usr/local/bin/auditwheel

ENV CC=gcc-15
ENV CXX=g++-15
ENV AUDITWHEEL_PLAT=manylinux_2_39_x86_64
ENV LC_ALL=C.UTF-8
ENV LANG=C.UTF-8
