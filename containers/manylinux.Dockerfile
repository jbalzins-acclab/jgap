FROM quay.io/pypa/manylinux_2_28_x86_64 AS manylinux

FROM ubuntu:24.04

ENV DEBIAN_FRONTEND=noninteractive

# Install base dependencies, GCC 15 toolchain, and libcrypt1 required by manylinux python binaries
RUN apt-get update && apt-get install -y --no-install-recommends \
    software-properties-common \
    ca-certificates \
    curl \
    git \
    patchelf \
    libcrypt1 \
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

# Copy pre-built Python environments and internal support tools (including auditwheel) from PyPA manylinux
COPY --from=manylinux /opt/_internal /opt/_internal
COPY --from=manylinux /opt/python /opt/python
COPY --from=manylinux /usr/local/bin/auditwheel /usr/local/bin/auditwheel

# Copy runtime shared libraries required by manylinux Python dynamic extensions (OpenSSL 1.1, libffi, libbz2)
COPY --from=manylinux /usr/lib64/libssl.so.1.1* /usr/local/lib/
COPY --from=manylinux /usr/lib64/libcrypto.so.1.1* /usr/local/lib/
COPY --from=manylinux /usr/lib64/libffi.so* /usr/local/lib/
COPY --from=manylinux /usr/lib64/libbz2.so* /usr/local/lib/

# Setup CA certificate paths and dynamic linker cache for OpenSSL 1.1 runtime
RUN mkdir -p /etc/pki/tls/certs /etc/pki/ca-trust/extracted/pem \
    && ln -sf /etc/ssl/certs/ca-certificates.crt /etc/pki/tls/certs/ca-bundle.crt \
    && ln -sf /etc/ssl/certs/ca-certificates.crt /etc/pki/tls/cert.pem \
    && ln -sf /etc/ssl/certs/ca-certificates.crt /etc/pki/ca-trust/extracted/pem/tls-ca-bundle.pem \
    && ldconfig

ENV CC=gcc-15
ENV CXX=g++-15
ENV AUDITWHEEL_PLAT=manylinux_2_39_x86_64
ENV SSL_CERT_FILE=/etc/ssl/certs/ca-certificates.crt
ENV LC_ALL=C.UTF-8
ENV LANG=C.UTF-8
