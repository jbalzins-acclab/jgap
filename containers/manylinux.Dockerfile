FROM quay.io/pypa/manylinux_2_28_x86_64

# Install EPEL repository and system development libraries (HDF5, OpenBLAS, Ninja)
RUN dnf install -y epel-release \
    && dnf install -y \
       hdf5-devel \
       openblas-devel \
       ninja-build \
    && dnf clean all

# Install GCC 15 toolchain from conda-forge (built against glibc 2.28 sysroot)
ENV MAMBA_ROOT_PREFIX=/opt/conda
RUN curl -Ls https://micro.mamba.pm/api/micromamba/linux-64/latest | tar -xj -C /usr/local/bin/ --strip-components=1 bin/micromamba \
    && micromamba create -y -r /opt/conda -n gcc15 -c conda-forge gcc=15 gxx=15 \
    && micromamba clean --all -y

# Symlink GCC 15 binaries into /usr/local/bin
RUN ln -sf /opt/conda/envs/gcc15/bin/x86_64-conda-linux-gnu-gcc /usr/local/bin/gcc-15 \
    && ln -sf /opt/conda/envs/gcc15/bin/x86_64-conda-linux-gnu-g++ /usr/local/bin/g++-15 \
    && ln -sf /usr/local/bin/gcc-15 /usr/local/bin/gcc \
    && ln -sf /usr/local/bin/g++-15 /usr/local/bin/g++

ENV PATH="/opt/conda/envs/gcc15/bin:${PATH}"
ENV CC=gcc-15
ENV CXX=g++-15
ENV CMAKE_INCLUDE_PATH="/usr/include"
ENV CMAKE_LIBRARY_PATH="/usr/lib64"
ENV AUDITWHEEL_PLAT=manylinux_2_28_x86_64
ENV LC_ALL=C.UTF-8
ENV LANG=C.UTF-8
