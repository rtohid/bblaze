#  Copyright (c) 2019 R. Tohid
#
#  Distributed under the Boost Software License, Version 1.0. (See
#  accompanying file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

FROM fedora

# Prerequisite packages
RUN dnf install -y sudo
RUN dnf install -y bash-completion
RUN dnf install -y vim-enhanced vim-nerdtree
RUN dnf install -y powerline powerline-fonts
RUN dnf install -y make cmake git bzip2 which zlib-devel
RUN dnf install -y hwloc-devel blas blas-devel lapack lapack-devel boost-devel libatomic
RUN dnf install -y gcc gcc-c++ gcc-gfortran openmpi-devel redhat-rpm-config clang
RUN dnf install -y llvm compiler-rt clang-devel llvm-devel
RUN dnf install -y findutils
RUN dnf install -y gdb wget curl lcov boost-devel hwloc-devel
RUN dnf install -y environment-modules 
RUN dnf install -y hdf5-openmpi-devel
RUN dnf install -y python3-devel python3 python3-pip python3-flake8

# Build arguments
ARG SRC
ARG CORES
ARG PREFIX
ARG MALLOC
ARG BUILD_TYPE

RUN pip3 install --trusted-host pypi.org --trusted-host files.pythonhosted.org numpy pytest
RUN pip3 install numpy pytest
WORKDIR /

# Install HPX
RUN git clone --depth 1 https://github.com/STEllAR-GROUP/hpx.git
RUN mkdir -p /hpx/build
WORKDIR /hpx/build
RUN pwd
RUN cmake \ 
      -DHPX_WITH_TESTS=OFF                                         \
      -DHPX_WITH_EXAMPLES=oFF                                      \
      -DHPX_WITH_MALLOC=${MALLOC}                                  \
      -DHPX_WITH_MAX_CPU_COUNT=80                                  \
      -DCMAKE_BUILD_TYPE=${BUILD_TYPE}                             \
      -DHPX_WITH_MORE_THAN_64_THREADS=ON                           \
      -DCMAKE_INSTALL_PREFIX=${PREFIX}/hpx/${BUILD_TYPE} \
      ..
RUN make -j ${CORES} install

# Install Blaze
RUN git clone https://bitbucket.org/blaze-lib/blaze.git ${SRC}
RUN mkdir -p ${SRC}/build 
WORKDIR /${SRC}/build
RUN cmake \
      -DCMAKE_BUILD_TYPE=${BUILD_TYPE}                             \
      -DBLAZE_SMP_THREADS=HPX                                      \
      -DHPX_DIR=${PREFIX}/hpx/Debug/lib64/cmake/HPX/               \
      -DCMAKE_INSTALL_PREFIX=${PREFIX}/bblaze/${BUILD_TYPE}        \
      ..

RUN make -j ${CORES} install

RUN echo "ALL            ALL = (ALL) NOPASSWD: ALL" >> /etc/sudoers
RUN useradd -m blaze
USER blaze
WORKDIR /home/blaze
CMD ["sleep", "infinity"]
