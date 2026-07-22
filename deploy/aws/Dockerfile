# GutIBM — AWS Batch GPU image (draft; see docs/AWS_BATCH.md)
#
# Build (from repo root):
#   docker build -f deploy/aws/Dockerfile -t gutibm:cuda \
#     --build-arg CUDA_ARCHS=75;86;89 .
#
# Runtime needs an NVIDIA Container Toolkit host (Batch GPU AMI / ECS GPU).

ARG CUDA_VERSION=12.4.1
ARG UBUNTU_VERSION=22.04

FROM nvidia/cuda:${CUDA_VERSION}-devel-ubuntu${UBUNTU_VERSION} AS build

ENV DEBIAN_FRONTEND=noninteractive
# The code is C++20 and uses <format>, which needs GCC 13+. The ubuntu 22.04 base
# ships GCC 11, so pull gcc-13/g++-13 from the toolchain PPA and build with them
# (CUDA 12.4 supports GCC 13 as host compiler). CI passes because its runner
# already ships GCC 13.
RUN apt-get update && apt-get install -y --no-install-recommends \
      build-essential cmake \
      software-properties-common \
 && add-apt-repository -y ppa:ubuntu-toolchain-r/test \
 && apt-get update && apt-get install -y --no-install-recommends \
      ca-certificates \
      g++-13 \
      gcc-13 \
      libhdf5-mpi-dev \
      libopenmpi-dev \
      openmpi-bin \
 && rm -rf /var/lib/apt/lists/*

ARG CUDA_ARCHS=75;86;89
WORKDIR /src
COPY . .

RUN CC=gcc-13 CXX=g++-13 cmake -B build \
      -DGUTIBM_USE_MPI=ON \
      -DGUTIBM_USE_HDF5=ON \
      -DGUTIBM_USE_CUDA=ON \
      -DCMAKE_BUILD_TYPE=Release \
      -DCMAKE_CUDA_HOST_COMPILER=g++-13 \
      -DCMAKE_CUDA_ARCHITECTURES="${CUDA_ARCHS}" \
 && cmake --build build -j"$(nproc)" --target gut_ibm

FROM nvidia/cuda:${CUDA_VERSION}-runtime-ubuntu${UBUNTU_VERSION}

ENV DEBIAN_FRONTEND=noninteractive
# The binary is built with GCC 13, so the runtime needs GCC 13's libstdc++6
# (for the GLIBCXX_3.4.3x symbols pulled in by <format>); the 22.04 base ships
# GCC 11's. Upgrade it from the same toolchain PPA.
RUN apt-get update && apt-get install -y --no-install-recommends \
      software-properties-common \
 && add-apt-repository -y ppa:ubuntu-toolchain-r/test \
 && apt-get update && apt-get install -y --no-install-recommends \
      awscli \
      ca-certificates \
      libhdf5-openmpi-103-1 \
      libstdc++6 \
      openmpi-bin \
      python3 \
 && apt-get purge -y software-properties-common \
 && apt-get autoremove -y \
 && rm -rf /var/lib/apt/lists/* \
 && useradd --create-home --uid 1000 --shell /bin/bash gutibm

WORKDIR /opt/gutibm
COPY --from=build --chown=gutibm:gutibm /src/build/gut_ibm /opt/gutibm/gut_ibm
COPY --chown=gutibm:gutibm deploy/aws/entry.sh /opt/gutibm/entry.sh
# Strip CR in case entry.sh was checked out with Windows CRLF endings: a `\r` in
# the shebang makes the container exit 127 ("/usr/bin/env: 'bash\r': not found").
RUN sed -i 's/\r$//' /opt/gutibm/entry.sh \
 && chmod +x /opt/gutibm/gut_ibm /opt/gutibm/entry.sh

# Non-root runtime (docker:S6471). Work dir defaults to /tmp which is writable.
USER gutibm
ENV GUTIBM_BINARY=/opt/gutibm/gut_ibm
ENTRYPOINT ["/opt/gutibm/entry.sh"]
