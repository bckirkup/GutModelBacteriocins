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
RUN apt-get update && apt-get install -y --no-install-recommends \
      build-essential cmake g++ \
      openmpi-bin libopenmpi-dev \
      libhdf5-mpi-dev \
      ca-certificates \
 && rm -rf /var/lib/apt/lists/*

ARG CUDA_ARCHS=75;86;89
WORKDIR /src
COPY . .

RUN CC=gcc CXX=g++ cmake -B build \
      -DGUTIBM_USE_MPI=ON \
      -DGUTIBM_USE_HDF5=ON \
      -DGUTIBM_USE_CUDA=ON \
      -DCMAKE_BUILD_TYPE=Release \
      -DCMAKE_CUDA_ARCHITECTURES="${CUDA_ARCHS}" \
 && cmake --build build -j"$(nproc)" --target gut_ibm

FROM nvidia/cuda:${CUDA_VERSION}-runtime-ubuntu${UBUNTU_VERSION}

ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update && apt-get install -y --no-install-recommends \
      openmpi-bin \
      libhdf5-openmpi-103-1 \
      python3 \
      awscli \
      ca-certificates \
 && rm -rf /var/lib/apt/lists/*

WORKDIR /opt/gutibm
COPY --from=build /src/build/gut_ibm /opt/gutibm/gut_ibm
COPY deploy/aws/entry.sh /opt/gutibm/entry.sh
RUN chmod +x /opt/gutibm/gut_ibm /opt/gutibm/entry.sh

ENV GUTIBM_BINARY=/opt/gutibm/gut_ibm
ENTRYPOINT ["/opt/gutibm/entry.sh"]
