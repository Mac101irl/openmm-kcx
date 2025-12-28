# --- Stage 1: Build Environment ---
FROM mambaorg/micromamba:1.5-focal AS builder
USER root
RUN apt-get update && apt-get install -y --no-install-recommends build-essential && rm -rf /var/lib/apt/lists/*
USER mamba

# Install A100-compatible OpenMM and AmberTools
RUN micromamba install -y -n base -c conda-forge -c nvidia \
    python=3.10 openmm cuda-version=12.1 cuda-nvcc cuda-toolkit \
    numpy scipy pandas parmed mdtraj ambertools openmmforcefields \
    && micromamba clean -afy

# --- Stage 2: Final Runtime Image ---
FROM nvidia/cuda:12.1.1-runtime-ubuntu22.04

# Fix for the "debconf" and "man-db" errors in your logs
ENV DEBIAN_FRONTEND=noninteractive
RUN mkdir -p /usr/share/man/man1 /usr/share/man/man7

RUN apt-get update && apt-get install -y --no-install-recommends \
    libgomp1 libxrender1 libxext6 libgl1 procps \
    && apt-get clean && rm -rf /var/lib/apt/lists/*

COPY --from=builder /opt/conda /opt/conda

ENV PATH="/opt/conda/bin:$PATH"
ENV AMBERHOME="/opt/conda"
ENV LD_LIBRARY_PATH="/opt/conda/lib:$LD_LIBRARY_PATH"
ENV OPENMM_CUDA_COMPILER="/opt/conda/bin/nvcc"

WORKDIR /app
# IMPORTANT: Ensure this filename matches your python script exactly
COPY run_openmm_kcx_v5.py /app/run_openmm_kcx_v5.py
RUN chmod +x /app/run_openmm_kcx_v5.py

ENTRYPOINT ["python", "/app/run_openmm_kcx_v5.py"]
