# --- Stage 1: Build Environment ---
FROM mambaorg/micromamba:1.5-focal AS builder

USER root
RUN apt-get update && apt-get install -y build-essential && rm -rf /var/lib/apt/lists/*
USER mamba

# Install OpenMM with CUDA support and the NVCC compiler
# Specifying cuda-version=12.1 ensures compatibility with A100/H100 drivers
RUN micromamba install -y -n base -c conda-forge -c nvidia \
    python=3.10 \
    openmm \
    cuda-version=12.1 \
    cuda-nvcc \
    cuda-toolkit \
    cuda-runtime \
    numpy scipy pandas parmed mdtraj \
    ambertools openmmforcefields \
    && micromamba clean -afy

# --- Stage 2: Final Runtime Image ---
# Use NVIDIA's official CUDA runtime for maximum performance
FROM nvidia/cuda:12.1.1-runtime-ubuntu22.04

# Avoid prompts from apt
ENV DEBIAN_FRONTEND=noninteractive

# Install system dependencies
RUN apt-get update && apt-get install -y \
    libgomp1 \
    libxrender1 \
    procps \
    && rm -rf /var/lib/apt/lists/*

# Copy the environment
COPY --from=builder /opt/conda /opt/conda

# Environment variables for CUDA and OpenMM
ENV PATH="/opt/conda/bin:$PATH"
ENV AMBERHOME="/opt/conda"
ENV LD_LIBRARY_PATH="/opt/conda/lib:$LD_LIBRARY_PATH"
# Tell OpenMM where to find the CUDA compiler for A100 optimization
ENV OPENMM_CUDA_COMPILER="/opt/conda/bin/nvcc"

WORKDIR /app

COPY run_openmm_kcx_v5.py /app/run_openmm_kcx_v5.py
RUN chmod +x /app/run_openmm_kcx_v5.py

# Verification: Ensure CUDA is detected during build
RUN python -c "import openmm; p = openmm.Platform.getPlatformByName('CUDA'); print('CUDA Platform found:', p.getName())"

ENTRYPOINT ["python", "/app/run_openmm_kcx_v5.py"]
