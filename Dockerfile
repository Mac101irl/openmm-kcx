# --- Stage 1: Build Environment ---
FROM mambaorg/micromamba:1.5-focal AS builder

# Switch to root to install system build tools
USER root
RUN apt-get update && apt-get install -y build-essential && rm -rf /var/lib/apt/lists/*
USER mamba

# Install Python and scientific packages
RUN micromamba install -y -n base -c conda-forge \
    python=3.10 \
    numpy \
    scipy \
    pandas \
    openmm \
    parmed \
    mdtraj \
    ambertools \
    openmmforcefields \
    && micromamba clean -afy

# --- Stage 2: Final Runtime Image ---
# Use Ubuntu 20.04 (Focal) to match the builder image and prevent GLIBC errors
FROM ubuntu:20.04

# Avoid prompts from apt
ENV DEBIAN_FRONTEND=noninteractive

# Install system runtime dependencies required by OpenMM and AmberTools
# libgomp1 is critical for OpenMP support in AmberTools
RUN apt-get update && apt-get install -y \
    libgomp1 \
    libxrender1 \
    libxext6 \
    libgl1 \
    procps \
    && rm -rf /var/lib/apt/lists/*

# Copy the Conda environment from the builder
COPY --from=builder /opt/conda /opt/conda

# Set environment variables
ENV PATH="/opt/conda/bin:$PATH"
ENV AMBERHOME="/opt/conda"

WORKDIR /app

# Copy the Python script (ensure filename matches exactly)
COPY run_openmm_kcx_v5.py /app/run_openmm_kcx_v5.py
RUN chmod +x /app/run_openmm_kcx_v5.py

# Verify OpenMM installation
RUN python -c "import openmm; print('OpenMM version:', openmm.__version__)"

# Set the entrypoint
ENTRYPOINT ["python", "/app/run_openmm_kcx_v5.py"]