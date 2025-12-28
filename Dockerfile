# --- Stage 1: Build Environment ---
FROM mambaorg/micromamba:1.5-jammy AS builder
USER root

RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    wget \
    && rm -rf /var/lib/apt/lists/*

USER $MAMBA_USER

# Install packages in a more reliable order
RUN micromamba install -y -n base -c conda-forge \
    python=3.10 \
    numpy \
    scipy \
    pandas \
    && micromamba clean -afy

RUN micromamba install -y -n base -c conda-forge \
    openmm \
    cudatoolkit=11.8 \
    parmed \
    mdtraj \
    ambertools \
    openmmforcefields \
    && micromamba clean -afy


# --- Stage 2: Runtime Image ---
FROM nvidia/cuda:11.8.0-runtime-ubuntu22.04

ENV DEBIAN_FRONTEND=noninteractive
RUN mkdir -p /usr/share/man/man1 /usr/share/man/man7

RUN apt-get update && apt-get install -y --no-install-recommends \
    libgomp1 \
    libxrender1 \
    libxext6 \
    libgl1 \
    libglib2.0-0 \
    procps \
    ca-certificates \
    && apt-get clean && rm -rf /var/lib/apt/lists/*

COPY --from=builder /opt/conda /opt/conda

ENV PATH="/opt/conda/bin: $PATH"
ENV AMBERHOME="/opt/conda"
ENV LD_LIBRARY_PATH="/opt/conda/lib:$LD_LIBRARY_PATH"
ENV PYTHONUNBUFFERED=1

WORKDIR /app
COPY run_openmm_kcx_v5.py /app/run_openmm_kcx_v5.py
RUN chmod +x /app/run_openmm_kcx_v5.py

ENTRYPOINT ["python", "/app/run_openmm_kcx_v5.py"]