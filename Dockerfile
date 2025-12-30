# --- Stage 1: Build Environment ---
FROM mambaorg/micromamba:1.5-jammy AS builder
USER root

RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    wget \
    && rm -rf /var/lib/apt/lists/*

USER $MAMBA_USER

RUN micromamba install -y -n base -c conda-forge \
    python=3.10 \
    numpy \
    scipy \
    pandas \
    matplotlib \
    openmm \
    cudatoolkit=11.8 \
    parmed \
    mdtraj \
    ambertools \
    openmmforcefields \
    && micromamba clean -afy \
    && find /opt/conda -type f -name '*.a' -delete \
    && find /opt/conda -type f -name '*.pyc' -delete \
    && find /opt/conda -type d -name '__pycache__' -exec rm -rf {} + 2>/dev/null || true \
    && find /opt/conda -type d -name 'tests' -exec rm -rf {} + 2>/dev/null || true \
    && find /opt/conda -type d -name 'test' -exec rm -rf {} + 2>/dev/null || true \
    && rm -rf /opt/conda/pkgs/* \
    && rm -rf /opt/conda/share/doc \
    && rm -rf /opt/conda/share/man \
    && rm -rf /opt/conda/share/gtk-doc \
    && rm -rf /opt/conda/include


# --- Stage 2: Runtime Image ---
FROM nvidia/cuda:11.8.0-runtime-ubuntu22.04

ENV DEBIAN_FRONTEND=noninteractive

RUN mkdir -p /usr/share/man/man1 /usr/share/man/man7 \
    && apt-get update && apt-get install -y --no-install-recommends \
    libgomp1 \
    libxrender1 \
    libxext6 \
    libgl1 \
    libglib2.0-0 \
    procps \
    ca-certificates \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* \
    && rm -rf /var/cache/apt/archives/*

COPY --from=builder /opt/conda /opt/conda

ENV PATH="/opt/conda/bin:$PATH" \
    AMBERHOME="/opt/conda" \
    LD_LIBRARY_PATH="/opt/conda/lib:$LD_LIBRARY_PATH" \
    PYTHONUNBUFFERED=1 \
    OPENMM_DEFAULT_PLATFORM="CUDA" \
    CUDA_PRECISION="mixed" \
    CUDA_CACHE_DISABLE=0 \
    CUDA_CACHE_MAXSIZE=2147483648 \
    CUDA_LAUNCH_BLOCKING=0 \
    CUDA_DEVICE_ORDER="PCI_BUS_ID" \
    CUDA_MPS_PIPE_DIRECTORY=/tmp/nvidia-mps \
    CUDA_MPS_LOG_DIRECTORY=/tmp/nvidia-log \
    NVIDIA_VISIBLE_DEVICES=all \
    NVIDIA_DRIVER_CAPABILITIES=compute,utility

WORKDIR /app

RUN mkdir -p /app/kcx_params /app/inputs /app/out

# Copy KCX parameter files
COPY kcx.lib kcx.frcmod /app/kcx_params/

# Copy the main simulation script
COPY run_openmm_kcx_v5.py /app/run_openmm_kcx_v5.py

RUN chmod +x /app/run_openmm_kcx_v5.py

ENTRYPOINT ["python", "/app/run_openmm_kcx_v5.py"]