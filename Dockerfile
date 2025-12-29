# --- Stage 1: Build Environment ---
FROM mambaorg/micromamba:1.5-jammy AS builder
USER root

RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    wget \
    && rm -rf /var/lib/apt/lists/*

USER $MAMBA_USER

# Combine all micromamba installs into a single layer to reduce image size
# and avoid duplicate dependency resolution
RUN micromamba install -y -n base -c conda-forge \
    python=3.10 \
    numpy \
    scipy \
    pandas \
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

ENV PATH="/opt/conda/bin: $PATH" \
    AMBERHOME="/opt/conda" \
    LD_LIBRARY_PATH="/opt/conda/lib: $LD_LIBRARY_PATH" \
    PYTHONUNBUFFERED=1

WORKDIR /app

# Create directory and copy files in fewer layers
COPY kcx. lib kcx.frcmod /app/kcx_params/
COPY run_openmm_kcx_v5.py /app/run_openmm_kcx_v5.py

ENTRYPOINT ["python", "/app/run_openmm_kcx_v5.py"]