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

# === GPU PERFORMANCE OPTIMIZATIONS ===

# CUDA environment variables for optimal GPU performance
ENV PATH="/opt/conda/bin:$PATH" \
    AMBERHOME="/opt/conda" \
    LD_LIBRARY_PATH="/opt/conda/lib: $LD_LIBRARY_PATH" \
    PYTHONUNBUFFERED=1 \
    # Force OpenMM to use CUDA platform
    OPENMM_DEFAULT_PLATFORM="CUDA" \
    # CUDA performance tuning
    CUDA_CACHE_DISABLE=0 \
    CUDA_CACHE_MAXSIZE=2147483648 \
    # Optimize CUDA kernel caching
    CUDA_FORCE_PTX_JIT=0 \
    # Disable CUDA debugging for performance
    CUDA_LAUNCH_BLOCKING=0 \
    # Enable TensorFloat-32 for supported GPUs (Ampere+)
    NVIDIA_TF32_OVERRIDE=1 \
    # Memory management - use default allocator for better performance
    PYTORCH_CUDA_ALLOC_CONF="" \
    # GPU device scheduling - yield for better multi-process performance
    CUDA_DEVICE_ORDER="PCI_BUS_ID"

# Set NVIDIA container runtime environment variables
ENV NVIDIA_VISIBLE_DEVICES=all \
    NVIDIA_DRIVER_CAPABILITIES=compute,utility

WORKDIR /app

# Copy parameter files
COPY kcx.lib kcx.frcmod /app/kcx_params/

# Copy main script
COPY run_openmm_kcx_v5.py /app/run_openmm_kcx_v5.py

# Health check to verify GPU is accessible
HEALTHCHECK --interval=30s --timeout=10s --start-period=5s --retries=3 \
    CMD python -c "import openmm; print(openmm. Platform.getPlatformByName('CUDA').getName())" || exit 1

ENTRYPOINT ["python", "/app/run_openmm_kcx_v5.py"]