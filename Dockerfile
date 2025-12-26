
# --- Stage 1: Build Environment ---
# Using micromamba image as it is much faster and lighter than mambaforge
FROM mambaorg/micromamba:1.5-focal AS builder

# Set up the environment
# We install all scientific dependencies here. Micromamba handles 
# dependency solving much faster than standard conda, preventing timeouts.
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
# We switch to a slim Debian image to keep the final size small for Tamarind.
FROM debian:bookworm-slim

LABEL maintainer="D-Hydantoinase Project"
LABEL description="OpenMM with KCX (N6-carboxylysine) force field support"

# Copy the environment from the builder stage
COPY --from=builder /opt/conda /opt/conda

# Set environment paths so Python can find OpenMM and AmberTools
ENV PATH="/opt/conda/bin:$PATH"
ENV AMBERHOME="/opt/conda"

WORKDIR /app

# Copy your specific parameter files and script
# These must exist in your GitHub repository root
COPY kcx.frcmod /app/kcx.frcmod
COPY kcx.lib /app/kcx.lib
COPY run_openmm_kcx_v5.py /app/run_openmm_kcx.py

# Make the script executable
RUN chmod +x /app/run_openmm_kcx.py

# Final verification step: If this fails, the GitHub Action will show a Red X
# so you know immediately if the installation is broken.
RUN python -c "import openmm; print('OpenMM version:', openmm.__version__)"

# Entrypoint for Tamarind to execute your custom simulation
ENTRYPOINT ["python", "/app/run_openmm_kcx.py"]
