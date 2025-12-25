FROM docker.io/condaforge/mambaforge:latest

LABEL maintainer="D-Hydantoinase Project"
LABEL description="OpenMM with KCX (N6-carboxylysine) force field support"

WORKDIR /app

# Install all packages in ONE step to avoid multiple dependency solves (reduces memory)
# Use --no-update-deps to prevent re-solving existing packages
RUN mamba install -y -c conda-forge --no-update-deps \
    python=3.10 \
    numpy \
    scipy \
    pandas \
    openmm \
    parmed \
    mdtraj \
    ambertools \
    openmmforcefields \
    && mamba clean -afy

# Copy force field files
COPY kcx.frcmod /app/kcx.frcmod
COPY kcx.lib /app/kcx.lib

# Copy the OpenMM KCX simulation script
COPY run_openmm_kcx_v5.py /app/run_openmm_kcx.py
RUN chmod +x /app/run_openmm_kcx.py

# Set environment
ENV AMBERHOME=/opt/conda
ENV PATH=/opt/conda/bin:$PATH

# Test installation
RUN python -c "import openmm; print('OpenMM version:', openmm.__version__)"

# Set entrypoint to the simulation script
ENTRYPOINT ["python", "/app/run_openmm_kcx.py"]
