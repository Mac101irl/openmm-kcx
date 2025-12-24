FROM continuumio/miniconda3:latest

LABEL maintainer="D-Hydantoinase Project"
LABEL description="OpenMM with KCX (N6-carboxylysine) force field support"

WORKDIR /app

# Install OpenMM and dependencies via conda
RUN conda install -y -c conda-forge \
    python=3.10 \
    openmm>=8.0 \
    parmed \
    ambertools>=22 \
    mdtraj \
    numpy \
    scipy \
    pandas \
    openmmforcefields \
    && conda clean -afy

# Copy force field files
COPY kcx.frcmod /app/kcx.frcmod
COPY kcx.lib /app/kcx.lib

# Set environment
ENV AMBERHOME=/opt/conda
ENV PATH=/opt/conda/bin:$PATH

# Test installation
RUN python -c "import openmm; print('OpenMM version:', openmm.__version__)"

CMD ["python", "-c", "import openmm; print('OpenMM-KCX ready')"]
