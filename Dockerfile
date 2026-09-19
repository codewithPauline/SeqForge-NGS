FROM mambaorg/micromamba:2.0.5
COPY --chown=$MAMBA_USER:$MAMBA_USER environment-linux-64.lock /tmp/environment-linux-64.lock
RUN micromamba install --yes --name base --file /tmp/environment-linux-64.lock && \
    micromamba clean --all --yes
# Nextflow starts a non-login shell and overrides the image entrypoint.
ENV PATH=/opt/conda/bin:$PATH
USER root
# Fail during image construction if Nextflow's metrics dependency is missing.
RUN ps --version
WORKDIR /data
LABEL org.opencontainers.image.title="SeqForge-NGS" \
      org.opencontainers.image.source="https://github.com/codewithPauline/SeqForge-NGS" \
      org.opencontainers.image.description="Tools for the SeqForge germline workflow"
