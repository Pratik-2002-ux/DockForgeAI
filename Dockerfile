FROM mambaorg/micromamba:2.3.2

USER root
RUN apt-get update && apt-get install -y --no-install-recommends curl tini \
    && rm -rf /var/lib/apt/lists/*
USER $MAMBA_USER

WORKDIR /app
COPY --chown=$MAMBA_USER:$MAMBA_USER environment.yml /tmp/environment.yml
RUN micromamba install --yes --name base --file /tmp/environment.yml \
    && micromamba clean --all --yes

COPY --chown=$MAMBA_USER:$MAMBA_USER app app
COPY --chown=$MAMBA_USER:$MAMBA_USER src src
COPY --chown=$MAMBA_USER:$MAMBA_USER scripts scripts
COPY --chown=$MAMBA_USER:$MAMBA_USER config config
COPY --chown=$MAMBA_USER:$MAMBA_USER data data
COPY --chown=$MAMBA_USER:$MAMBA_USER main.py README.md ./
RUN mkdir -p reports logs data/processed/final_reports data/processed/report_assets

EXPOSE 8501
HEALTHCHECK --interval=30s --timeout=5s --start-period=30s --retries=3 \
  CMD python -c "import urllib.request; urllib.request.urlopen('http://127.0.0.1:8501/_stcore/health', timeout=3)"

ENTRYPOINT ["/usr/bin/tini", "--", "/usr/local/bin/_entrypoint.sh"]
CMD ["streamlit", "run", "app/streamlit_app.py", "--server.address=0.0.0.0", "--server.port=8501"]
