FROM ensemblorg/ensembl-vep:latest

USER root
RUN apt-get update -qq && \
    apt-get install -y --no-install-recommends python3 samtools && \
    rm -rf /var/lib/apt/lists/*

COPY server.py /server.py

ENV VEP_CACHE_DIR=/vep_cache \
    VEP_ASSEMBLY=GRCh38 \
    VEP_PORT=8080

EXPOSE 8080

CMD ["python3", "/server.py"]
