FROM bioconductor/bioconductor_docker:RELEASE_3_20

LABEL org.opencontainers.image.source="https://github.com/DaleAnnear/openSesame"

RUN R -q -e "BiocManager::install(version='3.20', ask=FALSE); BiocManager::install(c('sesame','sesameData','SummarizedExperiment','limma','DMRcate','EPICv2manifest','IlluminaHumanMethylationEPICv2anno.20a1.hg38'), ask=FALSE, update=FALSE)" \
 && R -q -e "install.packages(c('optparse','jsonlite'), repos='https://cloud.r-project.org')" \
 && mkdir -p /home/rstudio/.cache/R/ExperimentHub \
 && chown -R rstudio:rstudio /home/rstudio/.cache

# Nextflow runs containers with the host UID (normally 1000/rstudio). Cache
# SeSAMe resources for that user, not only for the Docker build's root user.
USER rstudio
RUN R -q -e "library(sesameData); sesameDataCacheAll()"

USER root
WORKDIR /work
