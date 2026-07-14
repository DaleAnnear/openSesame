FROM bioconductor/bioconductor_docker:RELEASE_3_18

# Install required R packages: sesame, sesameData, and optparse
RUN R -e "BiocManager::install(c('sesame', 'sesameData', 'optparse', 'ExperimentHub', 'AnnotationHub'))"

# Cache SeSAMe data for EPICv2 arrays to speed up the pipeline
RUN R -e "library(sesameData); sesameDataCacheAll()"

# Create directories for data and scripts
RUN mkdir -p /data /scripts

# Set working directory
WORKDIR /data
