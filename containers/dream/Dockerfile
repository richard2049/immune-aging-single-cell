FROM bioconductor/bioconductor_docker:RELEASE_3_23@sha256:1d871e1ca9cca76b220eb16e22677e728f4352f81a9ee91aaf29e24aea43e624

ARG EXPECTED_R_VERSION=4.6
ARG EXPECTED_BIOCONDUCTOR_VERSION=3.23
ARG EXPECTED_EDGER_VERSION=4.10.1
ARG EXPECTED_VARIANCE_PARTITION_VERSION=1.42.0
ARG EDGER_SOURCE_URL=https://bioconductor.org/packages/3.23/bioc/src/contrib/Archive/edgeR/edgeR_4.10.1.tar.gz

RUN R --vanilla -e "\
    BiocManager::install(c('edgeR', 'variancePartition'), ask = FALSE, update = FALSE); \
    install.packages('${EDGER_SOURCE_URL}', repos = NULL, type = 'source'); \
    stopifnot(startsWith(as.character(getRversion()), '${EXPECTED_R_VERSION}')); \
    stopifnot(as.character(BiocManager::version()) == '${EXPECTED_BIOCONDUCTOR_VERSION}'); \
    stopifnot(as.character(packageVersion('edgeR')) == '${EXPECTED_EDGER_VERSION}'); \
    stopifnot(as.character(packageVersion('variancePartition')) == '${EXPECTED_VARIANCE_PARTITION_VERSION}')"

LABEL org.opencontainers.image.title="immune-aging repeated-measures runtime"
LABEL org.opencontainers.image.description="Pinned dream pseudobulk DE runtime"
LABEL org.opencontainers.image.version="bioc-3.23-dream-1.42.0"

WORKDIR /work
