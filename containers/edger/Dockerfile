FROM bioconductor/bioconductor_docker:RELEASE_3_23@sha256:1d871e1ca9cca76b220eb16e22677e728f4352f81a9ee91aaf29e24aea43e624

ARG EXPECTED_R_VERSION=4.6
ARG EXPECTED_BIOCONDUCTOR_VERSION=3.23
ARG EXPECTED_EDGER_VERSION=4.10.1

RUN R --vanilla -e "\
    BiocManager::install('edgeR', ask = FALSE, update = FALSE); \
    stopifnot(startsWith(as.character(getRversion()), '${EXPECTED_R_VERSION}')); \
    stopifnot(as.character(BiocManager::version()) == '${EXPECTED_BIOCONDUCTOR_VERSION}'); \
    stopifnot(as.character(packageVersion('edgeR')) == '${EXPECTED_EDGER_VERSION}')"

LABEL org.opencontainers.image.title="immune-aging edgeR runtime"
LABEL org.opencontainers.image.description="Pinned pseudobulk DE runtime"
LABEL org.opencontainers.image.version="bioc-3.23-edger-4.10.1"

WORKDIR /work
