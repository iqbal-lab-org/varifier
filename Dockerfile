FROM continuumio/miniconda3:4.9.2
ENV MINIMAP_VERSION 2.17
ENV MUMMER_VERSION "4.0.0rc1"
ENV BIOPY_VERSION 1.77
ENV CLUSTER_VERSION 0.13.2
ENV PANDAS_VERSION 1.1.0
ENV PYFASTAQ_VERSION 3.17.0
ENV PYMUMMER_VERSION 0.11.0
ENV PYSAM_VERSION 0.16
ENV SEABORN_VERSION 0.10.1
ENV PYTEST_VERSION 6.0
ENV BCFTOOLS_VERSION 1.10.2
ENV PROJECT "varifier"

#RUN apt update && apt install -y procps pigz && apt-get clean && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

RUN conda config --add channels conda-forge && \
    conda config --add channels bioconda && \
    conda config --add channels default

RUN conda install \
    bioconda::minimap2="$MINIMAP_VERSION" \
    bioconda::mummer4="$MUMMER_VERSION" \
    conda-forge::biopython="$BIOPY_VERSION" \
    bioconda::cluster_vcf_records="$CLUSTER_VERSION" \
    bioconda::mappy="$MINIMAP_VERSION" \
    conda-forge::pandas="$PANDAS_VERSION" \
    bioconda::pyfastaq="$PYFASTAQ_VERSION" \
    bioconda::pymummer="$PYMUMMER_VERSION" \
    bioconda::bcftools="$BCFTOOLS_VERSION" \
    bioconda::pysam="$PYSAM_VERSION" \
    conda-forge::seaborn="$SEABORN_VERSION" \
    conda-forge::pytest="$PYTEST_VERSION" \
    && conda clean -a

COPY . "/${PROJECT}"
WORKDIR "/${PROJECT}"
RUN python -m pip install --no-deps --ignore-installed .
RUN "$PROJECT" --help
RUN pytest
