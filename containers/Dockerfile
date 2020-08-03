# To build: docker build . -t leandroishilima/varifier
# However, prefer to pull it from dockerhub
FROM ubuntu:18.04
RUN apt update && apt install -y python3 python3-pip mummer curl git zlib1g-dev

# install minimap2
RUN curl -L https://github.com/lh3/minimap2/releases/download/v2.17/minimap2-2.17_x64-linux.tar.bz2 | tar -jxvf -

# install k8
RUN curl -L https://github.com/attractivechaos/k8/releases/download/v0.2.4/k8-0.2.4.tar.bz2 | tar -jxf -  && cp k8-0.2.4/k8-Linux /usr/bin

# install varifier
RUN git clone --recursive https://github.com/leoisl/varifier && cd varifier && git checkout pandora && pip3 install .

ENV PATH="/minimap2-2.17_x64-linux:${PATH}"