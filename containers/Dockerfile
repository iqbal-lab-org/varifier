FROM ubuntu:18.04

# install basic dependencies
RUN apt update && apt install -y python3 python3-pip mummer curl git zlib1g-dev

# install minimap2
RUN curl -L https://github.com/lh3/minimap2/releases/download/v2.17/minimap2-2.17_x64-linux.tar.bz2 | tar -jxvf -

# install k8
RUN curl -L https://github.com/attractivechaos/k8/releases/download/v0.2.4/k8-0.2.4.tar.bz2 | tar -jxf -  && cp k8-0.2.4/k8-Linux /usr/bin

# install varifier
ARG varifier_repo_URL
ARG commit_id_or_branch
RUN git clone --recursive $varifier_repo_URL && cd varifier && git checkout $commit_id_or_branch && pip3 install .

# alias python -> python3
RUN ln -s /usr/bin/python3 /usr/bin/python

# add minimap2 to the path
ENV PATH="/minimap2-2.17_x64-linux:${PATH}"