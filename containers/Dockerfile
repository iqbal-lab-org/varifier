FROM ubuntu:18.04

# install basic dependencies
RUN apt update && apt install -y python3 python3-pip mummer curl git zlib1g-dev

# alias python -> python3
RUN ln -s /usr/bin/python3 /usr/bin/python

# install minimap2
RUN curl -L https://github.com/lh3/minimap2/releases/download/v2.17/minimap2-2.17_x64-linux.tar.bz2 | tar -jxvf -

# add minimap2, paftools.js and k8 to the path (all these tools come bundled with minimap2 and are in the minimap2 dir)
ENV PATH="/minimap2-2.17_x64-linux:${PATH}"

# install varifier
ARG varifier_repo_URL
ARG commit_id_or_branch
RUN git clone --recursive $varifier_repo_URL && cd varifier && git checkout $commit_id_or_branch && pip3 install . && tox
