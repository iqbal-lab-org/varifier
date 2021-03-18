#!/usr/bin/env bash
set -vexu

install_root=$1

apt-get install -y software-properties-common
apt-add-repository universe
apt-get update

apt-get install -y \
  build-essential \
  cmake \
  curl \
  git \
  libbz2-dev \
  libcurl4-gnutls-dev \
  liblzma-dev \
  libssl-dev \
  python3-pip \
  python3-setuptools \
  tabix \
  libvcflib-tools \
  wget \
  zlib1g-dev



if [ ! -d $install_root ]; then
  mkdir $install_root
fi
cd $install_root


#_________________________ bcftools _______________________#
cd $install_root
wget https://github.com/samtools/bcftools/releases/download/1.10.2/bcftools-1.10.2.tar.bz2
tar xf bcftools-1.10.2.tar.bz2
cd bcftools-1.10.2/
make
cd ..
cp -s bcftools-1.10.2/bcftools .

#________________________ vt __________________________________#
cd $install_root
git clone https://github.com/atks/vt.git vt-git
cd vt-git
git checkout 2187ff6347086e38f71bd9f8ca622cd7dcfbb40c
make
cd ..
cp -s vt-git/vt .

#________________________ mummer ____________________________#
cd $install_root
wget https://github.com/mummer4/mummer/releases/download/v4.0.0rc1/mummer-4.0.0rc1.tar.gz
tar xf mummer-4.0.0rc1.tar.gz
rm mummer-4.0.0rc1.tar.gz
cd mummer-4.0.0rc1
./configure
make
make install

#________________________ minimap2 etc ______________________#
cd $install_root
git clone https://github.com/lh3/minimap2.git minimap2_git
cd minimap2_git
git checkout ccb0f7b05df3a17011f0d0f4388ddb301198871b
make
cd ..
cp -s minimap2_git/minimap2 .
curl -L https://github.com/attractivechaos/k8/releases/download/v0.2.4/k8-0.2.4.tar.bz2 | tar -jxf -
cp -s k8-0.2.4/k8-Linux k8
cp -s minimap2_git/misc/paftools.js .

#______________________ python ________________________________#
pip3 install tox "six>=1.14.0"

