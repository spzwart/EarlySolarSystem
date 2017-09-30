#! /bin/bash

# Download and extract AMUSEa
export GITHASH=79ec73d1ac2409c7564b927a5fa70b323d8362d0
wget https://github.com/rieder/amuse/archive/${GITHASH}.zip && \
  unzip ${GITHASH}.zip && \
  mv amuse-${GITHASH} amuse && \
  rm ${GITHASH}.zip 

export SCRIPT_DIR="${PWD}"
export AMUSE_DIR="${PWD}/amuse"

# Build AMUSE
cd ${AMUSE_DIR}
./configure && \
make framework

cd ${SCRIPT_DIR}
