#! /bin/bash
set -e
yum install -y gmp-devel gcc gcc-c++ glpk-devel || ( apt-get update && apt-get -y install g++ libgmp-dev libglpk-dev autoconf automake libtool make)
mkdir -p /build
cd /build
/src/configure --srcdir=/src  || ( echo '#### Contents of config.log: ####'; cat config.log; exit 1)
make -j2 && make check
