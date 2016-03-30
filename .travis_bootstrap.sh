#!/bin/bash

# does not work with clang yet
if [ "$CXX" = "clang++" ]; then exit; fi

git clone https://github.com/open-source-parsers/jsoncpp.git
cd jsoncpp

mkdir -p build/debug
cd build/debug

cmake -DCMAKE_INSTALL_PREFIX=/home/travis/jsoncpp/install -DCMAKE_BUILD_TYPE=debug -DJSONCPP_WITH_PKGCONFIG_SUPPORT=OFF -G "Unix Makefiles" ../..

make
make install
