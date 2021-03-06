#!/bin/bash

# make TARGET, can be overridden with env
: ${TARGET:=$HOME/miniconda}

function install_miniconda {
	if [ -d $TARGET ]; then echo "file exists"; return; fi
	echo "installing miniconda to $TARGET"
        if [[ "$TRAVIS_OS_NAME" == "linux" ]]; then
             platform="Linux"
        elif [[ "$TRAVIS_OS_NAME" == "osx" ]]; then
             platform="MacOSX"
        fi
	wget http://repo.continuum.io/miniconda/Miniconda3-latest-$platform-x86_64.sh -O mc.sh -o /dev/null
	bash mc.sh -b -f -p $TARGET
}

install_miniconda
export PATH=$TARGET/bin:$PATH

conda config --set always_yes true
#conda update --all
conda install -q conda-build=3.16.1
