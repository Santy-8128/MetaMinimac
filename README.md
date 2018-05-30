THIS IS CURRENTLY BEING UPDATED. PLEASE DO NOT USE IT.


# MetaMinimac

MetaMinimac is a tool to merge GWAS data imputed against different reference panels (using minimac4)

<<< SEE https://genome.sph.umich.edu/wiki/MetaMinimac FOR DOCUMENTATION >>>

Users should follow the following steps to compile MetaMinimac 

## Prerequisites

Automatic installation of MetaMinimac requires [cget](http://cget.readthedocs.io/en/latest/src/intro.html#installing-cget) and cmake v3.2. These prerequisites can be installed as follows:

Ubuntu 16.04
```
sudo apt-get install cmake python-pip python-dev
pip install cget
```
Ubuntu 14.04
```
sudo apt-get install software-properties-common
sudo add-apt-repository ppa:george-edison55/cmake-3.x
sudo apt-get update
sudo apt-get install cmake python-pip python-dev
pip install cget
```
MacOS
```
brew install cmake
sudo easy-install pip
pip install --user cget --ignore-installed six
```

## Installation
The easiest way to install MetaMinimac and its dependencies is to use the install.sh file provided.

```
cd MetaMinimac
bash install.sh
```

Alternatively, you can setup a dev environment cmake directly.
```bash
cd MetaMinimac
cget install -f ./requirements.txt                      # Install dependencies locally.
mkdir build && cd build                                 # Create out of source build directory.
cmake -DCMAKE_TOOLCHAIN_FILE=../cget/cget/cget.cmake .. # Configure project with dependency paths.
make                                                    # Build.
```

## Usage
A typical MetaMinimac command line for imputation is as follows
```bash
MetaMinimac -i ImputedPrefix1:ImputedPrefix2 \
            -o testRunPrefix 
```
Here ImputedPrefix1 and ImputedPrefix1 are the prefixes of the same GWAS data imputed against different reference panels (using minimac4)
