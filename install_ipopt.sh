wget https://www.coin-or.org/download/source/Ipopt/Ipopt-3.11.10.zip
unzip Ipopt-3.11.10.zip
cd Ipopt-3.11.10/ThirdParty/Mumps
./get.Mumps
cd ../../
mkdir build && cd build
../configure --prefix=/usr/local
make
make test
make install
export IPOPT_DIR=`pwd`
