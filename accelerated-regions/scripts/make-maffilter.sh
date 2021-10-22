#!/bin/bash

source config/settings.cfg

wd=$(pwd)/src

mkdir -p $wd

cd $wd

if [ ! -d "${wd}/bpp-core" ]; then
git clone https://github.com/BioPP/bpp-core.git
fi

if [ ! -d "${wd}/bpp-seq" ]; then
git clone https://github.com/BioPP/bpp-seq.git
fi

if [ ! -d "${wd}/bpp-seq-omics" ]; then
git clone https://github.com/kchiou/bpp-seq-omics.git
fi

if [ ! -d "${wd}/bpp-phyl" ]; then
git clone https://github.com/BioPP/bpp-phyl.git
fi

if [ ! -d "${wd}/bpp-phyl-omics" ]; then
git clone https://github.com/BioPP/bpp-phyl-omics.git
fi

if [ ! -d "${wd}/bpp-popgen" ]; then
git clone https://github.com/BioPP/bpp-popgen.git
fi

if [ ! -d "${wd}/maffilter" ]; then
git clone https://github.com/jydu/maffilter.git
fi

for i in bpp-core bpp-seq bpp-seq-omics bpp-phyl bpp-phyl-omics bpp-popgen maffilter; do
echo 'Now installing '${i}$'\n'
cd ${wd}/${i}
cmake -DCMAKE_INSTALL_PREFIX=${HOME}/${bpp_name}
make install
done

echo "To load maffilter, copy and paste the following:"
echo "export LD_LIBRARY_PATH=${HOME}/${bpp_name}/lib64:\$LD_LIBRARY_PATH"
echo "export PATH=${HOME}/${bpp_name}/bin:\$PATH"

exit