sudo yum install -y zsh tree wget screen git nfs-utils make gcc

wget https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-Linux-x86_64.sh
sh Mambaforge-Linux-x86_64.sh -b -p $HOME/mambaforge
rm -f Mambaforge-Linux-x86_64.sh
./mambaforge/bin/mamba init zsh
./mambaforge/bin/mamba init bash

mamba install -y gxx
exec /bin/zsh
