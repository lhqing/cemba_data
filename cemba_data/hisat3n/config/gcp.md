# Setup GCP image for mapping

## Create base image

```bash
# init install system tools
sudo yum install -y zsh tree wget screen git nfs-utils make gcc

# install mambaforge
wget https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-Linux-x86_64.sh
sh Mambaforge-Linux-x86_64.sh -b -p $HOME/mambaforge
rm -f Mambaforge-Linux-x86_64.sh
./mambaforge/bin/mamba init zsh
./mambaforge/bin/mamba init bash
exec /bin/zsh
mamba install -y gxx

# Create mapping env hisat3n_env.yml
wget https://raw.githubusercontent.com/lhqing/cemba_data/master/hisat3n_env.yml
mamba env update -f hisat3n_env.yml  # this should install things in the base env

# Install packages
mkdir -p ~/pkg

# install hisat-3n
cd ~/pkg
git clone https://github.com/DaehwanKimLab/hisat2.git hisat-3n
cd hisat-3n
git checkout hisat-3n-dev-directional-mapping-reverse
make
# put hisat-3n in the PATH
echo 'export PATH=$HOME/pkg/hisat-3n:$PATH' >> ~/.bashrc
source ~/.bashrc 
echo 'export PATH=$HOME/pkg/hisat-3n:$PATH' >> ~/.zshrc 
source ~/.zshrc

# make sure allcools and yap is upto date
cd ~/pkg
git clone https://github.com/lhqing/cemba_data.git
cd cemba_data
pip install -e .

cd ~/pkg
git clone https://github.com/lhqing/ALLCools.git
cd ALLCoools
pip install -e .

## Create genome reference

# add genome reference file
# prepare and copy specific genome reference file to $HOME

# prepare a $HOME/mapping.yaml file the records the path of required genome reference files

# clean unnecessary cache files
mamba clean -y -a
```

## Actual mapping

```bash
mkdir -p ~/mapping
cd ~/mapping
gsutil cp gs://PATH/TO/FASTQ_DIR/fastq ./
cp ~/pkg/cemba_data/hisat3n/snakefile/SNAKEFILE_YOU_WANT_TO_USE ./Snakefile

# run snakemake
snakemake --configfile ~/mapping.yaml -j
```

## Build hisat-3n index
```bash
# non-repeat index
hisat-3n-build --base-change C,T genome.fa genome
# repeat index
hisat-3n-build --base-change T,C --repeat-index genome.fa genome
# Build the repeat HISAT-3N integrated index with splice site information
hisat-3n-build --base-change C,T --repeat-index --ss genome.ss --exon genome.exon genome.fa genome 
```
