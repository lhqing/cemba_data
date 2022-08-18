# Setup GCP image for mapping

## Create base image

```bash
sudo yum install -y zsh tree wget screen git nfs-utils make gcc

wget https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-Linux-x86_64.sh
sh Mambaforge-Linux-x86_64.sh -b -p $HOME/mambaforge
rm -f Mambaforge-Linux-x86_64.sh
./mambaforge/bin/mamba init zsh
./mambaforge/bin/mamba init bash
exec /bin/zsh
mamba install -y gxx

# Create mapping env
mkdir pkg
cd pkg

# use hisat3n_mapping_env.yaml
wget https://raw.githubusercontent.com/lhqing/cemba_data/master/hisat3n_env.yml
mamba env update -y hisat3n_env.yml

# install hisat-3n
git clone https://github.com/DaehwanKimLab/hisat2.git hisat-3n
cd hisat-3n
git checkout hisat-3n-dev-directional-mapping-reverse
make

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

## Notes
yap-gcp run

1. pipeline snakefile
2. config yaml
3. source path
4. target path


yap-gcp make_file_list source path
1. mkdir, update software
2. gsutil cp filelist, create a flag if success
3. validate copy flag, if success, run yap snakemake

if mapping_summary exist, create and cp flag file to target and to source

# on daemon, after cloud job gone, check
yap-gcp validate_success source_path target_path
check flag on both side, remove source if completed, archive target file, delete non-active data (FASTQ / BAM) after archive


## Build hisat-3n index
```bash
# non-repeat index
hisat-3n-build --base-change C,T genome.fa genome
# repeat index
hisat-3n-build --base-change T,C --repeat-index genome.fa genome
# Build the repeat HISAT-3N integrated index with splice site information
hisat-3n-build --base-change C,T --repeat-index --ss genome.ss --exon genome.exon genome.fa genome 
```
