# Setup GCP for mapping

## Create base image

```bash
sudo yum install -y zsh tree wget screen git nfs-utils make gcc

wget https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-Linux-x86_64.sh
sh Mambaforge-Linux-x86_64.sh -b -p $HOME/mambaforge
rm -f Mambaforge-Linux-x86_64.sh
./mambaforge/bin/mamba init zsh
mamba install -y gxx
exec /bin/zsh

# Create mapping env
mkdir pkg
cd pkg
# use hisat3n_mapping_env.yaml

# install hisat-3n

# install 
```

## Create genome reference image

```bash
# add genome reference file
# build reference index
# save image to gcp
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



