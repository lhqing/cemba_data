# normal index
hisat-3n-build --base-change C,T -p THREAD \
~/ref/hg38/fasta/with_chrl/hg38_with_chrl.fa \
~/ref/hg38/fasta/with_chrl/hg38_with_chrl

# repeat index
hisat-3n-build --base-change C,T -p THREAD --repeat-index \
~/ref/hg38/fasta/with_chrl/hg38_with_chrl.fa \
~/ref/hg38/fasta/with_chrl/hg38_with_chrl.repeat

