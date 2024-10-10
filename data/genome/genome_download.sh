# Download data
datasets download genome accession GCF_023101765.2 --filename spodoptera_dataset.zip
datasets download genome accession GCA_027928905.1 --filename towne_dataset.zip
datasets download genome accession GCF_000001405.40 --filename human_GRCh38_dataset.zip
datasets download genome accession GCA_027926625.1 --filename  tb40E_dataset.zip 

# Unzip
unzip human_GRCh38_dataset.zip -d human_GRCh38_dataset
unzip towne_dataset.zip -d towne_dataset
unzip spodoptera_dataset.zip -d spodoptera_dataset
unzip tb40E_dataset.zip -d tb40E_dataset

# Simplify CMV genome names
sed -i .bak 's/^>.*/>KF297339.1/' tb40E_dataset/ncbi_dataset/data/GCA_027926625.1/GCA_027926625.1_ASM2792662v1_genomic.fna

sed -i .bak 's/^>.*/>FJ616285.1/' towne_dataset/ncbi_dataset/data/GCA_027928905.1/GCA_027928905.1_ASM2792890v1_genomic.fna

# Combine genomes 
cat human_GRCh38_dataset/ncbi_dataset/data/*/*.fna towne_dataset/ncbi_dataset/data/*/*.fna spodoptera_dataset/ncbi_dataset/data/*/*.fna | sed '/^>/ s/[[:space:]]/\_/g' > human_towne_spodoptera_combined_genome.fa

cat human_GRCh38_dataset/ncbi_dataset/data/*/*.fna tb40E_dataset/ncbi_dataset/data/*/*.fna spodoptera_dataset/ncbi_dataset/data/*/*.fna | sed '/^>/ s/[[:space:]]/\_/g' > human_tb40E_spodoptera_combined_genome.fa

# Build bowtie indices 
bowtie-build --threads 10 human_tb40E_spodoptera_combined_genome.fa human_tb40E_spodoptera_combined_genome

bowtie-build --threads 10 human_towne_spodoptera_combined_genome.fa human_tb40E_spodoptera_combined_genome