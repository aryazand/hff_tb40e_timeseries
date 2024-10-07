datasets download genome accession GCF_023101765.2 --filename spodoptera_dataset.zip
datasets download genome accession GCA_027928905.1 --filename towne_dataset.zip
datasets download genome accession GCF_000001405.40 --filename human_GRCh38_dataset.zip
datasets download genome accession GCA_027926625.1 --filename  tb40E_dataset.zip 

unzip human_GRCh38_dataset.zip -d human_GRCh38_dataset
unzip towne_dataset.zip -d towne_dataset
unzip spodoptera_dataset.zip -d spodoptera_dataset
unzip tb40E_dataset.zip -d tb40E_dataset

cat human_GRCh38_dataset/ncbi_dataset/data/*/*.fna towne_dataset/ncbi_dataset/data/*/*.fna spodoptera_dataset/ncbi_dataset/data/*/*.fna > human_towne_spodoptera_combined_genome.fa

cat human_GRCh38_dataset/ncbi_dataset/data/*/*.fna tb40E_dataset/ncbi_dataset/data/*/*.fna spodoptera_dataset/ncbi_dataset/data/*/*.fna > human_tb40E_spodoptera_combined_genome.fa