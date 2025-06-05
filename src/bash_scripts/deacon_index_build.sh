# Usage: deacon_index_build.sh
#!/bin/bash

deacon index build GCF_000001735.4_TAIR10.1_genomic.fna >GCF_000001735.4_TAIR10.1_genomic.fna.index

makeblastdb \
    -in /home/maurice/resources/host_genomes/arabidopsis_thaliana/whole/GCF_000001735.4_TAIR10.1_genomic.fna \
    -dbtype nucl \
    -out /home/maurice/resources/host_genomes/arabidopsis_thaliana/whole/blast/GCF_000001735.4_TAIR10.1_genomic

makeblastdb \
    -in /home/maurice/resources/host_genomes/lactuca_sativa/whole/GCF_002870075.4_Lsat_Salinas_v11_genomic.fna \
    -dbtype nucl \
    -out /home/maurice/resources/host_genomes/lactuca_sativa/whole/blast/GCF_002870075.4_Lsat_Salinas_v11_genomic

makeblastdb \
    -in /home/maurice/resources/host_genomes/eruca_sativa/whole/GCA_932364175.1_Eruca_sativa_genomic.fna \
    -dbtype nucl \
    -out /home/maurice/resources/host_genomes/eruca_sativa/whole/blast/GCA_932364175.1_Eruca_sativa_genomic

cat /home/maurice/resources/host_genomes/lactuca_sativa/whole/GCF_002870075.4_Lsat_Salinas_v11_genomic.fna \
    /home/maurice/resources/host_genomes/eruca_sativa/whole/GCA_932364175.1_Eruca_sativa_genomic.fna \
    /home/maurice/resources/host_genomes/arabidopsis_thaliana/whole/GCF_000001735.4_TAIR10.1_genomic.fna > \
    /home/maurice/resources/host_genomes/pan_crop/whole/pan_crop_genomic.fna

mkdir -p /home/maurice/resources/host_genomes/pan_crop/whole/blast

makeblastdb \
    -in /home/maurice/resources/host_genomes/pan_crop/whole/pan_crop_genomic.fna \
    -dbtype nucl \
    -out /home/maurice/resources/host_genomes/pan_crop/whole/blast/pan_crop_genomic

# lettuce
# cucumber
# tomato
# bell_pepper
# strawberry
# spinah
# beetroot
# rose
# Rocket/Arugula
# Pak choi; Tatsoi; Mizuna; Mustard red frills; Mustard red zest
# Fennel
# Kale
# Radish
# Chard
# Bok choy
# Sorrel
# Watercress

# Brassica juncea
# Brassica rapa
