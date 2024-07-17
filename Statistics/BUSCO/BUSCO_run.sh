#busco work example:

#!/bin/sh
busco -i DBG.fasta -c 20 -o DBG --augustus -m genome -l busco_downloads/lineages/aves_odb10 --offline

# Related software/library versions:
# BUSCO version is: 5.5.0 
# The lineage dataset is: aves_odb10 (Creation date: 2021-02-19, number of genomes: 62, number of BUSCOs: 8338)
# Dependencies and versions:
#  hmmsearch: 3.1
#  bbtools: 39.01
#  makeblastdb: 2.13.0+
#  tblastn: 2.13.0+
#  augustus: 3.5.0
#  busco: 5.5.0
