PROT_DB_DIR=/group/ctbrowngrp/gtdb/databases
NUCL_DB_DIR=/group/ctbrowngrp/gtdb/databases/ctb
DEST_DIR=test-data/databases


# protein alpha
#for k in 10 #7 10 11
#  do
#    sourmash sig extract --protein -k ${k} --no-dna --picklist gtdb-nine.picklist.csv:name:identprefix ${PROT_DB_DIR}/gtdb-rs202.protein.k${k}.sbt.zip | sourmash sig cat - -o ${DEST_DIR}/gtdb-nine.protein-k${k}.zip
#    sourmash index ${DEST_DIR}/gtdb-nine.protein-k${k}.sbt.zip ${DEST_DIR}/gtdb-nine.protein-k${k}.zip
#  done

# nucl alpha
for k in 31 51 #21
  do
    sourmash sig extract -k ${k} --dna --picklist gtdb-nine.picklist.csv:name:identprefix ${NUCL_DB_DIR}/gtdb-rs202.genomic.k${k}.sbt.zip | sourmash sig cat - -o ${DEST_DIR}/gtdb-nine.genomic-k${k}.zip
    sourmash index ${DEST_DIR}/gtdb-nine.nucleotide-k${k}.sbt.zip ${DEST_DIR}/gtdb-nine.genomic-k${k}.zip
  done
