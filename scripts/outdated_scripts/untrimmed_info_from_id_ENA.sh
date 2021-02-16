for study in `cat $1`; do
curl "https://www.ebi.ac.uk/ena/data/warehouse/filereport?accession=$study&result=read_run" >> $2
done
