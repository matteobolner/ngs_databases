for study in `cat $1`; do
	curl "https://www.ebi.ac.uk/ena/data/warehouse/filereport?accession=$study&result=read_run&fields=study_accession,sample_accession,secondary_sample_accession,experiment_accession,run_accession,tax_id,scientific_name,instrument_model,library_layout,fastq_ftp,fastq_galaxy,submitted_ftp,submitted_galaxy,sra_ftp,sra_galaxy,library_strategy,library_source,library_selection,sample_title,sample_alias,read_count" >> $2
done
