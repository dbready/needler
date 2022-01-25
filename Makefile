.PHONY: help build_env download fit proteome liver kinase

PYTHON=poetry run python

help:
	@echo "Makefile to build the peptide database"

build_env :
	poetry install

download: build/raw/uniprot_sprot-only2020_03.tar.gz build/raw/tissueatlas/MSB-15-e8503-s007.zip

munge : build/pepdb_rts_60min_gradient.csv build/liver_rts_60min_gradient.csv build/kinase_rts_60min_gradient.csv

# run algorithm on each of the subsets. advisable to do this with a -j parameter to run them in parallel
# this is a huge number of jobs, so consider the chosen parameter set before running
# the definitions for these jobs are contained in the generated makefile build/needler_jobs_list.mk
fit : proteome liver kinase

# This is some magic to generate a makefile of the individual fit jobs
# to customize the jobs configuration, may need to edit the python generation file
build/needler_jobs_list.mk : src/build_needle_jobs_makefile.py
	mkdir --parents build/
	${PYTHON} src/build_needle_jobs_makefile.py --rts=1,3,7,15,30,45,60,90,120,180,240,300 --targets=1,2,3,4,5,10,20,40,50,100,200,500,1000,2000,5000,10000 --replicates=3 --solve_seconds=86400 build/needler_jobs_list.mk

# `-include` will not fail if file does not exist, make will automatically build this missing file as needed from the above target
-include build/needler_jobs_list.mk

build/raw/uniprot_sprot-only2020_03.tar.gz :
	mkdir --parents build/raw/
	## Download uniprot release: 2020_03
	wget --no-clobber --directory-prefix=build/raw ftp://ftp.uniprot.org/pub/databases/uniprot/previous_releases/release-2020_03/knowledgebase/uniprot_sprot-only2020_03.tar.gz

build/raw/tissueatlas/MSB-15-e8503-s007.zip :
	mkdir --parents build/raw/
	# Download proteomics data table from publication: A deep proteome and transcriptome abundance atlas of 29 healthy human tissues
	# Proteome data table https://doi.org/10.15252/msb.20188503
	wget --no-clobber --directory-prefix=build/raw/tissueatlas https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6379049/bin/MSB-15-e8503-s007.zip
	touch build/raw/tissueatlas/MSB-15-e8503-s007.zip

build/tissueatlas/Table_EV5.xlsx : build/raw/tissueatlas/MSB-15-e8503-s007.zip
	unzip -o -d build/tissueatlas build/raw/tissueatlas/MSB-15-e8503-s007.zip
	touch build/tissueatlas/Table_EV5.xlsx

build/tissueatlas/liver_accessions.csv : src/extract_tissue_accessions.py build/tissueatlas/Table_EV5.xlsx
	# extract the identified liver proteins and convert identifers to uniprot accessions with web service
	${PYTHON} src/extract_tissue_accessions.py --tissuedb=build/tissueatlas/Table_EV5.xlsx build/tissueatlas/liver_accessions.csv

build/db/uniprot.fasta : build/raw/uniprot_sprot-only2020_03.tar.gz
	# extract uniprot artifact
	mkdir -p build/db
	tar -xvf build/raw/uniprot_sprot-only2020_03.tar.gz --directory=build/db uniprot_sprot.fasta.gz
	gunzip build/db/uniprot_sprot.fasta.gz
	mv build/db/uniprot_sprot.fasta build/db/uniprot.fasta
	touch build/db/uniprot.fasta

build/db/human_peptides.csv : build/db/uniprot.fasta src/extract_proteins.py src/digest_proteins.py
	# extract human proteins and generate distinct tryptic peptides
	${PYTHON} src/extract_proteins.py --fasta=build/db/uniprot.fasta build/db/human_proteins.csv
	${PYTHON} src/digest_proteins.py --proteins=build/db/human_proteins.csv --min_length=5 --max_length=30 --min_peptides=2 --unfiltered_filename=build/db/human_peptides_withmet.csv build/db/human_peptides.csv

build/pepdb_rts_60min_gradient.csv : src/predict_rt.py build/db/human_peptides.csv
	# assign iRT and RTs to the chosen peptide library
	${PYTHON} src/predict_rt.py --peptides=build/db/human_peptides.csv --procal=data/procal_supplementaltable_S1.csv --prosit=data/prosit_predicted_irt.csv build/pepdb_rts_60min_gradient.csv

build/liver_rts_60min_gradient.csv : src/filter_accessions.py build/pepdb_rts_60min_gradient.csv build/tissueatlas/liver_accessions.csv
	${PYTHON} src/filter_accessions.py --pepdb=build/pepdb_rts_60min_gradient.csv --filtersrc=build/tissueatlas/liver_accessions.csv build/liver_rts_60min_gradient.csv

build/kinase_rts_60min_gradient.csv : src/filter_accessions.py build/pepdb_rts_60min_gradient.csv
	${PYTHON} src/filter_accessions.py --pepdb=build/pepdb_rts_60min_gradient.csv --filtersrc=data/uniprot_pkinfam_202004.csv build/kinase_rts_60min_gradient.csv
