#!/bin/bash

annovar="$HOME/setup/software/annovar"

### For Avinput file as input ###
java -jar "$HOME/program/VIC/target/VIC-1.0.1.jar" \
	-b hg19 \
	-i $HOME/program/VIC/example/vic_test1.avinput \
	-o $HOME/program/VIC/example/vic_test1 \
	-input_type VCF \
	-db $HOME/program/VIC/vicdb \
	-table_annovar ${annovar}/table_annovar.pl \
	-convert2annovar ${annovar}/convert2annovar.pl \
	-annotate_variation ${annovar}/annotate_variation.pl \
	-d ${annovar}/humandb \
	-otherinfo true \
	-s $HOME/program/VIC/example/vic_test1.evidence_file \
	-skip_annovar

### For VCF file as input ###
#java -jar "$HOME/program/VIC/target/VIC-1.0.1.jar" \
#	-b hg19 \
#	-i $HOME/program/VIC/example/vic_test4.vcf \
#	-o $HOME/program/VIC/example/vic_test4 \
#	-input_type vcf \
#	-db $HOME/program/VIC/vicdb \
#	-table_annovar ${annovar}/table_annovar.pl \
#	-convert2annovar ${annovar}/convert2annovar.pl \
#	-annotate_variation ${annovar}/annotate_variation.pl \
#	-d ${annovar}/humandb \
#	-otherinfo true
