###README file for B.rapa (Chiifu) database files 

##### raw_data folder- commented out if I do not have these files
<!--(
<!-- 1. Chromosome file (fasta) - Brapa_sequence_v1.2.fa -->

<!-- 2. CDS (Coding sequence) file (fasta) - Brassica_rapa_final.cds -->

3. Annotation file (bed) - Brassica_rapa_final.bed

4. Gene description annotation file (txt) - Brassica_rapa_final.annotation.txt

5. GO annotation file (wego) - Brassica_rapa_final.annot.wego 

<!-- 6. SNP (between R500 & IMB211) file (vcf) - R500_IMB211.2014-03-25.vcf -->

Notes:
2015_01_22

# Genome DB
awk '/^>/ {OUT=substr($0,2) ".fa"}; OUT {print >OUT}' Brapa_sequence_v1.5.fa

# Promotor analysis
#installed genometools with homebrew and ran the following code to make a gff3 file from a bed file
gt bed_to_gff3 -featuretype gene -thicktype mRNA -blocktype CDS Brassica_rapa_v1.5_final.bed > Brassica_rapa_v1.5_final_correct.gff3

one of the regex searches for gene name in the extract-utr.pl was not working for my gff3 file
I changed line 126 to:
        my ($gene) = $attributes =~ /Name=(Bra.+)/;
        from:
        my ($gene) = $attributes =~ /Parent=(?:mRNA:)?([^;]+)/;

finally used these options on the perl script in order to get it to run

perl extract-utr.pl   --gff_file        /Users/Cody_2/git.repos/brassica_genome_db/raw_data/Brassica_rapa_v1.5_final_correct.gff3   --genome_fa_file  /Users/Cody_2/git.repos/brassica_genome_db/raw_data/Brapa_sequence_v1.5.fa   --cds_fa_file     /Users/Cody_2/git.repos/brassica_genome_db/raw_data/Brassica_rapa_v1.5_final.fa   --fiveprime   --utr_length      1000   --gene_length     0   --output_fa_file  Brapa_1000bp_upstream_3.fa

copied output file into brassica_genome_db directory )-->