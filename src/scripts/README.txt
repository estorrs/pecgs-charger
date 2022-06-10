##################################################################
## Germline workflow for PECGS pipeline 			##
## Steps from TinJasmine to CharGer 				##
## Contact: Fernanda Martins Rodrigues (fernanda@wustl.edu)	##
## Date modified: 05122022					##
##################################################################


## Steps taken from TinJasmine output to CharGer

The correct TinJasmine output file to use is the following: call-canonical_filter/execution/output/HotspotFiltered.vcf.gz
An example file for sample C3L-00081 from CPTAC LSCC cohort is placed here: /diskmnt/Projects/Users/fernanda/Projects/PECGS/PECGS_pipeline_CharGer/Data/TinJasmine_outputs/CPTAC_LSCC/C3L-00081/call-canonical_filter/execution/output/HotspotFiltered.vcf.gz

This VCF has already been through the following filters:
	- variants with allelic depth < 5 for the alternative allele have been removed
	- indels longer than 100bp are removed
	- only variants in regions of interest (ROI) are present (i.e. variants in coding regions, as per bed file; contains exonic variants plus variants affect splice sites)

# Step 1: prepare TinJasmine VCF file for CharGer by removing INFO fields so that CharGer will consider gnomAD AF and not the AF field from variant caller)

	- Script:  Scripts/format_vcf_for_CharGer.py
	- Command:

	$ python Scripts/format_vcf_for_CharGer.py -i Anlysis/test.vcf.gz -O Analysis/1.Preprocess_TinJasmine

	- Output: Analysis/1.Preprocess_TinJasmine/test.infofixed.vcf


# Step 2: run CharGer; using PanCan gene lists and other files for CharGer run

	- Commands:

	$ conda activate CharGer

	# pancan
	$ nohup charger --include-vcf-details -f Analysis/1.Preprocess_TinJasmine/test.infoFixed.vcf -o Analysis/2.CharGer/test.charged.tsv -O -D --inheritanceGeneList Data/CharGer_dependencyFiles/PanCan/cancer_pred_genes_160genes_011321_curated_forCharGer.txt --PP2GeneList Data/CharGer_dependencyFiles/PanCan/160cpgs.txt -z Data/CharGer_dependencyFiles/PanCan/emptyRemoved_20160428_pathogenic_variants_HGVSg_VEP_grch38lifOver.vcf -H Data/CharGer_dependencyFiles/PanCan/cptac_mc3_combined_noHypers_sorted.maf.3D_Proximity.pairwise.recurrence.l0.ad10.r10.clusters -l --mac-clinvar-tsv Data/CharGer_dependencyFiles/clinvar_alleles.single.b38.tsv.gz --rare-threshold 0.0005 1>Analysis/Logs/test.charger.out 2>Analysis/Logs/test.charger.err &


# Step 3: process CharGer output: fix HGVSc and HGVSp annotations and gather genotype information from VCF. Note that the script assumes that your CharGer output file is named as follows: sampleID.charged.tsv . It uses the sample ID within the code.

	- Script: /diskmnt/Projects/Users/fernanda/Projects/PECGS/PECGS_pipeline_CharGer/Scripts/post_CharGer.py
	- Command:

	$ nohup python Scripts/post_CharGer.py -i Analysis/1.Preprocess_TinJasmine/test.infoFixed.vcf -c Analysis/2.CharGer/test.charged.tsv -s C3L-00081 -O Analysis/3.PostCharGer/ 1>Analysis/Logs/test.postCharGer.out 2>Analysis/Logs/test.postCharGer.err &


# Step 4: filter out synonymous variants and variants with low CharGer score (< 4); also create a file containing only rare variants (MAF ≤ 0.05%)

	- Script: /diskmnt/Projects/Users/fernanda/Projects/PECGS/PECGS_pipeline_CharGer/Scripts/filter_CharGer.py
	- Command:

	$ nohup python Scripts/filter_CharGer.py -c Analysis/3.PostCharGer/C3L-00081.charged2vcf.tsv -a 0.0005 -o C3L-00081 -O Analysis/4.FilterCharGer/ 1>Analysis/Logs/test.filterCharGer.out 2>Analysis/Logs/test.filterCharGer.err &


# Step 5: generate file for igv batch mode to be used in manual review steps

	STILL WORKING ON TUTORIAL FOR THIS; NOT NECESSARY FOR PIPELINE IMPLEMENTATION








