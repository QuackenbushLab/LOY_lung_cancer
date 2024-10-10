
This is the folder containing the data necessary to perform the main analysis of the present study. Due to size or privacy restrictions, not all data is contained within this repository. We provide further information on how to obtain this data below.




- The human reference annotation (v32 for GRCh38.p13), can be obtained from GENCODE (https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_32/gencode.v32.annotation.gff3.gz)

- The TF motif prior can be downloaded through the GRAND database (https://granddb.s3.us-east-2.amazonaws.com/tissues/motif/tissues_motif.txt)

- The WGS (.bam) files fall under access restrictions as they contain sensitive information about individuals. An interested user has to separately apply for access to these files through GDC and the TCGA consortium. Once .bam files have been retrieved, the coverage counts per chromosome used in this study can be computed through samtools coverage (see here for reference http://www.htslib.org/doc/samtools-coverage.html).

- The bulk LUAD data from the Lung Cancer Consortium Singapour (LCCS) can be obtained through their GIS app with identifier GIS031 (https://src.gisapps.org/OncoSG_public/study/summary?id=GIS031)

- For sc data from the Human Tumor Atlas Network (HTAN), please follow the instructions in the sc_data folder to download Seurat .rds files (all cells) from the CELLxCELL portal.

- The validation sc data of Nayoung et al. can be obtained through GEO with identifier GSE131907, we here use files GSE131907_Lung_Cancer_cell_annotation.txt, GSE131907_Lung_Cancer_Feature_Summary.csv, and GSE131907_Lung_Cancer_normalized_log2TPM_matrix.rds
