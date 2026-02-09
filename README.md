# LOY_lung_cancer
This is the companion repository for our study on loss of Y chromosome (LOY) in lung cancer.

## Organisation of the repository

Besides this Readme, the repository contains the Supplementary Tables ("LOY_LUAD_supplementary_tables.zip") accompanying our manuscript (see Reference below).
For reproducibility, we share R code for recomputing the our analysis ("codebase/analysis/loy_lung_cancer.Rmd") and the routine for active subnetwork discovery ("codebase/analysis/bionet_loy.R").
All shareable data is deposited under "data/", containing intermediate results for easy use in our analysis script, as well as links to sources where data can be obtained from in the related READMEs in case data is too large to share or requires additional access grants due to licensing or safety.

## Reference

[BiorXiv](https://www.biorxiv.org/content/10.1101/2024.09.19.613876v2)

@article {Fischer:LOYLUAD,
	author = {Fischer, Jonas and Shutta, Katherine H. and Chen, Chen and Fanfani, Viola and Saha, Enakshi and Mandros, Panagiotis and Guebila, Marouen Ben and Xiu, Joanne and Nieva, Jorge and Liu, Stephen and Uprety, Dipesh and Spetzler, David and Lopes-Ramos, Camila M. and DeMeo, Dawn L. and Quackenbush, John},
	title = {Selective loss of Y chromosomes in lung adenocarcinoma modulates the tumor immune environment through cancer/testis antigens},
	year = {2025},
	doi = {10.1101/2024.09.19.613876},
	publisher = {Cold Spring Harbor Laboratory},
	URL = {https://www.biorxiv.org/content/early/2025/11/05/2024.09.19.613876},
	journal = {bioRxiv}
}

