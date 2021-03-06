Data and code for Simon et al (2020) Developmental Cell
DOI https://doi.org/10.1016/j.devcel.2020.09.030
https://github.com/therealkatlab/Simon2020

Included in this repository are:

1. CellProfiler pipelines for the segmentation and measurement of cytoplasmic and nuclear fluorescence intensities of ERK-KTR-mClover in confocal time-lapse images of hemizygousÊHprtERK-KTRÊand homozygousÊR26NLS-mKate2Êembryos during pre-implantation development. Pipeline A was used by default and Pipeline B and C were used for cells where segmentation was challenging. 
	CellProfiler/EmbryoPipeline_A.cppipe
	CellProfiler/EmbryoPipeline_B.cppipe
	CellProfiler/EmbryoPipeline_C.cppipe

2. Matlab script to calculate and output images of the ERK-KTR cytoplasmic : nuclear (C:N) ratio using the CellProfiler outputs.
	Matlab/CN_Ratio_mean.m

3. Data of ERK-KTR C:N values from tracked cells in Imaris in Control, MEKi and FGF treated embryos and untreated embryos that were fixed and immunostained to determine end-point lineage identities. 
	R_and_data/FGF.csv
	R_and_data/PD03.csv
	R_and_data/Lineage_2hr.csv
	R_and_data/Lineage_12hr.csv

4. R scripts for the analysis of the cell tracking data and ERK-KTR C:N used to plot figures and perform statistical tests.
	R_and_data/FGF_Analysis.R
	R_and_data/MEKi_Analysis.R
	R_and_data/Lineage_2hr_Analysis.R
	R_and_data/Lineage_12hr_Analysis.R

5. R and Matlab scripts used for setting ERK activity threshold and to plot figures for 2 hour movies.
	Matlab/Threshold_finder.m
	R_and_data/FGF_Threshold.R
	R_and_data/MEKi_Threshold.R
	R_and_data/Lineage_2hr_Threshold.R

6. R and Matlab scripts for identifying peaks, plot figures and perform statistical tests for 2 hour movies.
	Matlab/PeakCall_smooth_FGF.m
	Matlab/PeakCall_smooth_PD03.m
	Matlab/PeakCall_smooth.m
	R_and_data/FGF_Peak_Analysis.R
	R_and_data/MEKi_Peak_Analysis.R
	R_and_data/Lineage_2hr_Analysis.R

These scripts were compiled in CellProfiler 2.1 (Lamprecht, et al., 2007), Matlab R2019a, R Studio 1.0.143, and R version 3.5.2  

For raw imaging data see Figshare repository DOI dx.doi.org/10.6084/m9.figshare.12794810
