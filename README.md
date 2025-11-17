# Active-Microorganisms-during-N2O-production
Using BONCAT-FACS to probe the active microbial community during N2O production

These scripts accompany manuscript submitted by Gray, Kaye and Couradeau in 2025.  Each script or R markdown file has specific annotations throughout the file to explain the code used. Below is a list of the scripts and what each script contains.

16S_rRNA_data_processing
- This script contains the code to process our 16S rRNA amplicon sequences from the raw reads to a phyloseq object used in downstream analyses.

ASV_finder
- This script was specifically used to find specific ASVs and their associated 16S rRNA amplicon sequences when needed.

N2O_compiled_data
- Each N2O sampling time point was processed separately.  Here, they are compiled into a finalized data sheet used in downstream analysis.
- Raw data is converted into relevant units and linear models were run in this file.

NO3_NH4_analysis
- NO3 and NH4 were extracted at the ends of each incubation period. Here, the N concentrations are converted to relevant units and analyzed for outliers.
- Linear models were created

activity_data_analysis
- BONCAT-FACS activity data was processed initially using Kaluza software.  Data sheets were exported from that software and analyzed here.
- Cell counts, percent activity and percent inactivity were analyzed here.
- Fluorescence data was calculated by day
- cell count plots figures were made here

meqo_analysis_cleaner
- Code used in mEQO analysis

figures
- code used to make figures in manuscript
- compiles output files from other scripts
- includes necessary data manipulation if data is not ready to be plotted (phyloseq objects needed preprocessing steps before plotting)

T1_1
- individual time point N2O data processed here
- "T1" means time point 1, 0-3 hours,
- "_1" means it was the first sample in time point 1
