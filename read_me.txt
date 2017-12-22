In order to go from raw sequencing reads to estimates of genetic interactions, 
there are 3 main steps:
1) Parse raw fastq reads to obtain counts for possible BC/SDS combination.
2) Estimate and remove chimeric PCR products from the count data using R scripts.
3) Estimate fitness and interaction scores using R scripts.
The following documents walks through how each step is performed. Downstream analysis 
is then performed using additional R markdown files.

*****************************************************************************************
*****************************************************************************************
Part 1: Parsing fastq files using python scripts
*****************************************************************************************
*****************************************************************************************
***Overview***
This process is done in two sequential steps using the scripts named 
"pull_tags_double_crispri_v2.py" and "parse_discarded.py". The overarching strategy is
to loop through each read, identify where in the read the "anchor" primer sequence is,
then use that position to identify the location of the multiplexing tag, DNA barcode
sequence, and the site-directed sequence. The first script only allows perfect matches
to the known sequence. While the second script parses through unassigned reads from the
first script and allows one mismatch for each identifier (and 2 mismatches for the 
primer sequence, one in the first position, and an additional mismatch anywhere in the 
remaining sequence).

***Setting parameters***
Before running "pull_tags_double_crispri_v2.py", one must set the following experiment
specific parameters in the first part of the script:
no_lines: number of lines in the fastq files (use wc -l command in bash to get this)
read1_name and read2_name: the file names for the paired end reads
primer1 and primer2: the primer sequences used for library prep that are used for searching and anchoring data from each read
tag_key: a dictionary where each key is a multiplexing tag and each value is the sample id
BC_key: a dictionary where each key is a DNA barcode sequence, and each value is a query plasmid name
SDS_key: a dictionary where each key is a site directed sequence and each value is a starting pool guide name
row_key: a dictionary where each key is a potential species in the pool (based on all pairwise combinations of BC_key and SDS_key values), and each value is the row index where that species data will be stored.
col_key: a dictionary where each key is a sample id (from values in tag_key), and each value is the column index where the data for that sample will be stored.
col_names: a list of all samples in the sequencing data (i.e. all keys in the col_key, in order of their values)
row_names: a list of all possible species in the sequencing data (i.e. all keys in the row_key, in order of their values)

Before running "parse_discarded.py", one must set the following experimental parameter:
primer1, primer2, tag_key, BC_key, SDS_key, Row_key, Col_key, col_names, row_names : same as above (can be copied and pasted)
p1_regex and p2_regex: a regular expression containing all potential primer sequences with either 1 mismatch anywhere, or 1 mismatch in the first position plus an additional mismatch anywhere
 
***Output***
The script "pull_tags_double_crispri_v2.py" outputs the following files (for all files an index of -1 means no match for the given identifier was identified):
"non_matching_reads_primer.txt": a file containing information for reads where a primer sequence was not identified in either or both of the forward and reverse reads (columns:line number,read1 sequence,read2 sequence,primer1index,primer2index)
"non_matching_reads_tag.txt": a file containing information for reads where both primers were identified, but the multiplexing tag identified was not a perfect match to a known tag (columns:line,read1,read2,primer1index,primer2index,tag1,tag2,BC,SDS)
"non_matching_reads_BC.txt": a file containing information for reads where both primers and multiplexing tag were identified, but the DNA barcode identified was not a perfect match to a known BC (columns:line,read1,read2,primer1index,primer2index,tag1,tag2,BC,SDS)
"non_matching_reads_SDS.txt": a file containing information for reads where both primers, multiplexing tag, and DNA barcode were identified, but the site directed sequence identified was not a perfect match to a known SDS (columns:line,read1,read2,primer1index,primer2index,tag1,tag2,BC,SDS)
"all_counts.txt": matrix of counts identified from sequencing reads for each species (row) in each sample/library (column)
The script also prints the following information to the terminal:
Number of reads parsed
Number of reads where neither primer sequence was identified in the read
Number of reads where primer 1 wasn't identified
Number of reads where primer 2 wasn't identified
Number of reads where both primers were identified, but multiplexing tag identified was not a perfect match to a known multiplexing tag
Number of reads where both primers and multi-tag were identified, but DNA barcode identified was not a perfect match to a known DNA barcode
Number of reads where both primers and multi-tag, DNA barcode were identified, but site directed sequence identified was not a perfect match to a known SDS

The script "parse_discarded.py" outputs the following files:
"recovered_counts.txt": matrix of counts recovered from the "non_matching_reads.xxx.txt" files for each species (row) in each sample/library (column)
The script also prints the following to the terminal:
Total reads parsed
Estimated number of PhiX (based on reads where neither primer sequence was determined
Number of reads where primer 1 was now identifiable
Number of reads where primer 2 was now identifiable
Number of reads where multiplexing tag was now identifiable
Number of reads where DNA BC was now identifiable
Number of reads where site directed sequence was not identifiable
Total number of reads recovered and recorded to the count matrix

*****************************************************************************************
*****************************************************************************************
Part 2: Removing PCR chimeric reads from count matrix using R
*****************************************************************************************
*****************************************************************************************

***Overview***
During PCR, at some unknown rate, chimeric molecules are formed that contain two molecular
identifiers sourced from two distinct starting molecules. We developed the following 
procedure to estimate the rate of these events on a per sequencing library basis, and 
then correct count data for each library by subtracting out the estimated number of events.
This procedure requires that not all possible combinations of identifiers are present in the
experimental pool. Here, we call the two identifiers BC1 and BC2 (barcode 1 and 2), and 
each unique combination is a DBC (double barcode).

***R markdown file to run the analysis***
*Input*
The R code to run the analysis reads in the "all_counts.txt" and "recovered_counts.txt"
matrices generated from the python scripts in Part 1. It also requires a file 
"data_table_row_key.txt" which has the following columns:
row: row number of count matrix
strain: what combination of query and array guide identifiers the count data for this row represents (i.e. the DBC)
query: name of query guide for this row (i.e. BC1)
array: name of starting pool guide for this row (i.e. BC2)
present: whether the combination of query and array guide for this row was present in the experimental pool (not all possible DBC combinations were assayed in this experiment)
category: whether the strain is a WT (2 non-targeting guides), single (1 gene targeting and 1 non-targeting guide), or double (2 gene targeting guides) mutant
One must also enter the number of sequencing libraries represented in the count matrix.

*QC plotting*
For each library sample plot distributions of counts for rows/species with identifier combinations that are present (green) vs absent (blue) in the experimental pool
On one plot, look at coverage across library samples
For each library, look at representation of each query guide

*Details of chimeric read estimation and subtraction*
This analysis is performed on each library sample separately, as rates of chimeras will vary
depending on the prep. The process consists of the following steps:
1. Get observed frequencies of each individual BC1 by looping through each unique BC1, and summing all counts observed for any species carrying that BC1.
2. Get observed frequencies of each individual BC2 by repeating Step 1 for each unique BC2.
3. Calculate an expectation for DBC frequency by looping through each possible DBC and calculating the product of the observed frequencies of the individual BC1 and BC2 making up that DBC.
4. Model the rate of chimeras by using data from DBCs NOT PRESENT in the experimental pool. Fit a linear regression where the independent variable is expectation calculated in step 3, and dependent variable is the observed number of counts in the sequencing data.
5. Subtract the estimate of chimeric reads from each DBC PRESENT in the pool by using the expectation calculated in step 3 (x value) and the model generated in step 4 to estimate the number of chimeric reads (y value) for the given value of x. Subtract this value from the observed number of reads for this DBC.
6. Plot the data before and after subtraction of chimeric reads.

*Output*
Two RData objects are saved from this analysis:
"seqlibXX_raw_counts.RData": dataframe carrying sequencing counts before subtracting estimated PCR chimeras
"seqlibXX_chimera_normalized_counts.RData": dataframe carrying normalized counts after subtracting estimated PCR chimeras
in both objects, rows are experimental strains and columns are library samples


*****************************************************************************************
*****************************************************************************************
Part 3: Calculating fitness and genetic interaction scores using R
*****************************************************************************************
*****************************************************************************************

***Overview***
For this experiment we tested 17,069 strains in 5 growth conditions. 100 of our strains
carried 2 non-targeting guides, and were used as "WT" controls during normalization. 
We had 3 replicate cultures per condition, and performed amplicon sequencing at 3 time 
points per replicate culture. Before estimating fitness, we first applied starting and 
ending frequency thresholds, to remove strains present at non-observable starting 
frequencies, or strains with such severe fitness defects, they dropped from the pool. 
For each strain, we estimated fitness by its change in relative frequency in log space, 
using generations (as estimated by cell concentration measurements during experimental 
growth) as the x-axis. To estimate interaction scores, we first used the multiplicative 
model to calculate expected double mutant fitness, but subsequently moved to an empirical 
method based on the assumption that most gene pairs to do exhibit a genetic interaction. 
The genetic interaction score is then defined as the deviation of the observed double 
mutant fitness from expectation. 

***R markdown file to run the analysis***
*Input*
This file loads in RData files saved after step 2
The following inputs can also be saved in the first chunk of the markdown file
threshold1: frequency required at time point 1 for analysis of a particular strain in a replicate culture
threshold3: frequency required at time point 1 for analysis of a particular strain in a replicate culture
XXX_gen: for each condition, record the number of generations the pool went through during the fitness assay

*Output*
For each condition, the R script saves two tab-delimited text files (where "XXX" is the condition name):
"XXX_data.txt": contains results for all calculations performed during the analysis
"XXX_short.txt": contains just the necessary strain, fitness, and gi score information necessary for downstream analysis

