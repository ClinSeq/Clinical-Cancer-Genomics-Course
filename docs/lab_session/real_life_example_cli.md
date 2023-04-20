# Real life example of how to use command line features - vt 2023

This is an example of how to use various UNIX/bash/command line features for finding files, assigning variables, parsing and modifying text strings etc. The purpose in the example is to find files containing quality control metrics, and creating a table of certain variables from these files. 

This is not to be run as part of the lab, just a demonstration of a real life example when multiple various command line features are utilized.

First, move to the directory of interest and find the files:

```bash
cd ~/Courses/cancer\ bioinfo\ course\ v2/cmd_line_demo  # move to project directory
pwd  # check present working directory
ls -lhtr  # list files with detailed info (-l), human readable size (-h), in time order (-t), reveresed order with latest last (-r)
find . -name "*hsmetrics.txt"  # find the files in this directory that ends with "txt"
find . -name "*hsmetrics.txt" | wc -l  # count the number of found files
files=($(find . -name "*hsmetrics.txt" | sort))  # assign a sorted list of the files to an array
echo ${files[@]}  # print all items of the array
echo ${#files[@]}  # print number of items of the array
```
Output : 

```
## /Users/rebber/Courses/cancer bioinfo course v2/cmd_line_demo
## total 136
## -rw-r--r--  1 rebber  staff   7.0K Apr 14 12:01 PB-P-HD1-CFDNA-1811-KH20190925-PB20190927.PB.hsmetrics.slurm.log
## -rw-r--r--  1 rebber  staff   4.1K Apr 14 12:01 PB-P-HD1-CFDNA-1811-KH20190925-PB20190927.PB.hsmetrics.txt
## -rw-r--r--  1 rebber  staff   7.8K Apr 14 12:01 PB-P-HD2-CFDNA-1811-KH20190925-PB20190927.PB.hsmetrics.slurm.log
## -rw-r--r--  1 rebber  staff   4.1K Apr 14 12:01 PB-P-HD2-CFDNA-1811-KH20190925-PB20190927.PB.hsmetrics.txt
## -rw-r--r--  1 rebber  staff   8.1K Apr 14 12:01 PB-P-HD3-CFDNA-1811-KH20190925-PB20190927.PB.hsmetrics.slurm.log
## -rw-r--r--  1 rebber  staff   4.1K Apr 14 12:01 PB-P-HD3-CFDNA-1811-KH20190925-PB20190927.PB.hsmetrics.txt
## -rw-r--r--  1 rebber  staff   5.3K Apr 14 12:01 PB-P-HD6-CFDNA-1811-KH20190925-PB20190927.PB.hsmetrics.slurm.log
## -rw-r--r--  1 rebber  staff   4.1K Apr 14 12:01 PB-P-HD6-CFDNA-1811-KH20190925-PB20190927.PB.hsmetrics.txt
## ./PB-P-HD6-CFDNA-1811-KH20190925-PB20190927.PB.hsmetrics.txt
## ./PB-P-HD3-CFDNA-1811-KH20190925-PB20190927.PB.hsmetrics.txt
## ./PB-P-HD2-CFDNA-1811-KH20190925-PB20190927.PB.hsmetrics.txt
## ./PB-P-HD1-CFDNA-1811-KH20190925-PB20190927.PB.hsmetrics.txt
##        4
## ./PB-P-HD1-CFDNA-1811-KH20190925-PB20190927.PB.hsmetrics.txt ./PB-P-HD2-CFDNA-1811-KH20190925-PB20190927.PB.hsmetrics.txt ./PB-P-HD3-CFDNA-1811-KH20190925-PB20190927.PB.hsmetrics.txt ./PB-P-HD6-CFDNA-1811-KH20190925-PB20190927.PB.hsmetrics.txt
## 4
```

Now we have the files of interest in the array `$files`.
We can use this array for running further operations on the files.

Inspect the files:
```bash
ls -lh ${files[@]}  # list files with detailed view and human readable sizes
wc -l ${files[@]}  # count number of rows of all files
less -SN ${files[@]}  # open files in simple file viewer, with lines chopped (-S) and numbered (-N), next file can be accessed with :n, more options available in help page by pressing "h"
head -n10 ${files[0]}  # view 10 first lines of first file
tail -n3 ${files[0]}  # view 3 last lines of first file
```

Output :
```
## -rw-r--r--  1 rebber  staff   4.1K Apr 14 12:01 ./PB-P-HD1-CFDNA-1811-KH20190925-PB20190927.PB.hsmetrics.txt
## -rw-r--r--  1 rebber  staff   4.1K Apr 14 12:01 ./PB-P-HD2-CFDNA-1811-KH20190925-PB20190927.PB.hsmetrics.txt
## -rw-r--r--  1 rebber  staff   4.1K Apr 14 12:01 ./PB-P-HD3-CFDNA-1811-KH20190925-PB20190927.PB.hsmetrics.txt
## -rw-r--r--  1 rebber  staff   4.1K Apr 14 12:01 ./PB-P-HD6-CFDNA-1811-KH20190925-PB20190927.PB.hsmetrics.txt
##      213 ./PB-P-HD1-CFDNA-1811-KH20190925-PB20190927.PB.hsmetrics.txt
##      213 ./PB-P-HD2-CFDNA-1811-KH20190925-PB20190927.PB.hsmetrics.txt
##      213 ./PB-P-HD3-CFDNA-1811-KH20190925-PB20190927.PB.hsmetrics.txt
##      213 ./PB-P-HD6-CFDNA-1811-KH20190925-PB20190927.PB.hsmetrics.txt
##      852 total
## ## htsjdk.samtools.metrics.StringHeader
## # CollectHsMetrics BAIT_INTERVALS=[/nfs/PROBIO/new_design_eval/targets/all_baits_twist.interval_list] BAIT_SET_NAME=PB TARGET_INTERVALS=[/nfs/PROBIO/new_design_eval/targets/all_customized_targets_twist.interval_list] INPUT=/dev/stdin OUTPUT=metrics/realignment_cust_targets/HsMetrics/MarkDuplicates/PB/PB-P-HD1-CFDNA-1811-KH20190925-PB20190927.PB.hsmetrics.txt METRIC_ACCUMULATION_LEVEL=[ALL_READS] TMP_DIR=[/scratch/tmp/rebber] REFERENCE_SEQUENCE=/nfs/PROBIO/autoseq-genome/genome/human_g1k_v37_decoy.fasta    NEAR_DISTANCE=250 MINIMUM_MAPPING_QUALITY=20 MINIMUM_BASE_QUALITY=20 CLIP_OVERLAPPING_READS=true INCLUDE_INDELS=false COVERAGE_CAP=200 SAMPLE_SIZE=10000 ALLELE_FRACTION=[0.001, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.3, 0.5] VERBOSITY=INFO QUIET=false VALIDATION_STRINGENCY=STRICT COMPRESSION_LEVEL=5 MAX_RECORDS_IN_RAM=500000 CREATE_INDEX=false CREATE_MD5_FILE=false GA4GH_CLIENT_SECRETS=client_secrets.json USE_JDK_DEFLATER=false USE_JDK_INFLATER=false
## ## htsjdk.samtools.metrics.StringHeader
## # Started on: Tue May 26 17:12:11 CEST 2020
## 
## ## METRICS CLASS picard.analysis.directed.HsMetrics
## BAIT_SET BAIT_TERRITORY  BAIT_DESIGN_EFFICIENCY  ON_BAIT_BASES   NEAR_BAIT_BASES OFF_BAIT_BASES  PCT_SELECTED_BASES  PCT_OFF_BAIT    ON_BAIT_VS_SELECTED MEAN_BAIT_COVERAGE  PCT_USABLE_BASES_ON_BAIT    PCT_USABLE_BASES_ON_TARGET  FOLD_ENRICHMENT HS_LIBRARY_SIZE HS_PENALTY_10X  HS_PENALTY_20X  HS_PENALTY_30X  HS_PENALTY_40X  HS_PENALTY_50X  HS_PENALTY_100X TARGET_TERRITORY    GENOME_SIZE TOTAL_READS PF_READS    PF_BASES    PF_UNIQUE_READS PF_UQ_READS_ALIGNED PF_BASES_ALIGNED    PF_UQ_BASES_ALIGNED ON_TARGET_BASES PCT_PF_READS    PCT_PF_UQ_READS PCT_PF_UQ_READS_ALIGNED MEAN_TARGET_COVERAGE    MEDIAN_TARGET_COVERAGE  MAX_TARGET_COVERAGE ZERO_CVG_TARGETS_PCT    PCT_EXC_DUPE    PCT_EXC_ADAPTER PCT_EXC_MAPQ    PCT_EXC_BASEQ   PCT_EXC_OVERLAP PCT_EXC_OFF_TARGET  FOLD_80_BASE_PENALTY    PCT_TARGET_BASES_1X PCT_TARGET_BASES_2X PCT_TARGET_BASES_10X    PCT_TARGET_BASES_20X    PCT_TARGET_BASES_30X    PCT_TARGET_BASES_40X    PCT_TARGET_BASES_50X    PCT_TARGET_BASES_100X   AT_DROPOUT  GC_DROPOUT  HET_SNP_SENSITIVITY HET_SNP_Q   SAMPLE  LIBRARY READ_GROUP
## PB   290831  0.91208 682291515   195630134   1507472292  0.368041    0.631959    0.777167    2346.006839 0.281894    0.169875    3085.649544     0   0   0   0   0   0   265261  3137454505  16795512    16795512    2420380851  16795512    16679404    2385393941  2385393941  411162960   1   1   0.993087    1550.031705 200 7219    0   0   0   0.101408    0.004983    0.298104    0.423367    7.750159    0.999612    0.999585    0.999582    0.999529    0.999457    0.999416    0.999359    0.99928 7.694876    0.179254    0.999597    34          
## 
## ## HISTOGRAM java.lang.Integer
## 199  0   0
## 200  265028  0

```

Extract sample names from file names:
```bash
ls -lh ${files[@]}

# Method 1: all at once
# print content of array with line break as separator to get one item on each line
# use cut to get the second field of each line when separator is "/"
# use cut to get the first field of each line when separator is "."
printf '%s\n' ${files[@]} | cut -f 2 -d "/" | cut -f 1 -d "."

# Method 2: one at a time in for loop
# loop over each of the items in the $files array
# get the basename for the file, i.e. just the file name not full path 
# use sed to remove everything matching pattern ".PB.hsmetrics.txt" (actually replace it with nothing)
for f in ${files[@]}; do
    basename $f | sed 's/.PB.hsmetrics.txt//g'
done
```
Output:
```
## -rw-r--r--  1 rebber  staff   4.1K Apr 14 12:01 ./PB-P-HD1-CFDNA-1811-KH20190925-PB20190927.PB.hsmetrics.txt
## -rw-r--r--  1 rebber  staff   4.1K Apr 14 12:01 ./PB-P-HD2-CFDNA-1811-KH20190925-PB20190927.PB.hsmetrics.txt
## -rw-r--r--  1 rebber  staff   4.1K Apr 14 12:01 ./PB-P-HD3-CFDNA-1811-KH20190925-PB20190927.PB.hsmetrics.txt
## -rw-r--r--  1 rebber  staff   4.1K Apr 14 12:01 ./PB-P-HD6-CFDNA-1811-KH20190925-PB20190927.PB.hsmetrics.txt
## PB-P-HD1-CFDNA-1811-KH20190925-PB20190927
## PB-P-HD2-CFDNA-1811-KH20190925-PB20190927
## PB-P-HD3-CFDNA-1811-KH20190925-PB20190927
## PB-P-HD6-CFDNA-1811-KH20190925-PB20190927
## PB-P-HD1-CFDNA-1811-KH20190925-PB20190927
## PB-P-HD2-CFDNA-1811-KH20190925-PB20190927
## PB-P-HD3-CFDNA-1811-KH20190925-PB20190927
## PB-P-HD6-CFDNA-1811-KH20190925-PB20190927

```
Other combinations of the tools in the two methods can also be used to get the same results. 
For most tasks there are many possible ways to solve them on the command line. Just use whatever tools you are familiar with, and when there is something you can't figure out how to do, google offers tons of knowledge and ideas.

Some more things that can be done with a for loop:
```bash
for f in ${files[@]}; do
    # assign output of commands to a variable with the $() construction:
    samp=$(basename $f | sed 's/.PB.hsmetrics.txt//g')
    # check conditions with if, here if $samp varaible matches some regex:
    if [[ $samp =~ HD1 ]]; then
        echo "This is HD1 sample:"
    elif [[ $samp =~ HD2 ]]; then
        echo "This is HD2 sample:"
    else 
        echo "This another sample:"
    fi
    echo $samp
    # search for the sample name in the file content with grep, -c gives the count of hits
    hits=$(grep -c "$samp" $f)
    echo "File $f contains sample name $samp $hits times"
    echo  # empty line
done
```
Output:
```
## This is HD1 sample:
## PB-P-HD1-CFDNA-1811-KH20190925-PB20190927
## File ./PB-P-HD1-CFDNA-1811-KH20190925-PB20190927.PB.hsmetrics.txt contains sample name PB-P-HD1-CFDNA-1811-KH20190925-PB20190927 1 times
## 
## This is HD2 sample:
## PB-P-HD2-CFDNA-1811-KH20190925-PB20190927
## File ./PB-P-HD2-CFDNA-1811-KH20190925-PB20190927.PB.hsmetrics.txt contains sample name PB-P-HD2-CFDNA-1811-KH20190925-PB20190927 1 times
## 
## This another sample:
## PB-P-HD3-CFDNA-1811-KH20190925-PB20190927
## File ./PB-P-HD3-CFDNA-1811-KH20190925-PB20190927.PB.hsmetrics.txt contains sample name PB-P-HD3-CFDNA-1811-KH20190925-PB20190927 1 times
## 
## This another sample:
## PB-P-HD6-CFDNA-1811-KH20190925-PB20190927
## File ./PB-P-HD6-CFDNA-1811-KH20190925-PB20190927.PB.hsmetrics.txt contains sample name PB-P-HD6-CFDNA-1811-KH20190925-PB20190927 1 times
```
To create a table with the info from specific rows and columns we can use the awk program. In this case we will take columns MEAN_TARGET_COVERAGE (number 34) and FOLD_ENRICHMENT (number 13). The relevant data is on line 8, with the header on line 7. 
NR: row number in total; FNR: row number in file.
">" writes the output to a new file
column -t aligns the output columns nicely
```bash
awk -F "\t" -v OFS="\t" 'NR==7 {print "file", $34, $13}; FNR==8 {print FILENAME, $34, $13}' \
    ${files[@]} > qc_summary.txt
ls -lhtr 
cat qc_summary.txt
cat qc_summary.txt | column -t
```
Output:
```
## total 144
## -rw-r--r--  1 rebber  staff   7.0K Apr 14 12:01 PB-P-HD1-CFDNA-1811-KH20190925-PB20190927.PB.hsmetrics.slurm.log
## -rw-r--r--  1 rebber  staff   4.1K Apr 14 12:01 PB-P-HD1-CFDNA-1811-KH20190925-PB20190927.PB.hsmetrics.txt
## -rw-r--r--  1 rebber  staff   7.8K Apr 14 12:01 PB-P-HD2-CFDNA-1811-KH20190925-PB20190927.PB.hsmetrics.slurm.log
## -rw-r--r--  1 rebber  staff   4.1K Apr 14 12:01 PB-P-HD2-CFDNA-1811-KH20190925-PB20190927.PB.hsmetrics.txt
## -rw-r--r--  1 rebber  staff   8.1K Apr 14 12:01 PB-P-HD3-CFDNA-1811-KH20190925-PB20190927.PB.hsmetrics.slurm.log
## -rw-r--r--  1 rebber  staff   4.1K Apr 14 12:01 PB-P-HD3-CFDNA-1811-KH20190925-PB20190927.PB.hsmetrics.txt
## -rw-r--r--  1 rebber  staff   5.3K Apr 14 12:01 PB-P-HD6-CFDNA-1811-KH20190925-PB20190927.PB.hsmetrics.slurm.log
## -rw-r--r--  1 rebber  staff   4.1K Apr 14 12:01 PB-P-HD6-CFDNA-1811-KH20190925-PB20190927.PB.hsmetrics.txt
## -rw-r--r--  1 rebber  staff   382B Apr 14 17:56 qc_summary.txt
## file MEAN_TARGET_COVERAGE    FOLD_ENRICHMENT
## ./PB-P-HD1-CFDNA-1811-KH20190925-PB20190927.PB.hsmetrics.txt 1550.031705 3085.649544
## ./PB-P-HD2-CFDNA-1811-KH20190925-PB20190927.PB.hsmetrics.txt 2283.624973 3389.760665
## ./PB-P-HD3-CFDNA-1811-KH20190925-PB20190927.PB.hsmetrics.txt 2599.327089 3688.749254
## ./PB-P-HD6-CFDNA-1811-KH20190925-PB20190927.PB.hsmetrics.txt 1173.933741 4831.030898
## file                                                          MEAN_TARGET_COVERAGE  FOLD_ENRICHMENT
## ./PB-P-HD1-CFDNA-1811-KH20190925-PB20190927.PB.hsmetrics.txt  1550.031705           3085.649544
## ./PB-P-HD2-CFDNA-1811-KH20190925-PB20190927.PB.hsmetrics.txt  2283.624973           3389.760665
## ./PB-P-HD3-CFDNA-1811-KH20190925-PB20190927.PB.hsmetrics.txt  2599.327089           3688.749254
## ./PB-P-HD6-CFDNA-1811-KH20190925-PB20190927.PB.hsmetrics.txt  1173.933741           4831.030898
```