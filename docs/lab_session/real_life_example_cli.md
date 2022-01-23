# Real life example of how to use command line features

This is an example of how to use various UNIX/bash/command line features for finding files, assigning variables, parsing and modifying text strings etc. The purpose in the example is to run programs that creates bam files from fastq files and then calculates some quality control metrics. This is done as a loop over multiple samples, where the same commands are applied to every sample. At the end, a table of the resulting metrics of interest is created.

This is not to be run as part of the lab, just a demonstration of a real life example when multiple various command line features are utilized.

```bash
cd /path/to/project_dir  # move to project directory
ref=/path/to/human_g1k_v37_decoy.fasta  # assign path to file as variable
bwaRef=/path/to/bwa/human_g1k_v37_decoy.fasta  # assign path to file as variable
tmpDir=/scratch/tmp/rebber  # assign path to directory as variable
inbox=/path/to/fastq_dir  # assign path to directory as variable
# find fastq files of interest, 
# extract the sample names from found paths with basename, 
# sort and assign results to array (one sample per array entry):
samples=($(find $inbox -maxdepth 1 -mindepth 1 -newermt 2021-01-11 \
    ! -newermt 2021-01-12 -name "PB-P-HD*" | xargs -I {} basename {} | sort))
echo ${samples[@]}  # display all entries of array
echo ${#samples[@]}  # display number of entries in array
mkdir -p bams  # create output directory

# for loop over the samples:
for samp in ${samples[@]}; do
    # get correct design files:
    # extract info from sample name (7th field, everything before "2020"):
    design=$(echo $samp | cut -f 7 -d "-" | sed 's/2020.*//g')
    # --> e.g. PB-P-HD181117-N-7minfrag-KH20201002-PS20201005 --> PS
    # assign variables depending on the value of the $design variable:
    if [[ $design == "PB" ]]; then
        baits=/path/to/PB.baits.interval_list
        targets=/path/to/PB.targets.interval_list
    elif [[ $design == "PS" ]]; then
        baits=/path/to/PS.baits.interval_list
        targets=/path/to/PS.targets.interval_list
    else
        echo "Unknown design $design, skipping sample $samp"  # print message
        continue  # skip the rest of the loop for this value of $samp
    fi
    
    # skewer (trimming of fastq files)
    threads=6  # set variable to use for this tool
    dep_jids=""  # reset list of jobs to depend on (from previous iterations)
    mkdir -p skewer/$samp  #create output dir
    # find the "first in pair" fastq files for this sample:
    fq1s=($(find $inbox/$samp -name "*_1.fastq.gz" | sort))
    
    for fq1 in ${fq1s[@]}; do  # loop over the found fastq files
        # get the basename of the current fastq file,
        # sed replaces the pattern (to nothing), 
        # basename removes the path to the file:
        fqbase=$(echo $fq1 | sed 's/_1.fastq.gz//g'| xargs -I {} basename {})
        # find the corresponding second in pair fastq file, 
        # for the current first in pair fastq file:
        fq2=$(find $inbox/$samp -name ${fqbase}_2.fastq.gz)
        # assign output path and base name:
        outPrefix=skewer/$samp/skewer.$fqbase
        # use sbatch for submitting a job to the slurm queue, set options for it, 
        # then use "here document" to create script on the fly, 
        # that is submitted to slurm
        sbatch -o $outPrefix.slurm.log -n $threads -J skewer_$fqbase <(cat << EOF
#!/bin/bash  # initiate bash script with shebang
# activate a bash profile which activates access to e.g. a tool set:
. /path/to/.bash_profile
skewer -z -t $threads --quiet -o $outPrefix $fq1 $fq2  # run the desired program 
EOF
)

        # list jobs with squeue command, 
        # get job ID of previously submitted job:
        jid=$(squeue -o "%j %i" | grep skewer_$fqbase | cut -f 2 -d " ")
        dep_jids="${dep_jids}:$jid"  # append to list of job IDs from this sample
    done  # end of loop over fastq files
    
    # assign variable with job IDs that the next job is going to be dependent on:
    dep="--dependency=afterok$dep_jids"
    
    # mapping
    threads=12
    fastqDir=skewer/$samp
    outBam=bams/raw/$samp.mapped.bam
    # submit to slurm queue, with dependency on the previous jobs,
    # now a pre-written script which takes positional arguments
    sbatch -o $outBam.slurm.log -n $threads -J ${samp}_fastq2mappedBam $dep \
        /path/to/scripts/fastq2mappedBam.sh $fastqDir $outBam $ref $bwaRef \
        $threads $tmpDir
    
    # get jobid of submitted job
    jid=$(squeue -o "%j %i" | grep ${samp}_fastq2mappedBam | cut -f 2 -d " ")
    dep="--dependency=afterok:$jid"  # dependency
    
    # HsMetrics
    threads=1
    inbam=bams/raw/${samp}.mapped.bam
    outdir1=metrics/HsMetrics/$design
    mkdir -p $outdir1
    sbatch -o $outdir1/$samp.$design.hsmetrics.slurm.log -J ${samp}_HsMetrics \
        -n $threads $dep /path/to/scripts/runHsmetrics.sh $inbam $outdir1 \
        $targets $baits $design
    
done  # end of loop over samples

squeue | less -SN  # check the slurm queue
ls -lhtr skewer/*  # list content of skewer output dir
# find the log files and view them with less:
find skewer/ -newermt 2021-01-11 -name "*slurm.log" | xargs less -SN
# grep for any error message in the log files, 
# remove irrelevant hits, and view the remaining results:
grep -iP "(error|fail|exception|traceback|cancel|slurm|retry|No such file or directory)" \
    $(find skewer/* -newermt 2021-01-11 -name "*.slurm.log" | sort) | \
    grep -v "error ratio allowed" | less -SN
# more checks of output files can be done...

# make a table of columns 34 (MEAN_TARGET_COVERAGE) and 13 (FOLD_ENRICHMENT)
# for all the samples
# the relevant data is on line 8, with the header on line 7
# use awk to extract the relevant lines and columns
# NR: row number in total; FNR: row number in file
# column -t aligns the columns nicely
awk -F "\t" -v OFS="\t" 'NR==7 || FNR==8 {print FILENAME, $34, $13}' \
    metrics/HsMetrics/*/*txt | column -t
```
Resulting QC metrics table:
```
metrics/HsMetrics/PB/PB-P-HD181110-CFDNA-benchmarking-KH20201002-PB20201005.PB.hsmetrics.txt  MEAN_TARGET_COVERAGE  FOLD_ENRICHMENT
metrics/HsMetrics/PB/PB-P-HD181110-CFDNA-benchmarking-KH20201002-PB20201005.PB.hsmetrics.txt  36390.640776          6267.968261
metrics/HsMetrics/PB/PB-P-HD181115-CFDNA-benchmarking-KH20201002-PB20201005.PB.hsmetrics.txt  34498.125895          6251.155829
metrics/HsMetrics/PS/PB-P-HD181110-N-4minfrag-KH20201002-PS20201005.PS.hsmetrics.txt          64569.379242          16489.157949
metrics/HsMetrics/PS/PB-P-HD181115-N-4minfrag-KH20201002-PS20201005.PS.hsmetrics.txt          67690.957504          17337.593473
```
