---
title: Real life example of how to use command line features - vt 2025

---

# Real life example of how to use command line features - vt 2025

This is an example of how to use various UNIX/bash/command line features for finding files, assigning variables, parsing and modifying text strings etc. The purpose in the example is to find files containing quality control metrics, and creating a table of certain variables from these files. 


First, move to the directory of interest and find the files:

```bash
cd /nfs/course/inputs/cmd_line_demo  # move to project directory
pwd  # check present working directory
ls -lhtr  # list files with detailed info (-l), human readable size (-h), in time order (-t), reveresed order with latest last (-r)
find . -name "*hsmetrics.txt"  # find the files in this directory that ends with "txt"
find . -name "*hsmetrics.txt" | wc -l  # count the number of found files
files=($(find . -name "*hsmetrics.txt" | sort))  # assign a sorted list of the files to an array
echo ${files[@]}  # print all items of the array
echo ${#files[@]}  # print number of items of the array
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
Other combinations of the tools in the two methods can also be used to get the same results. 
For most tasks there are many possible ways to solve them on the command line. Just use whatever tools you are familiar with, and when there is something you can't figure out how to do, ChatGPT and google offers tons of knowledge and ideas.

Some more things that can be done with a for loop:
```bash
for f in ${files[@]}; do
    # assign output of commands to a variable with the $() construction:
    samp=$(basename $f | sed 's/.PB.hsmetrics.txt//g')  # get the base name and remove suffix
    # check conditions with if, here if $samp varaible matches some regex:
    if [[ $samp =~ HD1 ]]; then
        echo "This is HD1 sample:"
    elif [[ $samp =~ HD2 ]]; then
        echo "This is HD2 sample:"
    else 
        echo "This is another sample:"
    fi
    echo $samp
    # search for the sample name in the file content with grep, -c gives the count of hits
    hits=$(grep -c "$samp" $f)
    echo "File $f contains sample name $samp $hits times"
    echo  # empty line
done
```

To create a table with the info from specific rows and columns we can use the awk program. In this case we will take columns MEAN_TARGET_COVERAGE (number 34) and FOLD_ENRICHMENT (number 13). The relevant data is on line 8, with the header on line 7. 
Change the "my_user" to your user name.
NR: row number in total; FNR: row number in file.
">" writes the output to a new file
column -t aligns the output columns nicely
```bash
awk -F "\t" -v OFS="\t" 'NR==7 {print "file", $34, $13}; FNR==8 {print FILENAME, $34, $13}' \
    ${files[@]} > /nfs/course/students/my_user/qc_summary.txt
ls -lhtr 
cat /nfs/course/students/my_user/qc_summary.txt
cat /nfs/course/students/my_user/qc_summary.txt | column -t
```

The shell is not the best for doing maths, but some basic things can be acheived with e.g. awk.
Calculate mean MEAN_TARGET_COVERAGE:
```bash
# for all row numbers > 1, take the value of the second column and add to the total
# then after all lines have been parsed, print the total divided by number of rows
awk -F "\t" 'NR>1 {tot=tot+$2; rows=rows+1}; END {print tot, rows, tot/rows}' \
    /nfs/course/students/my_user/qc_summary.txt
```

Can you come up with a way to get the median MEAN_TARGET_COVERAGE, and explain what each part of it does? 


Now, let's look at the log files.
```bash
ls -lh *log
less -SN *log

# get the mean run time for the tool
grep "done. Elapsed time" *log  # find the lines with the total time info
grep "done. Elapsed time" *log | cut -f 11 -d " " # get the 11th space-separated field, where the time is
grep "done. Elapsed time" *log | cut -f 11 -d " " | \
    awk '{tot=tot+$0; rows=rows+1}; END {print tot/rows}' # calculate the mean

# get the total number of processed records for each sample
for file in *.log; do
  count=$(grep "CollectHsMetrics.*Processed" "$file" | tail -n 1 | awk '{print $6}')
  echo "$file: $count"
done
```

Can you now get the TOTAL_READS for each sample from the txt files?
Tip: `tr "\t" "\n"` can be used to convert the tabs in the column header line to linebreaks, to easier assess which column number the TOTAL_READS column has.
