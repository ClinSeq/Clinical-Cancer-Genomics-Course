# AWK crash course for Bioinformatics

## An introduction on how to analyze table data without MS Excel

## What is AWK?

AWK is an interpreted programming language designed for text processing and typically used as a data extraction and reporting tool.

The AWK language is a data-driven scripting language consisting of a set of actions to be taken against streams of textual data - either run directly on files or used as part of a pipeline - for purposes of extracting or transforming text, such as producing formatted reports. The language extensively uses the string datatype, associative arrays (that is, arrays indexed by key strings), and regular expressions.

AWK has a limited intended application domain, and was especially designed to support one-liner programs.

It is a standard feature of most Unix-like operating systems.

source: [Wikipedia](https://en.wikipedia.org/wiki/AWK)

## Why awk? 

You can replace a pipeline of 'stuff | grep | sed | cut...' with a single call to awk. For a simple script, most of the timelag is in loading these apps into memory, and it's much faster to do it all with one. This is ideal for something like an openbox pipe menu where you want to generate something on the fly. You can use awk to make a neat one-liner for some quick job in the terminal, or build an awk section into a shell script. 

### Simple AWK commands

An AWK program consists of a sequence of pattern-action statements and optional function definitions. It processes text files. AWK is a line oriented language. It divides a file into lines called records. Each line is broken up into a sequence of fields. The fields are accessed by special variables: $1 reads the first field, $2 the second and so on. The $0 variable refers to the whole record.

The structure of an AWK program has the following form:
```bash
pattern { action }
```
The pattern is a test that is performed on each of the records. If the condition is met then the action is performed. Either pattern or action can be omitted, but not both. The default pattern matches each line and the default action is to print the record.

```bash
awk -f program-file [file-list]
awk program [file-list]
```

An AWK program can be run in two basic ways: a) the program is read from a separate file; the name of the program follows the -f option, b) the program is specified on the command line enclosed by quote characters.

### Print all the lines from a file

By default, awk prints all lines of a file , so to print every line of above created file use below command 

```bash
awk '{print}' file
```
Note: In awk command ‘{print}’ is used print all fields along with their values.

### Print only specific field like 2nd & 3rd

In awk command, we use $ (dollar) symbol followed by field number to prints field values. 

```bash
awk -F "," '{print $2, $3}' file
```
In the above command we have used the option  -F “,”  which specifies that comma (,) is the field separator in the file. This is usually good practice when dealing with tables where the separators between each column could be single or multiple white-space (" ") or tab ("\t") or a colon (":") or a semicolon (";").

### AWK one-liners

AWK one-liners are simple one-shot programs run from the command line. Let us have the following text file: words.txt

We want to print all words included in the words.txt file that are longer than five characters.

```bash
wget https://course-cg-5534.s3.amazonaws.com/awk_exercise/words.txt -O words.txt
awk 'length($1) > 5 {print $0}' words.txt

storeroom
existence
ministerial
falcon
bookworm
bookcase
```
The AWK program is placed between two single quote characters. The first is the pattern; we specify that the length of the record is greater that five. The length function returns the length of the string. The $1 variable refers to the first field of the record; in our case there is only one field per record. The action is placed between curly brackets.

```bash
awk 'length($1) > 5' words.txt
storeroom
existence
ministerial
falcon
bookworm
bookcase
```
As we have specified earlier, the action can be omitted. In such a case a default action is performed — printing of the whole record.

```bash
awk 'length($1) == 3' words.txt
cup
sky
top
war
```
We print all words that have three characters.

```bash
awk '!(length($1) == 3)' words.txt
storeroom
tree
store
book
cloud
existence
ministerial
falcon
town
bookworm
bookcase
```
With the ! operator, we can negate the condition; we print all lines that do not have three characters.

```bash
awk '(length($1) == 3) || (length($1) == 4)' words.txt
tree
cup
book
town
sky
top
war
```
Next we apply conditions on numbers.

We have a file with scores of students - scores.txt.
```bash
wget https://course-cg-5534.s3.amazonaws.com/awk_exercise/scores.txt -O scores.txt
awk '$2 >= 90 { print $0 }' scores.txt
Lucia 95
Joe 92
Sophia 90
```
We print all students with scores 90+.

```bash
awk '$2 >= 90 { print }' scores.txt
Lucia 95
Joe 92
Sophia 90
```
If we omit an argument for the print function, the $0 is assumed.

```bash
awk '$2 >= 90' scores.txt
Lucia 95
Joe 92
Sophia 90
```
A missing { action } means print the matching line.

```bash
awk '{ if ($2 >= 90) print }' scores.txt
Lucia 95
Joe 92
Sophia 90
```
Instead of a pattern, we can also use an if condition in the action.

```bash
awk '{sum += $2} END { printf("The average score is %.2f\n", sum/NR) }' scores.txt
```
The average score is 77.56

This command calculates the average score. In the action block, we calculate the sum of scores. In the END block, we print the average score. We format the output with the built-in printf function. The %.2f is a format specifier; each specifier begins with the % character. The .2 is the precision -- the number of digits after the decimal point. The f expects a floating point value. The \n is not a part of the specifier; it is a newline character. It prints a newline after the string is shown on the terminal.

## AWK working with pipes

AWK can receive input and send output to other commands via the pipe.

```bash
echo -e "1 2 3 5\n2 2 3 8" | awk '{print $(NF)}'
5
8
```

In this case, AWK receives output from the echo command. It prints the values of last column.

```bash
awk -F: '$7 ~ /bash/ {print $1}' /etc/passwd | wc -l
3
```

Here, the AWK program sends data to the wc command via the pipe. In the AWK program, we find out those users who use bash. Their names are passed to the wc command which counts them. In our case, there are three users using bash.

### AWK for Bioinformatics - looking at transcriptome data

You can find a lot of online tutorials, but here we will try out a few steps which show how a bioinformatician analyses a **GTF file** using awk.

GTF is a special file format that contains information about different regions of the genome and their associated annotations. More on that here - [Ensembl File Formats-GTF](https://www.ensembl.org/info/website/upload/gff.html).

```bash
wget https://course-cg-5534.s3.amazonaws.com/awk_exercise/transcriptome.gtf -O transcriptome.gtf

head transcriptome.gtf | less -S

##description: evidence-based annotation of the human genome (GRCh37), version 18 (Ensembl 73)
##provider: GENCODE
##contact: gencode@sanger.ac.uk
##format: gtf
##date: 2013-09-02
chr1    HAVANA  exon    173753  173862  .   -   .   gene_id "ENSG00000241860.2"; transcript_id "ENST00000466557.2"; gene_type "processed_transcript"; gene_status "NOVEL"; gene_name "RP11-34P13.13"; transcript_type "lincRNA"; transcript_status "KNOWN"; transcript_name "RP11-34P13.13-001"; exon_number 1;  exon_id "ENSE00001947154.2";  level 2; tag "not_best_in_genome_evidence"; havana_gene "OTTHUMG00000002480.3"; havana_transcript "OTTHUMT00000007037.2";
chr1    HAVANA  transcript  1246986 1250550 .   -   .   gene_id "ENSG00000127054.14"; transcript_id "ENST00000478641.1"; gene_type "protein_coding"; gene_status "KNOWN"; gene_name "CPSF3L"; transcript_type "retained_intron"; transcript_status "KNOWN"; transcript_name "CPSF3L-006"; level 2; havana_gene "OTTHUMG00000003330.11"; havana_transcript "OTTHUMT00000009365.1";
chr1    HAVANA  CDS 1461841 1461911 .   +   0   gene_id "ENSG00000197785.9"; transcript_id "ENST00000378755.5"; gene_type "protein_coding"; gene_status "KNOWN"; gene_name "ATAD3A"; transcript_type "protein_coding"; transcript_status "KNOWN"; transcript_name "ATAD3A-003"; exon_number 13;  exon_id "ENSE00001664426.1";  level 2; tag "basic"; tag "CCDS"; ccdsid "CCDS31.1"; havana_gene "OTTHUMG00000000575.6"; havana_transcript "OTTHUMT00000001365.1";
chr1    HAVANA  exon    1693391 1693474 .   -   .   gene_id "ENSG00000008130.11"; transcript_id "ENST00000341991.3"; gene_type "protein_coding"; gene_status "KNOWN"; gene_name "NADK"; transcript_type "protein_coding"; transcript_status "KNOWN"; transcript_name "NADK-002"; exon_number 3;  exon_id "ENSE00003487616.1";  level 2; tag "basic"; tag "CCDS"; ccdsid "CCDS30565.1"; havana_gene "OTTHUMG00000000942.5"; havana_transcript "OTTHUMT00000002768.1";
chr1    HAVANA  CDS 1688280 1688321 .   -   0   gene_id "ENSG00000008130.11"; transcript_id "ENST00000497186.1"; gene_type "protein_coding"; gene_status "KNOWN"; gene_name "NADK"; transcript_type "nonsense_mediated_decay"; transcript_status "KNOWN"; transcript_name "NADK-008"; exon_number 2;  exon_id "ENSE00001856899.1";  level 2; tag "mRNA_start_NF"; tag "cds_start_NF"; havana_gene "OTTHUMG00000000942.5"; havana_transcript "OTTHUMT00000002774.3";
```
> The transcriptome has 9 columns. The first 8 are separated by tabs and look reasonable (chromosome, annotation source, feature type, start, end, score, strand, and phase), the last one is kind of hairy: it is made up of key-value pairs separated by semicolons, some fields are mandatory and others are optional, and the values are surrounded in double quotes. That’s no way to live a decent life. (*text copied from the source*)

let's get only the lines that have `gene` in the 3^th^ column.

``` bash
$ awk -F '$3 == "gene"' transcriptome.gtf | head | less -S

chr1    HAVANA  gene    11869   14412   .   +   .   gene_id "ENSG00000223972.4"; transcript_id "ENSG00000223972.4"; gene_type "pseudogene"; gene_status "KNOWN"; gene_name "DDX11L1"; transcript_type "pseudogene"; transcript_status "KNOWN"; transcript_name "DDX11L1"; level 2; havana_gene "OTTHUMG00000000961.2";
chr1    HAVANA  gene    14363   29806   .   -   .   gene_id "ENSG00000227232.4"; transcript_id "ENSG00000227232.4"; gene_type "pseudogene"; gene_status "KNOWN"; gene_name "WASH7P"; transcript_type "pseudogene"; transcript_status "KNOWN"; transcript_name "WASH7P"; level 2; havana_gene "OTTHUMG00000000958.1";
chr1    HAVANA  gene    29554   31109   .   +   .   gene_id "ENSG00000243485.2"; transcript_id "ENSG00000243485.2"; gene_type "lincRNA"; gene_status "NOVEL"; gene_name "MIR1302-11"; transcript_type "lincRNA"; transcript_status "NOVEL"; transcript_name "MIR1302-11"; level 2; havana_gene "OTTHUMG00000000959.2";
chr1    HAVANA  gene    34554   36081   .   -   .   gene_id "ENSG00000237613.2"; transcript_id "ENSG00000237613.2"; gene_type "lincRNA"; gene_status "KNOWN"; gene_name "FAM138A"; transcript_type "lincRNA"; transcript_status "KNOWN"; transcript_name "FAM138A"; level 2; havana_gene "OTTHUMG00000000960.1";
chr1    HAVANA  gene    52473   54936   .   +   .   gene_id "ENSG00000268020.2"; transcript_id "ENSG00000268020.2"; gene_type "pseudogene"; gene_status "KNOWN"; gene_name "OR4G4P"; transcript_type "pseudogene"; transcript_status "KNOWN"; transcript_name "OR4G4P"; level 2; havana_gene "OTTHUMG00000185779.1";
chr1    HAVANA  gene    62948   63887   .   +   .   gene_id "ENSG00000240361.1"; transcript_id "ENSG00000240361.1"; gene_type "pseudogene"; gene_status "KNOWN"; gene_name "OR4G11P"; transcript_type "pseudogene"; transcript_status "KNOWN"; transcript_name "OR4G11P"; level 2; havana_gene "OTTHUMG00000001095.2";
chr1    HAVANA  gene    69091   70008   .   +   .   gene_id "ENSG00000186092.4"; transcript_id "ENSG00000186092.4"; gene_type "protein_coding"; gene_status "KNOWN"; gene_name "OR4F5"; transcript_type "protein_coding"; transcript_status "KNOWN"; transcript_name "OR4F5"; level 2; havana_gene "OTTHUMG00000001094.1";
chr1    HAVANA  gene    89295   133566  .   -   .   gene_id "ENSG00000238009.2"; transcript_id "ENSG00000238009.2"; gene_type "lincRNA"; gene_status "NOVEL"; gene_name "RP11-34P13.7"; transcript_type "lincRNA"; transcript_status "NOVEL"; transcript_name "RP11-34P13.7"; level 2; havana_gene "OTTHUMG00000001096.2";
chr1    HAVANA  gene    89551   91105   .   -   .   gene_id "ENSG00000239945.1"; transcript_id "ENSG00000239945.1"; gene_type "lincRNA"; gene_status "NOVEL"; gene_name "RP11-34P13.8"; transcript_type "lincRNA"; transcript_status "NOVEL"; transcript_name "RP11-34P13.8"; level 2; havana_gene "OTTHUMG00000001097.2";
chr1    HAVANA  gene    131025  134836  .   +   .   gene_id "ENSG00000233750.3"; transcript_id "ENSG00000233750.3"; gene_type "pseudogene"; gene_status "KNOWN"; gene_name "CICP27"; transcript_type "pseudogene"; transcript_status "KNOWN"; transcript_name "CICP27"; level 1; tag "pseudo_consens"; havana_gene "OTTHUMG00000001257.3";
```

Perhaps filter a bit more and print the content of the 9^th^ column in the file.

``` bash
$ awk -F "\t" '$3 == "gene" { print $9 }' transcriptome.gtf | head | less -S

gene_id "ENSG00000223972.4"; transcript_id "ENSG00000223972.4"; gene_type "pseudogene"; gene_status "KNOWN"; gene_name "DDX11L1"; transcript_type "pseudogene"; transcript_status "KNOWN"; transcript_name "DDX11L1"; level 2; havana_gene "OTTHUMG00000000961.2";
gene_id "ENSG00000227232.4"; transcript_id "ENSG00000227232.4"; gene_type "pseudogene"; gene_status "KNOWN"; gene_name "WASH7P"; transcript_type "pseudogene"; transcript_status "KNOWN"; transcript_name "WASH7P"; level 2; havana_gene "OTTHUMG00000000958.1";
gene_id "ENSG00000243485.2"; transcript_id "ENSG00000243485.2"; gene_type "lincRNA"; gene_status "NOVEL"; gene_name "MIR1302-11"; transcript_type "lincRNA"; transcript_status "NOVEL"; transcript_name "MIR1302-11"; level 2; havana_gene "OTTHUMG00000000959.2";
gene_id "ENSG00000237613.2"; transcript_id "ENSG00000237613.2"; gene_type "lincRNA"; gene_status "KNOWN"; gene_name "FAM138A"; transcript_type "lincRNA"; transcript_status "KNOWN"; transcript_name "FAM138A"; level 2; havana_gene "OTTHUMG00000000960.1";
gene_id "ENSG00000268020.2"; transcript_id "ENSG00000268020.2"; gene_type "pseudogene"; gene_status "KNOWN"; gene_name "OR4G4P"; transcript_type "pseudogene"; transcript_status "KNOWN"; transcript_name "OR4G4P"; level 2; havana_gene "OTTHUMG00000185779.1";
gene_id "ENSG00000240361.1"; transcript_id "ENSG00000240361.1"; gene_type "pseudogene"; gene_status "KNOWN"; gene_name "OR4G11P"; transcript_type "pseudogene"; transcript_status "KNOWN"; transcript_name "OR4G11P"; level 2; havana_gene "OTTHUMG00000001095.2";
gene_id "ENSG00000186092.4"; transcript_id "ENSG00000186092.4"; gene_type "protein_coding"; gene_status "KNOWN"; gene_name "OR4F5"; transcript_type "protein_coding"; transcript_status "KNOWN"; transcript_name "OR4F5"; level 2; havana_gene "OTTHUMG00000001094.1";
gene_id "ENSG00000238009.2"; transcript_id "ENSG00000238009.2"; gene_type "lincRNA"; gene_status "NOVEL"; gene_name "RP11-34P13.7"; transcript_type "lincRNA"; transcript_status "NOVEL"; transcript_name "RP11-34P13.7"; level 2; havana_gene "OTTHUMG00000001096.2";
gene_id "ENSG00000239945.1"; transcript_id "ENSG00000239945.1"; gene_type "lincRNA"; gene_status "NOVEL"; gene_name "RP11-34P13.8"; transcript_type "lincRNA"; transcript_status "NOVEL"; transcript_name "RP11-34P13.8"; level 2; havana_gene "OTTHUMG00000001097.2";
gene_id "ENSG00000233750.3"; transcript_id "ENSG00000233750.3"; gene_type "pseudogene"; gene_status "KNOWN"; gene_name "CICP27"; transcript_type "pseudogene"; transcript_status "KNOWN"; transcript_name "CICP27"; level 1; tag "pseudo_consens"; havana_gene "OTTHUMG00000001257.3";
```

What about if we want just a specific piece from this information?
We can `|` the output from the first awk script in to a second one. Note that we will use different field separator `"; "`.

``` bash
$ awk -F "\t" '$3 == "gene" { print $9 }' transcriptome.gtf | awk -F "; " '{ print $3 }' | head

gene_type "pseudogene"
gene_type "pseudogene"
gene_type "lincRNA"
gene_type "lincRNA"
gene_type "pseudogene"
gene_type "pseudogene"
gene_type "protein_coding"
gene_type "lincRNA"
gene_type "lincRNA"
gene_type "pseudogene"
```
## Chaining AWK calls

We will start with the AWK call that we were using before, and we will append a pipe | so it can be used as input for the next AWK call, this time using a space and a semicolon as the delimiter to define what a column is:

```bash
awk -F "\t" '$3 == "gene" { print $9 }' transcriptome.gtf | awk -F "; " '{ print $3 }' | head | less -S

gene_type "pseudogene"
gene_type "pseudogene"
gene_type "lincRNA"
gene_type "lincRNA"
gene_type "pseudogene"
gene_type "pseudogene"
gene_type "protein_coding"
gene_type "lincRNA"
gene_type "lincRNA"
gene_type "pseudogene"
```
Now that we see what the third column looks like, we can filter for protein-coding genes

```bash
awk -F "\t" '$3 == "gene" { print $9 }' transcriptome.gtf | \
awk -F "; " '$3 == "gene_type \"protein_coding\""' | \
head | less -S

gene_id "ENSG00000186092.4"; transcript_id "ENSG00000186092.4"; gene_type "protein_coding"; gene_status "KNOWN"; gene_name "OR4F5"; transcript_type "protein_coding"; transcript_status "KNOWN"; transcript_name "OR4F5"; level 2; havana_gene "OTTHUMG00000001094.1";
gene_id "ENSG00000237683.5"; transcript_id "ENSG00000237683.5"; gene_type "protein_coding"; gene_status "KNOWN"; gene_name "AL627309.1"; transcript_type "protein_coding"; transcript_status "KNOWN"; transcript_name "AL627309.1"; level 3;
gene_id "ENSG00000235249.1"; transcript_id "ENSG00000235249.1"; gene_type "protein_coding"; gene_status "KNOWN"; gene_name "OR4F29"; transcript_type "protein_coding"; transcript_status "KNOWN"; transcript_name "OR4F29"; level 2; havana_gene "OTTHUMG00000002860.1";
gene_id "ENSG00000185097.2"; transcript_id "ENSG00000185097.2"; gene_type "protein_coding"; gene_status "KNOWN"; gene_name "OR4F16"; transcript_type "protein_coding"; transcript_status "KNOWN"; transcript_name "OR4F16"; level 2; havana_gene "OTTHUMG00000002581.1";
gene_id "ENSG00000269831.1"; transcript_id "ENSG00000269831.1"; gene_type "protein_coding"; gene_status "NOVEL"; gene_name "AL669831.1"; transcript_type "protein_coding"; transcript_status "NOVEL"; transcript_name "AL669831.1"; level 3;
gene_id "ENSG00000269308.1"; transcript_id "ENSG00000269308.1"; gene_type "protein_coding"; gene_status "NOVEL"; gene_name "AL645608.2"; transcript_type "protein_coding"; transcript_status "NOVEL"; transcript_name "AL645608.2"; level 3;
gene_id "ENSG00000187634.6"; transcript_id "ENSG00000187634.6"; gene_type "protein_coding"; gene_status "KNOWN"; gene_name "SAMD11"; transcript_type "protein_coding"; transcript_status "KNOWN"; transcript_name "SAMD11"; level 2; havana_gene "OTTHUMG00000040719.8";
gene_id "ENSG00000268179.1"; transcript_id "ENSG00000268179.1"; gene_type "protein_coding"; gene_status "NOVEL"; gene_name "AL645608.1"; transcript_type "protein_coding"; transcript_status "NOVEL"; transcript_name "AL645608.1"; level 3;
gene_id "ENSG00000188976.6"; transcript_id "ENSG00000188976.6"; gene_type "protein_coding"; gene_status "KNOWN"; gene_name "NOC2L"; transcript_type "protein_coding"; transcript_status "KNOWN"; transcript_name "NOC2L"; level 2; havana_gene "OTTHUMG00000040720.1";
gene_id "ENSG00000187961.9"; transcript_id "ENSG00000187961.9"; gene_type "protein_coding"; gene_status "KNOWN"; gene_name "KLHL17"; transcript_type "protein_coding"; transcript_status "KNOWN"; transcript_name "KLHL17"; level 2; havana_gene "OTTHUMG00000040721.6";
```
I added a space and a backslash \ (not to be confused with the regular slash /) after the first and second pipes to split the code into two lines; this makes it easier to read and it highlights that we are taking two separate steps.

The double quotes around protein_coding are escaped (also with a backslash \") because they are already contained inside double quotes. To avoid the backslashing drama we can use the partial matching operator ~ instead of the total equality operator ==.

```bash
awk -F "\t" '$3 == "gene" { print $9 }' transcriptome.gtf | \
awk -F "; " '$3 ~ "protein_coding"' | \
head | less -S
```
The output is the same as before: those lines that contain a protein_coding somewhere in their third column make the partial matching rule true, and they get printed (which is the default behavior when there are no curly braces).

Now we have all the protein-coding genes, but how do we get to the genes that only have one exon? Well, we have to revisit our initial AWK call: we selected lines that corresponded to genes, but we actually want lines that correspond to exons. That’s an easy fix, we just change the word “gene” for the word “exon”. Everything else stays the same.

```bash
awk -F "\t" '$3 == "exon" { print $9 }' transcriptome.gtf | \
awk -F "; " '$3 ~ "protein_coding"' | \
head | less -S

gene_id "ENSG00000186092.4"; transcript_id "ENST00000335137.3"; gene_type "protein_coding"; gene_status "KNOWN"; gene_name "OR4F5"; transcript_type "protein_coding"; transcript_status "KNOWN"; transcript_name "OR4F5-001"; exon_number 1;  exon_id "ENSE00002319515.1";  level 2; tag "basic"; tag "CCDS"; ccdsid "CCDS30547.1"; havana_gene "OTTHUMG00000001094.1"; havana_transcript "OTTHUMT00000003223.1";
gene_id "ENSG00000237683.5"; transcript_id "ENST00000423372.3"; gene_type "protein_coding"; gene_status "KNOWN"; gene_name "AL627309.1"; transcript_type "protein_coding"; transcript_status "KNOWN"; transcript_name "AL627309.1-201"; exon_number 1;  exon_id "ENSE00002221580.1";  level 3; tag "basic";
gene_id "ENSG00000237683.5"; transcript_id "ENST00000423372.3"; gene_type "protein_coding"; gene_status "KNOWN"; gene_name "AL627309.1"; transcript_type "protein_coding"; transcript_status "KNOWN"; transcript_name "AL627309.1-201"; exon_number 2;  exon_id "ENSE00002314092.1";  level 3; tag "basic";
gene_id "ENSG00000235249.1"; transcript_id "ENST00000426406.1"; gene_type "protein_coding"; gene_status "KNOWN"; gene_name "OR4F29"; transcript_type "protein_coding"; transcript_status "KNOWN"; transcript_name "OR4F29-001"; exon_number 1;  exon_id "ENSE00002316283.1";  level 2; tag "basic"; tag "CCDS"; ccdsid "CCDS41220.1"; havana_gene "OTTHUMG00000002860.1"; havana_transcript "OTTHUMT00000007999.1";
gene_id "ENSG00000185097.2"; transcript_id "ENST00000332831.2"; gene_type "protein_coding"; gene_status "KNOWN"; gene_name "OR4F16"; transcript_type "protein_coding"; transcript_status "KNOWN"; transcript_name "OR4F16-001"; exon_number 1;  exon_id "ENSE00002324228.1";  level 2; tag "basic"; tag "CCDS"; ccdsid "CCDS41221.1"; havana_gene "OTTHUMG00000002581.1"; havana_transcript "OTTHUMT00000007334.1";
gene_id "ENSG00000269831.1"; transcript_id "ENST00000599533.1"; gene_type "protein_coding"; gene_status "NOVEL"; gene_name "AL669831.1"; transcript_type "protein_coding"; transcript_status "NOVEL"; transcript_name "AL669831.1-201"; exon_number 1;  exon_id "ENSE00003063549.1";  level 3; tag "basic";
gene_id "ENSG00000269831.1"; transcript_id "ENST00000599533.1"; gene_type "protein_coding"; gene_status "NOVEL"; gene_name "AL669831.1"; transcript_type "protein_coding"; transcript_status "NOVEL"; transcript_name "AL669831.1-201"; exon_number 2;  exon_id "ENSE00003084653.1";  level 3; tag "basic";
gene_id "ENSG00000269831.1"; transcript_id "ENST00000599533.1"; gene_type "protein_coding"; gene_status "NOVEL"; gene_name "AL669831.1"; transcript_type "protein_coding"; transcript_status "NOVEL"; transcript_name "AL669831.1-201"; exon_number 3;  exon_id "ENSE00003138540.1";  level 3; tag "basic";
gene_id "ENSG00000269308.1"; transcript_id "ENST00000594233.1"; gene_type "protein_coding"; gene_status "NOVEL"; gene_name "AL645608.2"; transcript_type "protein_coding"; transcript_status "NOVEL"; transcript_name "AL645608.2-201"; exon_number 1;  exon_id "ENSE00003079649.1";  level 3; tag "basic";
gene_id "ENSG00000269308.1"; transcript_id "ENST00000594233.1"; gene_type "protein_coding"; gene_status "NOVEL"; gene_name "AL645608.2"; transcript_type "protein_coding"; transcript_status "NOVEL"; transcript_name "AL645608.2-201"; exon_number 2;  exon_id "ENSE00003048391.1";  level 3; tag "basic";
```
## Cleaning up the output

Before we try to count how many exons belong to the same protein-coding gene, let’s simplify the output so we only get the gene names (which are in column 5).

```bash
awk -F "\t" '$3 == "exon" { print $9 }' transcriptome.gtf | \
awk -F "; " '$3 ~ "protein_coding" {print $5}' | \
head

gene_name "OR4F5"
gene_name "AL627309.1"
gene_name "AL627309.1"
gene_name "OR4F29"
gene_name "OR4F16"
gene_name "AL669831.1"
gene_name "AL669831.1"
gene_name "AL669831.1"
gene_name "AL645608.2"
gene_name "AL645608.2"
```
This is sort of what we want. We could chain another AWK call using -F " ", and pick the second column (which would get rid of the gene_name). Feel free to try that approach if you are curious.

We can also take a shortcut by using the tr -d command, which deletes whatever characters appear in double quotes. For example, to remove every vowel from a sentence:
```bash
echo "This unix thing is cool" | tr -d "aeiou" # Ths nx thng s cl
```
Let’s try deleting all the semicolons and quotes before the second AWK call:

```bash
awk -F "\t" '$3 == "exon" { print $9 }' transcriptome.gtf | \
tr -d ";\"" | \
awk -F " " '$6 == "protein_coding" {print $10}' | \
head

OR4F5
AL627309.1
AL627309.1
OR4F29
OR4F16
AL669831.1
AL669831.1
AL669831.1
AL645608.2
AL645608.2
```
Run ```bash awk -F "\t" '$3 == "exon" { print $9 }' transcriptome.gtf | tr -d ";\"" | head ``` to understand what the input to the second AWK call looks like. It’s just words separated by spaces; the sixth word corresponds to the gene type, and the tenth word to the gene name.

## Counting genes

There is one more concept we need to introduce before we start counting. AWK uses a special rule called END, which is only true once the input is over. See an example:

```bash
echo -e "a\na\nb\nb\nb\nc" | \
awk '
    { print $1 }

END { print "Done with letters!" }
'

a
a
b
b
b
c
```
Done with letters!
The -e option tells echo to convert each \n into a new line, which is a convenient way of printing multiple lines from a single character string.

In AWK, any amount of whitespace is allowed between the initial and the final quote '. I separated the first rule from the END rule to make them easier to read.

Now we are ready for counting.

```bash
echo -e "a\na\nb\nb\nb\nc" | \
awk '
    { counter[$1] += 1 }

END {
    for (letter in counter){
        print letter, counter[letter]
    }
}
'

a 2
b 3
c 1
```
Wow, what is all that madness?

Instead of printing each letter, we manipulate a variable that we called counter. This variable is special because it is followed by brackets [ ], which makes it an associative array, a fancy way of calling a variable that stores key-value pairs.

In this case we chose the values of the first column $1 to be the keys of the counter variable, which means there are 3 keys (“a”, “b” and “c”). The values are initialized to 0. For every line in the input, we add a 1 to the value in the array whose key is equal to $1. We use the addition operator +=, a shortcut for counter[$1] = counter[$1] + 1.

When all the lines are read, the END rule becomes true, and the code between the curly braces { } is executed. The structure for (key in associate_array) { some_code } is called a for loop, and it executes some_code as many times as there are keys in the array. letter is the name that we chose for the variable that cycles through all the keys in counter, and counter[letter] gives the value stored in counter for each letter (which we we calculated in the previous curly brace chunk).

Now we can apply this to the real example:

```bash
awk -F "\t" '$3 == "exon" { print $9 }' transcriptome.gtf | \
tr -d ";\"" | \
awk -F " " '
$6 == "protein_coding" {
    gene_counter[$10] += 1
}

END {
    for (gene_name in gene_counter){
        print gene_name, gene_counter[gene_name]
    }
}' > number_of_exons_by_gene.txt

head number_of_exons_by_gene.txt
CAPN6 24
ARL14EPL 3
DACH1 38
IFNA13 1
HSP90AB1 36
CAPN7 52
DACH2 84
IFNA14 1
LARS 188
CAPN8 78
```
If you are using the real transcriptome, it takes less than a minute to count up one million exons. Pretty impressive.

We saved the output to a file, so now we can use AWK to see how many genes are made up of a single exon.

```bash
awk '$2 == 1' number_of_exons_by_gene.txt | wc -l # 1362
```


I will suggest to follow the original tutorial if you need to refer to these steps later on for your own data:  
[AWK GTF! How to Analyze a Transcriptome Like a Pro](http://reasoniamhere.com/2013/09/16/awk-gtf-how-to-analyze-a-transcriptome-like-a-pro-part-1/)