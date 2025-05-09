
Solution suggestion for median MEAN_TARGET_COVERAGE:
```bash
cut -f2 qc_summary.txt | grep -v 'MEAN_TARGET_COVERAGE' | sort -n | awk '{
    a[NR] = $1
}
END {
    if (NR % 2) {
        print a[(NR + 1) / 2]
    } else {
        print (a[NR / 2] + a[NR / 2 + 1]) / 2
    }
}'
```

Find TOTAL_READS values:
```bash
# find column number
grep "TOTAL_READS" PB-P-HD1-CFDNA-1811-KH20190925-PB20190927.PB.hsmetrics.txt | tr "\t" "\n" | less -SN
#--> 23
awk -F "\t" -v OFS="\t" 'NR==7 {print "file", $23}; FNR==8 {print FILENAME, $23}' \
    ${files[@]} | column -t
```
