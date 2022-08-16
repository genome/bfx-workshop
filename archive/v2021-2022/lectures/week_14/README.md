## Week 14: Exploring data from the command line, pt 2 + Office hours

- [Lecture Recording]()


### Bash history

- Use CTRL-R 
- `history | grep foo`


```
#Keep your bash history around for (nearly) ever
export HISTFILESIZE=10000000000
export HISTSIZE=1000000000
#don't save certain commands in your history
export HISTIGNORE="ll:clear:exit:ll -rt"
```

### Using cut/awk/sed to transform data
```
#cut uses tabs by default
cut -f 3 file
#but you can pass other delimiters
cut -f 3 -d - file
cut -f 7 -d / file

#awk splits on spaces by default, and allows for simple manipulation
awk '{print $1,$3}' file >output

#make a bedfile into igv-style coords:
awk '{print $1":"$2"-"$3}' file.bed >output

#convert a 1-based VCF into a bed file:
grep -v "^#" my.vcf | awk '{print $1,$2-1,$2}' >out.bed
```

### bash aliases to make life simpler - put these in your ~/.bashrc 

```
#make less use larger spaces by default
alias less='less -cSx 20'

#often better, "columnify" your data
cat file | head | column -t 

#a function to do this (that's a tab character after -s, not multiple spaces)
#also note that multiple empty tabs in a row may screw with this
function clf {
    column -t -s '      ' $1 | less
}

#use awk to quickly sum a column of data (ints or floats)
alias sumcol='awk '\''{ SUM += $1} END { print SUM}'\'
alias sumcolfloat='awk '\''{ SUM += $1} END { OFMT="%4.8f"; print SUM}'\'

#unique count - get a list of unique items, sorted by frequency
alias uc='sort | uniq -c | sort -nrk 1'

#sort a bed file
alias bedsort='sort -k1,1 -k2,2n'

#Sort human-readable sizes correctly. Save this as `hsort` somewhere in your path and make it executable:
perl -e 'sub h{%h=(K=>10,M=>20,G=>30,T=>40);($n,$u)=shift=~/([0-9.]+)(\D)/;return $n*2**$h{$u}}print sort{h($b)<=>h($a)}<>;'
