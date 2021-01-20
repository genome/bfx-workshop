## Week 7: Exploring data from the command line, pt 2 + Office hours

- Description

- [Lecture Recording](https://wustl.box.com/s/ivtmdn677mooltqz1irgelc2m4aiwjno)

### Exploring data with bash - part 2

Two ways to search your history:

- Use CTRL-R 
- `history | grep foo`

```
#Keep your bash history around for (nearly) ever
export HISTFILESIZE=10000000000
export HISTSIZE=1000000000
#don't save certain commands in your history
export HISTIGNORE="ll:clear:exit:ll -rt"


#have less use larger spaces by default
alias less='less -cSx 20'

#even better, "columnify" your data
function clf {
    column -t -s '  ' $1 | less
}
#that's a tab character after -s
#also note that multiple empty tabs in a row may screw with this

#use awk to quickly sum a column of data (ints or floats)
alias sumcol='awk '\''{ SUM += $1} END { print SUM}'\'
alias sumcolfloat='awk '\''{ SUM += $1} END { OFMT="%4.8f"; print SUM}'\'

#unique count - get a list of unique items, sorted by frequency
alias uc='sort | uniq -c | sort -nrk 1'

#sort a bed file
alias bedsort='sort -k1,1 -k2,2n'
```

Sort human-readable sizes correctly. Save this as `hsort` somewhere in your path and make it executable:

```
perl -e 'sub h{%h=(K=>10,M=>20,G=>30,T=>40);($n,$u)=shift=~/([0-9.]+)(\D)/;return $n*2**$h{$u}}print sort{h($b)<=>h($a)}<>;'
```