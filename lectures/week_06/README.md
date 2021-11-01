## Week 6: Exploring data from the command line + Office hours

This week's lecture was a mini demonstration of some bash utilities (cut, sort, uniq).

- [Lecture Recording]()

### Exploring data with bash - part 1

```
#there are two paths coming out of every command - stdout and stderr

#by default both write to the screen
date #output to stdout
date --asdf #error to stderr

#you can redirect stdout to a file
date >file.txt

#you can also redirect stderr to a file
date --asdf 2>err.log

#you can get fancy and send both to the same file
date 2>&1

#use wc -l to count the number of lines in a file
wc -l genes1

#use sort to reorder a file
sort genes1

#set operations
#find shared genes
cat genes1 genes2 | sort | uniq -d 

#find genes1 specific
cat genes1 genes2 genes2 | sort | uniq -u

#find genes2 specific
cat genes1 genes1 genes2 | sort | uniq -u

#get a certain column from a file
cut -f 10 tcga.aml.tsv | head

#find the most frequent items in a certain column
cut -f 10 tcga.aml.tsv | sort | uniq -c | sort -nk 1 | head
```
