## Week 3: Command line skills + Office hours

This week, we did a short talk about some unix command line utilities and how they can often get you to quick answers, followed by office hours to help with homework or answer general questions.

- [Lecture Recording](https://wustl.app.box.com/s/sxndl8ohr6q6eit5gmivfti3efh5oihv)

## Notes from the talk
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

# Note: you can grab the [genes1](genes1.txt), [genes2](genes2.txt), and [tcga.tsv](tcga.tsv) files with wget, but any files will do!

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

#find lines in one file that match genes from another
grep -wf genes1 tcga.aml.tsv | less

#what if we want to find lines in genes1 that match the genes column in our complicated file?  
# Nested expressions! <(command) treats the output of that command as a file
grep -wf <(cut -f 10 tcga.aml.tsv | sort | uniq) genes1
```
