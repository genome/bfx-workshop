## Manipulating files at the command line

Plain text files are fundamental to bioinformatics and data processing in general. Your favorite programming language will have lots of ways to read/write/manipulate these files. Sometimes, though, you need to manipulate or reformat files up front, or maybe just need a quick answer. That's where using command-line utilities can help. 

Open up your terminal, then navigate to the place where you checked out the course git directory last week.  If you can't remember, or have deleted it already, something like this might work to check it out again:

```
mkdir ~/workshop
cd ~/workshop
git clone https://github.com/genome/bfx-workshop.git
cd bfx-workshop/lectures/week_02
```
(what does the `~` in those commands mean?)

You can always use the `pwd -P` command to see where you are in the filesystem.

```
pwd -P 
```

### STDOUT and STDERR 
Let's start with some basics. Every time you run a command, there are two paths via which information can come out: `stdout` and `stderr`. You can think of them as "the path your data takes" and "the path that errors or warnings take".  Let's explore them a little using the "date" utility:

By default, date outputs to stdout, and by default, that goes to your screen.

```
date 
```
If you give it an invalid parameter, the resulting error will output to stderr

```
date --asdf #error to stderr
```

You can redirect stdout to a file

```
date >file.txt
```

Then you can check it's contents by using "cat", which takes a file and spits it out to stdout

```
cat file.txt 
```

You can also redirect stderr to a file using `2>`

```
date --asdf 2>err.log
cat err.log
```

Note that if you redirect the output of stderr, stdout still goes to your screen (and vice versa):

```
date 2>err.log
```

Using double arrows results in appending to a file (wait a second or two before running each cmd)

```
date >file.txt
date >>file.txt
date >>file.txt
cat file.txt
```



### sort and uniq

In this `week_02` directory, there are a few files with some gene names and mutation data. Use `ls` to look at them, or `ls -l` for more information. 

A common thing we'd like to do is count the number of lines in a file. Use ls to take a look Let's take a look at the counts of the genes1 file:

```
wc -l genes1.txt
```
or get line counts for all of the genes files:

```
wc -l genes*
```
(what is the * doing here?) 

Having files in an consistent order can be really important for many applications. The `sort` command handles that, alphabetically by default:

```
sort genes1.txt
```

We also have the `uniq` command, which removes identical lines that are _adjacent in the list_.  Compare the outputs of these two commands:

```
cat genes1.txt | uniq
cat genes1.txt | uniq | wc -l

sort genes1.txt | uniq
sort genes1.txt | uniq | wc -l
```

Here we've also introduced the pipe `|` operator.  It takes the `stdout` from one command and pushes it to the `stdin` of the next.  It's an incredibly useful way to chain commands together. 

By default `uniq` spits out only the unique entries, but you can also ask it to give you just the duplicate entries (appearing multiple times in the list (`uniq -d`) or just the unique entries (appearing exactly once in the list `uniq -u`). Try these out on your genes files.



### Set operations

Let's suppose that the genes1 and genes2 files contain the outputs of two different analyses, and I'd like to compare the lists.  These boil down to set operations, and we can do some of them right from the command line. 

If I want to find shared genes between the two lists, Let's first deduplicate each one:

```
sort genes1.txt | uniq >genes1.uniq.txt 
sort genes2.txt | uniq >genes2.uniq.txt
```

Great, now if I combine the two lists, I know that any gene appearing twice in the list must have been found in both experiments:

```
cat genes1.uniq.txt genes2.uniq.txt | sort | uniq -d
```

We can also use this same concept to find genes specific to one list or the other.

```
cat genes1.uniq.txt genes2.uniq.txt genes2.uniq.txt | sort | uniq -u
```
What happened there?  By including genes2 in my list twice, I have guaranteed that every single gene in genes2 will be duplicated and thus removed. Anything that is only in genes1 will be output.

Try composing a command that gives you the opposite - genes unique to the genes2 list. 


### Tab-delimited tables

Lots of data in bioinformatics uses tab-delimited data. Let's use the `head` command to take a look at this file, which contains mutations found in a large cohort of AML patients:

```
head tcga.tsv
```
Woah - that's kind of messy.  Instead of `head`, we can use the `less` command instead, which allows us to view the file:

```
less tcga.tsv
```
Once inside `less`, use the arrow keys to navigate, or pgup/pgdn.  Type `-S` to toggle wrapping of lines and immediately the data is easier to look at, but we can do even better.  Type `-x`, then `20` to set the tab stops to 20 characters. Note how it makes many things line up better.  Type `q` to exit less.

Our `sort` commands from above are useful here too, especially if we pipe them to less so that it doesn't scroll forever.

```
sort tcga.tsv | less
```

Or better yet, sort by a particular column or columns. What do these do? 

```
sort -k 2 tcga.tsv | less
sort -k 2,2 -k 3,3n tcga.tsv | less
```
 (`man sort` will bring up the manual, scroll with the arrows, again use `q` to exit)
 
Often, we don't want to work with a whole table, but with a certain portion of a file.  Let's grab just the gene names with the `cut` command:

```
cut -f 8 tcga.tsv | head
```

Now, let's use the commands we learned above to find the most frequently mutated genes in this column.  Note: `uniq -c` returns a count of how many times an item appears in a list:

```
cut -f 8 tcga.tsv | sort | uniq -c | head -n 20
```
Well, that's close - it gave us some genes and their counts, but how do I find out which ones are the most frequent from the whole long list?  Well, the output of uniq -c is itself a list with two columns, so we can sort it again by that first column showing the frequency! (`-n` for numeric and `-r` for reversed)

```
cut -f 8 tcga.tsv | sort | uniq -c | sort -nrk 1 | head -n 20
```

It takes time to get fluent in these commands, but these are great examples of ways to really quickly get simple answers from your data with only a few keystrokes. 

The assignment ends here, and if your head is already spinning, that's fine!  If you're eager to learn more, read on for some more advanced tips and tricks. (or come back to this page in a few weeks when you're more comfortable at the command line and ready to level up your skills!

### More complex concepts 
These can serve as a reference and ideas to play with as you become more advanced in your command line skills:

Find lines in a file that match a certain gene name:

```
grep "TP53" tcga.tsv | less
```

Find lines in one file that match genes from another file:

```
grep -wf genes1.txt tcga.tsv | less
```

Nested expressions turn the outputs of an expression into a "fake" file that can  be read by the surrounding command.  For example, what if we want to find which members of our `genes1.txt` output have entries in the genes column in our mutation file?  

```
grep -wf <(cut -f 8 tcga.tsv | sort | uniq) genes1.txt
```

Using loops are really useful as well to avoid typing the same thing over and over:

```
for i in 1 2 3 4;do echo $i;done
```

Here we do more complicated things in the loop

```
grep IDH genes1.txt | while read i;do echo "my favorite gene is $i";done
grep IDH genes1.txt | while read i;do grep $i tcga.tsv;done
```

We can mix types of loops (`for` and `while`) and also nest loops:

```
for num in 1 2 3 4;do 
  cat genes1.txt | while read gene;do 
    echo "$gene $num";
  done;
done
```

The `xargs` command can also be used instead of loops to split things out into separate commands, individually or in chunks

```
cat genes1.txt | xargs -n 1 echo
cat genes1.txt | xargs -n 5 echo
```

Instead of cutting one column, we can grab multiple columns from a file

```
cut -f 1,2 tcga.tsv | head
```

Note that doing it this way, doesn't change the order of the columns:

```
cut -f 2,1 tcga.tsv | head
```

But an awk command can help with this!

```
awk '{print $2,$1}' tcga.tsv | head
```

Some extra parameters can be used to keep things tab-delimited, instead of whitespace-delimited, which is the default;

```
awk -F"\t" -v OFS="\t" '{print $2,$1}' tcga.tsv | head
```

### Other tips

Use the up arrow to scroll through bash history, and the tab key to autocomplete. Use CTRL-R to search backwards through your history (keep hitting CTRL-R to scroll through matches) or type something like `history | less` or `history | grep python`

`cd` with no arguments always takes you to your home directory
`cd -` takes you back to the previous directory

`pushd` and `popd` can be used to navigate through directories saving your place, then jumping back to where you were

```
pushd /tmp
pushd ~/workshop
popd
popd
```

Advanced tip: You can get fancy and send both stdout and stderr to the same file

```
date >file.txt 2>&1
date --asdf >>file.txt 2>&1
```

Chain things together! Here's an example of change a gzipped VCF file from having "chr" prefixes to having no prefixes using only one line of code. First block gets the header, second block removes the chr prefixes

```
cat <(gunzip -c asdf.vcf.gz | grep "^#") <(gunzip -c asdf.vcf.gz | grep -v "^#" | sed 's/^chr//') | gzip >asdf.fixed.vcf.gz
```