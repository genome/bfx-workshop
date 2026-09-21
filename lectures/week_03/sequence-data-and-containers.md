# Working with Sequence Data: Downloading, Unpacking, and QC with Containers

Last week we learned how to log in to Compute2 and run jobs with Slurm. This week we'll put
that to work on some real sequencing data. By the end you should be able to:

- Download data from a collaborator or public repository straight onto the cluster
- Unpack `.tar` archives and `.gz` compressed files, and understand the difference between them
- Poke around in FASTQ seqeuence files and answer basic questions about them from the command line
- Explain what a container (Docker image) is and why we use them
- Run tools that aren't installed on the cluster (FastQC) by launching a Slurm job that runs a container

> [!NOTE]
> Just like last week, we strongly suggest typing the commands out yourself rather than copying and pasting. It really does help it stick.

Throughout this lesson, questions are followed by expandable **Hint** and **Solution** boxes. Click through to get help or check your work, but take a stab at doing it yourself first!

---

## 1. You've got data!

A collaborator emails you: *"Hey, here's that tumor exome data we talked about. Let me know what you find"* Along with it is a link:

```
https://storage.googleapis.com/bfx_workshop_tmp/Exome_Tumor.tar
```

This is a very common way to receive data, whether it's from a collaborator, a sequencing
core, or a public repository like the [SRA](https://www.ncbi.nlm.nih.gov/sra) or
[ENA](https://www.ebi.ac.uk/ena). You *could* click that link, download it to your laptop,
and then figure out how to move it to the cluster, but it's much faster (and for big datasets,
the only realistic option) to download it directly onto the cluster.

Log in to Compute2, and make sure your "workshop" directory still exists.  Replace \<washukey\> with your own username.

```bash
ssh <washukey>@c2-login-002.ris.wustl.edu
ls -l
```
(If not, refer back to [last's week's exercise](https://github.com/genome/bfx-workshop/blob/master/lectures/week_02/compute2-slurm-intro.md#3-where-is-my-data) to get that folder/shortcut set up)

Then, move into your workshop directory, make a folder for this week, and move into it.

```bash
cd ~/workshop
mkdir week03
cd week03
```


Next, you'll need to download the file. A common tool for doing this is `curl`.  Just like the other command line tool we've been using, `curl` spits it's output to STDOUT, so we'll redirect it to a file:

```bash
curl https://storage.googleapis.com/bfx_workshop_tmp/Exome_Tumor.tar >Exome_Tumor.tar
```

(Alternately, you can pass the `-O` flag, which tells `curl` to save the file using its original name. Either approach is fine)

```bash
curl -O https://storage.googleapis.com/bfx_workshop_tmp/Exome_Tumor.tar 
```

Check that it arrived, using `-h` to show "human readable" file sizes:

```bash
ls -lh
```

> [!NOTE]
> Downloading and moving data around is a perfectly fine thing to do on the login node.
> This file is tiny (2 MB). If you're pulling down hundreds of gigabytes, though, it's worth
> doing inside a Slurm job so you're not tying up the login node for hours.

---

## 2. Unpacking the data: `tar` and `gzip`

Our file ends in `.tar`, which tells us it's a **tar archive**. 

> ### An aside: What is a tar file?
> 
> If you've used `.zip` files, you know they do two jobs at once: they **bundle** many files
> and folders into a single file, and they **compress** that file to make it smaller.
> 
> In unix, we commonly split those two jobs between two different tools:
> 
> - **`tar`** (short for "tape archive") only does the *bundling*. It glues files and
> directories together into one big file that's easy to move around, keeping the folder structure intact. By itself, a `.tar` file is no smaller than the files inside it.
> 
> - **`gzip`** only does the *compressing*. It shrinks a single file
> and adds `.gz` to the name.
> 
> You'll often see files named `something.tar.gz` (or `.tgz`), which means a bunch of files were gathered into a bundle, and then compressed. 

Like the other unix tools we've ben using, tar can take flags.  The ones we're interested in today are `-xvf`

| Flag | Meaning |
|---|---|
| `-x` | e**x**tract files from an archive |
| `-v` | **v**erbose: print each file name as it goes |
| `-f` | the **f**ile name to extract (so `f` should always come last) |


> Note - here are a few other handy recipes:
>
> ```bash
> tar -tvf archive.tar             # peek inside without extracting
> tar -xzvf archive.tar.gz         # extract a gzipped tarball
> tar -czvf mydir.tar.gz mydir/    # bundle AND compress a directory
> ```

To extract our file:

```bash
tar -xvf Exome_Tumor.tar
```
Use `ls` to see what happened and notice that it created a directory. Move into it and look around:

```bash
cd Exome_Tumor
ls -lh
```

You'll see `2891351066_1.fastq.gz` and `2891351066_2.fastq.gz`

These are gzip-compressed FASTQ files. Many tools that work with fastq files can unzip them on the fly (because decompressing 100Gb of data takes a lot of disk!).  These are small, though, and we want to look more closely at what's inside. Decompress them with `gunzip`:

```bash
gunzip *.fastq.gz
ls -lh
```

```
-rw-r--r-- 1 washukey domain users 6.1M Sep 27  2020 2891351066_1.fastq
-rw-r--r-- 1 washukey domain users 6.1M Sep 27  2020 2891351066_2.fastq
```

Notice that each file got about 6 times bigger. Sequence data compresses *very* well, which is
why you'll almost always receive and store it compressed.

---

## 3. Exploring the FASTQ files

As a quick refresher from the lecture, a FASTQ file is made of **records**, and every record is
exactly **four lines**:

```
@read_name            <- 1. starts with @, then the read's name/ID
GATTTGGGGTTCAAAGCAG   <- 2. the sequence
+                     <- 3. a separator (sometimes repeats the name)
!''*((((***+))%%%++   <- 4. quality scores, one character per base
```

Let's explore these data and answer some questions about it.

### Question 1

Look at the first three records (not first three lines!) of each fastq file. Take a close look at the read names - do they match up across files?  Are there any differences?

<details>
<summary>Hint</summary>

Use the `head` command. By default it shows 10 lines, but you can change that with `-n`. How
many lines make up three records?

</details>

<details>
<summary>Solution</summary>

Three records × four lines per record = 12 lines:

```bash
head -n 12 2891351066_1.fastq
head -n 12 2891351066_2.fastq
```

The read names are **identical** between the two files except for the `/1` and `/2` at the end. This is **paired-end** data: each DNA fragment was sequenced once from each end. Read 1 of
a pair is in the `_1` file, its mate is in the `_2` file, and the two files are in the
**same order**. Aligners rely on that ordering, so never sort or filter one file without doing
the same to the other!

----
</details>


### Question 2

How many paired-end sequences do these files contain?

<details>
<summary>Hint</summary>

Last week we used `wc -l` to count the lines in a file. 

</details>

<details>
<summary>Solution</summary>

```bash
wc -l *.fastq
```

```
 100000 2891351066_1.fastq
 100000 2891351066_2.fastq
 200000 total
```

100,000 lines ÷ 4 lines per record = **25,000 reads in each file**. Since each pair has one read in each file, that's **25,000 read pairs** (or 50,000 individual reads).

</details>

It's worth checking that both files have the same number of records. If they didn't, something
went wrong (like a truncated download) and the pairs would no longer line up.

### Question 3

What is the read length? Is the read length consistent for every record?  This is going to require commands beyond the ones we've covered.  Use your favorite AI assistant to help you come up with a command that will help you figure it out.

<details>
<summary>Hint 1</summary>
Be sure to describe to the bot what you're doing in detail. Tell it what file type you're using, what you're trying to achieve, and have it talk you through how it approached the solution. Do you know a little python? or maybe some awk? If you do, have it use your favorite language.  

Keep in mind (and tell it!) that we're not working in a docker container and your ability to install software is limited. If it tries to write a python script that uses the `pysam` package, you might not be able to do in a straightforward way - you'd need a docker container with that package installed!
</details>

<details>
<summary>Hint 2</summary>
Try asking an LLM how to use `awk` to print the length of every 4th line, starting from line 2.

How can you do sanity checks by hand to determine whether the script is written correctly?  

In many tools, the invisible newline character at the end of each line wil count as a character.  How could you verify whether or not your solution takes this into account?  
</details>

<details>
<summary>Solution</summary>

One solution might look like this: 

Quick check of the first read:

```bash
head -n 2 2891351066_1.fastq | tail -n 1 | wc -c
```

```
101
```

That's 100 bases plus one newline character, so the reads look like they're **100 bp**.

To check every record:

```bash
awk 'NR % 4 == 2 {print length($0)}' 2891351066_1.fastq | sort | uniq -c
```

```
  25000 100
```

Breaking that down:

- `NR % 4 == 2` means "the line number divided by 4 has a remainder of 2", which is true for
  lines 2, 6, 10, 14... (the sequence lines)
- `{print length($0)}` prints the length of that whole line (`$0`)
- `sort | uniq -c` counts how many times each length appears (like the gene counting we did last week)

All 25,000 reads are 100 bp, so yes, the read length is **consistent**. Run the same command on
the `_2` file to confirm that it's true there too. (That won't always be the case! Reads that have
had adapters or low-quality bases trimmed off will have varying lengths.)

</details>

### Question 4

How many total nucleotides of sequence are contained in these two files?

<details>
<summary>Hint</summary>

Modify the `awk` command from the previous question to *add up* the lengths instead of printing them. Again, a chat with an AI assistant might be helpful!

Also, you already know how many reads there are and how long each one is. So you can verify your solution by doing the math! 

</details>

<details>
<summary>Solution</summary>

25,000 reads × 100 bp × 2 files = **5,000,000 nucleotides** (5 Mb).

You can verify it directly:

```bash
cat *.fastq | awk 'NR % 4 == 2 {total += length($0)} END {print total}'
```

```
5000000
```
Ask your LLM partner to explain the awk block to you so that you understand how it works.

For perspective, the human genome is about 3.1 billion bases, so this is a tiny subsample. A real exome would typically have tens of millions of read pairs.

</details>

### Question 5

Recompress these two fastq files to save space.

<details>
<summary>Hint</summary>

`gzip` is the reverse of `gunzip`. 
</details>

<details>
<summary>Solution</summary>

You could type a command for each file, but instead, we can use the wildcard "*" character to match all of the fastq files:

```bash
gzip *.fastq
ls -lh
```

```
-rw-r--r-- 1 washukey domain users 1020K Sep 27  2020 2891351066_1.fastq.gz
-rw-r--r-- 1 washukey domain users 1.1M Sep 27  2020 2891351066_2.fastq.gz
```

Back to ~1 MB each.

</details>

---

## 4. Quality control with FastQC

Now that we have a feel for the data, let's run a full quality check. The standard first-pass
tool for this is [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/), which
produces a report of a dozen different quality metrics.

Let's try running it:

```bash
fastqc --version
```

```
-bash: fastqc: command not found
```

Hmm. FastQC isn't installed on the cluster. This is going to happen *a lot*: there are
thousands of bioinformatics tools, each with many versions, and no cluster has all of them.
So how do we get access to the software we need?

---

## 5. Software environments

When we log in to the cluster, all of the tools that we might want to use aren't going to be there. Take samtools, for instance:

```
$ samtools
-bash: samtools: command not found
$ R
-bash: R: command not found
```

There are three major ways to get access to software (and different versions of software).

### A. Compile software yourself

This is often tricky, and sometimes lands you in 'dependency hell'. It's sometimes necessary, but shouldn't be your first stop.

### B. Load environments/software with `ml` modules

The cluster has certain tools installed, but they are not all on your `PATH` by
default — that would cause version conflicts. Instead, software is exposed through **modules**, which you manage with the `ml` command (short for Lmod's `module`).

See what is available:

```bash
ml avail
```

This prints a long list grouped by location. You can see even more by doing something like

```
module load ris
module --show_hidden avail
```

These trend heavily towards general-purpose packages, though, and there aren't a ton of
bioinformatics tools available here. Also, in order to install new software as a module, RIS has to be involved, and that takes time.

### C. Load software/environments with Docker

Docker is a containerization technology. Much like modules, Docker allows you to "load up" the software you need, in an environment that's pinned to a certain version, etc.

How does it work, at a high level?

- **A container is like a computer within the cluster, pre-installed with any software you want.**
  Someone (maybe you, maybe a community project) has already set up a mini Linux system with
  the tool, the right version of it, and every library it depends on.
- **It runs independently of all other containers and of the base system.** A container with
  Python 2 and a container with Python 3 can run side by side on the same node without ever
  knowing about each other. Nothing you install or break inside one container affects any other
  container, or the cluster itself. When you exit, the container goes away.
- **An *image* is the blueprint; a *container* is a running copy of it.** Images are stored in
  online registries like [Docker Hub](https://hub.docker.com/) and [Quay.io](https://quay.io/),
  and anyone can download and run them.

Image names follow a pattern of `registry/owner/name:tag`. For example:

```
quay.io/biocontainers/fastqc:0.11.9--0
  |          |          |      |
registry   owner      tool   tag (version)
```

The [BioContainers](https://biocontainers.pro/) project automatically builds images for
thousands of bioinformatics tools, so for most common tools, someone has already done the work
for you.

Because the tag pins an exact version, running the same image gives you the same software
today, next year, or on a completely different computer. That's great for reproducibility, and for your sanity (what happens if you need to re-run an analysis 2 years down the line?)

---

## 6. Running FastQC in a container with Slurm

On Compute2, Slurm can start your job *inside* a container for you. You just add one more
option, `--container-image`, to the same `srun` commands you learned last week.

First, make sure you're in the directory with the data:

```bash
cd ~/workshop/week03/Exome_Tumor
```

Then request an interactive job running inside the FastQC image:

```bash
srun -A compute2-workshop -p workshop -c 1 --mem=4G --pty --container-image=quay.io/biocontainers/fastqc:0.11.9--0 /bin/bash
```

Most of this should look familiar:

| Option | Meaning |
|---|---|
| `-A compute2-workshop` | charge the job to the workshop account |
| `-p workshop` | submit to the workshop's interactive partition |
| `-c 1 --mem=4G` | one CPU and 4 GB of RAM |
| `--pty ... /bin/bash` | give me an interactive shell |
| `--container-image=...` | **new!** run that shell inside this image |

> [!NOTE]
> The first time you use an image, it has to be downloaded, so it may take a minute or two to
> start. After that it's usually cached and starts much faster.

Okay, we're in!  Let's try to run `fastqc` again:

```
$ fastqc --version
FastQC v0.11.9
```

Great! It worked.  When I run `ls`, though, my data isn't there!  

```
$ ls
bin      dev      etc      home     lib      lib64    linuxrc  media    mnt      opt      proc     root     run      sbin     sys      tmp      usr      var
```

Where'd those fastq files go!?

### Problem #1 - slurm moved us to the root of the filesystem

The `pwd -P`  (print working directory) command lets us see where we are:

```
$ pwd -P
/
```

It just returns a slash, which means we're at the very top of the filesystem.  Okay, easy enough - let's just navigate back to where our data lives:

```
$ cd ~/workshop/
bash: cd: /home/c.a.miller/workshop/: No such file or directory
```

### Problem #2 - our filesystems aren't mounted

We talked about how docker images are little-self contained operating systems, and by default, only your home directory gets shared inside.  Look - the shortcut is still there, but is highlighted in red giving us a warning that the thing it's pointing to isn't there - the link is broken:

```
$ ls -l ~
```

Let's exit out of the docker container and back to the head node by simply typing

```
exit
```

To fix our data problem, we'll need to tell slurm and docker to pass that directory throughusing a new option to slurm called `--container-mounts`.  It takes a quoted list of colon-separated source and destination directories. To make things simple, let's just make them the same on both sies.

Adding it to our command gets us the full command we need:

```
srun -A compute2-workshop -p workshop -c 1 --mem=4G --pty \
--container-image=quay.io/biocontainers/fastqc:0.11.9--0 \
--container-mounts="/storage1/fs1/c.a.miller/Active/bfx-workshop-scratch:/storage1/fs1/c.a.miller/Active/bfx-workshop-scratch" \
/bin/bash
```

## 7. Finally running FASTQC on the cluster

Now that we've run the above command, we finally have all the pieces in place: We have the tool we need (via the docker container) and the data is mounted. 

Let's move to the right directory:

```
cd ~/workshop/week03/
```

Now let's run FastQC on our data:

```bash
fastqc *.fastq.gz
```

You'll see progress messages as it works through each file. When it finishes, look at what
it produced:

```bash
ls
```

```
2891351066_1.fastq.gz  2891351066_1_fastqc.html  2891351066_1_fastqc.zip
2891351066_2.fastq.gz  2891351066_2_fastqc.html  2891351066_2_fastqc.zip
```

For each input file, FastQC creates:

- an **HTML report**, with a user-friendly summary and plots
- a **zip file** with the same results as plain text files, which you can dig into if you need more
  detail or want to parse the results with your own scripts

### Question 6

In the command `fastqc *.fastq.gz`, what is the asterisk doing?

<details>
<summary>Hint</summary>

Try running `ls *.fastq.gz` and see what gets printed.

</details>

<details>
<summary>Solution</summary>

The `*` is a **wildcard** (or "glob") that matches any string of characters. Importantly,
it's the **shell** that expands it, *before* the command ever runs. So `fastqc` never sees a
`*` at all. What actually gets run is:

```bash
fastqc 2891351066_1.fastq.gz 2891351066_2.fastq.gz
```

That's why `ls *.fastq.gz` prints the list of file names. It's a handy trick for checking
what a wildcard will match before you use it in a command that deletes or overwrites files!

</details>

Now exit the job:

```bash
exit
```

You're back on the login node, and the container is gone. Try `fastqc --version` again to
convince yourself (it'll say 'not found'). The output files, however, are still in your directory, because they were written to your storage directory, not somewhere inside the container.

---

## 8. Looking at the results

The HTML reports are meant to be opened in a web browser, and there's no browser on the cluster.
We need to copy them to our laptop.

The `scp` ("secure copy") command copies files over the same connection that `ssh` uses. Open a
**new terminal window on your laptop** (not logged in to the cluster) and run:

```bash
scp "<washukey>@c2-login-002.ris.wustl.edu:~/workshop/week03/Exome_Tumor/*.html" .
```
(don't forget to replace \<washukey\> with your own username!)

The format is `scp <from> <to>`, and the `.` means "the directory I'm in right now".
The quotes stop your laptop's shell from trying to expand the `*` itself, so the wildcard gets
expanded on the cluster instead.

Then open the reports:

- **Mac:** `open 2891351066_1_fastqc.html`
- **Windows:** `start 2891351066_1_fastqc.html`
- or just open your browser, choose File > Open File and navigate to where you copied the files


### Question 7

Look through both reports. Do you see any potential issues with the sequence data?

<details>
<summary>Hint</summary>

The summary panel on the left gives each module a green check (pass), orange exclamation point
(warning), or red X (fail).

</details>

<details>
<summary>Solution</summary>

Both files pass every module except **Per sequence GC content**, which gets a red X.

FastQC compares the GC content of your reads to a smooth, bell-shaped curve that you'd expect
from random fragments of a genome. But this is **exome** data: the library was enriched for
protein-coding regions, which are more GC-rich than the genome as a whole (the Basic Statistics
module shows 51% GC, compared to about 41% for the human genome), and the mix of captured
regions makes for a lumpier distribution than the theoretical curve. So this "failure" is
**expected** for exome data and isn't a red flag by itself. If you saw a strong second peak in
*whole-genome* data, though, you'd want to think about contamination (bacteria, for example).

A few other things worth noticing:

- **Per base sequence quality** is almost perfectly flat. Remember the quality lines full of `?`
  characters? `?` is a quality score of 30, and nearly every base has exactly that score. Many
  sequencers and pipelines "bin" quality scores into just a few values to save space, so the
  plot looks unnaturally smooth.
- **Adapter content** and **overrepresented sequences** both pass, so there's no sign that
  the reads need adapter trimming.
- **Sequence duplication** looks fine, but keep in mind that with only 25,000 reads, you're
  unlikely to see many duplicates. The same library sequenced deeply might tell a different story.

Overall: this looks like good quality data!
</details>


---
## 9. Cheat sheet

```bash
# --- download ---
curl -O <url>                     # same thing, with curl

# --- archives and compression ---
tar -tvf archive.tar              # list what's inside a tar archive
tar -xvf archive.tar              # extract
tar -czvf dir.tar.gz dir/         # bundle and compress a directory
gunzip file.gz                    # decompress
gzip file                         # compress
zcat file.gz | head               # peek inside without decompressing

# --- FASTQ ---
head -n 12 reads.fastq                                  # first 3 records
wc -l reads.fastq                                       # lines (÷ 4 = records)
awk 'NR % 4 == 2 {print length($0)}' reads.fastq        # length of every read

# --- get an interactive shell in a container ---
srun -A compute2-workshop -p workshop --pty \
  --container-image=<image> \
  --container-mounts="/storage1/fs1/c.a.miller/Active/bfx-workshop-scratch:/storage1/fs1/c.a.miller/Active/bfx-workshop-scratch" \
  /bin/bash                  

# --- copy files to your laptop (run ON your laptop, don't forget the final ".") ---
scp "<washukey>@c2-login-002.ris.wustl.edu:<path>" .   
``` 

---

