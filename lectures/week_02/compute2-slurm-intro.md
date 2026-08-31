# Getting Started on the RIS Compute2 Cluster: Logging In and Running Jobs with Slurm

This is a hands-on lesson for people new to using a compute cluster. By the
end you should be able to:

- Explain what a login node, a compute node, and a scheduler are, and why the distinction matters
- Connect to Compute2 over SSH
- understand the difference between STDERR and STDOUT
- Run a quick interactive job with `srun`
- Write and submit a batch job with `sbatch`
- Watch your jobs with `squeue`, `sinfo`, and cancel them with `scancel`
- Ask the scheduler for the CPUs, memory, time, and GPUs your work actually needs

> **Note:** We strongly suggest that you type out the commands yourself, rather than using cut and paste.  Doing so really helps enormously with retaining this information. 

---

## 1. Cluster architecture

When you use a cluster you are really interfacing with a few different things. New users sometimes get into trouble when they confuse them.

```             
     your laptop
         |      
         |  ssh 
         v   
  +--------------+        +------------------------------------------o---+
  |  LOGIN NODE  |        |            COMPUTE NODES                     |
  |  (c2-login-  | -----> |  c2-node-001, c2-node-002, ... c2-gpu-001 ...|
  |   001/2/3)   | Slurm  |  where your actual computation runs          |
  +--------------+        +----------------------------------------------+
     ^                                    ^
     |                                    |
  you type here                    the scheduler (Slurm)
  edit files, move data            decides which job runs
  submit jobs                      on which node, and when
```

**Login nodes** (`c2-login-001`, `c2-login-002`, `c2-login-003`) are the front door. They are
*not* for computation. You use them to edit scripts, move data around, and submit jobs. Usage
on a login node is capped - if you run a heavy program
directly on a login node, it will be killed, and you will be affecting everyone else logged
into that node.

**Compute nodes** are where real work happens. You don't SSH to them directly. Instead you ask the scheduler for a slice of one (or many), and it hands you the resources.

**The scheduler** is a program called **Slurm**. It keeps a queue of everyone's job requests,
matches each request against free resources, and runs jobs in a fair order. Every command in
the second half of this lesson (`srun`, `sbatch`, `squeue`, `sinfo`, `scancel`) is you talking
to Slurm.

> **Why do we need a scheduler?** We have hundreds of people sharing a few thousand CPU cores. If everyone
> just ran things whenever they wanted, work would collide, you'd have to spend a bunch of time looking for an open node, etc. 
> The scheduler lets you say "I need 4 cores, 16 GB of RAM, and 2 hours," and it finds
> a place for that, runs it, and moves on to the next person.

**Exercise 1.**  Where would each of these tasks belong — login node or
compute node?

1. Opening a 50-line Python script in `nano` to fix a typo
2. Aligning 200 million sequencing reads to a genome
3. Copying a 4 GB file from your home directory to a project directory
4. Training a neural network for six hours

---

## 2. Logging in with SSH


First, open a terminal program. For macs, that's probably `Terminal`.  On windows, you can open `Command Prompt`.  And if you're using Linux, you probably already know this.

You connect to a login node with `ssh`. Do so now, by typing

```bash
ssh <washukey>@c2-login-002.ris.wustl.edu
```

> **Note:** The three login nodes are interchangeable. If `c2-login-001` feels slow or you
> can't reach it, try `-002` or `-003`. Your files are the same on all three because your home
> directory and files are shared across the cluster.

You'll be greeted with a splash screen and then be sitting at a command prompt.  This should be familiar after last week's assignment.

![unix.gif](unix.gif)

Let's look around and try out a few commands:

```
ls
ls -al
echo "Hello world"
grep --help
```
Yep, looks very similar to our shell from last week.

## 3. Where is my data?

There's one more important piece to the cluster, and that's the disk storage array:

``` 
                                   ________________________
     your laptop                  |   Disk Storage Array   |
         |           _____________|________________________|
         |  ssh     /                    |
         v         /                     |
  +--------------+        +----------------------------------------------+
  |  LOGIN NODE  |        |            COMPUTE NODES                     |
  |  (c2-login-  | -----> |  c2-node-001, c2-node-002, ... c2-gpu-001 ...|
  |   001/2/3)   | Slurm  |  where your actual computation runs          |
  +--------------+        +----------------------------------------------+
     ^                                    ^
     |                                    |
  you type here                    the scheduler (Slurm)
  edit files, move data            decides which job runs
  submit jobs                      on which node, and when
```

If you get up from your laptop and move to someone else's laptop, the data on your hard drive isn't there any more. On the cluster, though, the primary storage array is independent of the individual nodes.  This is convenient, as the same data will be available, no matter which node you're on. 

Your **Home directory** is where you are now. Use the `pwd` command to see the full path (print working directory)

```
pwd -P
```

These home directories are very small on the cluster. They're meant to store config files, links, and maybe little scripts or programs, but NOT big data.

**Storage allocations** are assigned to and paid for individual labs.  These are where you'll do most of your work.  For this workshop, we're going to create a folder on an allocation we have set up. 

```
mkdir /storage1/fs1/c.a.miller/Active/bfx-workshop-scratch/USERNAME
```
(again replacing USERNAME)

Then, for convenience, we're going to create a softlink ("shortcut") to that spot in your home directory

```
ln -s /storage1/fs1/c.a.miller/Active/bfx-workshop-scratch/USERNAME workshop
```

Do an `ls` and you can see that there is now a "workshop" item listed. If you do `ls -l` to get more details, you can see how it is a link to the location we indicated (note the arrow ->).

**Unless otherwise specified, all of your course work should be done inside this workshop folder**.  (If you have access to your lab's storage allocation, you can also do it there)

---

## 4. STDERR, STDOUT, and pipes

Last week you learned to do things like run `grep` to search through a file or `head` and `tail` to look at the beginning and end of a file.  Let's level that up a bit and talk about the paths that data can take in and out of commands.

Start by moving into your workspace (remember, `~` is a handy shortcut for your home directory)

```
cd ~/workshop/
```

Next, make a copy of a data file that we've put in a shared location:

```
cp /storage1/fs1/c.a.miller/Active/bfx-workshop/2026-27/lectures/week02/tcga.tsv .
```
(what does the `.` stand for?)

### Viewing a file with the 'less' command
Let's start by looking at this file using the `less` command:

```
less tcga.tsv
```

It's tab delimited data, and you can see some things that make sense - genomic coordinates, gene names, etc. 

To make it easier to view, let's use two tricks:

1) type `-S` then hit enter - this wraps long lines
2) type `-x20` then hit enter - this makes the space between the columns wider

Finally, type `q` to exit the program.


> **Note:** There are shortcuts built in that help with the command line.  1) Use 
> the up and down arrows to scroll through your command history.  2) Use `Tab` 
> to autocomplete a filename or command.  These can save you a lot of typing and 
> reduce errors!

### Using 'grep' to extract parts of a file

Let's say we were interested in extracting NPM1 mutations from this file.  Using your knowledge from last week we could try:

```
grep NPM1 tcga.tsv
```
What happened?  That's way too much data to view comfortably and it all appeared on the screen at once.

To be more precise, the output of that grep command went to **STDOUT**.  Since we're working in a terminal, that means it gets spit out to the screen.  We can change that though!

To send that output to a file, we can use the `>` operator, which says _redirect this output to a file_.  

```
grep NPM1 tcga.tsv >npm1_variants.tsv
```
Now run `ls` to see the file you've created, then use `less` to view the contents (remember, `q` to exit). 

Let's find out how many lines are in our NPM1 file, using the `wc` (word count) command. Using the `-l` parameter tells it to count the lines in a file (as opposed to characters or bytes).

```
wc -l npm1_variants.tsv
```

Okay, it took a couple of steps, but we've got our answer - 55 NPM1 variants in the file.  There's an easier way, though, and it involves using pipes.

### Pipes redirect output

The `|` (pipe) operator lets you redirect the output of one command into the input of another. 

To demonstrate it, let's count the NPM1 variants again:

```
grep NPM1 tcga.tsv | wc -l
```
What happened there?  We took the output from our `grep` command, and redirected it so that instead of going to the screen or a file, it went into the next command in the chain. Then that command's output was left undirected, so it came back to the screen.  

**Pipes are an incredibly powerful concept**

For example, if we wanted to get fancy, we could `cut` out the list of gene names, count unique instances of each, sort them by incidence, and then output the top 20 most frequently mutated genes:

```
cut -f 8 tcga.tsv | sort | uniq -c | sort -nrk 1 | head -n 20
```

This takes time to get fluent in, and we'll continue working on it throughout the course.   

**Exercise:** Ask your favorite AI model to explain the above command to you. Using that information, modify it to show you the top 10 most frequent **types** of mutations from that file (e.g. `missense`, `frame_shift_ins`)

### Understanding STDERR
So far we've talked about STDOUT and pushing it to files or to other programs, but there's one more path we need to understand, called STDERR.  It's usually used for error messages from tools, which you might not want mixed in with your output (where they're hidden and may be missed!)

To demonstrate, let's use the `date` command


By default, date outputs to STDOUT, and by default, that goes to your screen.

```
date
```

You can also redirect that output (STDOUT) to a file

```
date >today.txt
```

If we give date an invalid parameter, the resulting error will output to STDERR:

```
date --asdf
```

Errors do not go to STDOUT:

```
date --asdf >today.txt
```
You still see the error — which is what you want here.

If you want to log errors, you can redirect them to a file using `2>`

```
date --asdf >today.txt 2>today.err
```
----

## 5. Slurm: let's try it out

Slurm is the software that provisions compute for you. You tell it what to run, how many resources you need, and it handles the rest.  

Before we submit our first job, let's introduce two concepts:


### Accounts - pay to play

Every job is **charged to an account**. If you don't specify one, your job will be rejected. You pass it with `-A` / `--account`:

```
-A, --account=<name>     charge job to specified account
```

You may belong to more than one account (for example, one per lab or per grant). 

To start the year, RIS has graciously added everyone to the `compute2-workshop` group and agreed to cover the costs. 

### Partitions (queues)

A **partition** is a named pool of compute nodes with its own rules (time limits, node types,
priority). You choose one with `-p` / `--partition`. See what exists:

```bash
sinfo
```

Example output:

```
PARTITION       AVAIL  TIMELIMIT  NODES  STATE  NODELIST
bigmem             up   infinite     82   idle  c2-bigmem-[001-002],c2-node-[001-080]
general-cpu        up   infinite     80   idle  c2-node-[001-080]
general-short      up      30:00     80   idle  c2-node-[001-080]
gpu                up   infinite     94   idle  c2-gpu-[001-016],c2-node-[001-080]
```

Reading this:

- **`general-cup`** — the default partition. General-purpose CPU work.
- **`general-short`** — same nodes, but a short 30-minute time limit; jobs here often start
  sooner because the scheduler knows they'll finish quickly.
- **`bigmem`** — nodes with much more RAM, for memory-hungry jobs.
- **`gpu`** — nodes that have GPUs. (GPUs are also reachable from the general and general-short
  partitions)

### Our first Slurm job

Okay, with that info, let's run a job:

```
srun -A compute2-workshop -p general-short echo "Hello World!"
```

That didn't look like anything happened abnormally, but there was a lot happening behind the scenes:

- `srun` shot off a request that said "find me a node in the cluster that's got free resources to run this job"
- Slurm found a free CPU on a compute node and gave it to you for the duration of the command.
- When `echo` exited, the reservation was released automatically.

So how did Slurm know how many processors we needed to assign, or how much memory?  

### Resource requirements

By default a job gets 1 CPU and 4GB of RAM.  That's the standard "Unit" of consumption of the cluster.

Most times, we'll want to specify exactly what resources to use, using the resource flags:

| Option | Long form | Meaning |
|---|---|---|
| `-c` | `--cpus-per-task=<n>` | CPU cores per task (default 1) |
| `-N` | `--nodes=<n>` | number of nodes to spread across |
| | `--mem=<mem>` | total real memory for the job |

Okay, let's try our command again, giving it double the resources:

```
srun -A compute2-workshop -p general-short -c 2 -N 1 --mem=8G echo "Hello World!"
```

## 6. Interactive jobs

Sometimes you don't want to run one command — you want a shell *on* a compute node to poke around, test things, or run a short analysis interactively. 

To do that, we can ask `srun` for a pseudo-terminal (`--pty`) running `/bin/bash` (a shell).

```bash
srun -A compute2-workshop -p general-short -c 1 -N 1 --mem=4G --pty /bin/bash
```

How do you know anything happened? Try running the `hostname` command. It should output something like `c2-node-001`.  Then let's get out of our interactive job and try it again:

```
$ hostname
c2-node-001

$ exit
exit

$ hostname
c2-login-003
```

Typing `exit` dropped us out of the interactive job and back to our shell. 


> **Note:** An interactive shell holds a compute allocation for as long as you leave it open,
> even while you're doing other things.  You're paying for that time! Don't forget to `exit`. For long unattended work, see the next section

---

## 7. Batch jobs: `sbatch`

Interactive jobs are great for exploration, but when we do big jobs, we want them to run unattended in the
background, so that you can log off and come back later. That's `sbatch`.

A batch job is just a shell script with special comment lines at the top. Slurm reads the
`#SBATCH` lines as if they were command-line options, then runs the rest of the script on a
compute node.

We'll talk another day about ways to get files back and forth from the cluster, but for now, we're going to use a text editor called nano.  Let's use nano to create a file called `testjob.slurm`:

``` 
nano testjob.slurm
```
Copy/paste the following lines into the file.  (Windows users, you may have to use `CTRL-Shift-V` instead of just `CTRL-V`)

```bash
#!/bin/bash
#SBATCH --output=test.log
#SBATCH -p general-short
#SBATCH -A compute2-workshop
#SBATCH -t 10
#SBATCH --cpus-per-task=1
#SBATCH --mem=1000

echo "Hello World!"
```
Then hit `CTRL-S` to save and `CTRL-X` to exit

These options look familiar - they're the same options that we gave to `srun` earlier, with two new ones:

- `-t 10` - tell the scheduler to kill your job after X minutes
- `--output=test.log` - provide somewhere for the output of your job to go

Okay, let's submit this job:

```
sbatch testjob.slurm
```

You'll see output that looks like `Submitted batch job 2947946`

Again, this is going to run very quickly, so after just a second or two, you should be able to run `ls` and see something like this:

```
$ ls
test.log  testjob.slurm
```

and then we can use `cat` or `less` to see that the content of that file is just what we expected.  

> **Note:** For very simple commands, you can specify the requirements and job straight from the command line, using the `--wrap` parameter: 
> 
> ```
> sbatch -A compute2-workshop -p general-interactive -t 10 \
> -c 1 -N 1 --mem=4G --output=test.log --wrap="echo \"Hello World\""
> ```
> Note that we have to escape those quotations to get them through

---

## 8. Watching the queue: `squeue` and `sinfo`

`squeue` shows jobs currently in the system. 

Let's launch a job so we have something running on the cluster
```
sbatch -A compute2-workshop -p general-interactive -t 10 \
-c 1 -N 1 --mem=4G --output=test.log --wrap="sleep 60"
```
This command will just pause (or sleep) for 60 seconds, then exit. 

Now run `squeue` to see the details

```bash
squeue 
```

```
             JOBID PARTITION     NAME     USER ST       TIME  NODES NODELIST(REASON)
           2948249 general-i     wrap c.a.mill  R       0:04      1 c2-node-001
```           
 
Key columns:

- **ST** — state. `R` = running, `PD` = pending (waiting for resources), `CG` = completing.
- **TIME** — how long it has been running.
- **NODELIST(REASON)** — which node(s) it's on, or *why* it's still pending (e.g. `Priority`,
  `Resources`).


When you get a few dozen jobs running, it can be helpful to give them a job name so that you can tell them apart.  The `-J` parameter takes care of that.

```
sbatch -A compute2-workshop -p general-interactive -t 10 \
-J myjob1 -c 1 -N 1 --mem=4G --output=test.log --wrap="sleep 60"
```

Now run `squeue` and note that `myjob1` shows up as the name. 

### I forgot where to submit jobs
`sinfo` is a companion command — `squeue` tells you about **jobs**, `sinfo`
tells you about **nodes and partitions** (what's idle, what's allocated, what's down).  It's a useful command to see what queues are available 

### I forgot which details go with which job
Again, submit an example job so that we have something to work with:

```
sbatch -A compute2-workshop -p general-interactive -t 10 -J myjob1 -c 1 -N 1 --mem=4G --output=test.log --wrap="sleep 60"
```

Now look up its ID using `squeue`

```
$ squeue
             JOBID PARTITION     NAME     USER ST       TIME  NODES NODELIST(REASON)
           2948291 general-i   myjob1 c.a.mill  R       0:07      1 c2-node-001
```

Finally, we can see all the gory details by using `scontrol`:

```
scontrol show job 2948291
JobId=2948291 JobName=myjob1
   UserId=c.a.miller(1600109) GroupId=domain users(1000070) MCS_label=N/A
   Priority=1 Nice=0 Account=compute2-workshop QOS=compute2-workshop
   JobState=RUNNING Reason=None Dependency=(null)
   Requeue=1 Restarts=0 BatchFlag=1 Reboot=0 ExitCode=0:0
   RunTime=00:00:20 TimeLimit=00:10:00 TimeMin=N/A
   SubmitTime=2026-08-30T21:26:33 EligibleTime=2026-08-30T21:26:33
   AccrueTime=2026-08-30T21:26:33
   StartTime=2026-08-30T21:26:33 EndTime=2026-08-30T21:36:33 Deadline=N/A
   SuspendTime=None SecsPreSuspend=0 LastSchedEval=2026-08-30T21:26:33 Scheduler=Main
   Partition=general-interactive AllocNode:Sid=c2-login-003:2955780
   ReqNodeList=(null) ExcNodeList=(null)
   NodeList=c2-node-001
   BatchHost=c2-node-001
   NumNodes=1 NumCPUs=1 NumTasks=1 CPUs/Task=1 ReqB:S:C:T=0:0:*:*
   ReqTRES=cpu=1,mem=4G,node=1,billing=1
   AllocTRES=cpu=1,mem=4G,node=1,billing=1
   Socks/Node=* NtasksPerN:B:S:C=0:0:*:* CoreSpec=*
   MinCPUsNode=1 MinMemoryNode=4G MinTmpDiskNode=0
   Features=(null) DelayBoot=00:00:00
   OverSubscribe=OK Contiguous=0 Licenses=(null) LicensesAlloc=(null) Network=(null)
   Command=(null)
   WorkDir=/rdcw/fs2/c.a.miller/Active/bfx-workshop-scratch/c.a.miller
   StdErr=/rdcw/fs2/c.a.miller/Active/bfx-workshop-scratch/c.a.miller/test.log
   StdIn=/dev/null
   StdOut=/rdcw/fs2/c.a.miller/Active/bfx-workshop-scratch/c.a.miller/test.log
   TresPerTask=cpu=1
```

---

## 9. Stopping jobs: `scancel`

To cancel a job, give `scancel` its job ID (from `sbatch` or `squeue`):

```bash
scancel 12345
```

To cancel **all** of your jobs at once:

```bash
scancel --me
```

You can also cancel by attribute:

```bash
scancel -n python-test          # all jobs named python-test
scancel -u <washukey>      # all jobs belonging to this user
scancel -p general              # all your jobs in the general partition
```

Cancelling works the same whether the job is still pending or already running.

---

## 10. Cheat sheet

```bash
# --- connect ---
ssh <washukey>@c2-login-001.ris.wustl.edu     # + password + Duo push

# --- look around ---
sinfo                             # partitions and node states
squeue --me                       # my jobs
squeue -p general                 # jobs in a partition

# --- run something now (interactive) ---
srun -A <acct> -p general python script.py
srun -A <acct> -p general --pty /bin/bash          # interactive shell on a node
srun -A <acct> -p general --gpus=1 python gpu.py   # with a GPU

# --- run something unattended (batch) ---
sbatch testjob.slurm              # submit; prints a job ID

# --- stop something ---
scancel <jobid>                   # cancel one job
scancel --me                      # cancel all my jobs
```

Minimal `testjob.slurm` template:

```bash
#!/bin/bash
#SBATCH --job-name=my-job
#SBATCH --account=compute2-workshop
#SBATCH --partition=general-cpu
#SBATCH --time=60
#SBATCH --cpus-per-task=1
#SBATCH --mem=4000
#SBATCH --output=my-job_%j.log
```

(`%j` in the output filename is replaced with the job ID, so multiple runs don't overwrite
each other.)

----
## 11. Help!  I'm stuck!

- If your terminal isn't responding at all, try `Ctrl-C` to interrupt/kill a running process.   Commonly occurs when you run a command without expected arguments.  `grep NPM1` will sit there waiting forever because you didn't specify which file to search for NPM1 in!

- many interactive commands can be quit by typing `q` (like `less`).

- if you're in an interactive node, typing `exit` will get you back to the head node.  If you're on the head node, typing `exit` will drop you out of SSH and back onto your laptop. 

- If you're still stuck, ask an LLM for help - they're really good companions for debugging. Be sure to ask them to help you understand _why_ the error happened, not just how to fix it blindly

- Ask in our #bfx-workshop slack channel (you should have all gotten an invite by now).  We hang out in there regularly and help folks out. 
 
 
----

## 12. Assignment for this week

1. Visit the official RIS documentation site at https://washu.atlassian.net/wiki/spaces/RUD/overview.  Find the Compute2 Quickstart page, which has more details on all of this Slurm stuff, as well as links to other useful information

2. Ask an LLM to help you come up with a command or small script that will help you find all mutations in tcga.tsv that appear on chromosome 5 and are missense. Do NOT load the file into the chat (what if this file contained sensitive information!?).  What info did you need to communicate (column headers or numbers, data types, etc)?  Can you force it to use tools that are familiar to you from this lesson and last (instead of writing a script in python or rust).  How can you double check that the results of the answer you arrived at was correct? 
 
3. Right now, the `tcga.tsv` file we used above is ordered by patient id (UPN).  Let's say we wanted to look for mutation hotspots in this file.  It would be useful to be able to sort it by chromosome, then position.  Use the manual page `man sort` and your favorite AI to help you figure out the right command to do that chromosome, then position sort, then save it into your directory as `tcga_possorted.tsv`

4. Start an interactive job with one of the commands above. Run `hostname` and `nproc`. What blade did you land on? how many CPUs do you see? Compare with the same commands run on the login node.

5. Submit a batch job to the cluster that accomplishes the same thing as #4. 

6.  What did the job actually use? After your batch job from #5 finishes, `scontrol show job <id>` stops working — it only sees running jobs. Run `

```
sacct -j <jobid> --format=JobID,JobName,State,Elapsed,ReqMem,MaxRSS,ExitCode
````
How much memory did the job really use versus the 4GB default? What would you request next time?


7. We talked about STDOUT and STDERR above.  By default Slurm mixes the two into your ouput file.  Try submitting a batch job that runs the command:

```
cat tcga.tsv tcga_missing.tsv
```
Note that `tcga_missing.tsv` doesn't exist!

Open the output of the job with `less` and scroll through the output (the space bar acts like page down and jumps down a full page).  Is it obvious that something went wrong?

8.  Slurm has the ability to separate those STDOUT and STDERR streams by adding the `--error <filename>` flag.  Re-run the command from #6, adding this flag to catch any warnings or errors from STDERR.  Save it to `tcga_output.err`. 


9. Since we're going to be using this directory every week, let's clean it up when you're done.  Make a new directory called `week02`.  Move all of the files we created today into that folder.  (hint - use the wildcard character `*` to make it faster)


**For-credit students:** Paste the following into an email to John as proof of completion:

- the command you used in #2 and the resulting count
- The batch command you used in #5
