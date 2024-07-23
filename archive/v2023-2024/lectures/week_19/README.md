# Week 19: Genome and Pangenome Assembly

- [Lecture Recording](https://wustl.box.com/s/ebyw0832qhmp7ck19vo0c84gdxk2yiti)
- [Slides](Genome-Pangenome-Assembly.pdf)

## Exercise

### Graph Construction

For this week's exercise, we will be constructing a pangenome graph using the [Minigraph-Cactus pipeline](https://github.com/ComparativeGenomicsToolkit/cactus/blob/master/doc/pangenome.md). We will focus on a small region of chromosome 8, with sequences pulled from the T2T-CHM13 linear reference and from the phased assembly of HG00621, one of the individuals used in the Human Pangenome Reference Consortium (HPRC) pangenome reference.  

First, make a working directory. Download the following 2 files and place them in this directory:
- [Sequence file](https://storage.googleapis.com/icts-precision-health-bfx-workshop-public-data/pangenome/seqfile)
- [Script](https://storage.googleapis.com/icts-precision-health-bfx-workshop-public-data/pangenome/mcgb.sh)

Next, make a directory called `data` within the working directory, and place the following 3 files in that directory:

- [CHM13 fasta](https://storage.googleapis.com/icts-precision-health-bfx-workshop-public-data/pangenome/data/CHM13chr8.fasta)
- [HG00621 maternal fasta](https://storage.googleapis.com/icts-precision-health-bfx-workshop-public-data/pangenome/data/HG00621chr8_mat.fasta)
- [HG00621 paternal fasta](https://storage.googleapis.com/icts-precision-health-bfx-workshop-public-data/pangenome/data/HG00621chr8_pat.fasta)

Now, from the top level of the working directory (not inside the `data` directory), we can launch our construction docker container:
``docker run -it --rm --cpus 1 --memory 8589934592 --memory-swap 8589934592 -v `pwd`:/data mgibio/cactus:2.5.0-focal /bin/bash``
Inside the container, start the alignment and graph construction by running the provided script:
`bash mcgb.sh`
This will take about 10 minutes to run. After the script completes, you can exit the container. All outputs will be found in the `out` directory. 

#### Script overview
The provided script launches the minigraph-cactus pipeline. There is a wrapper command that runs the entire pipeline (`cactus-pangenome`), but due to some occasional bugs, we find it more reliable to manually run each command in the pipeline. Let's go through a high-level overview of each step:
1. `cactus-minigraph`: this uses minigraph to progressively build an initial graph, starting from the reference assembly and aligning large syntenic chunks from each additional assembly in the order provided. The resulting graph will only contain large SVs.
2. `cactus-graphmap`: this uses minigraph to map each assembly back to the graph constructed in the previous step
3. `cactus-graphmap-split`: this splits the assemblies and the mappings from the previous step by chromosome. Optional, but reduces memory in the next steps, which is especially important for our purposes, since we are running this exercise locally with limited RAM. 
4. `cactus-align`: this combines the mappings from the previous step into a multiple genome alignment, then converts that into a cactus graph. See the [Minigraph-Cactus paper](https://doi.org/10.1038/s41587-023-01793-w) for more information about these structures.
5. `cactus-graphmap-join`: this step runs several post-processing steps, such as normalizing, clipping, and filtering, to produce the final output graph. Also produces indices. The exact processing steps performed and indices generated may be adjusted based on a variety of flags as appropriate for the desired downstream analyses.

Note that some of these commands are themselves wrappers around multiple steps. For more information, see the [Minigraph-Cactus documentation](https://github.com/ComparativeGenomicsToolkit/cactus/blob/master/doc/pangenome.md#pipeline) and the [Minigraph-Cactus paper](https://doi.org/10.1038/s41587-023-01793-w).

#### Question
Look at the outputs in the `out` directory. Which files are indices, and which are graph files? What are some of the differences between the various graph file types?
- Hint: Try looking for the file extensions in: 
        - [The Minigraph-Cactus documentation](https://github.com/ComparativeGenomicsToolkit/cactus/blob/master/doc/pangenome.md#output)
        - [VG file format list](https://github.com/vgteam/vg/wiki/File-Types)- a fairly comprehensive overview of graph-specific file types
                - The "File Formats" and "Index Types" pages linked at the top of this page may also be useful- they describe some of the common formats in more detail
        - Note that `.gz` and `.tgz` extensions mean that a file has been compressed; you can ignore them when trying to determine the type of a file from its extensions

### Visualization
Download [Bandage](https://rrwick.github.io/Bandage/).
Unzip the file we will be loading: `gunzip -c bschr8.gfa.gz > bschr8.gfa`
Open Bandage, navigate to `File > Load graph`, and select `bschr8.gfa`
Under `Graph drawing` on the left tool bar, select `Entire graph` from the drop own, then click `Draw graph`. This may take a few minutes to load. 
Once the graph has loaded, under `Find nodes` on the right tool bar, enter `8737` next to `Node(s)`, then click `Find node(s)`.
Use the `+` and `-` keys to zoom in/out, and the scroll bars in the window or the scroll on your mouse to move around. You can also click and drag nodes (the colored line segments, representing sequences) to change their position.

#### Question
What seems to be happenening at this location? What would this look like in a linear reference, and what effect could this have on alignments at this locus?

## For-credit student homework:
Please send in a screenshot from Bandage of the region surrounding node `8737`, as well a brief answer to the questions in the [Graph Construction](#graph-construction) and [Visualization](#visualization) sections.
