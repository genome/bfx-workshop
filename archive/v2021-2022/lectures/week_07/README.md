## Week 3: Docker fundamentals and Germline Variant Calling

- [Lecture Recording](https://wustl.box.com/s/23iyu9ntm24uvsjwhyyh4dzro13pa4l2)
- [Slides](bfx_workshop_06_docker_variant_calling_2021.pdf)


## Week 07 Assignment


### Building a Docker image


**1)** Create a directory called `gatk`, and create a file called `Dockerfile` within it. Make the contents look like this:

```
FROM broadinstitute/gatk:4.1.3.0

RUN apt-get update && apt-get install -y libnss-sss && apt-get clean all
```

As is, this just adds a package for identity management to the GATK container

**2)** Update the image to use the most recent version of GATK. You'll find that info on on the [gatk dockerhub page](https://hub.docker.com/r/broadinstitute/gatk/)

**3)** Maybe we also want to do a little filtering of our VCFs.  Add in the `vcfpy` package. Typically, you'd install it using `pip install vcfpy`.

**4)** grab the depth filter script from this location and include it in your image:
https://raw.githubusercontent.com/genome/docker-depth-filter/master/depth_filter.py

Can you think of two different ways to include this script?  (hint: RUN/COPY)

**5)** when your Dockerfile is all set up, you're ready to execute the instructions and build your container.  Make sure that a) Docker is running on your computer (you'll see it in the taskbar). b) you're in the directory just above the `gatk` folder, then run :

```
$ docker build gatk
```

The outputs of your steps will scroll by, and if all goes well, you'll see something like: 

```
Successfully built 746b8dc67a4e
```

That hash is the id of your container. You can now run your image using that hash:

```
$ docker run -it 746b8dc67a4e
```
Voila! You're running in the enviroment that you just set up.  Play around if you wish, then type `exit` to close the container.

That hash is pretty hard to remember, though. "Tags" allow you to pick a name for your container that's easier to remember

``` 
$ docker build -t gatk:4.1.9.0 gatk 
```
You can tag an image with anything, but it's often useful to use the version number. 

Now you can run the same image using:

```
$ docker run -it gatk:4.1.9.0

```

We'll cover how to deploy these images to a public repository (so that you can share them and use them on a compute cluster) in a later session.

### Detecting Germline Variants

The starting point for this exercise is the alignment files that you generated in last week's [alignment workflow exercise](https://github.com/genome/bfx-workshop/blob/master/archive/lectures/week_05/cromwell_alignment_walkthrough.md).  I'll assume that you've copied the `final.bam` and `final.bai` files to your working directory.

1) run GATK HaplotypeCaller on this data, to produce a VCF file of germline variant calls.  Some hints:

  - you only need the three required parameters: `--input`, `--output`, and `--reference`
  - use the latest tagged GATK container from the broadinstitute dockerhub
  - get your reference file from the inputs.yaml that you used last week

2) Examine the VCF by using `less` or your favorite text editor. 
 
   - Note how each filed in the VCF has a definition in the header. 
   - Find a variant that has at least 100 reads of support for the variant allele
   - Find a site that is an indel (insertion or deletion)

3) Use command line utilities that tells you how many variants were called (don't count the header line!)

4) Open IGV, make sure you have genome build hg38 set, then load the bam file and the VCF using "Open File".  go to the position of the first variant in your VCF.  Do you have high confidence that this variant is homozygous?   

5) Now jump to the location of your variant with >100 reads of support. Are you more confident in this site or the previous one?

6) Click on the name of your VCF file in the lefthand panel, then hit "CTRL-F".  Note how it jumps to the next variant in the file.  Jump through until you find the following:

- an insertion and a deletion. Note how each is represented in IGV.
- a heterozygous variant
- a variant that lies in the exon of a gene






