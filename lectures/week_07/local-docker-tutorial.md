# Building a Docker image

1. Create a directory called `gatk`, and create a file called `Dockerfile` within it. Make the contents look like this:

NOTE: `$VERSION` is a placeholder, ie. a variable, that we either need to set in a shell, find/replace the text, or type a reasonable value in the appropriate command or file.
```
FROM broadinstitute/gatk:$VERSION

RUN apt-get update && apt-get install -y libnss-sss && apt-get clean all
```

As is, this just adds a package for identity management to the GATK container

2. Update the `Dockerfile` to use the most recent version of GATK. You'll find that info on on the [gatk dockerhub page](https://hub.docker.com/r/broadinstitute/gatk/)

3. Maybe we also want to do a little filtering of our VCFs.  Add in the `vcfpy` package. Typically, you'd install it using `pip install vcfpy`.

4. Grab the depth filter script from this location and include it in your image:
https://raw.githubusercontent.com/genome/docker-depth-filter/master/depth_filter.py

Can you think of two different ways to include this script?  (hint: RUN/COPY)

5. When your Dockerfile is all set up, you're ready to execute the instructions and build your container.  Make sure:
- Docker is running on your computer (you'll see it in the taskbar).
- You are in the directory just above the `gatk` folder.
- Run the following:
```
$ docker build gatk
```

The outputs of your steps will scroll by, and if all goes well, you'll see something like:
NOTE: `$HASH` is another placeholder for the hash returned during your execution.
```
Successfully built $HASH
```

That $HASH value is the id of your container. You can now run your image using that value:

```
$ docker run -it $HASH
```
Voila! You're running in the enviroment that you just set up.  Play around if you wish, then type `exit` to close the container.

That hash is pretty hard to remember, though. "Tags" allow you to pick a name for your container that's easier to remember

```
$ docker build -t gatk:$VERSION gatk
```
You can tag an image with anything, but it's often useful to use the version number.

Now you can run the same image using:

```
$ docker run -it gatk:$VERSION

```

Later in the [Cloud Build](cloudbuild-docker-tutorial.md) tutorial you deploy the image to a Google Artifact Registry for use on Google Cloud or other compute clusters.

# Detecting Germline Variants

The starting point for this tutorial is the HCC1395 Normal BAM file that you generated in last week's [DNA Alignment Workflow Tutorial](../week_06/bfx_workshop_06_alignment.md).

1. run GATK HaplotypeCaller on this data, to produce a VCF file of germline variant calls.  Some hints:
 - you only need the three required parameters: `--input`, `--output`, and `--reference`
 - use the latest tagged GATK container from the broadinstitute dockerhub
 - get your reference file from the inputs.yaml that you used last week
2. Examine the VCF by using `less` or your favorite text editor.
 - Note how each field in the VCF has a definition in the header.
 - Find a variant that has at least 100 reads of support for the variant allele
 - Find a site that is an indel (insertion or deletion)
3. Use command line utilities that tells you how many variants were called (don't count the header line!)
4. Open IGV, make sure you have genome build hg38 set, then load the bam file and the VCF using "Open File".  go to the position of the first variant in your VCF.  Do you have high confidence that this variant is homozygous?
5. Now jump to the location of your variant with >100 reads of support. Are you more confident in this site or the previous one?
6. Click on the name of your VCF file in the lefthand panel, then hit "CTRL-F".  Note how it jumps to the next variant in the file.  Jump through until you find the following:
 - an insertion and a deletion. Note how each is represented in IGV.
 - a heterozygous variant
 - a variant that lies in the exon of a gene
