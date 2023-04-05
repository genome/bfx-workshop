# Instructions for running RStudio on the RIS cluster 

To run the docker image with metagenomics packages pre-installed on the RIS cluster, you'll need to set up a personal R library to which R can install packages, and edit your .Rprofile to tell R where to store that data. (If you already have this set up, then skip steps 1 and 2)


1) Create a directory to store your libraries, replacing paths as necessary to point to your allocation:

```mkdir /storage1/fs1/${STORAGE_ALLOCATION}/Active/R_libraries/```

2) edit your .Rprofile using the text editor of your choice (e.g. `vim ~/.Rprofile`) and add this code:

```
devlib <- paste('/home/c.a.miller/lib/R',paste(R.version$major,R.version$minor,sep="."),sep="")
if (!file.exists(devlib))
  dir.create(devlib)
x <- .libPaths()
.libPaths(c(devlib,x))
rm(x,devlib)
```

3) if you haven't already, pull down the week 23 data/scripts and unzip them:

```wget -v -O week_23.tar.gz -L https://wustl.box.com/shared/static/8g2spm9vn8agzlkp2h8l4we1ssfj381p
tar -xzvf week_23.tar.gz
```

4) Now you can start the rstudio docker image by setting a password of your choosing (don't use your wustl key pass or anything - just something that prevents others from dropping in to your container)

```export PASSWORD=password```

Then launch a job that starts up Rstudio, replacing username/group params as necessary:

```
LSF_DOCKER_VOLUMES='/storage1/fs1/${user}/Active:/storage1/fs1/${user}/Active /scratch1/fs1/${user}:/scratch1/fs1/${user} /home/${user}:/home/${user}' PATH=/home/${user}:$PATH LSF_DOCKER_PORTS='8080:8080' bsub  -Is -n 16 -M 16GB -G compute-${user} -q general-interactive -R 'select[mem>16G,port8080=1] rusage[mem=16G]' -a 'docker(brusconi/rstudio_ris)' /bin/bash
```


5) You willl see something like this appear on the terminal that tells you which blade your job has started on: 

```<<Starting on compute1-exec-187.ris.wustl.edu>>```

Use that info to open up a rstudio session in a browser (chrome, firefox, etc). Just substitute in that value to this URL: 
`https://compute1-exec-187.compute.ris.wustl.edu:8080/vnc.html`

The password is whatever you set it to above. 

6) in the terminal type `rstudio`. That will launch an rstudio session connected to your storage that you called in with the LSF_DOCKER_VOLUMES. 

7) Use RStudio to open the Rmd files from your directory (that you mounted at `/home/rstudio/data` and step through the details.

Interactive sessions only last for 24h so remember to save your work!
