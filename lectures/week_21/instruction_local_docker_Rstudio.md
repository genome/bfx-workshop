# Instructions for running Rstudio locally

To run the docker image with metagenomics packages pre-installed on your laptop/desktop:

1) First pull the docker image

```docker pull brusconi/rstudio_workshop```

2) Pull down the week 23 data/scripts and unzip them:

```
wget -v -O week_23.tar.gz -L https://wustl.box.com/shared/static/8g2spm9vn8agzlkp2h8l4we1ssfj381p
tar -xzvf week_23.tar.gz
```

3) Mount the path where your data is stored when launching your docker image. 

```
docker run --rm -p 8888:8787 -v /path/to/week_23:/home/rstudio/data -e PASSWORD=password rstudio_workshop
```

4) Now open a browser like chrome and navigate to the URL `http://localhost:8888`

This will open the rstudio session and you can then login with 

```
Username: rstudio
Password: password
```

5) Use RStudio to open the Rmd files from your directory (that you mounted at `/home/rstudio/data` and step through the details.
