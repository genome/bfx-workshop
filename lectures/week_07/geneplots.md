## Plotting Gene/Protein expression data


Let's work through a real-world example of how we might pull data from a supplementary table of a manuscript, import it into R, and make some useful plots. 

First, load a bunch of libraries that we'll need:

```
library(tidyr)
library(ggplot2)
library(dplyr)
library(gridExtra)
```

Using your browser, go to 
[https://proteomics.leylab.org/](https://proteomics.leylab.org/) then click on **"Download the TMT, LFQ, and mRNA abundance data"** to pull down the excel spreadsheet. 

Open it in excel (or import to google sheets), then save the **lfq** tab as a CSV file (File > Save As). 

Now, use R to read in that csv file.  Change the path to the directory and filename where you saved the csv.

```
lfq = read.csv("/tmp/protein_lfq.csv")
head(lfq)
```

Some things to notice:

1. There are a lot more fields than we need here, and that makes it kind of messy to look at.
2. R prefixed all of the sample names with X.  What's with that? 

> [!NOTE] R really doesn't like using numbers as column names, because it could cause confusion.  (If I type `data[,21]` am I referring to the 21st column or the column named "21"?)  In this case, we'll override it, since our patient identifiers are numbers.
> 

```
# read it in without converting the column names 
lfq = read.csv("/tmp/protein_lfq.csv",check.names=F)
# Get rid of a bunch of columns that we don't intend to use here:
lfq = lfq[,1:51]
```

Now take a look at the nicely cleaned data:

```
head(lfq)
```

To me, this data looks great, but ggplot _hates_ this kind of data structure. It really wants it to be in what it calls a "tidy" format.  We can get there by pivoting the table.  Note the use of tidyverse pipe syntax:

```
prot.long = lfq %>% pivot_longer(cols = -Protein, names_to = "Sample", values_to = "Expression")
```

Look at the data now:

```
prot.long
#also look at the dimensions
dim(prot.long)
```

300k+ rows certainly isn't very human-readable, but it's "clean" in other ways that matter.

Before plotting, let's add one more column containing annotations for the groups we want to compare - in this case Healthy Donor (normal) samples vs all other (cancer) samples

```
# put all samples in the group "Case"
prot.long$Group = "Case"
# change just a subset of them to be "Control"
prot.long[grepl("HD",prot.long$Sample),"Group"] = "Control"
head(prot.long)
```

And do a quick sanity check to make sure we got it right:

```
table(prot.long$Group)
```

Now that we've got our data wrangled, we're ready to plot. Create an empty list to store the multiple plots we're going to make:

```
p = list()
```

Create a list of genes that we're going to plot:

```
genelist = c("NPM1", "A1BG", "ABCA1", "GOT2", "MPO")

```
This list could be chosen in many ways - by p-value if we had run statistical tests, or by a pathway of interest, etc.

Now, loop through those genes, and create a plot for each one:

```
for(gene in genelist){
	df = prot.long[prot.long$Protein==gene,]

   p[[gene]] = ggplot(df, aes(x=Group, y=Expression, fill=Group)) +
        geom_boxplot(outlier.shape = NA) + 
        geom_jitter(width=0.2) 
}
 
```

View an individual plot in RStudio by selecting that item of the list:

```
p[["GOT2"]]
```

use marrangeGrob to tile plots 1x2 across multiple pages:

```
ggsave(filename = "boxplots.pdf",
       plot = marrangeGrob(p, nrow=1, ncol=2, top=NULL),
       width = 8, height = 4)
```

Open that PDF from the place you saved it (use `getwd()` in R if you're not sure). 

Oh no, we don't know which gene is which in these plots!

Let's use the `ggtitle()` function to add gene names as a title.  Replot all of them and regenerate the pdf:

```
p = list()

for(gene in genelist){
	df = prot.long[prot.long$Protein==gene,]

   p[[gene]] = ggplot(df, aes(x=Group, y=Expression, fill=Group)) +
        geom_boxplot(outlier.shape = NA) + 
        geom_jitter(width=0.2) +
        ggtitle(gene)
}
ggsave(filename = "boxplots.pdf",
       plot = marrangeGrob(p, nrow=1, ncol=2, top=NULL),
       width = 8, height = 4)
 
```

## Get an LLM assist.

Describe this plotting code to your favorite LLM.  Ask it to help you modify ggplot to create plot that looks more publication-ready: 

- Make the boxplots smaller
- put them all on one page in a single row
- add panel lettering (A through E)
- change the colors to green and blue
- make the jittered points triangles
- remove the grey background and use light grey lines across the y-axis marks

Verify that the code works. 

### Verify your understanding
**Go through the LLM-generated plotting code line by line, adding comments that explain what each part does.**

For example, if I was commenting my code above, it might look like:

```
# setting up the plot - x axis by group, y will be expression levels, and 
# we'll set the fill color to match the group
ggplot(df, aes(x=Group, y=Expression, fill=Group)) +
		 # create a boxplot, but don't add outlier points (since they'll show up 
		 # in the jitter below and we don't want to double-plot)
        geom_boxplot(outlier.shape = NA) + 
	 	 # add points for every individual sample, jittered sideways so 
	 	 # that the points don't overlap
        geom_jitter(width=0.2)

```

If you don't know what the LLM code does, ask it to explain, or tweak it until you do understand. The key here is, that you shouldn't just accept untrusted code - you have to verify that it's doing what you intended (what if the LLM decided to hide outliers, or only plotted half the data?)  

Ultimately, **you** are the one with responsibility for that figure in the paper, so you need to understand where it came from!

### For-credit students:
Turn in your publication-ready boxplot as well as the code that created it, commented with your explanations of what the code does. 