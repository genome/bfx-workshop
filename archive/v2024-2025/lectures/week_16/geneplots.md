## Plotting Gene/Protein expression data


Load a bunch of libraries that we'll need:

```
library(tidyr)
library(ggplot2)
library(dplyr)
library(gridExtra)
```
Using your browser, download some supplemental data that includes protein abundance values
https://leylab.shinyapps.io/Proteomic_and_Phosphoproteomic_Landscapes_of_AML/

Download the sheet, then save the LFQ abundance tab as a csv file.  

Use R to read in that csv file

```
lfq = read.csv("/tmp/protein_lfq.csv")
head(lfq)
```

Some things to notice:
1) There's a lot more fields than I need here, and that makes it kind of messy to look at.
2) R prefixed all of the sample names with X.  What's with that?

```
# read it in without converting names (be aware of limitations!)
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
head(lfq)
#also look at the dimensions
dim(lfq)
```

That's certainly not very human-readable, but it's "clean" in other ways that matter. 

Before plotting, let's add one more column containing annotations for the groups we want to compare - in this case Healthy Donor (normal) samples vs all other (cancer) samples

```
# put all samples in the group "Case"
prot.long$Group = "Case"
# change just a subset of them to be "Control"
prot.long[grepl("HD",prot.long$Sample),"Group"] = "Control"
head(prot.long)
```

Now that we've got our data wrangled, we're ready to plot. Create an empty list to store the multiple plots we're going to make:

```
p = list()
```

Create a list of genes that we're going to plot:

```
genelist = c("NPM1", "A1BG", "ABCA1","GOT2","MPO")

```
This list could be chosen in many ways - by p-value if we had run statistical tests, or by pathway, etc.

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

Oh no, we don't know which gene is which in these plots!


Let's use the `ggtitle()` function to add gene names as a title.  Replot all of them and regenerate the pdf:

```
p = list()

for(gene in c("NPM1", "A1BG", "ABCA1")){
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

