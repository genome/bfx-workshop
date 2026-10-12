# Data munging and cleaning with R

Adapted from the [NCEAS dataviz module](https://nceas.github.io/oss-lessons/dataviz-and-interactive-tools/module-1-ggplot2.html) and [R for Data Science](https://r4ds.had.co.nz/data-visualisation.html)

## Learning outcomes

Students should

- Learn what it means to examine ("QC") and clean their data before analyzing it
- Learn to recognize and fix common issues with datasets
- Learn to use an LLM as a helper for this kind of work, and to check its output critically

## A note on using LLMs (ChatGPT, Claude, Gemini, Copilot, etc.)

LLMs can be very useful when you're learning R, and we'll work in concert with one at the bottom of this exercise.  As you go through the top half, there shouldn't be much need for LLM use. If you do get stuck and turn to an LLM for help, remember the following:

- **Try it yourself first.** The point of this exercise is for *you* to learn how to look at
  data. Use the LLM to get unstuck, to explain things, or to check your thinking, not to skip
  the thinking.
- **Verify everything.** LLMs sometimes make up function arguments, miscount things, or
  confidently describe output that doesn't match what R actually printed. Run the code
  yourself and check that the output makes sense.
- **Be careful what data you share.** This exercise uses a made-up dataset, so it's fine to
  paste it in. Never paste patient data, controlled-access genomic data (e.g. from dbGaP), or
  anything else sensitive into an LLM unless your institution has approved that tool for
  that kind of data.

## Lesson

Let's look at an example dataset and go over some ways
we might QA/QC this dataset.

It is an arbitrary temperature study. The **experimental design** was:

- There are **6 sites**, labeled `A` through `F`
- Within each site, there are **6 plots**, numbered `1` through `6`
- We took **10 temperature samples** (in degrees C) **per plot**

So if everything went perfectly, we'd expect 6 sites × 6 plots × 10 samples = **360
measurements**, with **60 per site** and **10 per site/plot combination**.

Writing down what you *expect* to see is often the most important step in QC. Most data problems
are only visible when you compare what's actually in the file to what should be there.
Keep this design handy, because we'll keep coming back to it.

### Step 1: Import and prepare

The `read.csv()` function reads a comma-separated file (from your computer or from a URL)
into a **data frame**, which is R's version of a spreadsheet table. `head()` shows the
first six rows.

``` r
fundata <- read.csv("http://genomedata.org/seq-tec-workshop/misc/fundata.csv")
head(fundata) # Look at the data
```

    ##   siteplot       date temp_c
    ## 1      A.1 03-21-2017  31.50
    ## 2      B.3 04-21-2017  33.63
    ## 3      C.2 05-21-2017  32.56
    ## 4      D.4 06-21-2017  31.22
    ## 5      E.1 07-21-2017  25.44
    ## 6      F.6 08-21-2017  33.55

`str()` (short for "structure") tells you how many rows ("obs." = observations) and columns
("variables") there are, and what *type* of data is in each column. `chr` means character
(text) and `num` means numeric.

``` r
str(fundata)
```

    ## 'data.frame':    355 obs. of  3 variables:
    ##  $ siteplot: chr  "A.1" "B.3" "C.2" "D.4" ...
    ##  $ date    : chr  "03-21-2017" "04-21-2017" "05-21-2017" "06-21-2017" ...
    ##  $ temp_c  : num  31.5 33.6 32.6 31.2 25.4 ...

`summary()` gives a quick overview of each column. For numeric columns, you get the minimum,
maximum, mean, median, the 1st and 3rd quartiles (the values that 25% and 75% of the data
fall below), and a count of missing values (`NA's`). For character columns, it just tells you
how many values there are.

``` r
summary(fundata)
```

    ##    siteplot             date               temp_c        
    ##  Length:355         Length:355         Min.   :-9999.00  
    ##  Class :character   Class :character   1st Qu.:   26.75  
    ##  Mode  :character   Mode  :character   Median :   30.20  
    ##                                        Mean   : -111.25  
    ##                                        3rd Qu.:   33.36  
    ##                                        Max.   :  180.00  
    ##                                        NA's   :1


What's our assessment of this dataset from the above commands? Compare it to the
experimental design above and write down every problem you can spot before opening the hints.

<details><summary>Hints</summary>
<p>

- How many rows did we expect? How many do we have?
- The site and plot codes are smushed together in one column.
- The dates look a bit funky.
- The range on the temperatures surely can't be right.
</p>
</details>

Let's get to work cleaning up this dataset and getting it ready for analysis. 

### Step 2: Cleaning the data


#### Cleanup Task #1: Split siteplot into two columns

There's a pretty good chance that we might want to run stats by site or plot, but right now, it's complicated to subset those out.  

The `separate()` function from the `tidyr` package splits one column into several. `tidyr` is
part of the "tidyverse" collection of packages, so loading `tidyverse` gives you access to it.
(If tidyverse isn't installed on your system, run `install.packages("tidyverse")` first.)

To explain it, we'll first use some code that makes a data frame with ugly smashed-together values:

``` r
library(tidyverse)
lettersdf <- data.frame(letters = paste(LETTERS, rev(LETTERS), sep = "."))
head(lettersdf)
```

    ##   letters
    ## 1     A.Z
    ## 2     B.Y
    ## 3     C.X
    ## 4     D.W
    ## 5     E.V
    ## 6     F.U

And then here's some example code that splits them into separate columns:

```r
head(separate(lettersdf, letters, c("letter_one", "letter_two"), sep = "\\."))
```

    ##   letter_one letter_two
    ## 1          A          Z
    ## 2          B          Y
    ## 3          C          X
    ## 4          D          W
    ## 5          E          V
    ## 6          F          U

> Note: 
Why `"\\."` instead of just `"."`? The `sep` argument is a **regular expression** (a pattern
language for matching text), and in regular expressions a plain `.` means "any character".
The `\\` "escapes" the period so that it means a literal period.

**Exercise:** Using your `fundata` object, split the `siteplot` column into
two columns called `site` and `plot`. Save the result back into `fundata`, then use `head()`
to check that it worked.

Don't forget - you can use `?separate` to consult the help page for an explanation of the arguments.

<details><summary>Answer</summary>
<p>
    
``` r
fundata = separate(fundata, siteplot, c("site", "plot"), sep = "\\.")
head(fundata)
```

    ##   site plot       date temp_c
    ## 1    A    1 03-21-2017  31.50
    ## 2    B    3 04-21-2017  33.63
    ## 3    C    2 05-21-2017  32.56
    ## 4    D    4 06-21-2017  31.22
    ## 5    E    1 07-21-2017  25.44
    ## 6    F    6 08-21-2017  33.55

Note the `fundata = ...` at the beginning. Most R functions don't change your data in place.
They return a *new* data frame. If you don't save the result to a variable, it just prints
to the screen and is gone.

</p>
</details> 

#### Cleanup Task #2: What's up with those dates?

These dates use a format that isn't wrong, but isn't ideal for plotting
and such. What happens when we try to sort them? 

``` r
head(sort(fundata$date), 10)
```

    ## "03-21-2016" "03-21-2017" "03-21-2017" "03-21-2017" "03-21-2017" "03-21-2017" "03-21-2017"
    ## "03-21-2017" "03-21-2017" "03-21-2017"

That's not super useful, because there are lots of repeated values. Let's display only
non-duplicate values using the `unique()` function:

``` r
sort(unique(fundata$date))
```

    ## "03-21-2016" "03-21-2017" "04-21-2016" "04-21-2017" "05-21-2016" "05-21-2017" "06-21-2016"
    ## "06-21-2017" "07-21-2017" "08-21-2016" "08-21-2017"

Now we can see what's going on, but ugh. These are character values, so R sorts them
alphabetically, one character at a time, the way a dictionary would. That's why March 2016
and March 2017 end up next to each other, ahead of April. 

When you have dates in R, it's usually best to convert them to a `Date` object. The function we'll use is `as.Date()`. We can't just run it blindly, though, because there are lots of date formats, some of which are ambiguous: Does `03/01/2017` represent March 1st (MM/DD/YYYY) or January 3rd (DD/MM/YYYY)? 

Here's an example of using `as.Date()` it on a few dates:

``` r
datestrings <- c("2000-08-03", "2017-02-20", "1980-04-27")
class(datestrings)
```

    ## [1] "character"

``` r
mydate <- as.Date(datestrings, format = "%Y-%m-%d")
mydate
```

    ## [1] "2000-08-03" "2017-02-20" "1980-04-27"

``` r
class(mydate)
```

    ## [1] "Date" 

The `format` argument tells R how to read the text. Each `%` code stands for one part of the
date, and everything else (like the `-`) must match the text exactly:

| Code | Meaning                   | Example |
|------|---------------------------|---------|
| `%Y` | 4-digit year              | `2017`  |
| `%y` | 2-digit year              | `17`    |
| `%m` | month as a number (01-12) | `03`    |
| `%d` | day of the month (01-31)  | `21`    |
| `%b` | abbreviated month name    | `Mar`   |

So `"%Y-%m-%d"` means "year, dash, month, dash, day". If the format doesn't match your text,
`as.Date()` won't give you an error. It just quietly returns `NA`, so always check
your result!

**Exercise:** Convert the `date` column in `fundata` from a character vector to a Date
vector with `as.Date()`. Look carefully at what order the month, day, and year are in. Then
check your work with `str()` and by sorting the dates again.

<details><summary>Answer</summary>
<p>

Our dates look like `03-21-2017`, which is month-day-year, so the format is `"%m-%d-%Y"`:

``` r
fundata$date = as.Date(fundata$date, format = "%m-%d-%Y")
str(fundata)
```

    ## 'data.frame':    355 obs. of  4 variables:
    ##  $ site  : chr  "A" "B" "C" "D" ...
    ##  $ plot  : chr  "1" "3" "2" "4" ...
    ##  $ date  : Date, format: "2017-03-21" "2017-04-21" "2017-05-21" "2017-06-21" ...
    ##  $ temp_c: num  31.5 33.6 32.6 31.2 25.4 ...

The `date` column is now type `Date`, and R prints it in the standard year-month-day format.
Sorting now works in time order:

``` r
head(unique(sort(fundata$date)), 10)
```

    ##  [1] "2016-03-21" "2016-04-21" "2016-05-21" "2016-05-21" "2016-06-21" "2016-06-21"
    ##  [7] "2016-08-21" "2017-03-21" "2017-03-21" "2017-03-21"

It's also a good idea to make sure none of the dates failed to convert:
`sum(is.na(fundata$date))` should be `0`. (We'll learn more about `is.na()` below.)

</p>
</details> 


#### Cleanup Task #3: Checking the sites and plots: 

Now that our columns are cleaned up, we can check the data against the experimental design.

The `table()` function is a great way to tally the occurrences of each of
the unique values in a vector:

``` r
letters = c("A", "A", "B", "C")
table(letters) 
```

    ## A B C 
    ## 2 1 1

``` r
fish <- c("grouper", "bass", "drum", "grouper")
table(fish)
```

    ## fish
    ##    bass    drum grouper 
    ##       1       1       2

If you give `table()` *two* vectors, it makes a **cross-tabulation**: a grid that counts how
many times each combination of values occurs. Rows are values from the first vector and
columns are values from the second.

``` r
fish  <- c("grouper", "bass", "drum", "grouper")
river <- c("north",   "north", "south", "south")
table(fish, river)
```

    ##          river
    ## fish      north south
    ##   bass        1     0
    ##   drum        0     1
    ##   grouper     1     1


**Exercise**: Apply `table()` to the `site` column. Does the result match with our experimental design at the beginning?

<details><summary>Answer</summary>
<p>

``` r
table(fundata$site)
```

    ## 
    ##  A  B  C  D  E  F 
    ## 59 59 60 59 58 60 

No. Sites C and F have all 60, but A, B, and D are each missing one observation, and E is
missing two. That's the 5 missing rows (355 instead of 360) that we noticed at the start.
Depending on our needs, we may need to go back to our field notes to find out what happened.
Missing samples aren't always a problem, but you need to know about them before you start
comparing sites.

</p>
</details> 

We can also apply table in two dimensions.  What happens if you run 

``` r 
table(fundata$site, fundata$plot)
```

That's a pretty handy way to view this data!

    
#### Cleanup task #4: Temperature readings

Our temps looks funky earlier.  Start by gettin some basic stats about the `temp_c` column with `summary()`:

``` r
summary(fundata$temp_c)
```

    ##     Min.  1st Qu.   Median     Mean  3rd Qu.     Max.     NA's 
    ## -9999.00    26.75    30.20  -111.25    33.36   180.00        1

Or look directly at the range of values with `range()`:

``` r
range(fundata$temp_c)
```

    ## [1] NA NA

Yikes. That's not very useful. When a vector contains even one `NA` (missing value), many R
functions return `NA`, because R can't know what the true minimum or maximum would be. You
can tell these functions to ignore missing values with `na.rm = TRUE`:

``` r
range(fundata$temp_c, na.rm = TRUE)
```

    ## [1] -9999   180

180 degrees seems implausible, whether we're in Celcius or Farenheit!  And unless we've invented some new physics, -9999 is absolutely a mistake.

In addition to looking at text output, plotting your data is always a good idea and can often uncover issues.  For a single numeric variable, a box plot is a great way to look at the distribution of values. The thick line is the median, the box covers the
middle 50% of the data, and points drawn individually beyond the "whiskers" are potential
outliers. Try it with some randomly generated data. This function draws 100 random numbers from a
normal distribution with mean 50 and standard deviation 5:

``` r
normdata = rnorm(100, 50, 5)
```

Then we can plot it using the boxplot() function:

``` r
boxplot(normdata)
```

**Exercise:** Use the `boxplot()` function to make a boxplot of your temperature values. What do you see?

<details><summary>Answer</summary>
<p>

``` r
boxplot(fundata$temp_c)
```

The plot is mostly useless. The -9999 values stretch the y-axis so far that all of the real
data gets squashed into a flat line near zero. That's a useful thing to know, though: when
a plot looks like this, it's a strong sign that there are extreme values in your data.

</p>
</details> 

Let's take a closer look at the extreme values. `which()` returns the row numbers where a
test is TRUE, and the `|` symbol means "or":

``` r
fundata[which(fundata$temp_c < -100 | fundata$temp_c > 100), ]
```

    ##     site plot       date temp_c
    ## 116    D    5 2017-06-21  -9999
    ## 164    D    5 2017-06-21  -9999
    ## 258    C    3 2017-05-21  -9999
    ## 305    D    4 2017-06-21    180
    ## 312    E    4 2017-07-21  -9999
    ## 320    A    6 2017-03-21  -9999

There are five values of exactly -9999.  The 180 value might be a typo (e.g. `18.0` with a misplaced decimal) or a broken sensor. Either way, we can't trust it.

Assuming we can't go back to our original data and fill these in, the best solution is often to drop these rows. There is more than one way to do this, but we're going to use the dplyr `filter` function to drop those rows.

Here's an example of using filter on another column of our data set:

``` r
filter(fundata, site == "B")
```

**Exercise:** Use filter to remove the implausible temperature rows. Remember to check your command first, then overwrite your `fundata` object

<details><summary>Answer</summary>
<p>

First, check our filter command 

``` r
filter(fundata, temp_c < 100)
```

That's way too much output, so let's flip the sign and make sure it grabs just the rows we want:

``` r
filter(fundata, temp_c > 100)
```

Okay, so that row is what we'll be removing. Now let's save the result.

``` r
fundata = filter(fundata, temp_c < 100)
```

Then repeat the process to remove the low value:

``` r
fundata = filter(fundata, temp_c > -100)
```

Worth noting that we could have gotten fancy and done them both with one command, using the ampersand to mean "AND":

``` r
fundata = filter(fundata, temp_c > -100 & temp_c < 100)
```
</p>
</details> 


Let's try the boxplot and summary again:

``` r
boxplot(fundata$temp_c)
summary(fundata$temp_c)
```

    ##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    ##   17.53   26.91   30.25   29.98   33.37   41.42 

Now the temperatures range from about 17.5 to 41.4 °C, which is plausible. We removed 7 rows
(355 goes down to 348).

Notice that this also dropped our "NA" values, since they didn't match the filter criteria.  This may or may not be desired behavior, but keep it in mind!

### Step 3: Summarize and save

Finally, save your cleaned data so you don't have to repeat all of this every time you want to work with it.
`row.names = FALSE` stops R from adding an extra column of row numbers:

``` r
write.csv(fundata, "fundata_clean.csv", row.names = FALSE)
```
Where did it write that file? Check `getwd()` to see your working directory on the filesystem.  You can also provide a full path in the save command:  `/Users/cmiller/fundata_clean.csv`

This data is now ready for analyses - we could do stats on the sites, or plot temperatures over time, etc.  **If we hadn't carefully examined our data first, though, all of those analyses would be wrong!**


## LLM assisted cleanup. 

Now that you've cleaned this dataset by hand, let's see how an LLM does at the same task, and
practice reading and understanding the R code that it generates. Make a new R script file for
this exercise (e.g. `fundata_llm_qc.R`), because you'll be turning it in.

### Part 1: Ask the LLM to QC the data

1. Get the **original, uncleaned** data. Download it to your computer with:

   ``` r
   download.file("http://genomedata.org/seq-tec-workshop/misc/fundata.csv", "fundata.csv")
   ```

2. Start a **new** conversation with an LLM (so it doesn't know about this exercise), and
   upload `fundata.csv` or paste its contents. Many LLMs can't download files from URLs, so
   don't just give it the link.

3. Ask it to QC the data. Write your own prompt, thinking about what details you want to know and what context the LLM needs in order to understand what you want to do with the data.

<details><summary>Hints:</summary>
<p>

Things you might want to mention: 

- I'm still learning R
- describe the study design a little
- what kinds of errors should it be looking for?
 
</p>
</details> 

Don't just ask it to fix the data! You need a reproducible way to get from this original file to a new "cleaned" file. Ask it to generate an R script that does the cleanup you decide on after going back and forth with the LLM about potential issues. 

### Part 2: Compare and contrast

- Did the LLM find **all** of the problems we found above? Which did it miss? 
- Did it find anything we **didn't** find? If so, examine the finding yourself. Is it a real problem?
- Did it report anything that is **wrong**? For example, check every number it gives you
  (how many rows are missing, how many `-9999`s, etc.) Verify by running the R code yourself and comparing. LLMs that
  don't actually run code can easily get counts wrong while sounding very confident.
- Did it make different **decisions** than we did? Were those choices reasonable? Did it
  explain why?
- How did its approach compare to ours? Did it use the same kinds of checks (`table()`,
  `summary()`, plots), or different ones?

### Part 3: Make sure you understand the code

1. Copy the R code that the LLM generated into RStudo and **run it, one chunk at a time.** 

2. Look for functions you don't recognize. LLMs often use tidyverse packages like `dplyr`
   (`mutate()`, `group_by()`, `summarise()`, `%>%` or `|>` pipes) and `ggplot2`
   (`ggplot()`, `geom_boxplot()`, etc.). These are great tools that we may learn later, but for this exercise, ask the LLM to keep it simple. For example:

   > Please rewrite this code without using piping (%>%) and stick to simple functions that are understandable to an R beginner.

3. Make sure the final script runs without errors.  Does the output match what the LLM said it would show? If you get an error, paste the error message back to the LLM and ask it to find the problem, explain what went wrong, and suggest a fix to the code.  Saying "Just fix it" doesn't help you understand and learn.


### Part 4: Comment the code yourself

This is the most important part. **Delete all of the comments the LLM wrote** (every line
starting with `#`). Then go through the script and write your own comment above each chunk of
code, explaining:

- **what** the code does,
- **why** it's a useful QC step, and
- **what the output showed** when you ran it on this data.

Write the comments in your own words. Don't paste the LLM's explanations back in! The goal is
to prove to yourself that you understand every step. If you can't explain a chunk, that's a sign to go back and figure it out.

For example, this comment isn't very useful:

``` r
# make a table
table(fundata$site)
```

This one shows that you understand the code and what it found:

``` r
# Count how many samples there are from each site. 
# From the study design, we expect 60 in each. 
table(fundata$site)
```


### What to turn in (for-credit students)

1. Your commented R script (`fundata_llm_qc.R`)
2. A few sentences of reflection: What was the LLM good at? Did it find anything that we didn't?  Did it miss anything we found?  Where would you *not* trust it without checking? Would you have caught its mistakes (if any) if you hadn't done this exercise by hand first?

