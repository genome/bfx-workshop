# R visualization - ggplot2

Adapted from the [NCEAS dataviz module](https://nceas.github.io/oss-lessons/dataviz-and-interactive-tools/module-1-ggplot2.html) and [R for Data Science](https://r4ds.had.co.nz/data-visualisation.html)

``` r
library(ggplot2)
library(dplyr)
```

## Overview

ggplot2 is a popular package for visualizing data in R.

From the [home page](http://ggplot2.tidyverse.org/):

> ggplot2 is a system for declaratively creating graphics, based on The Grammar of Graphics. You provide the data, tell ggplot2 how to map variables to aesthetics, what graphical primitives to use, and it takes care of the details.

It's been around for years and has pretty good documentation and tons of example code around the web (like on [StackOverflow](https://stackoverflow.com/questions/tagged/ggplot2)). This lesson will introduce you to the basic components of working with ggplot2.



### ggplot vs base vs lattice vs XYZ...

R provides **many** ways to get your data into a plot. Three common ones are:

- "base graphics" (`plot`, `hist`, etc.)
- lattice
- ggplot2

All of them work! I use base graphics for simple, quick and dirty plots. I use ggplot2 for most everything else.

ggplot2 excels at making complicated plots easy and easy plots simple enough.

## Geoms / Aesthetics

### First steps

Let's use our first graph to answer a question: Do cars with big engines use more fuel than cars with small engines? You probably already have an answer, but try to make your answer precise. What does the relationship between engine size and fuel efficiency look like? Is it positive? Negative? Linear? Nonlinear?

### The `mpg` data frame

You can test your answer with the `mpg` **data frame** found in ggplot2 (a.k.a. `ggplot2::mpg`). A data frame is a rectangular collection of variables (in the columns) and observations (in the rows). `mpg` contains observations collected by the US Environmental Protection Agency on 38 car models.

``` r
mpg
```

```
## # A tibble: 234 × 11
##    manufacturer model      displ  year   cyl trans    drv     cty   hwy fl    class
##    <chr>        <chr>      <dbl> <int> <int> <chr>    <chr> <int> <int> <chr> <chr>
##  1 audi         a4           1.8  1999     4 auto(l5) f        18    29 p     comp…
##  2 audi         a4           1.8  1999     4 manual(… f        21    29 p     comp…
##  3 audi         a4           2    2008     4 manual(… f        20    31 p     comp…
##  4 audi         a4           2    2008     4 auto(av) f        21    30 p     comp…
##  5 audi         a4           2.8  1999     6 auto(l5) f        16    26 p     comp…
##  6 audi         a4           2.8  1999     6 manual(… f        18    26 p     comp…
##  7 audi         a4           3.1  2008     6 auto(av) f        18    27 p     comp…
##  8 audi         a4 quattro   1.8  1999     4 manual(… 4        18    26 p     comp…
##  9 audi         a4 quattro   1.8  1999     4 auto(l5) 4        16    25 p     comp…
## 10 audi         a4 quattro   2    2008     4 manual(… 4        20    28 p     comp…
## # … with 224 more rows
## # ℹ Use `print(n = ...)` to see more rows
```

Among the variables in `mpg` are:

1.  `displ`, a car's engine size, in liters.

2.  `hwy`, a car's fuel efficiency on the highway, in miles per gallon (mpg). A car with a low fuel efficiency consumes more fuel than a car with a high fuel efficiency when they travel the same distance.

To learn more about `mpg`, open its help page by running `?mpg`.

If we wanted it to show up in our "Environment" pane of RStudio, we could run `data(mpg)`. At that point, you can click on it, then scroll around to explore more of the table.

### Creating a ggplot - scatterplots

To plot `mpg`, run this code to put `displ` on the x-axis and `hwy` on the y-axis:

Every graphic you make in ggplot2 will have at least one aesthetic and at least one geom (layer). The aesthetic maps your data to your geometry (layer). Your geometry specifies the type of plot we're making (point, bar, etc.).

``` r
ggplot(data = mpg, aes(x = displ, y = hwy)) + 
  geom_point()
```

You'll frequently see people use positional arguments and just write this as a shortcut.

``` r
ggplot(mpg, aes(displ,hwy)) + 
  geom_point()
```

That's okay to do and produces the same plot, but if in doubt, explicitly setting arguments can help avoid confusion.

What makes ggplot really powerful is how quickly we can make this plot visualize more aspects of our data. Coloring each point by class (compact, van, pickup, etc.) is just a quick extra bit of code:

``` r
ggplot(mpg, aes(x = displ, y = hwy, color = class)) +
  geom_point()
```

So how did I know to write `color = class`? `aes` will pass its arguments on to any geoms you use and we can find out what aesthetic mappings geom_point takes with `?geom_point` (see section "Aesthetics")

**Challenge:** Find another aesthetic mapping `geom_point` can take and use it to visualize the data in the `drv` column.

What if we just wanted the color of the points to be blue? Maybe we'd try this:

``` r
ggplot(mpg, aes(x = displ, y = hwy, color = "blue")) +
  geom_point()
```

Well that's weird – why are the points red?

What happened here? This is the difference between **setting** and **mapping** in ggplot. The `aes` function **only** takes mappings from our data onto our geom. If we want to make all the points blue, we need to **set** it inside the geom:

``` r
ggplot(mpg, aes(x = displ, y = hwy)) +
  geom_point(color = "blue")
```

**Challenge:** Using the aesthetic you discovered and tried above, **set** another aesthetic onto our points.

Sizing each point by the number of cylinders is easy:

``` r
ggplot(mpg, aes(x = displ, y = hwy, color = class, size = cyl)) +
  geom_point()
```

So it's clear we can make scatter and bubble plots. What other kinds of plots can we make? (Hint: Tons)

#### Histograms

Let's make a histogram:

``` r
ggplot(mpg, aes(x = hwy)) + 
  geom_histogram()
```

You'll see a message (red text):

> `stat_bin()` using `bins = 30`. Pick better value with `binwidth`.

ggplot2 can calculate statistics on our data such as frequencies and, in this case, it's doing that on our `hwy` column with the `stat_bin` function. Binning data requires choosing a bin size, and the choice of bin size can completely change our histogram (often resulting in misleading conclusions). We could change the `binwidth` argument to 1 if we didn't want to hide any of our frequencies:

``` r
ggplot(mpg, aes(x = hwy)) + 
  geom_histogram(binwidth = 1)
```

**Challenge:** Find an aesthetic `geom_histogram` supports and try it out.

#### Box plots

I'm a big fan of box plots and ggplot2 can plot these too:

``` r
ggplot(mpg, aes(x = cyl, y = hwy)) + 
  geom_boxplot()
```

Oops, we got a warning, and a single box instead of one per cylinder count:

> Warning message: Continuous x aesthetic – did you forget aes(group=…)?

That's because we need to convert `cyl` to a factor:

``` r
ggplot(mpg, aes(x = factor(cyl), y = hwy)) + 
  geom_boxplot()
```
---- 
#### Aside: numbers vs. factors

ggplot decides how to draw a variable based on its *type*:

- **Numeric** variables (integers, decimals) are **continuous**. ggplot treats them as points on a number line (either 5 or 6.6667 would be fine). Continuous variables get a numeric axis, color gradients, and size ramps.
- **Factors** (and character strings) are **discrete**, or categorical. They're a fixed set of distinct labels with nothing in between them. Discrete variables get one slot on the axis per category, one distinct color per category, and so on.

`cyl` is stored as a number (<int>), so to ggplot, it's just a continuous range from 4 to 8, and there's no obvious way to split it into groups. Wrapping it in `factor()` tells ggplot "treat these as four separate categories", and we get one box per group

---- 
#### Violin plots
Another type of visualization I use a lot for seeing my distributions is the violin plot:

``` r
ggplot(mpg, aes(x = class, y = hwy)) + geom_violin()
```

How could we make this more visually appealing? How could we color the violins? One option would be to set an aesthetic on the violins themselves:

``` r
ggplot(mpg, aes(x = factor(class), y = hwy)) + geom_violin(fill="cornflowerblue")
```

We could also give each class its own color by mapping class to the fill aesthetic:

``` r
ggplot(mpg, aes(x = factor(class), y = hwy, fill = factor(class))) + geom_violin()
```

So far we've made really simple plots: One geometry per plot. Let's layer multiple geometries on top of one another to show the raw points on top of the violins:

``` r
ggplot(mpg, aes(x = factor(class), y = hwy)) + 
  geom_violin() +
  geom_point(shape = 1)
```

Our points got added, but it's impossible to see the "pileup" of points because they're all plotted on top of each other. Adding some "jitter" can help with that.

``` r
ggplot(mpg, aes(x = factor(class), y = hwy)) + 
  geom_violin() +
  geom_point(shape = 1, position = "jitter")
```

#### Fitting lines

Some geoms can do even more than just show us our data. ggplot2 also helps us do some quick-and-dirty modeling:

``` r
ggplot(mpg, aes(x = cty, y = hwy)) + 
  geom_point() +
  geom_smooth()
```

Notice the message in red text:

> `geom_smooth()` using method = 'loess'

geom_smooth defaulted here to using a LOESS smoother. But `geom_smooth()` is pretty configurable. Here we set the method to `lm` instead of the default `loess`:

``` r
ggplot(mpg, aes(cty, hwy)) + 
  geom_point() +
  geom_smooth(method = "lm")
```

More on geoms built in to ggplot here: <http://ggplot2.tidyverse.org/reference/index.html#section-layer-geoms> and more that you can add by installing packages here: <https://exts.ggplot2.tidyverse.org/gallery/>


### Plot customization

ggplot2 is highly customizabile, and the `theme()` function and pre-set themes make it easy to change  how a plot looks. Let's use a different theme than the default:

``` r
ggplot(mpg, aes(displ, hwy, color = class)) + 
  geom_point() +
  theme_classic()
```

**Challenge:** Find another theme and use it instead.

#### Color schemes

The default ggplot colors are fine, but they aren't always easy to tell apart, and they don't hold up well for colorblind readers or in black-and-white printing. Rather than picking colors by hand, you can use an existing, well-designed color scheme. One of the most popular is **viridis**, which is built into ggplot2:

``` r
ggplot(mpg, aes(displ, hwy, color = class)) + 
  geom_point() +
  theme_classic() +
  scale_color_viridis_d()
```

The `_d` stands for **discrete**, because `class` is a categorical variable. For a continuous (numeric) variable, use `scale_color_viridis_c()` instead (this is the same numbers vs. factors distinction from above!). The `option` argument picks from several related palettes, such as "magma", "plasma", or "cividis":

``` r
ggplot(mpg, aes(displ, hwy, color = cty)) + 
  geom_point() +
  theme_classic() +
  scale_color_viridis_c(option = "plasma")
```

Another commonly used set of palettes comes from [ColorBrewer](https://colorbrewer2.org/), available with `scale_color_brewer()` (discrete) and `scale_color_distiller()` (continuous). Let's try one:

``` r
ggplot(mpg, aes(displ, hwy, color = class)) + 
  geom_point() +
  theme_classic() +
  scale_color_brewer(palette = "Set2")
```

> [!NOTE]
> If you're coloring with `fill` instead of `color` (like our violins), the same functions exist as `scale_fill_viridis_d()`, `scale_fill_brewer()`, etc.

For a much more detailed look at color scales and how to choose them, see the [Colour scales and legends](https://ggplot2-book.org/scales-colour.html) chapter of the *ggplot2* book.

#### Themes and labels - getting fancy

The legend in ggplot2 is a thematic element. Let's change a bunch of things about the way the legend displays:

``` r
ggplot(mpg, aes(displ, hwy, color = class)) + 
  geom_point() +
  theme_classic() +
  theme(legend.position = "bottom",
        legend.background = element_rect(fill = "#EEEEEE", color = "black"),
        legend.title = element_blank(),
        axis.title = element_text(size = 16))
```

And then, let's make another plot where we adjust our axis labels and title:

``` r
ggplot(mpg, aes(displ, hwy, color = class)) + 
  geom_point() +
  theme_classic() +
  theme(legend.position = "inside",
        legend.position.inside = c(1, 1),
        legend.justification = c(1,1),
        legend.direction = "horizontal",
        legend.title = element_blank()) +
  xlab("Engine Displacement") +
  ylab("Highway Fuel Economy (miles / gallon)") +
  ggtitle("Highway fuel economy versus engine displacement",
          "or why do you need that big truck again? ")
```

**Challenge:** Look at the help for `?theme` and try changing something else about the above plot.

More themes are available in a user-contributed package called [ggthemes](https://github.com/jrnold/ggthemes).

### Saving plots

Let's save that great plot we just made. Saving plots in ggplot is done with the `ggsave()` function:

``` r
ggsave("hwy_vs_displ.png")
```

`ggsave` automatically:

1. Saves the last plot you made (unless you specify otherwise)
2. Chooses the format based on your file extension 
3. Guesses a default image size. 

We can customize the size with the `width` and `height` arguments:

``` r
ggsave("hwy_vs_displ.png", width = 6, height = 6)
```


### Setting plot limits

Plot limits can be controlled in one of two ways:

- Filter the data (because limits are auto-calculated from the data ranges)
- Set them explicitly with `xlim()` and `ylim()`

Let's show this with an example plot:

``` r
ggplot(economics, aes(x = date, y = unemploy)) + 
  geom_line()
```

Since we're plotting data where the zero point on the vertical axis means something, maybe we want to start the vertical axis at 0:

``` r
ggplot(economics, aes(x = date, y = unemploy)) + 
  geom_line() +
  ylim(0, 20000)
```

Rather than hard-coding the upper limit, we can also calculate it from the data:

``` r
ggplot(economics, aes(x = date, y = unemploy)) + 
  geom_line() +
  ylim(0, max(economics$unemploy))
```

Or maybe we want to zoom in on just the 2000s up to 2020. Since our x-axis holds dates, the limits need to be dates too:

``` r
ggplot(economics, aes(x = date, y = unemploy)) + 
  geom_line() +
  ylim(0, max(economics$unemploy)) +
  xlim(as.Date("2000-01-01"), as.Date("2019-12-31"))
```

Note the warning message we received:

> Warning message: Removed 390 rows containing missing values or values outside the scale range (`geom_line()`).

That's normal when data in your input data.frame are outside the range we're plotting.


### Scales

The usual use case is to do things like changing scale limits or change the **way** our data are mapped onto our geom. We'll use scales in ggplot2 very often!

For example, how do we override the default colors ggplot2 uses here?

``` r
ggplot(mpg, aes(x = displ, y = hwy, color = class)) +
  geom_point()
```

**Tip:** Most scales follow the format `scale_{aesthetic}_{method}` where `aesthetic` is our aesthetic mappings such as color, fill, shape and `method` is how the colors, fill colors, and shapes are chosen. (The `scale_color_viridis_d()` we used above is one example.)

``` r
ggplot(mpg, aes(displ, hwy, color = class)) +
  geom_point() + 
  scale_color_manual(values = c("red", "orange", "yellow", "green", "blue", "purple", "violet")) # ROYGBIV
```

I'm sure that was a ton of fun to type out but we can make things easier on ourselves:

``` r
ggplot(mpg, aes(displ, hwy, color = class)) +
  geom_point() + 
  scale_color_hue(h = c(270, 360)) # blue to red
```

Above we were using scales to scale the `color` aesthetic. We can also use scales to rescale our data. Here's some census data, unscaled:

``` r
ggplot(midwest, aes(x = area, y = poptotal)) + 
  geom_point()
```

And scaled (log10):

``` r
ggplot(midwest, aes(area, poptotal)) + 
  geom_point() + 
  scale_y_log10()
```

We could have also log-scaled the y-axis like this:

``` r
ggplot(midwest, aes(x = area, y = log10(poptotal))) + 
    geom_point()
```

What's the difference? (look closely!)

Scales can also be used to change our axes. For example, we can override the labels:

``` r
hwybyclass = summarize(group_by(mpg,class), maxhwy = max(hwy))
ggplot(hwybyclass, aes(x = class, y = maxhwy)) +
    geom_col() +
    scale_x_discrete(labels = paste("group",1:7))
```

Or change the breaks:

``` r
ggplot(mpg, aes(displ, hwy)) + 
  geom_point() +
  scale_y_continuous(breaks = seq(min(mpg$hwy), max(mpg$hwy)))
```

### Facets

Facets allow us to create a powerful visualization called a small multiple. 

I use small multiples all the time when I have a variable like a site or year and I want to quickly compare across years. Let's compare highway fuel economy versus engine displacement across our two samples:

``` r
ggplot(mpg, aes(x = displ, y = hwy)) +
  geom_point() +
  facet_wrap(~ year)
```

Or fuel economy versus engine displacement across manufacturer:

``` r
ggplot(mpg, aes(x = displ, y = hwy)) +
  geom_point() +
  facet_wrap(~ manufacturer)
```

## Resources

- Multiple graphs in a single graph:
  - <http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/>
  - cowplot: <https://cran.r-project.org/web/packages/cowplot/vignettes/introduction.html>
  - <http://www.cookbook-r.com/Graphs/>
- Book on ggplot: <https://ggplot2-book.org/>
