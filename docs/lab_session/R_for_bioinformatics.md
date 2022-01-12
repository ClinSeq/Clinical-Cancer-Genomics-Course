# Programming in R for Bioinformatics
> Tip: In order to take most out of this tutorial you should not miss reading any lines in this tutorial and follow the flow and write codes on your computer.

### What is R?
R is a language and environment for statistical computing and graphics developed in 1993. It provides a wide variety of statistical and graphical techniques (linear and nonlinear modeling, statistical tests, time series analysis, classification, clustering, …), and is highly extensible, meaning that the user community can write new R tools. It is a GNU project (Free and Open Source).

The R language has its roots in the S language and environment which was developed at Bell Laboratories (formerly AT&T, now Lucent Technologies) by John Chambers and colleagues. R was created by Ross Ihaka and Robert Gentleman at the University of Auckland, New Zealand, and now, R is developed by the R Development Core Team, of which Chambers is a member. R is named partly after the first names of the first two R authors (Robert Gentleman and Ross Ihaka), and partly as a play on the name of S. R can be considered as a different implementation of S. There are some important differences, but much code written for S runs unaltered under R.

### Some of R’s strengths:

1. The ease with which well-designed publication-quality plots can be produced, including mathematical symbols and formulae where needed. Great care has been taken over the defaults for the minor design choices in graphics, but the user retains full control.
1. It compiles and runs on a wide variety of UNIX platforms and similar systems (including FreeBSD and Linux), Windows and MacOS.
1. R can be extended (easily) via packages.
1. R has its own LaTeX-like documentation format, which is used to supply comprehensive documentation, both on-line in a number of formats and in hardcopy.
1. It has a vast community both in academia and in business.
1. It’s FREE!

### The R environment
R is an integrated suite of software facilities for data manipulation, calculation and graphical display. It includes

* an effective data handling and storage facility,
* a suite of operators for calculations on arrays, in particular matrices,
* a large, coherent, integrated collection of intermediate tools for data analysis,
* graphical facilities for data analysis and display either on-screen or on hardcopy, and
* a well-developed, and effective programming language which includes conditionals, loops, user-defined recursive functions and input and output facilities.
* The term “environment” is intended to characterize it as a fully planned and coherent system, rather than an incremental accretion of very specific and inflexible tools, as is frequently the case with other data analysis software.

R, like S, is designed around a true computer language, and it allows users to add additional functionality by defining new functions. Much of the system is itself written in the R dialect of S, which makes it easy for users to follow the algorithmic choices made. For computationally-intensive tasks, C, C++ and Fortran code can be linked and called at run time. Advanced users can write C code to manipulate R objects directly.

Many users think of R as a statistics system. The R group prefers to think of it of an environment within which statistical techniques are implemented.

### The R Homepage
The R homepage has a wealth of information on it,

R-project.org

On the homepage you can:

Learn more about R
* Download R
* Get Documentation (official and user supplied)
* Get access to CRAN ‘Comprehensive R archival network’
* RStudio
* RStudio started in 2010, to offer R a more full featured integrated development environment (IDE) and modeled after matlabs IDE.

### RStudio has many features:

* syntax highlighting
* code completion
* smart indentation
* “Projects”
* workspace browser and data viewer
* embedded plots
* Markdown notebooks, Sweave authoring and knitr with one click pdf or html
* runs on all platforms and over the web
* etc. etc. etc.


RStudio and its team have contributed to many R packages. These include:

**Tidyverse** – R packages for data science, including ggplot2, dplyr, tidyr, and purrr
**Shiny** – An interactive web technology
**RMarkdown** – Insert R code into markdown documents
**knitr** – Dynamic reports combining R, TeX, Markdown & HTML
**packrat** – Package dependency tool
**devtools** – Package development tool

#### 1. Getting started

Let’s start RStudio

![](https://i.imgur.com/Rg8hId9.jpg)


#### 2. Open a new RScript File

File -> New File -> R Script

RStudio_newfile

Then save the new empty file as Intro2R.R

File -> Save as -> Intro2R.R

#### 3. Basics of your environment

The R prompt is the ‘>’ , when R is expecting more (command is not complete) you see a ‘+’

![](https://i.imgur.com/hNyKj3h.png)


#### 4. Writing and running R commands

In the source editor (top left by default) type

`getwd()`
Then on the line Control + Enter (Linux/Windows), Command + Enter (Mac) to execute the line.

#### 5. The assignment operator ( <- ) vs equals ( = )

The assignment operator is used assign data to a variable

```
x <- 1:10
x
[1]  1  2  3  4  5  6  7  8  9 10
```
In this case, the equal sign works as well

```
x = 1:10
x
[1]  1  2  3  4  5  6  7  8  9 10
```
But you should NEVER EVER DO THIS

```
1:10 -> x
x
[1]  1  2  3  4  5  6  7  8  9 10
```
The two act the same in most cases. The difference in assignment operators is clearer when you use them to set an argument value in a function call. For example:

```
median(x = 1:10)
x 
Error: object 'x' not found
```
In this case, x is declared within the scope of the function, so it does not exist in the user workspace.

```
median(x <- 1:10)
x
[1]  1  2  3  4  5  6  7  8  9 10
```
In this case, x is declared in the user workspace, so you can use it after the function call has been completed. There is a general preference among the R community for using <- for assignment (other than in function signatures)

#### 6. The RStudio Cheat Sheets

[rstudio-ide.pdf](https://github.com/rstudio/cheatsheets/raw/master/rstudio-ide.pdf)

spend 15m getting to know RStudio a little


Import and export data in R
====================================================

R base function read.table() is a general funciton that can be used to read a file in table format. The data will be imported as a data frame.


```r
# To read a local file. If you have downloaded the raw_counts.txt file to your local machine, you may use the following command to read it in, by providing the full path for the file location. The way to specify the full path is the same as taught in the command line session. Here we assume raw_counts.txt is in our current working directory
data <- read.table(file="./raw_counts.txt", sep="\t", header=T, stringsAsFactors=F)

# There is a very convenient way to read files from the internet.
data <- read.table(file="https://raw_counts.txt", sep="\t", header=T, stringsAsFactors=F)
```

Take a look at the beginning part of the data frame.

```r
head(data)
```

```
##            C61  C62  C63  C64  C91  C92  C93 C94 I561 I562 I563 I564 I591 I592
## AT1G01010  322  346  256  396  372  506  361 342  638  488  440  479  770  430
## AT1G01020  149   87  162  144  189  169  147 108  163  141  119  147  182  156
## AT1G01030   15   32   35   22   24   33   21  35   18    8   54   35   23    8
## AT1G01040  687  469  568  651  885  978  794 862  799  769  725  715  811  567
## AT1G01046    1    1    5    4    5    3    0   2    4    3    1    0    2    8
## AT1G01050 1447 1032 1083 1204 1413 1484 1138 938 1247 1516  984 1044 1374 1355
##           I593 I594 I861 I862 I863 I864 I891 I892 I893 I894
## AT1G01010  656  467  143  453  429  206  567  458  520  474
## AT1G01020  153  177   43  144  114   50  161  195  157  144
## AT1G01030   16   24   42   17   22   39   26   28   39   30
## AT1G01040  831  694  345  575  605  404  735  651  725  591
## AT1G01046    8    1    0    4    0    3    5    7    0    5
## AT1G01050 1437 1577  412 1338 1051  621 1434 1552 1248 1186
```


Depending on the format of the file, several variants of read.table() are available to make reading a file easier.

read.csv(): for reading "comma separated value" files (.csv).

read.csv2(): variant used in countries that use a comma "," as decimal point and a semicolon ";" as field separators.

read.delim(): for reading "tab separated value" files (".txt"). By default, point(".") is used as decimal point.

read.delim2(): for reading "tab separated value" files (".txt"). By default, comma (",") is used as decimal point.

Choosing the correct function (or parameters) is important! If we use `read.csv()` to read our tab-delimited file, it becomes a mess.

```r
data2 <- read.csv(file="https://raw_counts.txt", stringsAsFactors=F)

head(data2)
```

```
##                                         C61.C62.C63.C64.C91.C92.C93.C94.I561.I562.I563.I564.I591.I592.I593.I594.I861.I862.I863.I864.I891.I892.I893.I894
## 1                     AT1G01010\t322\t346\t256\t396\t372\t506\t361\t342\t638\t488\t440\t479\t770\t430\t656\t467\t143\t453\t429\t206\t567\t458\t520\t474
## 2                        AT1G01020\t149\t87\t162\t144\t189\t169\t147\t108\t163\t141\t119\t147\t182\t156\t153\t177\t43\t144\t114\t50\t161\t195\t157\t144
## 3                                               AT1G01030\t15\t32\t35\t22\t24\t33\t21\t35\t18\t8\t54\t35\t23\t8\t16\t24\t42\t17\t22\t39\t26\t28\t39\t30
## 4                     AT1G01040\t687\t469\t568\t651\t885\t978\t794\t862\t799\t769\t725\t715\t811\t567\t831\t694\t345\t575\t605\t404\t735\t651\t725\t591
## 5                                                                     AT1G01046\t1\t1\t5\t4\t5\t3\t0\t2\t4\t3\t1\t0\t2\t8\t8\t1\t0\t4\t0\t3\t5\t7\t0\t5
## 6 AT1G01050\t1447\t1032\t1083\t1204\t1413\t1484\t1138\t938\t1247\t1516\t984\t1044\t1374\t1355\t1437\t1577\t412\t1338\t1051\t621\t1434\t1552\t1248\t1186
```

However, the `read.csv()` function is appropriate for a comma-delimited file.

```r
data3 <- read.csv(file="https://raw_counts.csv", stringsAsFactors=F)

head(data3)
```

```
##            C61  C62  C63  C64  C91  C92  C93 C94 I561 I562 I563 I564 I591 I592
## AT1G01010  322  346  256  396  372  506  361 342  638  488  440  479  770  430
## AT1G01020  149   87  162  144  189  169  147 108  163  141  119  147  182  156
## AT1G01030   15   32   35   22   24   33   21  35   18    8   54   35   23    8
## AT1G01040  687  469  568  651  885  978  794 862  799  769  725  715  811  567
## AT1G01046    1    1    5    4    5    3    0   2    4    3    1    0    2    8
## AT1G01050 1447 1032 1083 1204 1413 1484 1138 938 1247 1516  984 1044 1374 1355
##           I593 I594 I861 I862 I863 I864 I891 I892 I893 I894
## AT1G01010  656  467  143  453  429  206  567  458  520  474
## AT1G01020  153  177   43  144  114   50  161  195  157  144
## AT1G01030   16   24   42   17   22   39   26   28   39   30
## AT1G01040  831  694  345  575  605  404  735  651  725  591
## AT1G01046    8    1    0    4    0    3    5    7    0    5
## AT1G01050 1437 1577  412 1338 1051  621 1434 1552 1248 1186
```

Since the data contained in these files is the same, we don't need to keep three copies.

```r
identical(data, data3)
```

```
## [1] TRUE
```

```r
rm(data2, data3)
```


R base function write.table() can be used to export data to a file.


```r
# To write to a file called "output.txt" in your current working directory.
write.table(data[1:20,], file="output.txt", sep="\t", quote=F, row.names=T, col.names=T)
```

It is also possible to export data to a csv file.

write.csv()

write.csv2()

---

Basic statistics in R
====================================================

| Description                            |  R_function  |
|----------------------------------------|--------------|
|  Mean                                  |  mean()      |
|  Standard deviation                    |  sd()        |
|  Variance                              |  var()       |
|  Minimum                               |  min()       |
|  Maximum                               |  max()       |
|  Median                                |  median()    |
|  Range of values: minimum and maximum  |  range()     |
|  Sample quantiles                      |  quantile()  |
|  Generic function                      |  summary()   |
|  Interquartile range                   |  IQR()       |

Calculate the mean expression for each sample.


```r
apply(data, 2, mean)
```

```
##         V1         V2         V3         V4         V5         V6         V7 
##  0.1800418 -0.2534468 -0.1897400 -0.5426865 -0.1894511 -0.2531511 -0.3651648
```

Calculate the range of expression for each sample.


```r
apply(data, 2, range)
```

```
##              V1        V2         V3         V4        V5         V6        V7
## [1,] -0.6225434 -2.452988 -1.1885474 -1.5352739 -1.311441 -1.3231301 -3.108072
## [2,]  0.9157785  1.558738  0.5844574  0.5020093  1.986243  0.8357217  1.515365
```

Calculate the quantiles of each samples.


```r
apply(data, 2, quantile)
```

```
##               V1         V2         V3          V4          V5         V6
## 0%   -0.62254340 -2.4529883 -1.1885474 -1.53527385 -1.31144098 -1.3231301
## 25%  -0.24319891 -1.4878457 -0.4367767 -1.03930207 -0.83116242 -0.9900499
## 50%   0.08886876  0.2487017 -0.3765026 -0.58285066 -0.50073601 -0.4826310
## 75%   0.68229320  0.9235564  0.2629829 -0.05204306  0.08105049  0.5890407
## 100%  0.91577852  1.5587378  0.5844574  0.50200930  1.98624311  0.8357217
##               V7
## 0%   -3.10807211
## 25%  -0.84236920
## 50%  -0.05402889
## 75%   0.38766034
## 100%  1.51536487
```

---

Simple data visualization
====================================================

Scatter plot and line plot can be produced using the function plot().


```r
x <- c(1:50)
y <- 1 + sqrt(x)/2
plot(x,y)
```

![Plot1](https://course-5534.s3.amazonaws.com/Course_img/Rplot_01.png)

```r
plot(x,y, type="l")
```

![Plot2](https://course-5534.s3.amazonaws.com/Course_img/Rplot_02.png)

```r
# plot both the points and lines
## first plot points
plot(x,y)
lines(x,y, type="l")
```

![Plot2](https://course-5534.s3.amazonaws.com/Course_img/Rplot_03.png)

```r
## lines() can only be used to add information to a graph, while it cannot produce a graph on its own.
```


boxplot() can be used to summarize data.


```r
boxplot(data, xlab="Sample ID", ylab="Raw Counts")
```

![Plot4](plot_uploaded--HERE)<!-- -->

add more details to the plot.


```r
boxplot(data, xlab="Sample ID", ylab="Raw Counts", main="Expression levels", col="blue", border="black")
```

![Plot5](plot_uploaded--HERE)<!-- -->



```r
x <- rnorm(1000)
boxplot(x)
```

![Plot6](https://course-5534.s3.amazonaws.com/Course_img/Rplot_06.png)

hist() can be used to create histograms of data.

```r
hist(x)
```

![Plot7](https://course-5534.s3.amazonaws.com/Course_img/Rplot_07.png)

```r
# use user defined break points
hist(x, breaks=seq(range(x)[1]-1, range(x)[2]+1, by=0.5))
```

![Plot8](https://course-5534.s3.amazonaws.com/Course_img/Rplot_08.png)


```r
# clear plotting device/area
dev.off()
```

```
## null device 
##           1
```

---

ggplot2
====================================================

ggplot2 is an R package for making visualizations and is based on the ‘Grammar of Graphics’ . It is often referred as one of the best packages for visualizations . Actually , there are other tools for making graphs in R such as Lattice and base R graphics . However , ggplot2 is straightforwardly easy to make any kind of graphs even complex ones . It is easy to add complexity or take it away which makes ggplot2 superior among these three . Since you can’t go back and forth once you created a plot in base R or in lattice . That’s enough for comparison , let’s start learning it now.

> Tip: In order to take most out of this tutorial you should not miss reading any lines in this tutorial and follow the flow and write codes on your computer.

***Grammar of Graphics***

Just as the grammar of language helps us construct meaningful sentences out of words, the Grammar of Graphics helps us to construct graphical figures out of different visual elements. This grammar gives us a way to talk about parts of a plot: all the circles, lines, arrows, and words that are combined into a diagram for visualizing data. Originally developed by Leland Wilkinson, the Grammar of Graphics was adapted by Hadley Wickham to describe the components of a plot, including

* the data being plotted
* the geometric objects (circles, lines, etc.) that appear on the plot
* a set of mappings from variables in the data to the aesthetics (appearance) of the geometric objects
* a statistical transformation used to calculate the data values used in the plot
* a position adjustment for locating each geometric object on the plot
* a scale (e.g., range of values) for each aesthetic mapping used
* a coordinate system used to organize the geometric objects
* the facets or groups of data shown in different plots

Similar to how tidyr and dplyr provide efficient data transformation and manipulation, ggplot2 provides more efficient ways to create specific visual images. 

ggplot objects can be highly complex, but the basic order of layers will usually look like this:

* Begin with the baseline ggplot() command - this “opens” the ggplot and allow subsequent functions to be added with +. Typically the dataset is also specified in this command
* Add “geom” layers - these functions visualize the data as geometries (shapes), e.g. as a bar graph, line plot, scatter plot, histogram (or a combination!). These functions all start with geom_ as a prefix.
* Add design elements to the plot such as axis labels, title, fonts, sizes, color schemes, legends, or axes rotation
* A simple example of skeleton code is as follows. We will explain each component in the sections below.

```r
# plot data from my_data columns as red points
ggplot(data = my_data)+                   # use the dataset "my_data"
  geom_point(                             # add a layer of points (dots)
    mapping = aes(x = col1, y = col2),    # "map" data column to axes
    color = "red")+                       # other specification for the geom
  labs()+                                 # here you add titles, axes labels, etc.
  theme()                                 # here you adjust color, font, size etc of non-data plot elements (axes, title, etc.) 
```

## Starting up with ggplot

Firstly , we need to install ggplot2 package
```r
install.packages("ggplot2")
library("ggplot2")
```

### The Basics of plotting with ggplot

In order to create a plot, you:

1. Call the ggplot() function which creates a blank canvas
1. Specify aesthetic mappings, which specifies how you want to map variables to visual aspects. In this case we are simply mapping the displ and hwy variables to the x- and y-axes.
1. You then add new layers that are geometric objects which will show up on the plot. In this case we add geom_point to add a layer with points (dot) elements as the geometric shapes to represent the data.

```r
# create canvas
ggplot(data)

# variables of interest mapped
ggplot(data, aes(x = displ, y = hwy))

# data plotted
ggplot(data, aes(x = displ, y = hwy)) +
  geom_point()
```
![](https://i.imgur.com/g5rCVvp.png)


To try ggplot hands-on, we will follow the original ggplot basics chapter from [The Epidemiologist R Handbook](https://www.epirhandbook.com/en/index.html) (Batra, Neale, et al. The Epidemiologist R Handbook. 2021).

## Now let's go to [<b>ggplot hands-on exercise.</b>](https://www.epirhandbook.com/en/ggplot-basics.html)
(clicking on the link above will take you to the handbook's website, please feel free to bookmark the page for later reference.)
 
