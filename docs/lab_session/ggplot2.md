# ggplot2


ggplot2 is an R package for making visualizations and is based on the ‘Grammar of Graphics’. It is often referred as one of the best packages for visualizations. 

Actually, there are other tools for making graphs in R such as Lattice and base R graphics. However, ggplot2 is straightforwardly easy to make any kind of graphs even complex ones. 

It is easy to add complexity or take it away which makes ggplot2 superior among these three. Since you can’t go back and forth once you created a plot in base R or in lattice. That’s enough for comparison , let’s start learning it now.

> Tip: In order to take most out of this tutorial you should not miss reading any lines in this tutorial and follow the flow and write codes on your computer.


### Grammar of Graphics

Each language has a grammar and this grammar dictates the structure of sentences and the position of each element in a given sentence. For example in English grammar, we could take the sentence **"The little monkey hangs confidently by a branch"**

![](https://i.imgur.com/WamzQaZ.png)

Notice the need for these elements (words) in specific positions to make a meaningful sentence. This is governed by the grammar. Just as the grammar of language helps us construct meaningful sentences out of words, the Grammar of Graphics helps us to construct graphical figures out of different visual elements. This grammar gives us a way to talk about parts of a plot: all the circles, lines, arrows, and words that are combined into a diagram for visualizing data. 

Originally developed by Leland Wilkinson, the [Grammar of Graphics](https://www.worldcat.org/title/934513040){:target="_blank"} was adapted by Hadley Wickham to describe the components of a plot, including

* the data being plotted
* the geometric objects (circles, lines, etc.) that appear on the plot
* a set of mappings from variables in the data to the aesthetics (appearance) of the geometric objects
* a statistical transformation used to calculate the data values used in the plot
* a position adjustment for locating each geometric object on the plot
* a scale (e.g., range of values) for each aesthetic mapping used
* a coordinate system used to organize the geometric objects
* the facets or groups of data shown in different plots

Similar to how tidyr and dplyr provide efficient data transformation and manipulation, ggplot2 provides more efficient ways to create specific visual images.

### Working with ggplot2 is like playing with LEGO bricks and Tetris!
Imagine you have a big pile of LEGO bricks and you want to build a castle. You start by selecting the bricks you want to use and arranging them in a certain way to create the foundation of the castle. This is similar to how you start building a plot in ggplot2 by selecting your data and specifying the basic plot components.

Next, you begin to add **layers** to your LEGO castle, such as towers, walls, and a drawbridge. Each layer adds a different aspect to the castle and helps to make it more complete and visually interesting. In ggplot2, you can add layers to your plot to display different **elements** of your data, such as points, lines, and labels. Each layer can be customized separately to create a unique and informative plot.

Now imagine you are playing a game of Tetris. Each piece that falls from the top of the screen is like a layer in ggplot2. You have to rotate and position each piece carefully to create a complete row and prevent gaps. Similarly, in ggplot2, you have to add and position each layer carefully to create a complete and informative plot.

> ***Just like building a LEGO castle or playing Tetris, creating a plot in ggplot2 requires careful selection and arrangement of the different components to create a complete and visually appealing result.*** 

So with that inspiration, we can summarise that ggplot objects can be highly complex. The basic order of layers will usually look like this:

1. **Begin with the baseline ggplot() command** - this “opens” the ggplot and allow subsequent functions to be added with +. Typically the dataset is also specified in this command
1. **Add “geom” layers** - these functions visualize the data as geometries (shapes), e.g. as a bar graph, line plot, scatter plot, histogram (or a combination!). These functions all start with geom_ as a prefix.
1. **Add design elements to the plot** such as axis labels, title, fonts, sizes, color schemes, legends, or axes rotation.

Read about all **the building blocks available to you** from the ggplot2 Function reference documentation -- [ggplot2 function reference](https://ggplot2.tidyverse.org/reference/index.html). 

A simple example of skeleton code is as follows. We will explain each component in the sections below.

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
![](https://i.imgur.com/g5rCVvp.png){:target="_blank"}


## ggplot2 hands-on exercise

> **WE WILL DO THIS ON YOUR LAPTOP!**

For our hands-on exercise, we will use a simple, easy to understand step-by-step tutorial from Library Carpentry project. (Library Carpentry is a global community teaching software and data skills to people working in library- and information-related roles. Read more about the project [here](https://librarycarpentry.org/){:target="_blank"}).

In the exercise, we will gradually build up on plotting data using ggplot2 in R.

Before We Start lesson set up your directories and data

1. Under the File menu, click on New project, choose New directory, then New project
1. Enter the name **ki_cg_course** for this new folder (or “directory”). This will be your working directory for the rest of the day.
1. Click on Create project
1. Create a new file where we will type our scripts. Go to File > New File > R script. Click the save icon on your toolbar and save your script as “script.R”.

Copy and paste the below lines of code to create three new subdirectories and download the original and the reformatted books data:

```R
install.packages("fs")
library(fs)   
# https://fs.r-lib.org/.  
# fs is a cross-platform, uniform interface to file system operations via R. 

# let's create data directories
dir_create("data")
dir_create("data_output")
dir_create("fig_output")

# now lets download the files we will be using in our exercise

download.file("https://ndownloader.figshare.com/files/22031487",
              "data/books.csv", mode = "wb")

download.file("https://ndownloader.figshare.com/files/22051506",
              "data_output/books_reformatted.csv", mode = "wb")
```

### CLICK HERE  ---> [Data Visualisation with ggplot2](https://librarycarpentry.org/lc-r/04-data-viz-ggplot/index.html#load-the-tidyverse-and-data-frame-into-your-r-session){:target="_blank"}

(clicking on the link above will take you to the **Data Visualisation with ggplot2** chapter from **Introduction to R** series in Library Carpentry's website, please feel free to bookmark the page for later reference.)

--------
#### Extra tips

If you have completed the above tutorial and would like to explore more on using ggplot in epidemiological context, please feel free to use the EpiRHandbook resources. 


To try ggplot hands-on, you could follow the original ggplot basics chapter from [The Epidemiologist R Handbook](https://www.epirhandbook.com/en/index.html){:target="_blank"} (Batra, Neale, et al. The Epidemiologist R Handbook. 2021).