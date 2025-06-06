# Demonstration 2
Complete the following demonstration in RStudio. All the sections below match those in Xia et al. (2018), Chapter 4: Introduction to R, RStudio and ggplot2. Review them to interpret scripts and outcomes of your analyses

## 4.2. Introduction to the dplyr Package
dplyr is a grammar of data manipulation, providing a consistent set of verbs that help you solve the most common data manipulation challenges
You can read more about dplyr at https://cran.r-project.org/web/packages/dplyr/index.html

```r
setwd("your working directoty")
tab <- read.csv("hsb2demo.csv")
head(tab)

install.packages("dplyr")
library(dplyr)

#Select columns: id, read, write and math
head(select(tab, id, read, write, math))

# Same as above using %>% (pipe operator)
tab %>% 
  select(id, read, write, math) %>% 
  head 

#Select all columns between read and socst (inclusive)
head(select(tab, read:socst))

#Select all columns except female
head(select(tab, -female))

#Select all columns except those from female to prog (inclusive)
head(select(tab, -(female:prog )))

#Select all columns that start with the character string "s"
head(select(tab, starts_with("s")))

#Filter the rows for students with reading score greater than or equal 70.
filter(tab, read >= 70)

#Filter the rows for students with both reading and math scores greater than or equal 70
filter(tab, read >= 70, math >= 70)

#Delete rows 1 to 195
slice(tab, -c(1:195))

#Re-order by id, read and write
head(arrange(tab, id, read, write))

#Use desc() to order a column in descending order
head(arrange(tab, desc(read)))

#To re-order rows by a particular column(female)
tab %>% arrange(female) %>% head

#Select three columns id, gender, read from tab
#Arrange the rows by the gender and then by read
#Then return the head of the final data frame
tab%>%select(id, female, read) %>%
  arrange(female, read) %>% 
  head

#Filter the rows for read with score greater or equal to 70
tab %>% select(id, female, read) %>%
  arrange(female, read) %>% 
  filter(read >= 70)

#Arrange the rows for read in a descending order
tab %>% select(id, female, read) %>%
  arrange(female, desc(read)) %>% 
  filter(read >= 70)

#Create new columns using mutate()
#Calculate average read and write scores
head(mutate(tab, avg_read = read/n()))

#To keep only the new variables, use transmute()
head(transmute(tab,avg_read = read/n()))

#Create new columns using mutate() and pipe operator
tab %>% mutate(avg_read = read/n()) %>%
  head

#To collapses a data frame to a single row.
summarise(tab, avg_read = mean(read, na.rm = TRUE))


#Create summaries of the data frame using summarise() and pipe operator
tab %>% summarise(avg_read = mean(read), 
                  min_read = min(read),
                  max_read = max(read),
                  n = n())

#First group by gender, and then get the summary statistics of reading by gender
by_gender <- group_by(tab, female)
read_by_gender <- summarise(by_gender,
                            n = n(),
                            avg_read = mean(read, na.rm = TRUE),
                            min_read = min(read,na.rm = TRUE),
                            max_read = max(read,na.rm = TRUE))
read_by_gender

#Create summaries of the data frame using summarise() and pipe operator
tab %>% group_by(female) %>%
  summarise(n = n(),
            avg_read  = mean(read), 
            min_read = min(read),
            max_read = max(read))


#Use sample_n() to sample a fixed random number of rows
sample_n(tab, 5)

#Use sample_frac() to sample a fixed random fraction of rows
sample_frac(tab, 0.02)

#Use replace = TRUE to perform a bootstrap sampling
sample_n(tab, 5,replace = TRUE)
```
## 4.3 Introduction to the ggplot2 Package
ggplot2 is a system for declaratively creating graphics, based on The Grammar of Graphics. You provide the data, tell ggplot2 how to map variables to aesthetics, what graphical primitives to use, and it takes care of the details. You can read more about ggplot2 at https://ggplot2.tidyverse.org

## 4.3.1	ggplot2 and the Grammar of Graphics
```r
citation("ggplot2")
install.packages("ggplot2")
library(ggplot2)
```

## 4.3.3	Creating a plot Using ggplot
### 4.3.3.1 Creating a plot Layer by Layer with ggplot
```r
data(iris)
head(iris)
library(ggplot2)

# step by step
p <- ggplot(iris, aes(x=Sepal.Width, y=Sepal.Length)) 
# Sepal.Width and Sepal.Length are columns in iris dataframe
p   
summary(p)

#Add scatterplot geom (layer1)
p1 <- p + geom_point()  
summary(p1)

ggplot(iris, aes(x=Sepal.Width, y=Sepal.Length)) + geom_point()

#Add smoothing geom (layer2)
p2 <- p1 + geom_smooth(method="lm") 
p2  
summary(p2)

#set se = FALSE to turn off confidence bands
p1 + geom_smooth(method="lm", se = FALSE) 

p3 <- ggplot(iris, aes(x=Sepal.Width, y=Sepal.Length)) + 
  #Add scatterplot geom (layer1)
  geom_point(col="blue", size=3) + 
  #Add smoothing geom (layer2)
  geom_smooth(method="lm",col="red",size=2) 
p3 

p4 <- ggplot(iris, aes(x=Sepal.Width, y=Sepal.Length)) + 
  #Add scatterplot geom (layer1)
  geom_point(aes(col=Species), size=3) + 
  #Add smoothing geom (layer2)
  geom_smooth(method="lm",col="red",size=2) 
p4

p5 <- p4 + coord_cartesian(xlim=c(2.2,4.2), ylim=c(4, 7)) # zooms in
plot(p5)

#Add Title and Labels using labs()
p6 <- p5 + labs(title="Sepal width vs sepal length", subtitle="Using iris dataset",
                y="Length of Sepal", x="Width of Sepal")
print(p6)#Or plot(p6)

#Add Title and Labels using ggtitle(), xlab() and ylab()
p7 <-p5 +  ggtitle("Sepal width vs sepal length", subtitle="Using iris dataset") +
  ylab("Length of Sepal") + xlab("Width of Sepal")
print(p7)

# Create the full scatterplot call in one step
ggplot(iris, aes(x=Sepal.Width, y=Sepal.Length)) +
  geom_point(aes(col=Species), size=3) + 
  geom_smooth(method="lm",col="red",size=2) +
  coord_cartesian(xlim=c(2.2,4.2), ylim=c(4, 7)) +
  labs(title="Sepal width vs sepal length", subtitle="Using iris dataset", 
       y="Length of Sepal", x="Width of Sepal")
```
![Alt text](image1.png)
### 4.3.3.5 Using Faceting to Detect Patterns Across Conditions
```r
#Spliting plots by rows
ggplot(iris, aes(x=Sepal.Width, y=Sepal.Length)) +
  geom_point(aes(col=Species), size=3) + 
  geom_smooth(method="lm",col="red",size=2) +
  coord_cartesian(xlim=c(2.2,4.2), ylim=c(4, 7)) +
  # Add Facet Grid
  facet_grid(Species ~.) 


#Spliting plots by columns
ggplot(iris, aes(x=Sepal.Width, y=Sepal.Length)) +
  geom_point(aes(col=Species), size=3) + 
  geom_smooth(method="lm",col="red",size=2) +
  coord_cartesian(xlim=c(2.2,4.2), ylim=c(4, 7)) +
  # Add Facet Grid
  facet_grid(.~ Species) 


#Spliting plots by columns
ggplot(iris, aes(x=Sepal.Width, y=Sepal.Length)) +
  geom_point(aes(col=Species), size=3) + 
  geom_smooth(method="lm",col="red",size=2) +
  coord_cartesian(xlim=c(2.2,4.2), ylim=c(4, 7)) +
  # Add Facet Grid
  facet_grid(.~ Species, margin=TRUE) 

#Facet Wrap
#Spliting plots by columns
ggplot(iris, aes(x=Sepal.Width, y=Sepal.Length)) +
  geom_point(aes(col=Species), size=3) + 
  geom_smooth(method="lm",col="red",size=2) +
  coord_cartesian(xlim=c(2.2,4.2), ylim=c(4, 7)) +
  #Add Facet Wrap
  facet_wrap(~ Species, nrow=2) 
```
## Creating plots with ggpubr

### Box plots
```r

install.packages("ggpubr")
library(ggpubr)

data(iris)
df <- iris
head(df, 3)

# Use outline colors for groups: Species
# Use custom color palette
# Add jitter points and use different shapes for groups

p <- ggboxplot(df, x = "Species", y = "Sepal.Length",
               color = "Species", palette =c("#00AFBB", "#E7B800", "#FC4E07"),
               add = "jitter", shape = "Species") # Adds jittered points (individual data points) to the boxplot for visualization and their shape varies by species
p

# Compare species groups using statistics (Wilcoxon test and Kruskal-Wallis test)

# Specify the pairwise group comparisons
comps <- list( c("setosa", "versicolor"), c("setosa", "virginica"), c
               ("versicolor", "virginica") )
p + stat_compare_means(comparisons = comps, method = "wilcox.test") + # Add global p-value and p-values for pairwise comparisons
  stat_compare_means(label.y = 10) # range of y-axis

p + stat_compare_means(comparisons = comps,method = "wilcox.test", label = "p.signif") + # Show the significance levels for all
  stat_compare_means(label.y = 10) # range of y-axis

```
### Histogram plots
```r
# Change outline and fill different color for groups ("Species")
# Use custom color palette

gghistogram(df, x = "Sepal.Width",
            add = "mean",
            color = "Species", fill = "Species",
            bins = 20,
            palette = c("#00AFBB", "#E7B800","#FC4E07"))

```
