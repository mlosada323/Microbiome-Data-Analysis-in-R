# Demonstration 2
Complete the following demonstration in RStudio. Create a markdown file of your script. You can follow detail instructions in Xia et al. (2018), Chapter 4: Introduction to R, RStudio and ggplot2. All the sections below match the sections in the book

## 4.2. Introduction to the dplyr Package
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

#Re-order by read and write
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
head(mutate(tab, avg_read = sum(read)/n()))

#To keep only the new variables, use transmute()
head(transmute(tab,avg_read = sum(read)/n()))

#Create new columns using mutate() and pipe operator
tab %>% mutate(avg_read = sum(read/n())) %>%
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


#Use sample_n() to sample a fixed number
sample_n(tab, 5)

#Use sample_frac() to sample a fixed fraction
sample_frac(tab, 0.02)

#Use replace = TRUE to perform a bootstrap sampling
sample_n(tab, 5,replace = TRUE)
```
