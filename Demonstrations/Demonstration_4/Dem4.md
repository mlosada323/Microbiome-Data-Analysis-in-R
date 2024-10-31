# Demonstration 5

Complete the following demonstration in RStudio. Create a markdown file of your script. You can follow detail instructions in Xia et al. (2018), Chapter 5: Power and Sample Size Calculations for Microbiome Data. All the sections below match the sections in the book

# Power Analysis for Microbiome Data

## 5.2.2 Diversity Data for ALS Study

# upload the ALS data set including Shannon diversity estimates for 16 samples belonging to two groups
df_H_G93BUm3 <- read.csv("ALS.csv",row.names=1,check.names=FALSE) 
df_H_G93BUm3

library(ggplot2) 
# Create histogram plot of estimates for both groups 
p<-ggplot(df_H_G93BUm3, aes(x=value))+
  geom_histogram(color="black", fill="black")+
  facet_grid(Group ~ .)
p
#Calculate the mean of each group
#Calculate the average Shannon diversity of each group using the package plyr
library(plyr)
mu <- ddply(df_H_G93BUm3, "Group", summarise, grp.mean=mean(value))
head(mu)

#add mean lines
p+geom_vline(data=mu, aes(xintercept=grp.mean, color="red"),
             linetype="dashed")


## 5.2.3 Calculating Power or Sample Size Using R Function power.t.test


# calculate variance
mu <- ddply(df_H_G93BUm3, "Group", summarise, grp.mean=mean(value));mu

#         Group grp.mean
#1    BUm3to3.5    2.504
#2  NOBUm3to3.5    2.205

var <- ddply(df_H_G93BUm3, "Group", summarise, grp.var=var(value));var

#        Group   value
#1   BUm3to3.5 0.02892
#2 NOBUm3to3.5 0.04349

n1 <- 9
n2 <-7
s1<-sqrt(0.02892)
s1
#[1] 0.1701
s2<-sqrt(0.04349)
s2
#[1] 0.2085
s=sqrt((n1-1)*s1^2+(n2-1)*s2^2)/(n1+n2-2)
s
#[1] 0.05012

power.t.test(n=2:10,delta=2.504-2.205,sd=0.05012)
df_P <-data.frame(n,power)
df_P

n = c(2, 3, 4, 5, 6, 7, 8, 9, 10)  
power = c(0.8324, 0.9994, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000)

power <- sapply(n, function (x) power.t.test(n=x, delta=2.504-2.205,sd=0.05012)$power)
plot(n, power, xlab  = "Sample Size per group", ylab  = "Power to reject null",
     main="Power curve for\n t-test with delta = 0.05",
     lwd=2, col="red", type="l")

abline(h = 0.90, col="blue")


power.t.test(n=2:10,delta=2.504-2.205,sd=0.05012, type = "one.sample" )


## 5.3 Power Analysis for Comparing Diversity Across More than Two Groups Using ANOVA	

### 5.3.2 Calculating Power or Sample Size Using R Function pwr.avova.test	
df_H_G93WTm1N4 <- filter(df_H_G6,Group%in%c("G93m1","WTm1","G93m4","WTm4"))
df_H_G93WTm1N4 

fit = lm(formula = value~Group,data=df_H_G93WTm1N4)
anova (fit)

install.packages("pwr")
library(pwr)
pwr.anova.test(f= 0.23,k=4,n=45:55,sig.level=0.05)
