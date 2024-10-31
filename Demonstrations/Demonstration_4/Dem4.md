# Demonstration 5

Complete the following demonstration in RStudio. Create a markdown file of your script. You can follow detail instructions in Xia et al. (2018), Chapter 5: Power and Sample Size Calculations for Microbiome Data. All the sections below match the sections in the book

# Power Analysis for Microbiome Data

## 5.2.2 Diversity Data for ALS Study
```r
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

```
## 5.2.3 Calculating Power or Sample Size Using R Function power.t.test
Here, we focus on illustrating how to calculate the power or sample size using R software. In R, the function power.t.test() in basic R and the function pwr.t.test() in the pwr package can be used to conduct power analysis. We use the power.t.test. The usage of this function is shown below:
power.t.test (n = sample size, delta = effect size, sd = standard deviation, sig. level = 0.05, power = NULL, type = c(“two.sample”, “one.sample”, “paired”), alternative = (“two.sided”, “one.sided”))
where, n is the number of sample size per group, delta is true difference in means, sd is the standard deviation, sig.level is the significance level (Type I error probability), power is the power of test (1 minus Type II error probability), type is the type of t test, and alternative is one-or-two sided test.
Since the standard deviation of the mean difference is unknown, it needs to be estimated using formula in this chapter (5.5).

```r
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

  n power
1 2 0.8324
2 3 0.9994
3 4 1.0000
4 5 1.0000
5 6 1.0000
6 7 1.0000
7 8 1.0000
8 9 1.0000
9 10 1.0000
```
From above power analysis, we can see that a size sample of 2 G93A mice per group, randomly assigned to butyrate treatment or no treatment control, will provide 83% power to reject the null hypothesis of no difference in the Shannon diversity in the two groups. If the sample size increases to 3 per group, the power will increase to more than 99%. 
```r
# We can generate power and sample size graphs to visualize the power and sample size we need to reject the null hypothesis using following R codes

n = c(2, 3, 4, 5, 6, 7, 8, 9, 10)  
power = c(0.8324, 0.9994, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000)

power <- sapply(n, function (x) power.t.test(n=x, delta=2.504-2.205,sd=0.05012)$power)
plot(n, power, xlab  = "Sample Size per group", ylab  = "Power to reject null",
     main="Power curve for\n t-test with delta = 0.05",
     lwd=2, col="red", type="l")

abline(h = 0.90, col="blue")

power.t.test(n=2:10,delta=2.504-2.205,sd=0.05012, type = "one.sample" )
```

## 5.3 Power Analysis for Comparing Diversity Across More than Two Groups Using ANOVA	

### 5.3.2 Calculating Power or Sample Size Using R Function pwr.avova.test	
df_H_G93WTm1N4 <- filter(df_H_G6,Group%in%c("G93m1","WTm1","G93m4","WTm4"))
df_H_G93WTm1N4 

fit = lm(formula = value~Group,data=df_H_G93WTm1N4)
anova (fit)

install.packages("pwr")
library(pwr)
pwr.anova.test(f= 0.23,k=4,n=45:55,sig.level=0.05)
