# Linear mixed effects modelling in R

What is mixed effects modelling and why does it matter? 
Longitudinal data are often complex and messy. We can have different grouping factors like populations, species, sites where we collect the data, etc. Sample sizes might leave something to be desired too, especially if we are trying to fit complicated models with many parameters. On top of that, our data points might not be truly independent. For instance, we might be using quadrats within our sites to collect the data (and so there is structure to our data: quadrats are nested within the sites). This is why mixed models were developed, to deal with such messy data and to allow us to use all our data, even when we have low sample sizes, structured data and many covariates to fit. Oh, and on top of all that, mixed models allow us to save degrees of freedom compared to running standard linear models! Sounds good, doesn’t it? 

We are going to focus on a fictional study system, dragons, so that we don’t have to get too distracted with the specifics of this example. Imagine that we decided to train dragons and so we went out into the mountains and collected data on dragon intelligence ( testScore) as a prerequisite. We sampled individuals with a range of body lengths across three sites in eight different mountain ranges. 
```r
# load the data and have a look at it
load("dragons.RData")
head(dragons)

# Let's say we want to know how the body length affects test scores.
## Have a look at the data distribution:
hist(dragons$testScore)  

# seems close to normal distribution - good!
```
It is good practice to standardise your explanatory variables before proceeding so that they have a mean of zero (“centering”) and standard deviation of one (“scaling”). It ensures that the estimated coefficients are all on the same scale, making it easier to compare effect sizes. You can use scale() to do that:

scale() centers the data (the column mean is subtracted from the values in the column) and then scales it (the centered column values are divided by the column’s standard deviation)
```r
dragons$bodyLength2 <- scale(dragons$bodyLength)
head(dragons)
```
Back to our question: is test score affected by body length?
One way to analyse this data would be to try fitting a linear model to all our data, ignoring the sites and the mountain ranges for now.
Fit the model with testScore as the response and bodyLength2 as the predictor and have a look at the output
```r

library(lme4)

basic.lm <- lm(testScore ~ bodyLength2, data = dragons)

summary(basic.lm)

# Let's plot the data with ggplot2

library(ggplot2)

ggplot(dragons, aes(x = bodyLength, y = testScore)) +
  geom_point()+
  geom_smooth(method = "lm")

```
Okay, so both from the linear model and from the plot, it seems like bigger dragons do better in our intelligence test. That seems a bit odd: size shouldn’t really affect the test scores.

But are assumptions of the lm met?
```r
# Plot the residuals - the red line should be close to being flat, like the dashed grey line

plot(basic.lm, which = 1)  

# not perfect, but look alright

# Have a quick look at the  qqplot too - point should ideally fall onto the diagonal dashed line

plot(basic.lm, which = 2)  

# a bit off at the extremes, but that's often the case; again doesn't look too bad
```
However, what about observation independence? Are our data independent?
We collected multiple samples from eight mountain ranges
It's perfectly plausible that the data from within each mountain range are more similar to each other than the data from different mountain ranges, they are correlated.
```r
# Have a look at the data to see if above is true
boxplot(testScore ~ mountainRange, data = dragons)  

# certainly looks like something is going on here

# We could also plot it colouring points by mountain range
ggplot(dragons, aes(x = bodyLength, y = testScore, colour = mountainRange))+
  geom_point(size = 2)+
  theme_classic()+
    theme(legend.position = "none")
```
From the above plots it looks like our mountain ranges vary both in the dragon body length and in their test scores. This confirms that our observations from within each of the ranges aren't independent. We can't ignore that. 
So what do we do?

## Option 1: Run multiple analyses

We could run many separate analyses and fit a regression for each of the mountain ranges
```r
# Lets have a quick look at the data split by mountain range. We use the facet_wrap to do that

ggplot(aes(bodyLength, testScore), data = dragons) + geom_point() +
    facet_wrap(~ mountainRange) +
    xlab("length") + ylab("test score")
```
That’s eight analyses. Oh wait, we also have different sites in each mountain range, which similarly to mountain ranges aren’t independent. So we could run an analysis for each site in each range separately.
To do the above, we would have to estimate a slope and intercept parameter for each regression. That’s two parameters, three sites and eight mountain ranges, which means 48 parameter estimates (2 x 3 x 8 = 48)! Moreover, the sample size for each analysis would be only 20 (dragons per site).
This presents problems: not only are we hugely decreasing our sample size, but we are also increasing chances of a Type I Error (where you falsely reject the null hypothesis) by carrying out multiple comparisons. Not ideal!

## Option 2: Modify the model

We want to use all the data, but account for the data coming from different mountain ranges
```r
# let's add mountain range as a fixed effect to our basic.lm

mountain.lm <- lm(testScore ~ bodyLength2 + mountainRange, data = dragons)
summary(mountain.lm)
```
Now body length is not significant. But let’s think about what we are doing here for a second. The above model is estimating the difference in test scores between the mountain ranges - we can see all of them in the model output returned by summary(). But we are not interested in quantifying test scores for each specific mountain range: we just want to know whether body length affects test scores and we want to simply control for the variation coming from mountain ranges. This is what we refer to as “random factors” and so we arrive at mixed effects models. Ta-daa!

## Option 3: Linear mixed effects models (LMM)

A mixed model is a good choice here: it will allow us to use all the data we have (higher sample size) and account for the correlations between data coming from the sites and mountain ranges.  We will also estimate fewer parameters and avoid problems with multiple comparisons that we would encounter while using separate regressions

Fixed and random effects
Let’s talk a little about the difference between fixed and random effects first. It’s important to note that this difference has little to do with the variables themselves, and a lot to do with your research question! In many cases, the same variable could be considered either a random or a fixed effect (and sometimes even both at the same time!) so always refer to your questions and hypotheses to construct your models accordingly.

Should my variables be fixed or random effects?
In broad terms, fixed effects are variables that we expect will have an effect on the dependent/response variable: they’re what you call explanatory variables in a standard linear regression. In our case, we are interested in making conclusions about how dragon body length impacts the dragon’s test score. So body length is a fixed effect and test score is the dependent variable.
On the other hand, random effects are usually grouping factors for which we are trying to control. They are always categorical, as you can’t force R to treat a continuous variable as a random effect. A lot of the time we are not specifically interested in their impact on the response variable, but we know that they might be influencing the patterns we see. Additionally, the data for our random effect is just a sample of all the possibilities: with unlimited time and funding we might have sampled every mountain where dragons live, but we usually tend to generalise results to a whole population based on representative sampling.
In our particular case, we are looking to control for the effects of mountain range. We haven’t sampled all the mountain ranges in the world (we have eight) so our data are just a sample of all the existing mountain ranges. We are not really interested in the effect of each specific mountain range on the test score: we hope our model would also be generalisable to dragons from other mountain ranges! However, we know that the test scores from within the ranges might be correlated so we want to control for that. If we specifically chose eight particular mountain ranges a priori and we were interested in those ranges and wanted to make predictions about them, then mountain range would be fitted as a fixed effect
Note that the golden rule is that you generally want your random effect to have at least five levels. So, for instance, if we wanted to control for the effects of dragon’s sex on intelligence, we would fit sex (a two level factor: male or female) as a fixed, not random, effect.

Let’s fit our first mixed model
We have a response variable, the test score and we are attempting to explain part of the variation in test score through fitting body length as a
fixed effect. But the response variable has some residual variation (i.e. unexplained variation) associated with mountain ranges. By using random effects, we are modeling that unexplained variation through variance

```r
# Fit the LMM
library(lme4)

# LMM general formula: lmer(response variable ~ fixed effect + random effect, data)
  
mixed.lmer <- lmer(testScore ~ bodyLength2 + (1|mountainRange), data = dragons)
summary(mixed.lmer)
```
Linear mixed model fit by REML
Formula: testScore ~ bodyLength2 + (1 | mountainRange)
   Data: dragons

REML criterion at convergence: 3985.6

Scaled residuals: 
    Min      1Q  Median      3Q     Max 
-3.4815 -0.6513  0.0066  0.6685  2.9583 

Random effects:
 Groups        Name        Variance Std.Dev.
 mountainRange (Intercept) 339.7    18.43   
 Residual                  223.8    14.96   
Number of obs: 480, groups:  mountainRange, 8

Fixed effects:
            Estimate Std. Error t value
(Intercept)  50.3860     6.5517   7.690
bodyLength2   0.5377     1.2750   0.422

Correlation of Fixed Effects:
            (Intr)
bodyLength2 0.000 

We can see the variance for mountainRange = 339.7. Mountain ranges are clearly important: they explain a lot of variation. How do we know that? We can take the variance for the mountainRange and divide it by the total variance:

339.7/(339.7 + 223.8)=0.60 ~60 %

So the differences between mountain ranges explain ~60% of the variance that’s “left over” after the variance explained by our fixed effects

## compare LMM models

You should use maximum likelihood when comparing models with different fixed effects, so REML = FALSE
```r
full.lmer <- lmer(testScore ~ bodyLength2 + (1|mountainRange), data = dragons, REML = FALSE)

# a reduced LMM model in which we dropped our fixed effect
reduced.lmer <- lmer(testScore ~ 1 + (1|mountainRange), data = dragons, REML = FALSE)

anova(reduced.lmer, full.lmer)  
# the two models are not significantly different

```
