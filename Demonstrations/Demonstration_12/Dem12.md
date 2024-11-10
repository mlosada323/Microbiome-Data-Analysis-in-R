# Linear mixed effects modelling in R

# What is mixed effects modelling and why does it matter? 
# Longitudinal data are often complex and messy. We can have different grouping factors like populations, species, sites where we collect the data, etc. Sample sizes might leave something to be desired too, especially if we are trying to fit complicated models with many parameters. On top of that, our data points might not be truly independent. For instance, we might be using quadrats within our sites to collect the data (and so there is structure to our data: quadrats are nested within the sites). This is why mixed models were developed, to deal with such messy data and to allow us to use all our data, even when we have low sample sizes, structured data and many covariates to fit. Oh, and on top of all that, mixed models allow us to save degrees of freedom compared to running standard linear models! Sounds good, doesn’t it? 

# We are going to focus on a fictional study system, dragons, so that we don’t have to get too distracted with the specifics of this example. Imagine that we decided to train dragons and so we went out into the mountains and collected data on dragon intelligence ( testScore) as a prerequisite. We sampled individuals with a range of body lengths across three sites in eight different mountain ranges. 

# load the data and have a look at it
load("dragons.RData")
head(dragons)

# Let's say we want to know how the body length affects test scores.
## Have a look at the data distribution:
hist(dragons$testScore)  

# seems close to normal distribution - good!

# It is good practice to standardise your explanatory variables before proceeding so that they have a mean of zero (“centering”) and standard deviation of one (“scaling”). It ensures that the estimated coefficients are all on the same scale, making it easier to compare effect sizes. You can use scale() to do that:

dragons$bodyLength2 <- scale(dragons$bodyLength)
head(dragons)

# scale() centers the data (the column mean is subtracted from the values in the column) and then scales it (the centered column values are divided by the column’s standard deviation).

# Back to our question: is test score affected by body length?

# One way to analyse this data would be to try fitting a linear model to all our data, ignoring the sites and the mountain ranges for now.
# Fit the model with testScore as the response and bodyLength2 as the predictor and have a look at the output

library(lme4)

basic.lm <- lm(testScore ~ bodyLength2, data = dragons)

summary(basic.lm)

# Let's plot the data with ggplot2

library(ggplot2)

ggplot(dragons, aes(x = bodyLength, y = testScore)) +
  geom_point()+
  geom_smooth(method = "lm")

# Okay, so both from the linear model and from the plot, it seems like bigger dragons do better in our intelligence test. That seems a bit odd: size shouldn’t really affect the test scores.

### Are assumptions of the lm met?

## Plot the residuals - the red line should be close to being flat, like the dashed grey line

plot(basic.lm, which = 1)  

# not perfect, but look alright

## Have a quick look at the  qqplot too - point should ideally fall onto the diagonal dashed line

plot(basic.lm, which = 2)  

# a bit off at the extremes, but that's often the case; again doesn't look too bad

## However, what about observation independence? Are our data independent?
## We collected multiple samples from eight mountain ranges
## It's perfectly plausible that the data from within each mountain range are more similar to each other than the data from different mountain ranges, they are correlated.

## Have a look at the data to see if above is true
boxplot(testScore ~ mountainRange, data = dragons)  

# certainly looks like something is going on here

# We could also plot it colouring points by mountain range
ggplot(dragons, aes(x = bodyLength, y = testScore, colour = mountainRange))+
  geom_point(size = 2)+
  theme_classic()+
    theme(legend.position = "none")

# From the above plots it looks like our mountain ranges vary both in the dragon body length and in their test scores. This confirms that our observations from within each of the ranges aren't independent. We can't ignore that. 
# So what do we do?

## Option 1: Run multiple analyses

# We could run many separate analyses and fit a regression for each of the mountain ranges

# Lets have a quick look at the data split by mountain range. We use the facet_wrap to do that

ggplot(aes(bodyLength, testScore), data = dragons) + geom_point() +
    facet_wrap(~ mountainRange) +
    xlab("length") + ylab("test score")

# That’s eight analyses. Oh wait, we also have different sites in each mountain range, which similarly to mountain ranges aren’t independent. So we could run an analysis for each site in each range separately.
# To do the above, we would have to estimate a slope and intercept parameter for each regression. That’s two parameters, three sites and eight mountain ranges, which means 48 parameter estimates (2 x 3 x 8 = 48)! Moreover, the sample size for each analysis would be only 20 (dragons per site).
# This presents problems: not only are we hugely decreasing our sample size, but we are also increasing chances of a Type I Error (where you falsely reject the null hypothesis) by carrying out multiple comparisons. Not ideal!

## Option 2: Modify the model

# We want to use all the data, but account for the data coming from different mountain ranges
# let's add mountain range as a fixed effect to our basic.lm

mountain.lm <- lm(testScore ~ bodyLength2 + mountainRange, data = dragons)
summary(mountain.lm)

# Now body length is not significant. But let’s think about what we are doing here for a second. The above model is estimating the difference in test scores between the mountain ranges - we can see all of them in the model output returned by summary(). But we are not interested in quantifying test scores for each specific mountain range: we just want to know whether body length affects test scores and we want to simply control for the variation coming from mountain ranges. This is what we refer to as “random factors” and so we arrive at mixed effects models. Ta-daa!

## Option 3: Linear mixed effects models (LMM)

# A mixed model is a good choice here: it will allow us to use all the data we have (higher sample size) and account for the correlations between data coming from the sites and mountain ranges. We will also estimate fewer parameters and avoid problems with

## First LMM

library(lme4)

mixed.lmer <- lmer(testScore ~ bodyLength2 + (1|mountainRange), data = dragons)
summary(mixed.lmer)



### plots

plot(mixed.lmer)

qqnorm(resid(mixed.lmer))
qqline(resid(mixed.lmer))  # points fall nicely onto the line - good!

### summary

### variance accounted for by mountain ranges

339.7/(339.7 + 223.8)  # ~60 %

##-- implicit vs explicit nesting --##

head(dragons)  # we have site and mountainRange
str(dragons)  # we took samples from three sites per mountain range and eight mountain ranges in total

### create new "sample" variable

dragons <- within(dragons, sample <- factor(mountainRange:site))

##----- Second mixed model -----##

### model

mixed.WRONG <- lmer(testScore ~ bodyLength2 + (1|mountainRange) + (1|site), data = dragons)  # treats the two random effects as if they are crossed
summary(mixed.WRONG)

# syntax 1 for lmer

mixed.lmer2 <- lmer(testScore ~ bodyLength2 + (1|mountainRange) + (1|sample), data = dragons)  # the syntax stays the same, but now the nesting is taken into account
summary(mixed.lmer2)

# syntax 2 for lmer

mixed.lmer2 <- lmer(testScore ~ bodyLength2 + (1|mountainRange/site), data = dragons)  # the syntax stays the same, but now the nesting is taken into account
summary(mixed.lmer2)

# syntax 3 for lmer

mixed.lmer2 <- lmer(testScore ~ bodyLength2 + (1|mountainRange) + (1|mountainRange:site), data = dragons)  # the syntax stays the same, but now the nesting is taken into account
summary(mixed.lmer2)

### summary

### plot

(mm_plot <- ggplot(dragons, aes(x = bodyLength, y = testScore, colour = site)) +
   facet_wrap(~mountainRange, nrow=2) +   # a panel for each mountain range
   geom_point(alpha = 0.5) +
   theme_classic() +
   geom_line(data = cbind(dragons, pred = predict(mixed.lmer2)), aes(y = pred), size = 1) +  # adding predicted line from mixed model 
   theme(legend.position = "none",
         panel.spacing = unit(2, "lines"))  # adding space between panels
)

##----- Model selection for the keen -----##

### full model with random slopes and intercepts

# check intercepts

coef(mixed.lmer2)

# random slopes - it gives a singlular fit error

mixed.ranslope <- lmer(testScore ~ bodyLength2 + (1+bodyLength2|mountainRange/site), data = dragons) 

summary(mixed.slope)

(mm_plot <- ggplot(dragons, aes(x = bodyLength, y = testScore, colour = site)) +
   facet_wrap(~mountainRange, nrow=2) +   # a panel for each mountain range
   geom_point(alpha = 0.5) +
   theme_classic() +
   geom_line(data = cbind(dragons, pred = predict(mixed.ranslope)), aes(y = pred), size = 1) +  # adding predicted line from mixed model 
   theme(legend.position = "none",
         panel.spacing = unit(2, "lines"))  # adding space between panels
)

install.packages("ggeffects")
library(ggeffects)

# Extract the prediction data frame
pred.mm <- ggpredict(mixed.lmer2, terms = c("bodyLength2"))  # this gives overall predictions for the model

# Plot the predictions 

(ggplot(pred.mm) + 
   geom_line(aes(x = x, y = predicted)) +          # slope
   geom_ribbon(aes(x = x, ymin = predicted - std.error, ymax = predicted + std.error), 
               fill = "lightgrey", alpha = 0.5) +  # error band
   geom_point(data = dragons,                      # adding the raw data (scaled values)
              aes(x = bodyLength2, y = testScore, colour = mountainRange)) + 
   labs(x = "Body Length (indexed)", y = "Test Score", 
        title = "Body length does not affect intelligence in dragons") + 
   theme_minimal()
)

ggpredict(mixed.lmer2, terms = c("bodyLength2", "mountainRange"), type = "re") %>% 
  plot() +
  labs(x = "Body Length", y = "Test Score", title = "Effect of body size on intelligence in dragons") + 
  theme_minimal()

### reduced model

reduced.lmer <- lmer(testScore ~ 1 + (1|mountainRange) + (1|sample), 
                     data = dragons, REML = FALSE)


### comparison
# You should use maximum likelihood when comparing models with different fixed effects, so REML = FALSE

full.lmer <- lmer(testScore ~ bodyLength2 + (1|mountainRange) + (1|sample), 
                  data = dragons, REML = FALSE)

reduced.lmer <- lmer(testScore ~ 1 + (1|mountainRange) + (1|sample), 
                     data = dragons, REML = FALSE)

anova(reduced.lmer, full.lmer)  # the two models are not significantly different

full.lmer <- lmer(testScore ~ bodyLength2 + (1|mountainRange) + (1|sample), 
                  data = dragons, REML = FALSE)

reduced.lmer <- lmer(testScore ~ bodyLength2 + (1|mountainRange), 
                     data = dragons, REML = FALSE)

anova(reduced.lmer, full.lmer)  # the two models are not significantly different


