# make figure 3 from https://doi.org/10.1016/j.scitotenv.2023.164828

library(tidyverse)
library(plm)
library(modelbased)
library(patchwork)
library(fixest)
library(texreg)


# a function to calculate marginal effects of each fire severity category
margeff <- function(model, effect, interaction, varcov="default", conf=.95) {
  # Extract Variance Covariance matrix
  if (varcov == "default"){
    covMat = vcov(model)
  }else{
    covMat = varcov
  }
  
  # Get coefficients of variables
  beta_1 = model$coefficients[[effect]]
  beta_3 = model$coefficients[[interaction]]
  
  # Get coefficients of variables
  beta_1 = model$coefficients[[effect]]
  beta_3 = model$coefficients[[interaction]]
  
  # Create list of moderator values at which marginal effect is evaluated
  x_2 <- c(0, 1)
  
  # Compute marginal effects
  delta_1 = beta_1 + beta_3*x_2
  
  # Compute variances
  var_1 = covMat[effect,effect] + (x_2^2)*covMat[interaction, interaction] + 2*x_2*covMat[effect, interaction]
  
  # Standard errors
  se_1 = sqrt(var_1)
  
  # Upper and lower confidence bounds
  z_score = qnorm(1 - ((1 - conf)/2))
  upper_bound = delta_1 + z_score*se_1
  lower_bound = delta_1 - z_score*se_1
  
  df <- data.frame(delta_1, lower_bound, upper_bound)
  
  df
}

# make sure working directory is set to ~/data/samples
# read in our samples from dfout and dffit
# dfout contains the undifferenced Levels model data
# dffit contains the differenced SFD model data
# we will run four models 

# model 1: levels_manual_1_2345 (Levels with severity categories 1_2345 where 1 = unburnt and 2345 = all severity categories grouped)
# model 2: sfd_manual_1_2345 (SFD with severity categories 1_2345 where 1 = unburnt and 2345 = all severity categories grouped)
# model 3: levels_manual_12345 (Levels with severity categories 12345 where marginal effects of all severity categories are calculated individually)
# model 4: sfd_manual_12345 (SFD with severity categories 12345 where where marginal effects of all severity categories are calculated individually)


loc = 'C://Users//1dand//PycharmProjects//njf-fire-effects//data//samples'
setwd(loc)

dfout_1_2345 = read.csv("dfout_a50_g1000-2000.csv")
dfout_1_2345$f_prop = dfout_1_2345$f_prop * 100 # our response variable is the proportion (%) of flowering trees relative to non-flowering trees

dffit_1_2345 = read.csv("dffit_a50_g1000-2000.csv")
dffit_1_2345$f_prop = dffit_1_2345$f_prop * 100


dfout_12345 = read.csv("dfout_a50_g1000-2000.csv")
dfout_12345$f_prop = dfout_12345$f_prop * 100

dffit_12345 = read.csv("dffit_a50_g1000-2000.csv")
dffit_12345$f_prop = dffit_12345$f_prop * 100



# estimate the models
# burnt/unburnt (1_2345)
levels_manual_1_2345 <- feols(f_prop ~ cat_1  + cat_1_int + ys + factor(tpi_zone), cluster = ~se5, data = dfout_1_2345)
sfd_manual_1_2345 <-  feols(f_prop ~ cat_1  + cat_1_int + ys , cluster = ~se5, data = dffit_1_2345)

# all severities
levels_manual_12345 <- feols(f_prop ~ cat_1 + cat_2 + cat_3 + cat_4 + ys + 
                               cat_1_int + cat_2_int + cat_3_int + cat_4_int + factor(tpi_zone), cluster = ~se5, data = dfout_12345)
sfd_manual_12345 <- feols(f_prop ~ cat_1 + cat_2 + cat_3 + cat_4 + ys + 
                            cat_1_int + cat_2_int + cat_3_int + cat_4_int, cluster = ~se5, data = dffit_12345)

# calculate marginal effects
# 1-2345 levels and sfd
meL_1_2345_cat1 <- margeff(levels_manual_1_2345, "ys", "cat_1_int")[2,]
meL_1_2345_cat2345 <- margeff(levels_manual_1_2345, "ys", "cat_1_int")[1,]

meSFD_1_2345_cat1 <- margeff(sfd_manual_1_2345, "ys", "cat_1_int")[2,]
meSFD_1_2345_cat2345 <- margeff(sfd_manual_1_2345, "ys", "cat_1_int")[1,]

data1 = rbind(meL_1_2345_cat1, meL_1_2345_cat2345, meSFD_1_2345_cat1, meSFD_1_2345_cat2345)
data1$model <- c('Levels','Levels','SFD','SFD')
data1$sev <- c('1-Unburnt (within fire)','All severities > 1-Unburnt')
m1_obs1 = unlist(list(c(nrow(subset(dfout_1_2345 , in_out == 1 & cat_1 == 1)), 
                        nrow(subset(dfout_1_2345 , in_out == 1 & cat_2345 == 1)))))
m1_obs2 = m1_obs1 * 2
data1$obs = c(m1_obs2, m1_obs1 )

data1$morder = 1
data1$sorder <- c(1,2)

# levels 
meL_12345_cat1 <- margeff(levels_manual_12345, "ys", "cat_1_int")[2,]
meL_12345_cat2 <- margeff(levels_manual_12345, "ys", "cat_2_int")[2,]
meL_12345_cat3 <- margeff(levels_manual_12345, "ys", "cat_3_int")[2,]
meL_12345_cat4 <- margeff(levels_manual_12345, "ys", "cat_4_int")[2,]
meL_12345_cat5 <- margeff(levels_manual_12345, "ys", "cat_4_int")[1,]

# sfd
meSFD_12345_cat1 <- margeff(sfd_manual_12345, "ys", "cat_1_int")[2,]
meSFD_12345_cat2 <- margeff(sfd_manual_12345, "ys", "cat_2_int")[2,]
meSFD_12345_cat3 <- margeff(sfd_manual_12345, "ys", "cat_3_int")[2,]
meSFD_12345_cat4 <- margeff(sfd_manual_12345, "ys", "cat_4_int")[2,]
meSFD_12345_cat5 <- margeff(sfd_manual_12345, "ys", "cat_4_int")[1,]

# put into a table to plot
data3 = rbind( meL_12345_cat2, meL_12345_cat3, meL_12345_cat4, meL_12345_cat5,
               meSFD_12345_cat2, meSFD_12345_cat3, meSFD_12345_cat4, meSFD_12345_cat5)
data3$model <- c('Levels','Levels','Levels','Levels','SFD','SFD','SFD','SFD')
data3$sev <- c('2-Low canopy scorch','3-Mid canopy scorch','4-High canopy scorch','5-Canopy burnt')

m2_obs1 = unlist(list(c(#nrow(subset(dfout_12345 , in_out == 1 & hsev == 1)), 
  nrow(subset(dfout_12345 , in_out == 1 & hsev == 2)),
  nrow(subset(dfout_12345 , in_out == 1 & hsev == 3)),
  nrow(subset(dfout_12345 , in_out == 1 & hsev == 4)),
  nrow(subset(dfout_12345 , in_out == 1 & hsev == 5)))))


m2_obs2 = m2_obs1 * 2
data3$obs = c(m2_obs2, m2_obs1)

data1$morder = 1
data1$sorder <- c(1,2)

data3$morder = 3
data3$sorder <- c(2,3,4,5)

data = rbind(data1, data3)

cat1 = "#2c7bb6"
cat2 = "#abd9e9"
cat3 = "#ffffbf"
cat4 = "#fdae61"
cat5 = "#d7191c"

data$color = c(cat1, 'gray', cat1, 'gray',  
               cat2, cat3, cat4, cat5, cat2, cat3, cat4, cat5)

data <- data[with(data, order(morder, sorder)),]
data$x = seq(-nrow(data), -1)


# make the figure 

par(mgp=c(2,0.6,0), mar=c(3,0,1,0), lend=1)
windowsFonts(A = windowsFont("Arial"))  # Specify font
# 20 = number of rows

xlab = expression(paste("                                            ", Delta, " Flower % / ", " ", Delta, " Year"))
plot(1,type="n",xlim=c(-0.8,0.8),ylim=c(1, 13),axes=F,ylab='',xlab=xlab,cex.lab=0.8,  family = "A")
tp <- data$delta_1 
cilo <- data$lower_bound
cihi <- data$upper_bound 
xx <- seq(0,0.4,0.1)
colz = data$color#rep("black",length(tp))

n=length(tp):1
nn <- n
abline(v=xx,lwd=0.5,col="grey",lty=2)
abline(v=0,lwd=2)

lz <- c(10.5, 16.5)

# the points
segments(cilo,nn,cihi,nn,lwd =1.5)
points(tp,nn,pch=21,bg=colz,cex=1.86,lwd=1.8)

# adding the effect size and 95% CI
txt_effect <- paste0(formatC(round(tp,2),digits=2,format="f"), " (", formatC(round(cilo,2),digits=2,format="f"), ",", formatC(round(cihi,2),digits=2,format="f"), ")")
text(0.42,nn, txt_effect, cex=0.65,pos=4)
text(0.40, 13, 'Effect size (95% CI)', cex=0.75,pos=4)

# adding in Levels SFD column
model_type <- data$model
text(-0.16, nn, model_type, cex=0.65,pos=4)
text(-0.16, 13, 'Model', cex=0.75,pos=4)

# adding in observations count
obs_counts = data$obs
obs_count <- format(obs_counts,big.mark = ",")
text(-0.38, nn, obs_count, cex=0.65,pos=4)
text(-0.41, 13, 'Observations', cex=0.75,pos=4)

# adding in severity info
sevlist = c("1-Unburnt (within fire)", "All categories (> 1-unburnt)", "2-Low canopy scorch","3-Mid canopy scorch",  "4-High canopy scorch", "5-Canopy consumed")
sev_type <- sevlist#data$sev
text(-0.78, c(12, 10, 8, 6, 4, 2), sev_type, cex=0.65,pos=4)
text(-0.78, 13, 'Fire Severity', cex=0.75,pos=4)

axis(1,at=xx,labels=xx, cex.axis=0.75)







