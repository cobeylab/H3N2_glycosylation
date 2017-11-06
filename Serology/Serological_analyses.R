#############################################
# Analysis of VE Zost et al. (2017) 
#
#
# This script includes the code to generate 
# the figures and perform statistical
# analysis of the serological data.
#
#
# SECTIONS
# - Data import & preparation (must execute)
# - (Statistics) Correlations by age, 160K/T
# - (Figures) Fold change by vaccine group
# - (Figures) Vaccination history & boost
#
#############################################


# WT = glycosylated (T160)
# 160K = unglycosylated

# FZ = Fluzone, egg 
# FB = Flublok, recombinant 
# FCV = Flucelvax, MDCK 

require(dplyr)
require(ggplot2)

############################################
####### Data import & preparation ##########
############################################

# Import FRNT data (with two columns of Vx history)
titers_raw <- read.table("2017-04-19_FRNT_for_import.csv",header=TRUE,sep=",",na.strings="Unknown")

# Merge data and add basic statistics
titers <- titers_raw
titers$YOB <- 2017 - titers$Age

geo_mean <- function(t1,t2,t3) {
	tmp_df <- data.frame(t1,t2,t3)
	tmp_df <- log(tmp_df)
	tmp_df <- rowMeans(tmp_df)
	tmp_df <- exp(tmp_df)
	return(tmp_df)
}
titers$GMT_Pre_WT <- geo_mean(titers$Pre_WT_T1,titers$Pre_WT_T2,titers$Pre_WT_T3)
titers$GMT_Pre_160K <- geo_mean(titers$Pre_160K_T1,titers$Pre_160K_T2,titers$Pre_160K_T3)
titers$GMT_Post_WT <- geo_mean(titers$Post_WT_T1,titers$Post_WT_T2,titers$Post_WT_T3)
titers$GMT_Post_160K <- geo_mean(titers$Post_160K_T1,titers$Post_160K_T2,titers$Post_160K_T3)
titers$GM_FC_WT <- titers$GMT_Post_WT/titers$GMT_Pre_WT
titers$GM_FC_160K <- titers$GMT_Post_160K/titers$GMT_Pre_160K
titers$Pre_ratio_WT_160K <- titers$GMT_Pre_WT/titers$GMT_Pre_160K


# Define ordinal variables based on vaccination history
# and add probability of primary exposure to H3
titers$Vx_hx <- NA
for (p in 1:dim(titers)[1]) {
	if ( is.na(titers$Vx_2014[p])==FALSE & is.na(titers$Vx_2015[p])==FALSE ) {
		if ( (titers$Vx_2014[p] == 'Inactivated') & (titers$Vx_2015[p] =='Inactivated') ) {
			titers$Vx_hx[p] <- 2
		} else if ( (titers$Vx_2014[p] == 'Inactivated') & (titers$Vx_2015[p] =='None') ) {
				titers$Vx_hx[p] <- 1
		} else if ( (titers$Vx_2014[p] == 'None') & (titers$Vx_2015[p]=='Inactivated') ) {
			titers$Vx_hx[p] <- 1
		} else if ( (titers$Vx_2014[p] == 'None') & (titers$Vx_2015[p] =='None') ) {
			titers$Vx_hx[p] <- 0
		}
	}
}

#### Calculate people with any unknown status
any_unknown <- titers %>% filter(is.na(Vx_2014) | is.na(Vx_2015))
nrow(any_unknown)

#### Define partitions
titers.FB <- subset(titers,titers$Group=="FB")
titers.FCV <- subset(titers,titers$Group=="FCV")
titers.FZ <- subset(titers,titers$Group=="FZ")
titers.FB_FZ <- titers %>% filter(Group=="FB" | Group=="FZ")
titers.FB_FCV <- titers %>% filter(Group=="FB" | Group=="FCV")
No_vax <- titers %>% filter(Vx_2014=="None" & Vx_2015=="None")
Only_2014 <- titers %>% filter(Vx_2014=="Inactivated" & Vx_2015=="None")
Only_2015 <- titers %>% filter(Vx_2014=="None" & Vx_2015=="Inactivated")
Both <- titers %>% filter(Vx_2014=="Inactivated" & Vx_2015=="Inactivated")
titers.post1979 <- titers %>% filter(YOB >= 1979)
titers.pre1979 <- titers %>% filter(YOB < 1979)

both2015 <- data.frame(rbind(Both,Only_2015))
none2014 <- data.frame(rbind(Only_2014,No_vax))
nonboth_df <- data.frame(rbind(No_vax,Only_2014,Only_2015))
both_df <- data.frame(Both)
any_df <- data.frame(rbind(Both,Only_2014,Only_2015))
none_df <- data.frame(No_vax)

FZ_col <- rgb(237, 163, 28,100,maxColorValue=255)
FCV_col <- rgb(41,137,192,100,maxColorValue=255)
FB_col <- rgb(18,153,20,100,maxColorValue=255)
WT_col <- rgb(0,0,1,1/2)
K_col <- rgb(1,0,0,1/2)


############################################
#### (Statistics) Correlations by age ######
############################################

cor.test(titers$Age,titers$GMT_Pre_WT,method='spearman')
# Spearman's rank correlation rho

# data:  titers$Age and titers$GMT_Pre_WT
# S = 65015, p-value = 0.2563
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
       # rho 
# -0.1375183 

# Warning message:
# In cor.test.default(titers$Age, titers$GMT_Pre_WT, method = "spearman") :
  # Cannot compute exact p-value with ties


cor.test(titers$Age,titers$GMT_Pre_160K,method='spearman')
	# Spearman's rank correlation rho

# data:  titers$Age and titers$GMT_Pre_160K
# S = 75621, p-value = 0.006371
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
     # rho 
# -0.32309 

# Warning message:
# In cor.test.default(titers$Age, titers$GMT_Pre_160K, method = "spearman") :
  # Cannot compute exact p-value with ties 

cor.test(titers$Age,titers$GMT_Post_WT,method='spearman')
	# Spearman's rank correlation rho

# data:  titers$Age and titers$GMT_Post_WT
# S = 68533, p-value = 0.09849
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
      # rho 
# -0.199081 

# Warning message:
# In cor.test.default(titers$Age, titers$GMT_Post_WT, method = "spearman") :
  # Cannot compute exact p-value with ties

cor.test(titers$Age,titers$GMT_Post_160K,method='spearman')
# Spearman's rank correlation rho

# data:  titers$Age and titers$GMT_Post_160K
# S = 79850, p-value = 0.0006653
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
       # rho 
# -0.3970838 

# Warning message:
# In cor.test.default(titers$Age, titers$GMT_Post_160K, method = "spearman") :
  # Cannot compute exact p-value with ties

# 160K v. 160T, pre-vaccination
wilcox.test(titers$GMT_Pre_160K,titers$GMT_Pre_WT,"greater")
# Wilcoxon rank sum test with continuity correction
#
# data:  titers$GMT_Pre_160K and titers$GMT_Pre_WT
# W = 3976.5, p-value = 7.703e-11
# alternative hypothesis: true location shift is greater than 0

# 160K v. 160T, post-vaccination
wilcox.test(titers$GMT_Post_160K,titers$GMT_Post_WT,"greater")
# Wilcoxon rank sum test with continuity correction
#
# data:  titers$GMT_Post_160K and titers$GMT_Post_WT
# W = 3990, p-value = 6.458e-11
# alternative hypothesis: true location shift is greater than 0

  
############################################
## (Figures) Fold change by vaccine group ##
############################################

pdf("FC_group_combined.pdf",width=7.2,height=3.2)
par(mfrow=c(1,2),mar=c(4, 4, 1, 1) + 0.1)
plot(titers$Group,log(titers$GM_FC_WT,2),xlab="Vaccine type",ylab=expression("T160 HA fold change"),col=alpha('blue',0.5),ylim=c(-1,7),axes=FALSE,outpch=NA)
points(jitter(as.numeric(titers$Group),factor=0.4),log(titers$GM_FC_WT,2),axes=FALSE,pch=16,col=alpha('blue',0.5))
axis(1,at=c(1,2,3),labels=c("Flublok","Flucelvax","Fluzone"))
axis(2,at=c(0,2,4,6),labels=c(1,4,16,64))
plot(titers$Group,log(titers$GM_FC_160K,2),xlab="Vaccine type",ylab=expression("K160 HA fold change"),col=alpha('red',0.5),ylim=c(-1,7),axes=FALSE, pch=16,outpch=NA)
points(jitter(as.numeric(titers$Group),factor=0.4),log(titers$GM_FC_160K,2),axes=FALSE,pch=16,col=alpha('red',0.5))
axis(1,at=c(1,2,3),labels=c("Flublok","Flucelvax","Fluzone"))
axis(2,at=c(0,2,4,6),labels=c(1,4,16,64))
dev.off()


############################################
## (Figures) Vaccination history and boost #
############################################
pdf("Pre_post_FC_Vx_hx.pdf",width=8,height=6)
par (mfrow=c(2,3))

# 160K pre-vacc titer by vaccination history
plot(jitter(titers$Vx_hx),titers$GMT_Pre_160K,xlab="Vaccination history",ylab="Pre-vaccination titer",main="K160 HA",log="y",ylim=c(7,10245),axes=FALSE,col=alpha('red',0.5),pch=16)
axis(1,at=c(0,1,2),labels=c("None","One","Both"))
axis(2, at = c(10,80,640,5120))

# 160K post-vacc titer by vaccination history
plot(jitter(titers$Vx_hx),titers$GMT_Post_160K,xlab="Vaccination history",ylab="Post-vaccination titer",main="K160 HA",log="y",ylim=c(7,10245),axes=FALSE,col=alpha('red',0.5),pch=16)
axis(1,at=c(0,1,2),labels=c("None","One","Both"))
axis(2, at = c(10,80,640,5120))

# 160K change by vaccination history
plot(jitter(titers$Vx_hx),titers$GM_FC_160K,xlab="Vaccination history",ylab="Fold change",main="K160 HA",ylim=c(0.4,81),axes=FALSE,log="y",col=alpha('red',0.5),pch=16)
axis(1,at=c(0,1,2),labels=c("None","One","Both"))
axis(2,at=c(0.5,2.0,10.0,50.0))

# WT pre-vacc by vaccination history
plot(jitter(titers$Vx_hx),titers$GMT_Pre_WT,xlab="Vaccination history",ylab="Pre-vaccination titer",main="T160 HA", log="y",ylim=c(7,10245),axes=FALSE,col=alpha('blue',0.5),pch=16)
axis(1,at=c(0,1,2),labels=c("None","One","Both"))
axis(2, at = c(10,80,640,5120))

# WT post-vacc by vaccination history
plot(jitter(titers$Vx_hx),titers$GMT_Post_WT,xlab="Vaccination history",ylab="Post-vaccination titer",main="T160 HA",log="y",ylim=c(7,10245),axes=FALSE,col=alpha('blue',0.5),pch=16)
axis(1,at=c(0,1,2),labels=c("None","One","Both"))
axis(2, at = c(10,80,640,5120))

# WT change by vaccination history
plot(jitter(titers$Vx_hx),titers$GM_FC_WT,xlab="Vaccination history",ylab="Fold change",main="T160 HA", log="y",ylim=c(0.4,81),axes=FALSE,col=alpha('blue',0.5),pch=16)
axis(1,at=c(0,1,2),labels=c("None","One","Both"))
axis(2,at=c(0.5,2.0,10.0,50.0))
dev.off()  

############################################
########## Regression models ###############
############################################

############################################
## Prevacc titers to K160

fit <- lm( log(GMT_Pre_160K,2) ~ Age + factor(Vx_hx), data=titers)
summary(fit); logLik(fit); extractAIC(fit)

# Call:
# lm(formula = log(GMT_Pre_160K, 2) ~ Age + factor(Vx_hx), data = titers)

# Residuals:
    # Min      1Q  Median      3Q     Max 
# -3.9844 -1.0637  0.2044  1.1688  3.2474 

# Coefficients:
               # Estimate Std. Error t value Pr(>|t|)    
# (Intercept)     8.10997    0.88653   9.148 7.57e-13 ***
# Age            -0.07295    0.02448  -2.980   0.0042 ** 
# factor(Vx_hx)1  0.69526    0.71296   0.975   0.3335    
# factor(Vx_hx)2  0.65545    0.64459   1.017   0.3134    
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Residual standard error: 1.726 on 58 degrees of freedom
  # (8 observations deleted due to missingness)
# Multiple R-squared:  0.1424,	Adjusted R-squared:  0.09801 
# F-statistic: 3.209 on 3 and 58 DF,  p-value: 0.02952

# 'log Lik.' -119.739 (df=5)
# [1]  4.00000 71.52954


############################################
## Prevacc titers to WT

fit <- lm( log(GMT_Pre_WT,2) ~ Age + factor(Vx_hx), data=titers)
summary(fit); logLik(fit); extractAIC(fit)
# Call:
# lm(formula = log(GMT_Pre_WT, 2) ~ Age + factor(Vx_hx), data = titers)

# Residuals:
     # Min       1Q   Median       3Q      Max 
# -1.92020 -1.08809 -0.07095  0.83215  2.53211 

# Coefficients:
               # Estimate Std. Error t value Pr(>|t|)    
# (Intercept)     4.96254    0.61633   8.052 5.02e-11 ***
# Age            -0.02056    0.01702  -1.208    0.232    
# factor(Vx_hx)1  0.69079    0.49567   1.394    0.169    
# factor(Vx_hx)2  0.17487    0.44813   0.390    0.698    
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Residual standard error: 1.2 on 58 degrees of freedom
  # (8 observations deleted due to missingness)
# Multiple R-squared:  0.06539,	Adjusted R-squared:  0.01705 
# F-statistic: 1.353 on 3 and 58 DF,  p-value: 0.2663

# 'log Lik.' -97.20013 (df=5)
# [1]  4.00000 26.45189


############################################
## Postvacc titers to K160

fit <- lm( log(GMT_Post_160K,2) ~ Age + factor(Vx_hx) + Group + log(GMT_Pre_160K,2), data=titers)
summary(fit); logLik(fit); extractAIC(fit)

# Call:
# lm(formula = log(GMT_Post_160K, 2) ~ Age + factor(Vx_hx) + Group + 
    # log(GMT_Pre_160K, 2), data = titers)

# Residuals:
     # Min       1Q   Median       3Q      Max 
# -2.09618 -0.78038 -0.08893  0.77508  2.39687 

# Coefficients:
                     # Estimate Std. Error t value Pr(>|t|)    
# (Intercept)           7.74734    1.00506   7.708 2.61e-10 ***
# Age                  -0.03798    0.01804  -2.105 0.039838 *  
# factor(Vx_hx)1       -1.38345    0.48948  -2.826 0.006551 ** 
# factor(Vx_hx)2       -1.85282    0.45119  -4.107 0.000134 ***
# GroupFCV             -0.36871    0.36462  -1.011 0.316334    
# GroupFZ              -0.08527    0.38386  -0.222 0.825020    
# log(GMT_Pre_160K, 2)  0.46420    0.09288   4.998 6.24e-06 ***
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Residual standard error: 1.172 on 55 degrees of freedom
  # (8 observations deleted due to missingness)
# Multiple R-squared:  0.5352,	Adjusted R-squared:  0.4844 
# F-statistic: 10.55 on 6 and 55 DF,  p-value: 9.078e-08

# 'log Lik.' -94.12552 (df=8)
# [1]  7.00000 26.30267


############################################
## Postvacc titers to WT

fit <- lm( log(GMT_Post_WT,2) ~ Age + factor(Vx_hx) + Group + log(GMT_Pre_WT,2), data=titers)
summary(fit); logLik(fit); extractAIC(fit)

# Call:
# lm(formula = log(GMT_Post_WT, 2) ~ Age + factor(Vx_hx) + Group + 
    # log(GMT_Pre_WT, 2), data = titers)

# Residuals:
    # Min      1Q  Median      3Q     Max 
# -3.2809 -0.7337 -0.2436  0.8351  2.7716 

# Coefficients:
                   # Estimate Std. Error t value Pr(>|t|)    
# (Intercept)         4.94816    1.02500   4.827 1.14e-05 ***
# Age                -0.01837    0.01897  -0.968  0.33719    
# factor(Vx_hx)1     -1.66731    0.55470  -3.006  0.00399 ** 
# factor(Vx_hx)2     -2.20743    0.50569  -4.365 5.65e-05 ***
# GroupFCV           -1.04215    0.40471  -2.575  0.01274 *  
# GroupFZ            -0.89441    0.43076  -2.076  0.04255 *  
# log(GMT_Pre_WT, 2)  0.84208    0.14638   5.753 4.02e-07 ***
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Residual standard error: 1.318 on 55 degrees of freedom
  # (8 observations deleted due to missingness)
# Multiple R-squared:  0.5547,	Adjusted R-squared:  0.5062 
# F-statistic: 11.42 on 6 and 55 DF,  p-value: 2.974e-08

# 'log Lik.' -101.3966 (df=8)
# [1]  7.00000 40.84475
