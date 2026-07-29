

#rm(list = ls())

options(scipen = 999) 

project_dir <- "C:/Users/ulric/Documents/zcurve3"

source(file.path(project_dir, "zcurve3S_package.R"))
source(file.path(project_dir, "zcurve_plot.R"))

#Read in data, making "" = NA

dat <- read.csv("https://raw.githubusercontent.com/UlrichSchimmack/zcurve3.0/refs/heads/main/data.csv")

dat$Title[415] = "Why does the Sinner Act Prosocially"

colnames(dat)
table(dat$effectsize)

bad_titles <- which(
  is.na(iconv(dat$Title, from = "UTF-8", to = "UTF-8"))
)
dat$Title[bad_titles] = paste0("Title",bad_titles)

write.csv(dat,"temp.dat")
dat = read.csv("temp.dat")

dat$Title[bad_titles]

dat$se = sqrt(dat$mn)

dat$cluster.id = as.character(dat$uniquestudy)
tab = table(dat$cluster.id)
length(tab)

dat$row_id = 1:nrow(dat)
dim(dat)

colnames(dat)
table(!is.na(dat$p.control),!is.na(dat$p.prime))

dat = dat[is.na(dat$p.control),]

dim(dat)

dat$z = dat$d / dat$se
table(dat$z > 4)



###

# Fitzsimons & Bargh (2003), Study 1

# 5 extreme z > 4 results
# results could not be matched to article, 

suspect = dat$Title == dat$Title[dat$row_id == 240]
table(suspect)

dat = dat[!suspect,]
dim(dat)

table(dat$z > 4)

table(!is.na(dat$z))

table(dat$z > 4)
mean(dat$z > 4)

table(cut(dat$z,c(-10,-1.96,0,1.65,1.96,4,6,100)))


###


boot_iter          <- 500
est_method         <- "EM"
ncp                <- seq(0,6,.5)
z_sd               <- rep(1,length(ncp))
ncp_fixed          <- TRUE
z_sd_fixed         <- TRUE
directional        <- FALSE
folded             <- TRUE
bootstrap_parallel <- TRUE
show_plots         <- TRUE
int_beg            <- 2.10

### Figure 1b

z_res <- zcurve(
  zval = abs(dat$z),
  yi = abs(dat$d),
  sei = dat$se,
  cluster_id = dat$cluster.id,
  est_method = est_method,
  directional = directional,
  folded = folded,
  show_plot = show_plots,
  show_summary = FALSE,
  ncp = ncp,
  z_sd = z_sd,
  ncp_fixed = ncp_fixed,  
  z_sd_fixed = z_sd_fixed,
  int_beg = int_beg,  
  boot_iter = boot_iter,
  control = zcurve_control(
    parallel = bootstrap_parallel,
    cores = 16
  )
)

z_res_study = z_res$individual_effects
nrow(z_res_study)

tab = table(z_res_study$effect_adjusted_lb > .2)
tab
tab / nrow(dat)

tab = table(dat$z > 4, z_res_study$effect_adjusted_lb > .2)
tab
tab / nrow(dat)


### Figure 1a

z_res_sorted = z_res_study[order(z_res_study$effect_adjusted_lb,
  decreasing=TRUE),]

z_res_sorted$sort = 1:nrow(z_res_sorted)

plot(z_res_sorted$sort,z_res_sorted$effect_adjusted_lb,ylim=c(0,1),
  ylab="Effect Size (SMD)",xlab = "Effect sorted by lower bound",pch=18,cex=1.5)
lines(z_res_sorted$sort,z_res_sorted$effect_adjusted_lb,ylim=c(0,1),
  ylab="",xlab = "",type="l",lwd=3)
abline(h=c(.2,.5,.8),lty=2)


### Figure 1c

dat$study_id = 1:nrow(dat)

z.forest = zcurve_forest(z_res, 
   text_size = 1.5,
   es_lb_min = .20,
   k_max = 50,
   z_min = 4.0,
   sort_by = "lower_bound",
   labels = setNames(dat$Title, dat$study_id),
   show_id = TRUE
   )

dim(z.forest)

max = nchar("Money Cues Increase Agency and Decrease Prosociality ")
substring(dat$Title[z.forest$study_id],1,max)


###############################


library(ashr)

betahat = abs(dat$d) 
sebetahat = dat$se 

K = (length(ncp)-1)
ncp2 = seq(-1,1,1/K)
ncp2
K2 = length(ncp2)

ash_res = ash(   
   betahat = betahat,
   sebetahat = sebetahat,
   g = normalmix(
       pi          = rep(1/K2, K2),
       mean        = ncp2,
       sd          = rep(.01, K2)),  # ~point masses
       fixg        = FALSE,            # estimate the weights
       prior       = "uniform",        # <- penalty OFF, so it matches unpenalized ML
       mixcompdist = "normal",
       optmethod   = "mixEM",    
       control = list(tol = 1e-8,maxiter = 5000, trace = FALSE)
    )

   pm <- get_pm(ash_res)
   length(pm)


   ash_pe <- get_pm(ash_res)
   ci_ash <- ashci(
     ash_res,
     level = 0.95
   )

   dat$ash_lb <- ci_ash[, 1]
   dat$ash_ub <- ci_ash[, 2]


   plot(z_res_study$effect_adjusted_lb,ash_lb,xlim=c(0,2),ylim=c(0,2))
   abline(a = 0,b = 1)

   z_res_study$study_id
   z_res_study$effect_adjusted
   z_res_study$Title = dat$Title
   z_res_study$ash_lb = dat$ash_lb

   z_res_study$obs_lb = abs(dat$d) - 2*dat$se

   table(z_res_study$ash_lb > .33,z_res_study$effect_adjusted_lb > .2)
   table(z_res_study$obs_lb > .5,z_res_study$effect_adjusted_lb > .2)
   table(z_res_study$obs_lb > .5,z_res_study$ash_lb > .33)


   z_res_study[order(dat$ash_lb,decreasing=TRUE),c("study_id","ash_lb","effect_adjusted_lb","Title")][1:31,]
   
   

   plot(dat$d,z_res_study$effect_adjusted,xlim=c(0,2),ylim=c(0,2))
   par(new=TRUE)
   plot(dat$d,pm,xlim=c(0,2),ylim=c(0,2),col="blue")
   abline(a = 0,b = 1)



