
rm(list = ls())

options(scipen = 999)

install.packages("remotes")
remotes::install_github("UlrichSchimmack/zcurve3")

library(zcurve3)

library(dplyr)
library(osfr)


data2=read.csv("https://osf.io/xq4b2/?action=download")

osf_retrieve_file("cne28") %>% osf_download(conflicts="skip")          # download file CDSR.Rdata
load("CDSR.Rdata")                                             # load data.frame "data"
dim(data)
dim(data2)

dat = data.frame(data,data2)
dim(dat)

colnames(dat)

sel = dat$RCT == "yes"
table(sel)

dat = dat[sel,]

dat$outliers = 
   (dat$outcome.flag == "DICH" & dat$effect.es < .05) |
   (dat$outcome.flag == "DICH" & dat$effect.es > 20) |
   (dat$outcome.flag == "CONT" & dat$effect.es < -1) |
   (dat$outcome.flag == "CONT" & dat$effect.es > 1) 
table(dat$outliers)

dat$outliers[dat$effect.N < 10] = TRUE
dat$outliers[dat$effect.N > 10000] = TRUE

dat$sei = sqrt(1/dat$total1) + sqrt(1/dat$total2)

dat$outliers[dat$sei < .05] = TRUE
dat$outliers[dat$sei > .60] = TRUE

plot(dat$effect.t,dat$z)

dat$yi = dat$effect.t * dat$sei
dat$vi = dat$sei^2

head(dat$id)

dat$id2 <- paste0(dat$id, ".", dat$outcome.nr)
dat$cluster_id = as.numeric(factor(dat$id2))
tab = table(dat$cluster_id)
length(tab)

dat$outliers[abs(dat$yi) > 2] = TRUE

sel = dat$outcome.group == "efficacy"

n1 <- tapply(dat$outcome.nr, dat$id, function(x) sum(x %in% 1:5))
table(n1)

sel_k = dat$id %in% names(n1)[n1 >= 5]
table(sel,sel_k)

sel = sel & sel_k

table(sel,dat$outlier)

sel = sel & !dat$outlier


summary(dat$z[sel])
summary(dat$yi[sel])
summary(dat$sei[sel])

tab = table(dat$id)
length(tab[tab > 0])

tab = table(dat$id[sel])
length(tab[tab > 0])

tab = table(dat$id2[sel])
length(tab[tab > 0])

length(table(dat$id[sel]))


vanZwet_zero <- zcurve(
  zval = abs(dat$z)[sel],
#  yi = abs(dat$yi)[sel],
#  sei = dat$sei[sel],
  cluster_id = dat$cluster_id[sel],
  boot_iter = 500,
  ncp = c(0,0,0),
  z_sd = c(1,3,5),
  ncp_fixed = TRUE,
  z_sd_fixed = FALSE,
  control = zcurve_control(
    parallel = TRUE,
    cores = 16
  )
)

plot(vanZwet_zero)

summary(vanZwet_zero)


vanZwet_free <- zcurve(
  zval = abs(dat$z)[sel],
#  yi = abs(dat$yi)[sel],
#  sei = dat$sei[sel],
  cluster_id = dat$cluster_id[sel],
  boot_iter = 500,
  ncp = c(0,2,4),
  z_sd = c(1,3,5),
  ncp_fixed = FALSE,
  z_sd_fixed = FALSE,
  control = zcurve_control(
    parallel = TRUE,
    cores = 16
  )
)

plot(vanZwet_free)

summary(vanZwet_free)


discrete_fixed <- zcurve(
  zval = abs(dat$z)[sel],
#  yi = abs(dat$yi)[sel],
#  sei = dat$sei[sel],
  cluster_id = dat$cluster_id[sel],
  boot_iter = 500,
  control = zcurve_control(
    parallel = TRUE,
    cores = 16
  )
)

plot(discrete_fixed)

summary(discrete_fixed)


discrete_wide <- zcurve(
  zval = abs(dat$z)[sel],
#  yi = abs(dat$yi)[sel],
#  sei = dat$sei[sel],
  cluster_id = dat$cluster_id[sel],
  ncp = c(0,2,4,6),
  z_sd = c(1,1,1,1),
  boot_iter = 500,
  control = zcurve_control(
    parallel = TRUE,
    cores = 16
  )
)

summary(discrete_wide)

plot(discrete_wide)



