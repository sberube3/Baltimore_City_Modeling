library(R0)
library(dplyr)
library(ggplot2)
library(lubridate)
library(cowplot)
library(grid)
library(gridExtra) 
library(reshape2)

bc <- read.csv("~/Dropbox (UFL)/_active/COVID19/Baltimore_City_Modeling/balt_confirmed_cases_20200514.csv")
bc <- bc[which(bc$Date>=min(bc$Date[which(bc$Date>=20200315&bc$New.Cases!=0)])),]
bc_incid <- setNames(bc$New.Cases,1:nrow(bc))
bc_pop <- 620961
  # 620961 Baltimore city
  # 6045680 Maryland

md <- read.csv("~/Dropbox (UFL)/_active/COVID19/Baltimore_City_Modeling/md_confirmed_cases_20200514.csv")
md <- md[which(md$Date>=min(md$Date[which(md$Date>=20200315&md$New.Cases!=0)])),]
md_incid <- setNames(md$New.Cases,1:nrow(md))
md_pop <- 6045680

ga <- read.csv("~/Dropbox (UFL)/_active/COVID19/Baltimore_City_Modeling/ga_confirmed_cases_20200514.csv")
ga <- ga[which(ga$Date>=min(ga$Date[which(ga$Date>=20200315&ga$New.Cases!=0)])),]
ga_incid <- setNames(ga$New.Cases,1:nrow(ga))
ga_pop <- 10617423

mGT <- generation.time("gamma",c(6.5,4.5))

estBC <- estimate.R(bc_incid, mGT, begin=1, end=as.numeric(length(bc_incid)), 
                    methods=c("TD"), pop.size=bc_pop, nsim=1000)
estMD <- estimate.R(md_incid, mGT, begin=1, end=as.numeric(length(md_incid)), 
                    methods=c("TD"), pop.size=md_pop, nsim=1000)
estGA <- estimate.R(ga_incid, mGT, begin=1, end=as.numeric(length(ga_incid)), 
                    methods=c("TD"), pop.size=ga_pop, nsim=1000)

rtBC <- cbind(bc$Date,estBC$estimates$TD$R,estBC$estimates$TD$conf.int)
colnames(rtBC) <- c("date","estimate","lower","upper")
rtBC$lty <- ifelse(rtBC$date>=(max(rtBC$date)-7), "B", "A")
rtBC <- rtBC[which(rtBC$estimate!=0),]
rtBC$loc <- "Baltimore"
rtBC$date <- as.Date(as.character(rtBC$date), format="%Y%m%d")

rtMD <- cbind(md$Date,estMD$estimates$TD$R,estMD$estimates$TD$conf.int)
colnames(rtMD) <- c("date","estimate","lower","upper")
rtMD$lty <- ifelse(rtMD$date>=(max(rtMD$date)-7), "B", "A")
rtMD <- rtMD[which(rtMD$estimate!=0),]
rtMD$loc <- "Maryland"
rtMD$date <- as.Date(as.character(rtMD$date), format="%Y%m%d")

rtGA <- cbind(ga$Date,estGA$estimates$TD$R,estGA$estimates$TD$conf.int)
colnames(rtGA) <- c("date","estimate","lower","upper")
rtGA$lty <- ifelse(rtGA$date>=(max(rtGA$date)-7), "B", "A")
rtGA <- rtGA[which(rtGA$estimate!=0),]
rtGA$loc <- "Georgia"
rtGA$date <- as.Date(as.character(rtGA$date), format="%Y%m%d")

dat <- rbind(rtBC, rtMD, rtGA)
dat <- rbind(rtBC, rtMD)

ggplot(dat %>% filter(date>=as.Date("2020-03-20")), aes(x=date, y=estimate, ymin=lower, ymax=upper, color=loc, group=interaction(loc,lty), fill=loc)) +
  geom_line(aes(linetype=lty)) +
  geom_ribbon(alpha=0.25) +
  xlab("Date") +
  ylab("Effective Reproductive Number") +
  ylim(-0.5,3.5) +
  geom_hline(yintercept=1,size=0.63,color="gray28", linetype="dashed") +
  theme_minimal()

summary(rtBC %>% filter(date > as.Date("2020-03-20") & date <= as.Date("2020-03-29")))
summary(rtBC %>% filter(date > as.Date("2020-03-29") & date <= as.Date("2020-04-19")))
summary(rtBC %>% filter(date > as.Date("2020-04-19") & date <= as.Date("2020-04-26")))
summary(rtBC %>% filter(date > as.Date("2020-04-26") & date <= as.Date("2020-05-04")))

summary(rtGA %>% filter(date > as.Date("2020-04-30") & date <= as.Date("2020-05-10")))

