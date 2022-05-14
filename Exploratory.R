library(readxl)
library(writexl)
library(ggplot2)
library(gridExtra)
library(survminer)
library(magrittr)
library(dplyr)
library(survival)
directory <- "Data/positives VL data 9nov21_mjs.xlsx"

## Variables of interest
dat_var <- c("study_id","day","vl","ct","non_study_sample","viral_culture_pos")
demo_var <- c("Study ID","Vaccination Status","Vaccine Type","Days Since First Vaccine",
              "Days Since Second Vaccine","Variant Lineage", "Time from Symptoms to Positive Test",
              "Symptomatic","Participant Type")

## Read in data + demographics
dat <- read_excel(directory, sheet = 1) %>% select(all_of(dat_var))
demo <- read_excel(directory, sheet = 4) %>% select(all_of(demo_var))

## Data wrangling
demo %<>% rename(variant=`Variant Lineage`, days_fv=`Days Since Second Vaccine`, id=`Study ID`,
                 v_status=`Vaccination Status`, v_type=`Vaccine Type`,
                 time_symp = `Time from Symptoms to Positive Test`)
demo %<>% mutate(variant = ifelse(id <= 126, "Non-delta", "Delta"))
demo %<>% mutate(days_fv = ifelse(demo$days_fv == "N/A", `Days Since First Vaccine`, days_fv))
demo %<>% mutate(days_fv = as.double(days_fv), time_symp = as.double(time_symp)) %>% mutate_if(is.character,as.factor)
demo %<>% mutate(month_fv = days_fv/30, month_cat = ifelse(month_fv >= 3, ">3 months", "<3 months"))
dat %<>% rename(id=study_id)
dat %<>% mutate(ct=ifelse(ct==".",0,ct),ct=as.numeric(ct))

## Create merged dataset
fullDat <- dat %>% left_join(demo)

## Remove problem entries
# Unvaccinated
fullDat %<>% filter(v_status=="Y")
# Transplant ID 101 + false positive ID 109
fullDat %<>% filter(!id %in% c(101,109))
# Remove non-study samples
fullDat %<>% mutate(non_study_sample = ifelse(is.na(non_study_sample),0,1))
fullDat %<>% filter(!(non_study_sample==1 & day!=0))
# Remove inpatients
fullDat %<>% filter(`Participant Type`=="Outpatient")

# Sub-analysis of only symptomatic individuals
#fullDat %<>% filter(Symptomatic=="Y")
# Update demographics set
sDemo <- demo %>% filter(id %in% unique(fullDat$id))





## Calculate time to event outcomes
# Time to 0 VL
unique_id <- sDemo$id
vlDat <- tibble(id=unique_id, day=0, status=1)
for (i in unique_id){
   # Current longitudinal observation
   curr <- fullDat %>% filter(id==i)
   
   # Index for vlDat
   ind <- vlDat$id==i
   
   for (k in 2:nrow(curr)){
      if (curr[k,]$vl == 0){
         vlDat[ind,]$day <- curr[k,]$day
         break
      } else if (k==nrow(curr)){
         vlDat[ind,]$day <- curr[k,]$day
         vlDat[ind,]$status <- 0
      }
   }
}
vlDat %<>% left_join(sDemo)
vlDat %<>% within(variant <- relevel(variant, ref = "Non-delta"))
vlDat %<>% mutate(day = ifelse(is.na(time_symp), day, day+time_symp))

# Time to negative viral culture
vcDat <- tibble(id=unique_id, day=0, status=1)
for (i in unique_id){
   curr <- fullDat %>% filter(id==i)
   ind <- vcDat$id==i
   seenPos <- 0
   
   # If all 0 cultures, set to first date
   if (sum(curr$viral_culture_pos[2:nrow(curr)])==0){
      vcDat[ind,]$day <- curr$day[2]
   }
   
   for (k in 2:nrow(curr)){
      # Look for first positive value
      if (curr[k,]$viral_culture_pos>0 && seenPos==0){
         seenPos <- 1
      } else if (seenPos==1 && curr[k,]$viral_culture_pos==0){
         vcDat[ind,]$day <- curr[k,]$day
         break
      } else if (k==nrow(curr) && seenPos==1){
         vcDat[ind,]$day <- curr[k,]$day
         vcDat[ind,]$status <- 0
      }
   }
}
vcDat %<>% left_join(sDemo)
vcDat %<>% within(variant <- relevel(variant, ref = "Non-delta"))
vcDat %<>% mutate(day = ifelse(is.na(time_symp), day, day+time_symp))

# Time to >30 ct
ctDat <- tibble(id=unique_id, day=0, status=1)
for (i in unique_id){
   curr <- fullDat %>% filter(id==i)
   ind <- ctDat$id==i
   
   for (k in 1:nrow(curr)){
      if (curr[k,]$ct >= 30){
         ctDat[ind,]$day <- curr[k,]$day
         break
      } else if (k==nrow(curr)){
         ctDat[ind,]$day <- curr[k,]$day
         ctDat[ind,]$status <- 0
      }
   }
}
ctDat %<>% left_join(sDemo)
ctDat %<>% within(variant <- relevel(variant, ref = "Non-delta"))
ctDat %<>% mutate(day = ifelse(is.na(time_symp), day, day+time_symp)) %>% mutate(day = ifelse(day<0,0,day))




## Kaplan Meier plots
confintFlag <- TRUE
title_size <- 33
font_size <- 30
axis_size <- 28
rt_size <- 10


# Negative VL - variant
fit1 <- survfit(Surv(day, status) ~ variant, data = vlDat)
fit1.p <- ggsurvplot(fit1, data=vlDat, legend.labs =
           c("Non-delta", "Delta"), title="Days to Confirmed Undetectable Viral Load", 
           xlab="Days",
           ylab="Prob. of Detectable SARS-CoV-2 Viral Load",
           risk.table=TRUE,
           font.title=c(title_size,"bold"),
           font.x = c(font_size),
           font.y = c(font_size),
           font.legend = c(font_size),
           font.tickslab = c(axis_size),
           conf.int = confintFlag,
           fontsize=rt_size,
           break.x.by=5,
           xlim=c(0,20))

fit1.p$table$theme$axis.text.y$size <- font_size
fit1.p$table <- fit1.p$table +
   theme(plot.title = element_text(size = font_size),
         axis.title.x = element_text(size = font_size),
         axis.title.y = element_text(size = font_size),
         axis.text.x = element_text(size = font_size),
         axis.ticks.length=unit(.25, "cm"))

fit1.p$plot <- fit1.p$plot + geom_hline(yintercept = 0.5, linetype="dashed") + theme(axis.ticks.length=unit(.25, "cm"))
png(file="Plots/nvl_variant.png",width=3600, height=3780, res=300)
fit1.p
dev.off()


# Negative culture - variant
fit2 <- survfit(Surv(day, status) ~ variant, data = vcDat, conf.type="log-log")
fit2.p <- ggsurvplot(fit2, data=vcDat, legend.labs =
           c("Non-delta", "Delta"), title="Days to Confirmed Negative Viral Culture",
           xlab="Days",
           ylab="Prob. of Positive SARS-CoV-2 Viral Culture",
           risk.table=TRUE,
           font.title=c(title_size,"bold"),
           font.x = c(font_size),
           font.y = c(font_size),
           font.legend = c(font_size),
           font.tickslab = c(axis_size),
           conf.int = confintFlag,
           fontsize=rt_size,
           break.x.by=5,
           xlim=c(0,20))

fit2.p$table$theme$axis.text.y$size <- font_size
fit2.p$table <- fit2.p$table +
   theme(plot.title = element_text(size = font_size),
         axis.title.x = element_text(size = font_size),
         axis.title.y = element_text(size = font_size),
         axis.text.x = element_text(size = font_size),
         axis.ticks.length=unit(.25, "cm"))


fit2.p$plot <- fit2.p$plot + geom_hline(yintercept = 0.5, linetype="dashed") + theme(axis.ticks.length=unit(.25, "cm"))
png(file="Plots/nvc_variant.png",width=3600, height=3780, res=300)
fit2.p
dev.off()


# CT >30 - variant
fit3 <- survfit(Surv(day, status) ~ variant, data = ctDat)
fit3.p <- ggsurvplot(fit3, data=ctDat, legend.labs =
              c("Non-delta", "Delta"), title="Days to Confirmed CT value >30",
           xlab="Days",
           ylab="Prob. of SARS-CoV-2 CT <30",
           risk.table=TRUE,
           font.title=c(title_size,"bold"),
           font.x = c(font_size),
           font.y = c(font_size),
           font.legend = c(font_size),
           font.tickslab = c(axis_size),
           conf.int = confintFlag,
           fontsize=rt_size,
           break.x.by=5,
           xlim=c(0,20))

fit3.p$table$theme$axis.text.y$size <- font_size
fit3.p$table <- fit3.p$table +
   theme(plot.title = element_text(size = font_size),
         axis.title.x = element_text(size = font_size),
         axis.title.y = element_text(size = font_size),
         axis.text.x = element_text(size = font_size),
         axis.ticks.length=unit(.25, "cm"))

fit3.p$plot <- fit3.p$plot + geom_hline(yintercept = 0.5, linetype="dashed") + theme(axis.ticks.length=unit(.25, "cm"))
png(file="Plots/ct_variant.png",width=3600, height=3780, res=300)
fit3.p
dev.off()


# Negative VL - by month category
fit4 <- survfit(Surv(day, status) ~ month_cat, data = vlDat)
fit4.p <- ggsurvplot(fit4, data=vlDat, legend.labs =
                        c("Months <3", "Months >3"), title="Days to Confirmed Undetectable Viral Load", 
                     xlab="Days",
                     ylab="Prob. of Detectable SARS-CoV-2 Viral Load",
                     risk.table=TRUE,
                     font.title=c(title_size,"bold"),
                     font.x = c(font_size),
                     font.y = c(font_size),
                     font.legend = c(font_size),palette=c("purple","orange"),
                     font.tickslab = c(axis_size),
                     conf.int=confintFlag,
                     fontsize=rt_size,
                     break.x.by=5,
                     xlim=c(0,20))

fit4.p$table$theme$axis.text.y$size <- font_size
fit4.p$table <- fit4.p$table +
   theme(plot.title = element_text(size = font_size),
         axis.title.x = element_text(size = font_size),
         axis.title.y = element_text(size = font_size),
         axis.text.x = element_text(size = font_size),
         axis.ticks.length=unit(.25, "cm"))

fit4.p$plot <- fit4.p$plot + geom_hline(yintercept = 0.5, linetype="dashed") + theme(axis.ticks.length=unit(.25, "cm"))
png(file="Plots/nvl_time.png",width=3600, height=3780, res=300)
fit4.p
dev.off()

fit5 <- survfit(Surv(day, status) ~ month_cat, data = vcDat, conf.type="log-log")
fit5.p <- ggsurvplot(fit5, data=vcDat, legend.labs =
                        c("Months <3", "Months >3"), title="Days to Confirmed Negative Viral Culture",
                     xlab="Days",
                     ylab="Prob. of Positive SARS-CoV-2 Viral Culture",
                     risk.table=TRUE,
                     font.title=c(title_size,"bold"),
                     font.x = c(font_size),
                     font.y = c(font_size),
                     font.legend = c(font_size),palette=c("purple","orange"),
                     font.tickslab = c(axis_size),
                     conf.int=confintFlag,
                     fontsize=rt_size,
                     break.x.by=5,
                     xlim=c(0,20))

fit5.p$table$theme$axis.text.y$size <- font_size
fit5.p$table <- fit5.p$table +
   theme(plot.title = element_text(size = font_size),
         axis.title.x = element_text(size = font_size),
         axis.title.y = element_text(size = font_size),
         axis.text.x = element_text(size = font_size),
         axis.ticks.length=unit(.25, "cm"))

fit5.p$plot  <- fit5.p$plot + geom_hline(yintercept = 0.5, linetype="dashed") + theme(axis.ticks.length=unit(.25, "cm"))
png(file="Plots/nvc_time.png",width=3600, height=3780, res=300)
fit5.p
dev.off()





# Cox models
getHRvals <- function(fit){
   return(summary(fit)[[8]][c(1,3,4)])
}

cox_vl <- coxph(Surv(day,status) ~ variant, data=vlDat)
cox_vc <- coxph(Surv(day,status) ~ variant, data=vcDat)
cox_ct <- coxph(Surv(day,status) ~ variant, data=ctDat)
cox_vl2 <- coxph(Surv(day, status) ~ month_cat, data = vlDat)
cox_vc2 <- coxph(Surv(day, status) ~ month_cat, data = vcDat)

tibble(Label=c("HR","95%L","95%U"), 
                 VL=getHRvals(cox_vl), VC=getHRvals(cox_vc), CT=getHRvals(cox_ct),
       VL_time=getHRvals(cox_vl2), VC_time=getHRvals(cox_vc2))

# Median survival times
survfit(Surv(day,status) ~ variant, data=vlDat)
survfit(Surv(day,status) ~ variant, data=vcDat)
survfit(Surv(day,status) ~ variant, data=ctDat)
survfit(Surv(day,status) ~ month_cat, data=vlDat)
survfit(Surv(day,status) ~ month_cat, data=vcDat)





## Viral decay trajectory
# Non-delta
nz_nd <- fullDat %>% filter(variant=="Non-delta") %>% 
   group_by(id) %>% filter(vl[2]>0) %>% ungroup %>% select(id, day, vl, variant)
z_nd <- fullDat %>% anti_join(nz_nd, by="id") %>% filter(variant=="Non-delta")
nz_nd %<>% mutate(id=as.factor(id))
z_nd %<>% filter(day!=0) %>% mutate(id=as.factor(id)) %>% select(id,day,vl,variant)
z_nd$vl2 <- z_nd %>% group_by(id) %>% group_indices() * -0.1-0.2 

# Delta
nz_d <- fullDat %>% filter(variant=="Delta") %>% 
   group_by(id) %>% filter(vl[2]>0) %>% ungroup %>% select(id, day, vl, variant)
nz_d %<>% mutate(id=as.factor(id))

# Predictive plot
dp <- fullDat %>% filter(variant=="Delta") %>% select(id,day,vl,variant) %>% 
   mutate(day2=day^2, day3=day^3)
dp_fit <- lm(vl ~ day + day2 + day3, data=dp)
dp_temp <- tibble(day=0:20) %>% mutate(day2=day^2,day3=day^3)
dp_vals <- predict(dp_fit,dp_temp)
dp_plot <- tibble(xval=dp_temp$day, yval=dp_vals,id=0)

ndp <- fullDat %>% filter(variant=="Non-delta") %>% select(id,day,vl,variant) %>% 
   mutate(day2=day^2, day3=day^3)
ndp_fit <- lm(vl ~ day + day2, data=ndp)
ndp_temp <- tibble(day=0:20) %>% mutate(day2=day^2)
ndp_vals <- predict(ndp_fit,ndp_temp)
ndp_plot <- tibble(xval=ndp_temp$day, yval=ndp_vals,id=0)

allZero <- fullDat %>% filter(variant=="Non-delta") %>% group_by(id) %>% slice(2) %>% filter(vl<=3) %>% pull(id)
ndp_cond <- ndp %>% filter(!(id %in% allZero))
ndp_cond_fit <- lm(vl ~ day + day2, data=ndp_cond)
ndp_cond_vals <- predict(ndp_cond_fit,ndp_temp)
ndp_cond_plot <- tibble(xval=ndp_temp$day, yval=ndp_cond_vals,id=0)

# Trajectory plot code
p1 <- ggplot(aes(x=day,y=vl,group=id), data = nz_nd) + geom_line(na.rm = TRUE) + 
   geom_point(color="#F8766d", size=3, na.rm = TRUE) + 
   ylim(-1,8) + ylab("Viral load (log RNA copies/mL)") +
   xlab("Days Since Initial Positive PCR") +
   ggtitle("Trajectory of Non-Delta Viral Decay") +
   geom_line(data = z_nd, aes(x=day,y=vl2,group=id), na.rm = TRUE) +
   geom_point(data = z_nd, aes(x=day,y=vl2,group=id), size=2,
              color="#F8766d",na.rm = TRUE)

p2 <- ggplot(aes(x=day,y=vl,group=id), data = nz_d) + geom_line(na.rm = TRUE) + 
   geom_point(color="#00BFC4", shape=17, size=3, na.rm = TRUE) + 
   ylim(-1,8) + ylab("Viral load (log RNA copies/mL)") +
   xlab("Days Since Initial Positive PCR") +
   ggtitle("Trajectory of Delta Viral Decay")

p3 <- p2 + geom_line(data=dp_plot, aes(x=xval,y=yval), linetype="dashed", 
                     na.rm=TRUE, color="#00BFC4", size=2) +
   theme(legend.position = "none", axis.title = element_text(size=font_size),
         plot.title = element_text(size=title_size, face="bold"),
         axis.text=element_text(size=axis_size),
         axis.title.y=element_blank(),
         axis.ticks.length=unit(.25, "cm"))

p4 <- p1 + geom_line(data=ndp_plot, aes(x=xval,y=yval), linetype="dashed", 
                     na.rm=TRUE, color="#F8766d", size=2) + 
   theme(legend.position = "none", axis.title = element_text(size=font_size), 
         plot.title = element_text(size=title_size, face="bold"),
         axis.text = element_text(size=axis_size),
         axis.ticks.length=unit(.25, "cm")) + 
   geom_line(data=ndp_cond_plot, aes(x=xval,y=yval), linetype="dotted", 
             na.rm=TRUE, color="#F8766d", size=2)

png(file="Plots/vl_trajectory.png",width=7200, height=3780, res=300)
grid.arrange(p4,p3,ncol=2)
dev.off()





# CT trajectory
# Filter out those with only a single CT value
comb_bucket <- fullDat %>% 
  mutate(ct=ifelse(ct==0,NA,ct)) %>%
  filter(!is.na(ct)) %>%
  count(id) %>%
  filter(n>1) %>% 
  pull(id)

ct_bucket <- fullDat %>% filter(id %in% comb_bucket) %>% 
   mutate(ct=ifelse(ct==0,NA,ct)) %>%
   mutate(id=as.factor(id))

# CT trajectory plot code
shapes <- c(rep(16,4),rep(17,10))
colors <- c(rep("#F8766d",4),rep("#00BFC4",10))
p5 <- ggplot(aes(x=day,y=ct,group=id), data = ct_bucket) + geom_line(na.rm = TRUE) + 
   geom_point(aes(color=id,shape=id), size=3, na.rm = TRUE) + 
   ylab("CT value") +
   scale_color_manual(values = colors) +
   scale_shape_manual(values=shapes) +
   xlab("Days Since Initial Positive PCR") +
   ggtitle("Trajectory of CT values") +
   scale_y_reverse() + 
   theme(legend.position = "none", axis.title = element_text(size=font_size), 
         plot.title = element_text(size=title_size, face="bold"),
         axis.text=element_text(size=axis_size),
         axis.ticks.length=unit(.25, "cm")) + xlim(0,20)

png(file="Plots/ct_trajectory.png",width=3600, height=3780, res=300)
p5
dev.off()




# Export dataset
desc_text <- c("Study ID", "Days since initial postive PCR",
               "Viral load", "CT value",
               "Flag for non study sample",
               "Indicator for positive viral culture",
               "Vaccination status",
               "Vaccine type",
               "Days since first vaccine",
               "Days since second vaccine",
               "COVID-19 Variant",
               "Time from symptoms to positive test",
               "Flag for symptoms",
               "Participant type in/outpatient",
               "Months since final vaccine",
               "Month category")
description <- tibble(Variables = colnames(fullDat),
                      Description = desc_text)
write_xlsx(list(Data=fullDat, Description=description), "Plots/exported_data.xlsx")

