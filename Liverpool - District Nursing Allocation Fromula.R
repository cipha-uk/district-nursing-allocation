#############################################################################################
# THE DATA- Importing from CSV and appending
#############################################################################################

gc() #clean memory

library(dplyr)
library(tidyr)
require(ggplot2)
require(pscl)
require(boot)
library(data.table)
library(speedglm)
library(VGAM)
library(broom)
library(broom.mixed)
library(glmmTMB)
library(ggspatial)
require(pscl)
library(officer)
library(data.table)
library(sf) 
library(dplyr)

library(janitor)


library(ggthemes)
library(RColorBrewer)
library(colorspace)


library(tableone)
library(gdata)
library(gtools)


library(ggspatial)
library(flextable)
library(cowplot)

library("parameters")

#  prepare district nurse look up from FOI request
# source("./scripts/clean_nurses2.R")
rm(list=ls()) # remove everything in the environment
#  you will need to change this file path
load("./data/liverpool_anon_district.RData")
gps<-unique(CM[,.(gp_code)])

#  missing GP practices??? maybe not given consent
# ALBION SURGERY (N82095)

#  read it in
dist_lk<-fread("./data/dist_lk.csv")

# gp_ccg<-fread("./data/epcmem.csv")
# gp_ccg<-unique(gp_ccg[V2=="99A", .(gp_code=V1, ccgcode=V2)])
#  Liverpool=99A
CM[,gp_code:=as.character(gp_code)]
dist_lk[, dup:=.N, by=.(gp_code)]
CM<-merge(CM, dist_lk, by.x="gp_code", by.y="gp_code", all.x = T)
`%notin%` <- Negate(`%in%`)
check<-dist_lk[gp_code %notin% gps$gp_code]
CM<-CM[ccgcode=="99A"]
# CM<-merge(CM, gp_ccg, by="gp_code", all.x = T)

#  truncate contacts at 365
CM[contacts>365, contacts:=365]

inc<-unique(CM[is.na(dn_ppop)==F]$gp_code)

dist_lk[, notinc:=!(gp_code %in% inc)]
length(unique(CM[is.na(dn_ppop)==T]$gp_code))

notinc_liv<-unique(CM[is.na(dn_ppop)==T]$gp_code)
notinc_liv

CM[, notinc:=!(gp_code %in% notinc_liv)]

table(CM$notinc,CM$ccgcode, useNA = "ifany")
#  all in liverpool CCG merged
#  "N82052" "N82103" - both closed. 

CM[, c_imd:=(IMD-mean(IMD))/sd(IMD)]
CM[, c_age:=(age-mean(age))/sd(age)]
CM[, age_IMD:=c_age*c_imd]
CM[, agesq:=age^2]
CM[, agesq_imd:=agesq*c_imd]

CM[, age_gender:=c_age*as.numeric(gender==1)]
CM[, agesq_gender:=agesq*as.numeric(gender==1)]

CM$age_int<-NULL
CM$IMD_decile<-NULL
CM$gender_IMD<-NULL

# 
# fit_zipoisson <- glmmTMB(NCalls~(FT+ArrivalTime)*SexParent+
#                            offset(log(BroodSize))+(1|Nest),
#                          data=Owls,
#                          ziformula=~1,
#                          family=poisson)


# reg1<- glmmTMB(contacts ~ age+IMD+gender+EthnicGroup+dead+log_time_alive+NursingCareHomeFlag+Living_alone+
#                  living_with_over_64_flag+Palliative_EOL_Reg+Dementia+Cancer+Learning_Disability+Physical_Disability+
#                  CVD+CLD+Asthma+COPD+Neurological+has_mh_any+aae1a2_attend_12+Substance_or_alcohol+
#                  admissions_electives_12+admissions_emergency_12+gp_consultations_12+dn_ppop,
#                data = CM,    ziformula=~ .,
#                family=poisson)
# summary(reg1)

table(CM$EthnicGroup)
CM[, bame:=as.numeric(EthnicGroup!="White")]
CM[EthnicGroup=="Not known", bame:=NA]

CM[, bame1:=bame]
# CM[EthnicGroup=="Not known", bame:=2]
# CM[, bame:=as.factor(bame)]
prop.table(table(CM$contacts, useNA = "ifany"))


reg1<- zeroinfl(contacts ~ c_age+c_imd+gender+age_gender+age_IMD+
                  dead+bame+log_time_alive+NursingCareHomeFlag+Living_alone+
                  Palliative_EOL_Reg+Learning_Disability+
                  Dementia+Cancer+CVD+CLD+COPD+Neurological+has_mh_any+
                  aae1a2_attend_12+ admissions_electives_12+admissions_emergency_12+gp_consultations_12+
                  dn_ppop |c_age+c_imd+gender+age_gender+age_IMD+
                  dead+bame1+log_time_alive+NursingCareHomeFlag+Living_alone+
                  Palliative_EOL_Reg+Learning_Disability+
                  Dementia+Cancer+CVD+CLD+COPD+Neurological+has_mh_any+
                  aae1a2_attend_12+ admissions_electives_12+admissions_emergency_12+gp_consultations_12+
                  dn_ppop,
                data = CM)


summary(reg1)

zi_results<-as.data.table(model_parameters(reg1, component="zi"))
cond_results<-as.data.table(model_parameters(reg1, component="conditional"))


reg2<-  zeroinfl(contacts ~ c_age+c_imd+gender+age_gender+age_IMD+
                   dead+log_time_alive+NursingCareHomeFlag+Living_alone+
                   Palliative_EOL_Reg+Learning_Disability+
                   Dementia+Cancer+CVD+CLD+COPD+Neurological+has_mh_any+
                   aae1a2_attend_12+ admissions_electives_12+admissions_emergency_12+gp_consultations_12+
                   dn_ppop,
                 data = CM)
summary(reg2)



cond1_results<-cond_results[, .(Parameter, estimate=exp(Coefficient), conf.high=exp(CI_high), conf.low=exp(CI_low), p.value=p)]

zi1_results<-zi_results[, .(Parameter, estimate=exp(Coefficient), conf.high=exp(CI_high), conf.low=exp(CI_low), p.value=p)]

fwrite(cond1_results, file="./results/cond1_results.csv")
fwrite(zi1_results, file="./results/zi1_results.csv")


#############################################################################################
#  Predicted Values 
#############################################################################################

CM[, old_dn_ppop:= dn_ppop]
CM[, dn_ppopp:= mean(dn_ppop)]
#  set bame to 0 in 1st stage
CM[, bame1:= 0]
# getting the predicted values of the main equation
CM[, pred1_need:=predict(reg1,newdata=CM, type="response")]
CM[, pred1_need_noeth:=predict(reg2,newdata=CM, type="response")]

plot(CM$contacts, CM$pred1_need)

summary(lm(CM$contacts ~ CM$pred1_need))

#############################################################################################
#  caculate need weights
#############################################################################################

# agreggate sum / mean / needs weights [needs weighted mean for the GP divided by mean for the total population] by GP
gp_pred<-CM[, list(pred1_need=sum(pred1_need, na.rm = T),
                   pred1_need_noeth=sum( pred1_need_noeth, na.rm = T),
                   contacts=sum(contacts, na.rm = T),
                   m_pred1_need=mean(pred1_need, na.rm = T),
                   m_pred1_need_noeth=mean(pred1_need_noeth, na.rm = T),
                   m_contacts=mean(contacts, na.rm = T),
                   nw_pred1_need=mean(pred1_need, na.rm = T)/mean(CM$pred1_need, na.rm = T),
                   nw_pred1_need_noeth=mean(pred1_need_noeth, na.rm = T)/mean(CM$pred1_need_noeth, na.rm = T),
                   nw_contacts=mean(contacts, na.rm = T)/mean(CM$contacts, na.rm = T),
                   reg_pop=.N, 
                   IMD=mean(IMD), 
                   age=mean(age)), by=.(gp_code,dn_ppop, nat_nw)]

# # 
# nat_nw<-as.data.table(readxl::read_xlsx("./data/nat_nw.xlsx", sheet=5))
# 
# 
# nat_nw<-nat_nw[, c("...1", "...2", "...6","...24")]
# names(nat_nw)<-c("gp_code", "ccg", "dist_contact", "nw2021")
# nat_nw[, nw2021:=as.numeric(nw2021)]
# nat_nw[, dist_contact:=as.numeric(dist_contact)]
# nat_nw<-nat_nw[is.na(nw2021)==F]
# nat_nw<-nat_nw[is.na(gp_code)==F]
# gp_pred<-merge(gp_pred,nat_nw, by="gp_code", all.x=T)
plot(gp_pred$nat_nw, gp_pred$nw_pred1_need)
summary(lm(nat_nw ~ IMD+age, data=gp_pred))
summary(lm(nw_pred1_need ~ IMD+age, data=gp_pred))

plot(gp_pred$nw_contacts, gp_pred$nw_pred1_need)
plot(gp_pred$nw_pred1_need_noeth, gp_pred$nw_pred1_need)
# including ethnicity - doesnt make much difference
plot(gp_pred$pred1_need_noeth, gp_pred$pred1_need)
plot(gp_pred$dn_ppop, gp_pred$nw_pred1_need)
plot(gp_pred$dn_ppop, gp_pred$nw_pred1_need)
plot( gp_pred$IMD, gp_pred$dn_ppop)
plot(gp_pred$age, gp_pred$dn_ppop)

#  new nurse distribution
dist_nurses<-unique(dist_lk[, .(team, district_nurses)])
gp_pred[, dn_need:=pred1_need*sum(dist_nurses$district_nurses)/sum(gp_pred$pred1_need)]

gp_pred[, dn_need_pp:=dn_need*10000/reg_pop]

gp_pred[, diff_dn:=dn_need_pp-dn_ppop]
summary(lm(diff_dn ~ nw_pred1_need, data=gp_pred))

gp_pred[, ab_diff_nat_nw:=nat_nw-nw_pred1_need]

gp_pred[, rel_diff_nat_nw:=nat_nw/nw_pred1_need]

#  IMD comparing local index to national index
summary(lm(gp_pred$rel_diff_nat_nw~gp_pred$IMD+gp_pred$age))
summary(lm(gp_pred$ab_diff_nat_nw~gp_pred$IMD))

plot(gp_pred$IMD, gp_pred$ab_diff_nat_nw)

#agregated to lsoa
lsoa_pred<-CM[, list(pred1_need=sum(pred1_need, na.rm = T),
                     pred1_need_noeth=sum(pred1_need_noeth, na.rm = T),
                     contacts=sum(contacts, na.rm = T),
                     m_pred1_need=mean(pred1_need, na.rm = T),
                     m_pred1_need_noeth=mean(pred1_need_noeth, na.rm = T),
                     m_contacts=mean(contacts, na.rm = T),
                     nw_pred1_need=mean(pred1_need, na.rm = T)/mean(CM$pred1_need, na.rm = T),
                     nw_pred1_need_noeth=mean(pred1_need_noeth, na.rm = T)/mean(CM$pred1_need_noeth, na.rm = T),
                     nw_contacts=mean(contacts, na.rm = T)/mean(CM$contacts, na.rm = T),
                     reg_pop=.N,
                     dn_ppop=mean(dn_ppop),
                     age=mean(age)), by=.(LSOA, IMD)]


lsoa_pred[, dn_need:=pred1_need*sum(dist_nurses$district_nurses)/sum(lsoa_pred$pred1_need)]
lsoa_pred[, dn_need_pp:=dn_need*10000/reg_pop]

#aggregated to DN team
team_pred<-CM[, list(pred1_need=sum(pred1_need, na.rm = T),
                     pred1_need_noeth=sum(pred1_need_noeth, na.rm = T),
                     contacts=sum(contacts, na.rm = T),
                     m_pred1_need=mean(pred1_need, na.rm = T),
                     m_pred1_need_noeth=mean(pred1_need_noeth, na.rm = T),
                     m_contacts=mean(contacts, na.rm = T),
                     nw_pred1_need=mean(pred1_need, na.rm = T)/mean(CM$pred1_need, na.rm = T),
                     nw_pred1_need_noeth=mean(pred1_need_noeth, na.rm = T)/mean(CM$pred1_need_noeth, na.rm = T),
                     nw_contacts=mean(contacts, na.rm = T)/mean(CM$contacts, na.rm = T),
                     reg_pop=.N,
                     dn_ppop=mean(dn_ppop),
                     age=mean(age), 
                     IMD=mean(IMD)), by=.(district_nurses, team)]

sum(dist_nurses$district_nurses)
sum(team_pred$district_nurses)

team_pred[, dn_need:=pred1_need*sum(dist_nurses$district_nurses)/sum(team_pred$pred1_need)]
#  change if redistribute by need - i.e positive number means gain nurses
team_pred[, ab_diff:=dn_need-district_nurses]
team_pred[, dn_need_pp:=dn_need*10000/reg_pop]
team_pred[, curr_dn_pp:=district_nurses*10000/reg_pop]
#  change if redistribute by need - i.e positive number means gain nurses
team_pred[, pp_diff:=dn_need_pp-curr_dn_pp]






#deprivation age and gender
# 
CM[, quint_imd:=ntile(IMD, n=5)]
# Inverting quintiles so they match the IMD technical reports
CM$quint_imd <- 6 - CM$quint_imd
# age intervals
CM[, age_g:=cut(age, breaks = c(-Inf,30,50,70,90,Inf), # creating age intervals
                
                labels = c("18-29 years" , "30-49 years", "50-69 years",     # labeling age intervals
                           "70-89 years",  "90+ years"),
                right = TRUE)]

imd_age_pred<-CM[is.na(gender)==F, list(pred1_need=sum(pred1_need, na.rm = T),
                                        pred1_need_noeth=sum(pred1_need_noeth, na.rm = T),
                                        contacts=sum(contacts, na.rm = T),
                                        m_pred1_need=mean(pred1_need, na.rm = T),
                                        m_pred1_need_noeth=mean(pred1_need_noeth, na.rm = T),
                                        m_contacts=mean(contacts, na.rm = T),
                                        nw_pred1_need=mean(pred1_need, na.rm = T)/mean(CM$pred1_need, na.rm = T),
                                        nw_pred1_need_noeth=mean(pred1_need_noeth, na.rm = T)/mean(CM$pred1_need_noeth, na.rm = T),
                                        nw_contacts=mean(contacts, na.rm = T)/mean(CM$contacts, na.rm = T),
                                        reg_pop=.N), by=.(quint_imd,age_g, gender)]

imd_age_pred[, gender:=factor(gender, labels = c("Men", "Women"))]

imd_age_pred[, quint_imd:=factor(quint_imd, labels = c("Most deprived", "Quintile 2","Quintile 3","Quintile 4","Least deprived"))]


# extracting it to a csv table
fwrite(gp_pred,file="./results/gp_liv_need_index.csv", row.names = TRUE)

fwrite(lsoa_pred,file="./results/lsoa_liv_need_index.csv", row.names = TRUE)

fwrite(imd_age_pred,file="./results/imd_age_liv_need_index.csv", row.names = TRUE)

fwrite(team_pred,file="./results/team_liv_need_index.csv", row.names = TRUE)

# present results - age , sex and deprivation distribution

# imd_age_pred<-fread("./results/imd_age_liv_need_index.csv")

hmap1<-ggplot(data=imd_age_pred, aes(y=age_g, x=quint_imd,fill=nw_pred1_need ))+geom_tile()+facet_wrap(gender~.) +
  scale_fill_continuous_sequential(palette = "SunsetDark", name="District Nursing needs index")+
  theme(legend.text=element_text(size=12))+labs(y="Age Group", x="IMD Quintile")+
  theme_minimal()+
  theme(legend.position = "top")+
  theme(legend.direction = "horizontal")
hmap1

ggsave("./results/dn_need_age_dep_hm.pdf", height=5, width = 10)
ggsave("./results/dn_need_age_dep_hm.svg", height=5, width = 10)


lsoa_pred<-fread("./results/lsoa_liv_need_index.csv")
sp_lsoa <- st_read(normalizePath("./data/shp/Lower_Layer_Super_Output_Areas__December_2011__Boundaries_EW_BGC.shp"))
sp_lsoa <- sp_lsoa[sp_lsoa$LSOA11CD %in% lsoa_pred$LSOA, ]

#  map1
lsoa_pred[, val_need:=nw_pred1_need]
lsoa_pred[, pc_need:=ntile(val_need,n=5)]

rate_lsoa <- merge(sp_lsoa,lsoa_pred, by.x="LSOA11CD", by.y="LSOA", all.x=T)



label1<-paste0("<", round(max(lsoa_pred[pc_need==1]$val_need), 2))
label2<-paste0( round(min(lsoa_pred[pc_need==2]$val_need), 2), "-", round(max(lsoa_pred[pc_need==2]$val_need), 2))
label3<-paste0( round(min(lsoa_pred[pc_need==3]$val_need), 2), "-", round(max(lsoa_pred[pc_need==3]$val_need), 2))
label4<-paste0( round(min(lsoa_pred[pc_need==4]$val_need), 2), "-", round(max(lsoa_pred[pc_need==4]$val_need), 2))
label5<-paste0(">", round(min(lsoa_pred[pc_need==5]$val_need), 2))
map1<-ggplot() +
  geom_sf(data = rate_lsoa,aes(fill = as.factor(pc_need)), lwd=0) +
  scale_fill_discrete_sequential(palette = "SunsetDark", name="quintile",alpha = 1, labels=c(label1, label2, label3,label4, label5))+
  guides(fill = guide_legend(title = "")) +
  theme(legend.text=element_text(size=12))+
  theme_void()+
  annotation_scale(location = "bl", width_hint = 0.5) +
  annotation_north_arrow(location = "tr", which_north = "true",
                         pad_x = unit(0.75, "in"), pad_y = unit(0.5, "in"),
                         style = north_arrow_fancy_orienteering)+ggtitle("1. District Nursing needs weights")

map1




#  map2 - current staffing
lsoa_pred[, val_need:=dn_ppop]
lsoa_pred[, pc_need:=ntile(val_need,n=5)]

rate_lsoa <- merge(sp_lsoa,lsoa_pred, by.x="LSOA11CD", by.y="LSOA", all.x=T)

label1<-paste0("<", round(max(lsoa_pred[pc_need==1]$val_need), 2))
label2<-paste0( round(min(lsoa_pred[pc_need==2]$val_need), 2), "-", round(max(lsoa_pred[pc_need==2]$val_need), 2))
label3<-paste0( round(min(lsoa_pred[pc_need==3]$val_need), 2), "-", round(max(lsoa_pred[pc_need==3]$val_need), 2))
label4<-paste0( round(min(lsoa_pred[pc_need==4]$val_need), 2), "-", round(max(lsoa_pred[pc_need==4]$val_need), 2))
label5<-paste0(">", round(min(lsoa_pred[pc_need==5]$val_need), 2))
map2<-ggplot() +
  geom_sf(data = rate_lsoa,aes(fill = as.factor(pc_need)), lwd=0) +
  scale_fill_discrete_sequential(palette = "SunsetDark", name="quintile",alpha = 1, labels=c(label1, label2, label3,label4, label5))+
  guides(fill = guide_legend(title = "")) +
  theme(legend.text=element_text(size=12))+
  theme_void()+
  annotation_scale(location = "bl", width_hint = 0.5) +
  annotation_north_arrow(location = "tr", which_north = "true",
                         pad_x = unit(0.75, "in"), pad_y = unit(0.5, "in"),
                         style = north_arrow_fancy_orienteering)+ggtitle("2. Current District Nursing staffing FTE/10,000 population")

map2

map12<-plot_grid(map1, map2, ncol = 2)
map12
ggsave("./results/maps_dn_nw.pdf", height=5, width = 10)
ggsave("./results/maps_dn_nw.svg", height=5, width = 10)



team_pred[order(IMD), team_id:=1:.N]

plot<-ggplot(data=team_pred, aes(y=ab_diff, x=team_id ))+geom_col() +
  theme_minimal()+labs(y="Change in District Nursing per 10,000 population", x="DN teams (ranked by catchment IMD)")
plot

plot<-ggplot(data=team_pred, aes(y= pp_diff, x=team_id , fill=nw_pred1_need))+geom_col() +
  theme_minimal()+labs(y="Change in District Nursing per 10,000 population", x="District Nursing teams (ranked by catchment IMD)")+
  scale_fill_continuous_sequential(palette = "SunsetDark", name="District Nursing needs index")+
  scale_x_continuous(breaks = c(1:16))+scale_y_continuous(breaks = c(-7:7))
plot

ggsave("./results/change_DN_teamneed_imd.svg", height=7, width = 10)
ggsave("./results/change_DN_teamneed_imd.pdf", height=7, width = 10)


team_pred[order(age), team_id:=1:.N]

plot<-ggplot(data=team_pred, aes(y= pp_diff, x=team_id , fill=nw_pred1_need))+geom_col() +
  theme_minimal()+labs(y="Change in District Nursing per 10,000 population", x="District Nursing teams (ranked by catchment mean Age)")+
  scale_fill_continuous_sequential(palette = "SunsetDark", name="District Nursing needs index")+
  scale_x_continuous(breaks = c(1:16))+scale_y_continuous(breaks = c(-7:7))

plot

ggsave("./results/change_DN_teamneed_age.pdf", height=7, width = 10)
ggsave("./results/change_DN_teamneed_age.svg", height=7, width = 10)


team_pred[order(nw_pred1_need), team_id:=1:.N]
plot<-ggplot(data=team_pred, aes(y= pp_diff, x=team_id , fill=nw_pred1_need))+geom_col() +
  theme_minimal()+labs(y="Change in District Nursing per 10,000 population", x="District Nursing teams (ranked by catchment mean need)")+
  scale_fill_continuous_sequential(palette = "SunsetDark", name="District Nursing needs index")+
  scale_x_continuous(breaks = c(1:16))+scale_y_continuous(breaks = c(-7:7))
plot

ggsave("./results/change_DN_teamneed_need.pdf", height=7, width = 10)
ggsave("./results/change_DN_teamneed_need.svg", height=7, width = 10)


#discriptive stats

final <- CM[, c("contacts", "age_g", "IMD", "gender", "bame", "dead",
                "log_time_alive", "NursingCareHomeFlag", 
                "Living_alone", "Palliative_EOL_Reg", 
                "Learning_Disability", "Dementia", "Cancer", 
                "CVD", "CLD", "COPD", "Neurological", 
                "has_mh_any", "aae1a2_attend_12", 
                "admissions_electives_12", 
                "admissions_emergency_12", 
                "gp_consultations_12", "dn_ppop", "EthnicGroup")]


final <- na.omit(final)


summary_table <- CreateTableOne(data = final,
                                vars = c("contacts", "age_g", "IMD", "gender", "bame", "dead",
                                         "log_time_alive", "NursingCareHomeFlag", 
                                         "Living_alone", "Palliative_EOL_Reg", 
                                         "Learning_Disability", "Dementia", "Cancer", 
                                         "CVD", "CLD", "COPD", "Neurological", 
                                         "has_mh_any", "aae1a2_attend_12", 
                                         "admissions_electives_12", 
                                         "admissions_emergency_12", 
                                         "gp_consultations_12", "dn_ppop"))

# Print the summary table


print(summary_table)

positive_contacts <- final[final$contacts > 0, ]

# Create summary table for individuals with positive contacts
summary_table_positive <- CreateTableOne(data = positive_contacts,
                                         vars = c("contacts", "age_g", "IMD", "gender", "bame", "dead",
                                                  "log_time_alive", "NursingCareHomeFlag", 
                                                  "Living_alone", "Palliative_EOL_Reg", 
                                                  "Learning_Disability", "Dementia", "Cancer", 
                                                  "CVD", "CLD", "COPD", "Neurological", 
                                                  "has_mh_any", "aae1a2_attend_12", 
                                                  "admissions_electives_12", 
                                                  "admissions_emergency_12", 
                                                  "gp_consultations_12", "dn_ppop"))

