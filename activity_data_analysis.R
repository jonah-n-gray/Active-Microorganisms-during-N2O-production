#This file contains data organization and analysis of the cell sorting data collected in summer 2024.  This file contains 1 figure for cell count data


setwd("C:/Users/murra/OneDrive/Desktop/Research/Experiments/Act_N2O/Boncat data/R_analysis")

library(tidyverse)
library(ggplot2)
#read in packages
t=c(3,5.67,8.42,12,15)
t=rep(t, each = 5)


fluor_data=read_csv("complete_analysis-Info.csv")
sort_count_data=read_csv("complete_analysis-Statistics.csv")
#reading in the data
#type selects the correct data. S is for activity, c is for counts, 0 is for control.

glimpse(fluor_data)
head(fluor_data)
tail(fluor_data)

glimpse(sort_count_data)
head(sort_count_data)
tail(sort_count_data)


###
### DATA CLEANING
#this section cleans up the data by removing columns and separating data into more informative tibbles
###

#sort and count data
s_c_clean=sort_count_data %>%
  select('Data Set', Gate, Number, '%Total','%Gated',Logic,TYPE)
#cleaning out unnessary columns

active_sort=s_c_clean %>%
  select('Data Set', Gate, '%Gated', TYPE) %>%
  filter(Gate == "Active")%>%
  filter(TYPE == "s")%>%
  mutate(index=c(36,37,47,3,14,15,25,26,4,5,48,49,1,2,12,13,16,34,35,24,27,38,45,46,4,5,15,16,26,27,38,1,12,23,34,36,37,47,48,49,2,14,24,35,3,23,25,45,45,46))%>%
  mutate(date=c(717,717,717,718,718,718,718,718,719,719,719,719,722,722,722,722,722,722,722,723,723,723,723,723,813,813,813,813,813,813,813,815,815,815,815,815,815,815,815,815,816,816,816,816,819,819,819,819,819,819))%>%
  arrange(index)
#parsing active cell SORT data, added index column.  Will use for percent activity analysis

active_count=s_c_clean %>%
  select('Data Set', Gate, Number, '%Gated', TYPE) %>%
  filter(Gate == "Active")%>%
  filter(TYPE == "c")%>%
  mutate(index=c(36,37,47,3,14,15,25,26,4,5,48,49,1,2,12,13,16,34,35,24,27,38,45,46))%>%
  mutate(date=c(717,717,717,718,718,718,718,718,719,719,719,719,722,722,722,722,722,722,722,723,723,723,723,723))%>%
  arrange(index)
#parsing active cell COUNT data, added index column.  Will be used to estimate absolute abundance in future analysis
##NOTE: these are not pooled yet between replicates and still need to be pooled since the PCR procedure pooled active cell vials

inactive_sort=s_c_clean %>%
  select('Data Set', Gate, '%Gated', TYPE) %>%
  filter(Gate == "Inactive")%>%
  filter(TYPE == "s")%>%
  mutate(index=c(36,37,47,3,14,15,25,26,4,5,48,49,1,2,12,13,16,34,35,24,27,38,45,46,4,5,15,16,26,27,38,1,12,23,34,36,37,47,48,49,2,14,24,35,3,23,25,45,45,46))%>%
  mutate(date=c(717,717,717,718,718,718,718,718,719,719,719,719,722,722,722,722,722,722,722,723,723,723,723,723,813,813,813,813,813,813,813,815,815,815,815,815,815,815,815,815,816,816,816,816,819,819,819,819,819,819))%>%
  arrange(index)

inactive_count=s_c_clean %>%
  select('Data Set', Gate, Number, '%Gated', TYPE) %>%
  filter(Gate == "Inactive")%>%
  filter(TYPE == "c")%>%
  mutate(index=c(36,37,47,3,14,15,25,26,4,5,48,49,1,2,12,13,16,34,35,24,27,38,45,46))%>%
  mutate(date=c(717,717,717,718,718,718,718,718,719,719,719,719,722,722,722,722,722,722,722,723,723,723,723,723))%>%
  arrange(index)
#parsing inactive cell COUNT data, added index column.  Will be used to estimate absolute abundance in future analysis


total_count=s_c_clean %>%
  select('Data Set', Gate, Number, '%Gated', TYPE) %>%
  filter(Gate == "A")%>%
  filter(TYPE == "c")%>%
  mutate(index=c(36,37,47,3,14,15,25,26,4,5,48,49,1,2,12,13,16,34,35,24,27,38,45,46))%>%
  mutate(date=c(717,717,717,718,718,718,718,718,719,719,719,719,722,722,722,722,722,722,722,723,723,723,723,723))%>%
  arrange(index)
#parsing TOTAL cell COUNT data, added index column.  Will be used to estimate absolute abundance in future analysis

controls=s_c_clean%>%
  select('Data Set', Gate, Number, '%Gated', TYPE) %>%
  filter(TYPE=="0")

all_count_data=bind_rows(active_count,inactive_count,total_count)

#write.csv(all_count_data, "C:/Users/murra/OneDrive/Desktop/Research/Experiments/Act_N2O/Boncat data/R_analysis/all_count_data.csv")

##NOTE: inactive and active counts =/= total counts because a gap was left between active and inactive gates due to uncertainty
  ##SYBR control indicates that we are sure of inactivity, SYBR+af_647 control was limited to .05% positivity rate. 
  ## total - active = may be a better representatation of our assumed inactive population.

##Fluor data
active_fluor=fluor_data %>%
  select(-`$BTIM`,-`$ETIM`) %>%
  filter(TYPE == "s")%>%
  mutate(ind=c(36,37,47,3,14,15,25,26,4,5,48,49,1,2,12,13,16,34,35,24,27,38,45,46,4,5,15,16,26,27,38,1,12,23,34,36,37,47,48,49,2,14,24,35,3,23,25,45,45,46))%>%
  mutate(date=c(717,717,717,718,718,718,718,718,719,719,719,719,722,722,722,722,722,722,722,723,723,723,723,723,813,813,813,813,813,813,813,815,815,815,815,815,815,815,815,815,816,816,816,816,819,819,819,819,819,819))%>%
  mutate(`[Active] ACTIVE Alexa 647-A Median` = as.numeric(`[Active] ACTIVE Alexa 647-A Median`),
       `[Active] ACTIVE Alexa 647-A Geometric Mean` = as.numeric(`[Active] ACTIVE Alexa 647-A Geometric Mean`),
        `[Active] ACTIVE Alexa 647-A Standard Deviation` = as.numeric(gsub("[^0-9.]", "",`[Active] ACTIVE Alexa 647-A Standard Deviation`)))%>%
  group_by(ind)%>%
  mutate(af_median=mean(`[Active] ACTIVE Alexa 647-A Median`),af_geom_mean=mean(`[Active] ACTIVE Alexa 647-A Geometric Mean`), af_sd=mean(`[Active] ACTIVE Alexa 647-A Standard Deviation`))%>%
  arrange(ind)
#Fluor data w/o controls
#averaged technical replicates

control_fluor=fluor_data %>%
  select(-`$BTIM`,-`$ETIM`) %>%
  filter(TYPE == "0")%>%
  mutate(index=c(0,0,9,0,0,0,9,0,9,9,0,0,9,0,0,0,0,9,9,0,0,9,0,0,9,0,0,9,0,0,9,0,0,9))%>%
  filter(index=="9")%>%
  mutate(date=c(717,718,719,719,722,723,723,723,813,815,816,819))%>%
  mutate(`[Active] ACTIVE Alexa 647-A Median` = as.numeric(`[Active] ACTIVE Alexa 647-A Median`),
         `[Active] ACTIVE Alexa 647-A Geometric Mean` = as.numeric(`[Active] ACTIVE Alexa 647-A Geometric Mean`),
         `[Active] ACTIVE Alexa 647-A Standard Deviation` = as.numeric(gsub("[^0-9.]", "",`[Active] ACTIVE Alexa 647-A Standard Deviation`)))%>%
  select(-index)%>%
  group_by(date)%>%
  summarise(af_median=mean(`[Active] ACTIVE Alexa 647-A Median`),af_geom_mean=mean(`[Active] ACTIVE Alexa 647-A Geometric Mean`), af_sd=mean(`[Active] ACTIVE Alexa 647-A Standard Deviation`))
#fluor data from the af_647 controls only.  Will be used to calculate relative fluorescnce
#averaged to reflect 1 value per date

#for fluoresence data, first find RFU per date then average by time point

###
#Data Manipulation
###
#calculating averages for technical replicates and calculating RFU
#Some averages are already calculated - control fluor and sorted flour are done

#sort data averaging between technical replicates
active_sort_clean=active_sort%>%
  mutate(percent_active=as.numeric(`%Gated`))%>%
  group_by(index)%>%
  summarise(percent_active=mean(percent_active))%>%
  mutate(time=c(1,1,1,1,1,2,2,2,2,2,3,3,3,3,3,4,4,4,4,4,5,5,5,5,5))

inactive_sort_clean=inactive_sort%>%
  mutate(percent_inactive=as.numeric(`%Gated`))%>%
  group_by(index)%>%
  summarise(percent_inactive=mean(percent_inactive))%>%
  mutate(time=c(1,1,1,1,1,2,2,2,2,2,3,3,3,3,3,4,4,4,4,4,5,5,5,5,5))
all_sort_data=bind_rows(active_sort_clean,inactive_sort_clean)
#write.csv(all_sort_data, "C:/Users/murra/OneDrive/Desktop/Research/Experiments/Act_N2O/Boncat data/R_analysis/all_sort_data.csv")
#Count data

ul_used=64
#ul used in sort based of flow rate from Mitch.  32 uL per min x 2 mins = 64 ul

ul_total=1000
#total volume of sample tube - 500ul sample +500 pbs (dilution)


sort_dilution=2
#dilution used to sort.  doubled volume to sort (see above)
storage_dilution=1.67
#glycerol dilution to freeze samples 
extract_percent=.720/(8+.310+.132)
#is aliquot is x percent of extract volume 

gdw_soil=1.50 #from soil data calculations

inactive_count=inactive_count%>%
  mutate(inactive_per_gdw=inactive_count$Number/ul_used*ul_total*sort_dilution*storage_dilution/extract_percent/gdw_soil)
active_count=active_count%>%
  mutate(active_per_gdw=active_count$Number/ul_used*ul_total*sort_dilution*storage_dilution/extract_percent/gdw_soil)
total_count=total_count%>%
  mutate(total_per_gdw=(total_count$Number)/ul_used*ul_total*sort_dilution*storage_dilution/extract_percent/gdw_soil)



###
###Graphs and models
###
t=c(3,5.67,8.42,12,15)
t=rep(t, each = 5)
active_sort_clean$t= as.numeric(t)
#plot(active_sort_clean$percent_active~active_sort_clean$t)

plot(log(unique(active_fluor$af_median))~active_sort_clean$t)
plot((active_fluor_RFU_clean$af_med_RFU)~active_sort_clean$t)
#plot((active_fluor$af_geom_mean)~active_sort_clean$time)

fluor_lm_logmedian=lm(log(unique(active_fluor$af_median))~active_sort_clean$t)
summary(fluor_lm_logmedian)

fluor_lm_loggeom=lm(log(active_fluor$af_geom_mean)~active_sort_clean$time)
summary(fluor_lm_loggeom)

RFU_lm_logmedian=lm((log(active_fluor_RFU_clean$af_med_RFU))~active_sort_clean$t)
summary(RFU_lm_logmedian)

t.2=t[-12]
total_count$t=t.2
plot(total_count$total_per_gdw~t.2)
total_count_lm=lm(total_count$total_per_gdw~t.2)
summary(total_count_lm)

### anova test for difference in percent of active organisms by time.
one_way=aov(percent_active~as.factor(time), active_sort_clean)

anova=as.data.frame(summary(one_way))

tukey=TukeyHSD(one_way)
tukey=as.data.frame(tukey$`as.factor(time)`)
write.csv(tukey, "tukey_output.csv")


cell_count= ggplot(total_count, aes(x = t, y = total_per_gdw, group=t)) +
  #geom_smooth(method = "lm", colour = "blue", fill = "gray", size=1) +
  geom_boxplot(colour = "darkblue",fill = "gray")+
  geom_jitter(width = 0.2, color="blue", size=3)+
  # Set x and y limits
  scale_x_continuous(breaks = seq(0, 15, by =3),limits = c(0, 18)) +
  theme_classic()+xlab("Time (hr)") + ylab("Absolute cell abundance per gram soil") +  # Set axis labels
  #ggtitle("Level of activity vs N2O-N production rate") +  # Set plot title
  theme(title = element_text(size = 20),
        axis.text = element_text(size = 20))+
  theme(axis.title.x = element_text(size = 20), 
        axis.title.y = element_text(size = 20))
cell_count


ggsave("cell_count.pdf", cell_count, width = 3000, height = 3000, units = "px")
