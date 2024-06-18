#Purpose: QC/QA check for Kate's dust data 

#Load packages----
library(haven)
library(tidyverse)
library(Publish)
library(readxl)
library(broom)
library(RColorBrewer)
library(wesanderson)
library(haven)

#Set working dir----
setwd("C:\\Users\\afossa\\Dropbox (Brown)\\Data")

#Dust PFAS concentrations and loadings----

##Get area of floor sampled from original dust dataset----
cdc_dust<-read_sas("C:\\Users\\afossa\\Dropbox (Brown)\\Data\\exp_dust.sas7bdat")

dust_area_and_weight<-cdc_dust %>%
  dplyr::dplyr::filter(Visit=="24M") %>% 
  dplyr::select(Participant_ID,Sampling_Area_m2,SW_Sieved) %>% 
  mutate(Participant_ID_chr=as.character(Participant_ID)) %>% 
  dplyr::select(-Participant_ID) %>% 
  rename(Participant_ID=Participant_ID_chr)

##Get names for dust concentration sheet----
dust_pfas_names<-read_xlsx(
  file.path(getwd(),"PFAS Data HOME Dust Sept 2023_9.27.2023.xlsx"),
  sheet="Dust Concentrations",
  range="A2:GB2",
  col_names = F) %>% 
  as.character()

dust_pfas_names[2]<-"LOD"

##Import values----
dust_pfas<-read_xlsx(
  file.path(getwd(),"PFAS Data HOME Dust Sept 2023_9.27.2023.xlsx"),
  sheet="Dust Concentrations",
  range="A4:GB45",
  col_names=dust_pfas_names,
  na=c("N/F","N/A")
  )

##Melt to tidy dataset----
dust_pfas %>% 
  pivot_longer(
    cols = 3:184,
    names_to = "label", 
    values_to = "result"
  ) %>%
  mutate(
    batch=str_split_i(label,"_",1),
    Participant_ID=str_split_i(label,"_",2),
    batch_f=factor(
      case_when(
        batch=="B1"~1,
        batch=="B2"~2,
        batch=="B3"~3,
        batch=="B4"~4
      ),
      levels=c(1:4),
      labels=c("Batch 1","Batch 2","Batch 3","Batch 4")
    ),
    detect=factor(
      case_when(
        is.na(result)~1,
        !is.na(result) & (result<LOD)~2,
        !is.na(result) & (result>=LOD)~3
      ),
      levels=c(1:3),
      labels=c("Non-detect","Below LOD","Detect")
    ),
    result_corrected=case_when(
      detect %in% c("Detect","Below LOD")~result,
      detect=="Non-detect"~LOD/sqrt(2)
    ),
    class=case_when(
      `Sample (Batch _PID)` %in% c("PFBA","PFPeA","PFHxA","PFHpA","PFOA","PFNA","PFDA","PFUnA","PFDoA","PFTrDA","PFTeDA")~"PF carboxylic acids",
      `Sample (Batch _PID)` %in% c("PFBS","PFPeS","PFHxS","PFHpS","PFOS","PFNS","PFDS","PFDoS")~"PF sulfonic acids",
      `Sample (Batch _PID)` %in% c("4:2-FTS","6:2-FTS","8:2-FTS")~"FT sulfonic acids",
      `Sample (Batch _PID)` %in% c("PFOSA","NMeFOSA","NEtFOSA")~"PFO sulfonamides",
      `Sample (Batch _PID)` %in% c("NMeFOSAA","NEtFOSAA")~"PFO sulfonamidoacetic acids",
      `Sample (Batch _PID)` %in% c("NMeFOSE","NEtFOSE")~"PFO sulfonamide ethanols",
      `Sample (Batch _PID)` %in% c("HFPO-DA","ADONA","PFMPA","PFMBA","NFDHA")~"PFE carboxylic acids",
      `Sample (Batch _PID)` %in% c("9Cl-PF3ONS","11Cl-PF3OUdS","PFEESA")~"Ether sulfonic acids",
      `Sample (Batch _PID)` %in% c("3:3 FTCA","5:3 FTCA","7:3 FTCA")~"FT carboxylic acids",
      `Sample (Batch _PID)` %in% c("FBSA","FHxSA","NFDHA","HFPODA-GenX")~"Other"
      )
 ) %>% 
  rename(
    analyte=`Sample (Batch _PID)`
  )->dust_pfas_long
  
  
dust_pfas_long_w_weights_and_area<-reduce(
    list(dust_pfas_long,dust_area_and_weight),
    left_join,
    by="Participant_ID"
  ) %>% 
  mutate(
    Sampling_Area_m2=if_else(is.na(Sampling_Area_m2),mean(cdc_dust$Sampling_Area_m2,na.rm=T),Sampling_Area_m2),
    SW_Sieved=if_else(is.na(SW_Sieved),mean(cdc_dust$SW_Sieved,na.rm=T),SW_Sieved),
    loading=(result_corrected*SW_Sieved)/Sampling_Area_m2
  )

##Export to flat file----
dust_pfas_long_w_weights_and_area %>% write_csv("C:\\Users\\afossa\\Dropbox (Brown)\\Dust and EDC R21\\data\\dust_pfas_long.csv")

##Explore and visualize----
###Boxplots----
round2<-function(x){round(x,1)}

dust_pfas_long %>% 
  ggplot()+
  geom_boxplot(
    aes(y=analyte,x=result_corrected,fill=class)
  )+
  facet_wrap(~class,scales = "free")+
  scale_fill_manual(values=brewer.pal(10,"Paired"),guide=NULL,name="")+
  scale_x_continuous(trans="log10",n.breaks = 6,labels = round2)+
  theme_classic()+
  theme(
    axis.text = element_text(colour="black")
  )+
  labs(
    title="PFAS in HOME Participant House Dust",
    x="Analyte Concentrations (ug/kg)",
    y="",
    caption = "Note: X-axis is on log scale"
  )

ggsave("C:\\Users\\afossa\\Dropbox (Brown)\\Dust and EDC R21\\output\\dust_pfas_box.tiff",width=1920,height=1080,scale=2,units='px')

###Boxplot by batch
dust_pfas_long %>% 
  ggplot()+
  geom_boxplot(
    aes(x=batch_f,y=result_corrected,fill=class)
  )+
  facet_wrap(~analyte,scales = "free_y")+
  scale_y_continuous(trans="log",n.breaks = 5,labels = round2)+
  scale_fill_manual(values=brewer.pal(10,"Paired"),name="")+
  theme_classic()+
  theme(axis.text.x = element_text(angle=20))+
  labs(
    title="PFAS in House Dust by Batch",
    y="Analyte Concentrations (ug/kg)",
    x="",
    caption="Note: X-axis is on log scale"
  )

ggsave("C:\\Users\\afossa\\Dropbox (Brown)\\Dust and EDC R21\\output\\dust_pfas_box_by_batch.tiff",width=1920,height=1080,scale=2.5,units='px')

###Violin plots----
#Overall
dust_pfas_long %>% 
  ggplot()+
  geom_violin(
    aes(y=analyte,x=result_corrected,fill=class),
    scale="width",
    draw_quantiles = c(0.25, 0.5, 0.75)
  )+
  facet_wrap(~class,scales = "free")+
  scale_fill_manual(values=brewer.pal(10,"Paired"),guide=NULL,name="")+
  scale_x_continuous(trans="log",n.breaks = 6,labels = round2)+
  theme_classic()+
  theme(
    axis.text = element_text(colour="black")
  )+
  labs(
    title="PFAS in HOME Participant House Dust",
    x="Analyte Concentrations (ug/kg)",
    y="",
    caption = "Note: X-axis is on log scale\n*LOD"
  )+
  geom_point(aes(x=LOD,y=analyte),color='black',shape=8)

ggsave("C:\\Users\\afossa\\Dropbox (Brown)\\Dust and EDC R21\\output\\dust_pfas_violin.tiff",width=1920,height=1080,scale=2,units='px')

#By batch
dust_pfas_long %>% 
  ggplot()+
  geom_violin(
    aes(x=batch_f,y=result_corrected,fill=class),
    scale="width"
  )+
  facet_wrap(~analyte,scales = "free_y")+
  scale_y_continuous(trans="log",n.breaks = 5,labels = round2)+
  scale_fill_manual(values=brewer.pal(10,"Paired"),name="")+
  theme_classic()+
  theme(axis.text.x = element_text(angle=25))+
  labs(
    title="PFAS in House Dust by Batch",
    y="Analyte Concentrations (ug/kg)",
    x="",
    caption="Note: X-axis is on log scale\n------------ LOD"
  )+
  geom_hline(aes(yintercept=LOD),linetype='dashed')

ggsave("C:\\Users\\afossa\\Dropbox (Brown)\\Dust and EDC R21\\output\\dust_pfas_violin_by_batch.tiff",width=1920,height=1080,scale=2,units='px')

#By batch and intervention status
#Import core dataset
demo<-read_sas("C:\\Users\\afossa\\Dropbox (Brown)\\Data\\home_data_20230502.sas7bdat") %>% select(Participant_ID,Intervention) %>% 
  mutate(
    id_char=as.character(Participant_ID)
  ) %>% 
  select(
    -Participant_ID
  ) %>% 
  rename(
    Participant_ID=id_char
  ) #Need to do a type conversion routine here :(

dust_pfas_long_w_int<-left_join(dust_pfas_long,demo,by="Participant_ID")

dust_pfas_long_w_int %>% 
  dplyr::filter(Intervention !="No intervention") %>% 
  ggplot()+
  geom_violin(
    aes(x=batch_f,y=result_corrected,fill=Intervention),
    scale="width",
    position = position_dodge(width=0.75),
    alpha=0.75
  )+
  facet_wrap(~analyte,scales = "free_y")+
  scale_y_continuous(trans="log",n.breaks = 5,labels = round2)+
  scale_fill_manual(values=brewer.pal(3,"RdBu")[c(1,3)],name="")+
  theme_classic()+
  theme(axis.text.x = element_text(angle=25))+
  labs(
    title="PFAS in House Dust by Batch and Intervention Status",
    y="Analyte Concentrations (ug/kg)",
    x="",
    caption="Note: X-axis is on log scale\n------------ LOD"
  )+
  geom_hline(aes(yintercept=LOD),linetype='dashed')
  
ggsave("C:\\Users\\afossa\\Dropbox (Brown)\\Dust and EDC R21\\output\\dust_pfas_violin_by_batch_and_intervention.tiff",width=1920,height=1080,scale=2,units='px')

###Basic descriptives----
dust_pfas_long %>% 
  group_by(analyte) %>%
  summarise(
    `First Quartile`=quantile(result_corrected,probs = c(0.25),na.rm=T),
    Median=median(result_corrected,na.rm=T),
    `Third Quartile`=quantile(result_corrected,probs = c(0.75),na.rm=T),
    `Total n`=n(),
    `Detect n`=sum(detect=="Detect"),
    `Percent Non-detect`=mean(detect=="Non-detect"),
    `Percent <LOD`=mean(detect=="Below LOD"),
    `CV`=sqrt(var(result_corrected,na.rm=T))/mean(result_corrected,na.rm=T)
  ) %>% write_csv("C:\\Users\\afossa\\Dropbox (Brown)\\Dust and EDC R21\\output\\dust_pfas_desc.csv")

##By batch
dust_pfas_long %>% 
  group_by(analyte,batch_f) %>%
  summarise(
    `First Quartile`=quantile(result_corrected,probs = c(0.25),na.rm=T),
    Median=median(result_corrected,na.rm=T),
    `Third Quartile`=quantile(result_corrected,probs = c(0.75),na.rm=T),
    `Total n`=n(),
    `Detect n`=sum(detect=="Detect"),
    `Percent Non-detect`=mean(detect=="Non-detect"),
    `Percent <LOD`=mean(detect=="Below LOD"),
    `CV`=sqrt(var(result_corrected,na.rm=T))/mean(result_corrected,na.rm=T)
  ) %>% write_csv("C:\\Users\\afossa\\Dropbox (Brown)\\Dust and EDC R21\\output\\dust_pfas_desc_by_batch.csv")

###Check for batch effects----
dust_pfas_long %>% 
  dplyr::filter(!(analyte %in% c("3:3 FTCA","5:3 FTCA",""))) %>% #dplyr::filter out analytes that have <=1 batch with detects
  with(.,split(.,analyte))->dust_pfas_long_by_analyte

dust_batch_fits<-map_dfr(dust_pfas_long_by_analyte,~tidy(lm(log(result_corrected)~batch_f,data=.x)),.id="analyte")

dust_batch_fits %>% write_csv("C:\\Users\\afossa\\Dropbox (Brown)\\Dust and EDC R21\\output\\dust_batch_regs.csv")

#NSTF CRM----
##Get names from NIST sheet----
crm_names<-read_xlsx(
  file.path(getwd(),"PFAS Data HOME Dust Sept 2023_9.27.2023.xlsx"),
  sheet="NIST 2585 Dust",
  range="T3:AG3",
  col_names = F
) %>% as.character()

##Import value columns----
crm_values<-read_xlsx(
  file.path(getwd(),"PFAS Data HOME Dust Sept 2023_9.27.2023.xlsx"),
  sheet="NIST 2585 Dust",
  range="T5:AG46",
  col_names = crm_names,
  na="N/F"
)

##Import analyte and LOD columns----
crm_analytes<-read_xlsx(
  file.path(getwd(),"PFAS Data HOME Dust Sept 2023_9.27.2023.xlsx"),
  sheet="NIST 2585 Dust",
  range="A5:B46",
  col_names = c("analyte","LOD"),
  na="N/F"
)

##Merge analyte column and value columns----
crm<-bind_cols(crm_analytes,crm_values)

##Melt and tidy----
crm %>% 
  pivot_longer(
    cols=3:16,
    names_to="label",
    values_to = "result"
  ) %>% 
  mutate(
    batch=str_split_i(label,"_",1),
    sample=str_split_i(label,"_",2),
    batch_f=factor(
      case_when(
        batch=="B1"~1,
        batch=="B2"~2,
        batch=="B3"~3,
        batch=="B4"~4
      ),
      levels=c(1:4),
      labels=c("Batch 1","Batch 2","Batch 3","Batch 4")
    ),
    detect=factor(
      case_when(
        is.na(result)~1,
        !is.na(result) & (result<LOD)~2,
        !is.na(result) & (result>=LOD)~3
      ),
      levels=c(1:3),
      labels=c("Non-detect","Below LOD","Detect")
    ),
    result_corrected=case_when(
      detect %in% c("Detect","Below LOD")~result,
      detect=="Non-detect"~LOD/sqrt(2)
    ),
    class=case_when(
      analyte %in% c("PFBA","PFPeA","PFHxA","PFHpA","PFOA","PFNA","PFDA","PFUnA","PFDoA","PFTrDA","PFTeDA")~"PF carboxylic acids",
      analyte %in% c("PFBS","PFPeS","PFHxS","PFHpS","PFOS","PFNS","PFDS","PFDoS")~"PF sulfonic acids",
      analyte %in% c("4:2-FTS","6:2-FTS","8:2-FTS")~"FT sulfonic acids",
      analyte %in% c("PFOSA","NMeFOSA","NEtFOSA")~"PFO sulfonamides",
      analyte %in% c("NMeFOSAA","NEtFOSAA")~"PFO sulfonamidoacetic acids",
      analyte %in% c("NMeFOSE","NEtFOSE")~"PFO sulfonamide ethanols",
      analyte %in% c("HFPO-DA","ADONA","PFMPA","PFMBA","NFDHA")~"PFE carboxylic acids",
      analyte %in% c("9Cl-PF3ONS","11Cl-PF3OUdS","PFEESA")~"Ether sulfonic acids",
      analyte %in% c("3:3 FTCA","5:3 FTCA","7:3 FTCA")~"FT carboxylic acids",
      analyte %in% c("FBSA","FHxSA","NFDHA","HFPODA-GenX")~"Other"
    )
  )->crm_long

##Explore and visualize----
###Boxplots----
crm_long %>% 
  ggplot()+
  geom_boxplot(
    aes(y=analyte,x=log(result_corrected))
  )+
  facet_wrap(~class,scales ="free")+
  theme_classic()+
  labs(
    title="PFAS in NIST CRMs",
    x="Log-transformed Analyte Concentrations (ug/kg)",
    y=""
  )

###Basic descriptives----
crm_long %>% 
  group_by(analyte) %>%
  summarise(
    `First Quartile`=quantile(result_corrected,probs = c(0.25),na.rm=T),
    Median=median(result_corrected,na.rm=T),
    `Third Quartile`=quantile(result_corrected,probs = c(0.75),na.rm=T),
    `Total n`=n(),
    `Detect n`=sum(detect=="Detect"),
    `Percent Non-detect`=mean(detect=="Non-detect"),
    `Percent <LOD`=mean(detect=="Below LOD"),
    `CV`=sqrt(var(result_corrected,na.rm=T))/mean(result_corrected,na.rm=T)
  ) %>% write_csv("C:\\Users\\afossa\\Dropbox (Brown)\\Dust and EDC R21\\output\\nist_dust_desc.csv")

###Percent yield for NIST certified species
crm_cert_values<-read_xlsx(
  file.path(getwd(),"PFAS Data HOME Dust Sept 2023_9.27.2023.xlsx"),
  sheet="NIST 2585 Dust",
  range="AH5:AJ46",
  col_names = c("average","CV","ug/kg"),
  na="N/F"
)

crm_cert<-bind_cols(crm_analytes,crm_cert_values) %>% 
  dplyr::filter(!is.na(`ug/kg`)) %>% 
  mutate(
    pct_yeild=average/`ug/kg`
  ) %>% write_csv("C:\\Users\\afossa\\Dropbox (Brown)\\Dust and EDC R21\\output\\nist_pct_yield.csv")

###Look for batch effects----
crm_long %>% 
  dplyr::filter(analyte %in% c("PFBA","PFDoA","PFHpA","PFHxA","PFHxS","PFNA","PFOS","PFTrDA")) %>% #dplyr::filter to analytes that are NIST certified
  with(.,split(.,analyte))->crm_long_by_analyte

crm_batch_fits<-map(crm_long_by_analyte,~summary(lm(log(result_corrected)~batch_f,data=.x)))

#Joe's dust----
##Get names from Joe dust sheet----
joe_names<-read_xlsx(
  file.path(getwd(),"PFAS Data HOME Dust Sept 2023_9.27.2023.xlsx"),
  sheet="Joe Dust",
  range = "T3:AG3",
  col_names=F
) %>% as.character()

##Import value columns----
joe_values<-read_xlsx(
  file.path(getwd(),"PFAS Data HOME Dust Sept 2023_9.27.2023.xlsx"),
  sheet="Joe Dust",
  range = "T5:AG46",
  col_names=joe_names,
  na="N/F"
)

##Import analyte and LOD columns----
joe_analytes<-read_xlsx(
  file.path(getwd(),"PFAS Data HOME Dust Sept 2023_9.27.2023.xlsx"),
  sheet="Joe Dust",
  range = "A5:B46",
  col_names=c("analyte","LOD")
)

##Merge analyte column and value columns----
joe<-bind_cols(joe_analytes,joe_values)

##Melt and tidy----
joe %>% 
  pivot_longer(
    cols=3:16,
    names_to="label",
    values_to = "result"
  ) %>% 
  mutate(
    batch=str_split_i(label,"_",1),
    sample=str_split_i(label,"_",2),
    batch_f=factor(
      case_when(
        batch=="B1"~1,
        batch=="B2"~2,
        batch=="B3"~3,
        batch=="B4"~4
      ),
      levels=c(1:4),
      labels=c("Batch 1","Batch 2","Batch 3","Batch 4")
    ),
    detect=factor(
      case_when(
        is.na(result)~1,
        !is.na(result) & (result<LOD)~2,
        !is.na(result) & (result>=LOD)~3
      ),
      levels=c(1:3),
      labels=c("Non-detect","Below LOD","Detect")
    ),
    result_corrected=case_when(
      detect %in% c("Detect","Below LOD")~result,
      detect=="Non-detect"~LOD/sqrt(2)
    ),
    class=case_when(
      analyte %in% c("PFBA","PFPeA","PFHxA","PFHpA","PFOA","PFNA","PFDA","PFUnA","PFDoA","PFTrDA","PFTeDA")~"PF carboxylic acids",
      analyte %in% c("PFBS","PFPeS","PFHxS","PFHpS","PFOS","PFNS","PFDS","PFDoS")~"PF sulfonic acids",
      analyte %in% c("4:2-FTS","6:2-FTS","8:2-FTS")~"FT sulfonic acids",
      analyte %in% c("PFOSA","NMeFOSA","NEtFOSA")~"PFO sulfonamides",
      analyte %in% c("NMeFOSAA","NEtFOSAA")~"PFO sulfonamidoacetic acids",
      analyte %in% c("NMeFOSE","NEtFOSE")~"PFO sulfonamide ethanols",
      analyte %in% c("HFPO-DA","ADONA","PFMPA","PFMBA","NFDHA")~"PFE carboxylic acids",
      analyte %in% c("9Cl-PF3ONS","11Cl-PF3OUdS","PFEESA")~"Ether sulfonic acids",
      analyte %in% c("3:3 FTCA","5:3 FTCA","7:3 FTCA")~"FT carboxylic acids",
      analyte %in% c("FBSA","FHxSA","NFDHA","HFPODA-GenX")~"Other"
    )
  )->joe_long

##Explore and visualize----
###Boxplots----
joe_long %>% 
  ggplot()+
  geom_boxplot(
    aes(y=analyte,x=result_corrected,fill=class)
  )+
  facet_wrap(~class,scales = "free")+
  scale_fill_manual(values=brewer.pal(10,"Paired"),guide=NULL,name="")+
  scale_x_continuous(trans="log",n.breaks = 6,labels = round2)+
  theme_classic()+
  theme(
    axis.text = element_text(colour="black")
  )+
  labs(
    title="PFAS in Joe's House Dust",
    x="Analyte Concentrations (ug/kg)",
    y="",
    caption = "Note: X-axis is on log scale"
  )

ggsave("C:\\Users\\afossa\\Dropbox (Brown)\\Dust and EDC R21\\output\\joe_pfas_box.tiff",width=1920,height=1080,scale=2,units='px')

###Violin plots----
joe_long %>% 
  ggplot()+
  geom_violin(
    aes(y=analyte,x=result_corrected,fill=class),
    scale="width",
    draw_quantiles = c(0.25, 0.5, 0.75)
  )+
  facet_wrap(~class,scales = "free")+
  scale_fill_manual(values=brewer.pal(10,"Paired"),guide=NULL,name="")+
  scale_x_continuous(trans="log",n.breaks = 6,labels = round2)+
  theme_classic()+
  theme(
    axis.text = element_text(colour="black")
  )+
  labs(
    title="Joe's House Dust",
    x="Analyte Concentrations (ug/kg)",
    y="",
    caption = "Note: X-axis is on log scale\n*LOD"
  )+
  geom_point(aes(x=LOD,y=analyte),color='black',shape=8)

ggsave("C:\\Users\\afossa\\Dropbox (Brown)\\Dust and EDC R21\\output\\joe_pfas_violin.tiff",width=1920,height=1080,scale=2,units='px')

###Basic descriptives----
joe_long %>% 
  group_by(analyte) %>%
  summarise(
    `First Quartile`=quantile(result_corrected,probs = c(0.25),na.rm=T),
    Median=median(result_corrected,na.rm=T),
    `Third Quartile`=quantile(result_corrected,probs = c(0.75),na.rm=T),
    `Total n`=n(),
    `Detect n`=sum(detect=="Detect"),
    `Percent Non-detect`=mean(detect=="Non-detect"),
    `Percent <LOD`=mean(detect=="Below LOD"),
    `CV`=sqrt(var(result_corrected,na.rm=T))/mean(result_corrected,na.rm=T)
  ) %>% write_csv("C:\\Users\\afossa\\Dropbox (Brown)\\Dust and EDC R21\\output\\joe_dust_desc.csv")

#Blanks----
##Import data----
blanks<-read_xlsx(
  file.path(getwd(),"PFAS Data HOME Dust Sept 2023_9.27.2023.xlsx"),
  sheet="Blanks",
  range="A3:P45",
  col_names = T,
  na="N/F"
) %>% 
  rename(
    analyte=Sample,
    LOD=`LOD (ng/L)`
  )

blanks %>% 
  pivot_longer(
    cols=3:16,
    names_to="label",
    values_to = "result"
  ) %>% 
mutate(
    batch=str_split_i(label,"_",1),
    sample=str_split_i(label,"_",2),
    batch_f=factor(
      case_when(
        batch=="B1"~1,
        batch=="B2"~2,
        batch=="B3"~3,
        batch=="B4"~4
      ),
      levels=c(1:4),
      labels=c("Batch 1","Batch 2","Batch 3","Batch 4")
    ),
    detect=factor(
      case_when(
        is.na(result)~1,
        !is.na(result) & (result<LOD)~2,
        !is.na(result) & (result>=LOD)~3
      ),
      levels=c(1:3),
      labels=c("Non-detect","Below LOD","Detect")
    ),
    result_corrected=case_when(
      detect %in% c("Detect","Below LOD")~result,
      detect=="Non-detect"~LOD/sqrt(2)
    ),
    class=case_when(
      analyte %in% c("PFBA","PFPeA","PFHxA","PFHpA","PFOA","PFNA","PFDA","PFUnA","PFDoA","PFTrDA","PFTeDA")~"PF carboxylic acids",
      analyte %in% c("PFBS","PFPeS","PFHxS","PFHpS","PFOS","PFNS","PFDS","PFDoS")~"PF sulfonic acids",
      analyte %in% c("4:2-FTS","6:2-FTS","8:2-FTS")~"FT sulfonic acids",
      analyte %in% c("PFOSA","NMeFOSA","NEtFOSA")~"PFO sulfonamides",
      analyte %in% c("NMeFOSAA","NEtFOSAA")~"PFO sulfonamidoacetic acids",
      analyte %in% c("NMeFOSE","NEtFOSE")~"PFO sulfonamide ethanols",
      analyte %in% c("HFPO-DA","ADONA","PFMPA","PFMBA","NFDHA")~"PFE carboxylic acids",
      analyte %in% c("9Cl-PF3ONS","11Cl-PF3OUdS","PFEESA")~"Ether sulfonic acids",
      analyte %in% c("3:3 FTCA","5:3 FTCA","7:3 FTCA")~"FT carboxylic acids",
      analyte %in% c("FBSA","FHxSA","NFDHA","HFPODA-GenX")~"Other"
    )
) %>% dplyr::filter(!str_detect(label,"Blank8"))->blanks_long #dplyr::filtering out results from blank 8 because of an issue with residual PFAS in the machines. See Kate's email on July 10th 2023 in the chain "Meeting tomorrow". 

##Explore and visualize----
###Basic descriptive----
blanks_long %>% 
  group_by(analyte) %>%
  summarise(
    `First Quartile`=quantile(result_corrected,probs = c(0.25),na.rm=T),
    Median=median(result_corrected,na.rm=T),
    `Third Quartile`=quantile(result_corrected,probs = c(0.75),na.rm=T),
    `Total n`=n(),
    `Detect n`=sum(detect=="Detect"),
    `Percent Non-detect`=mean(detect=="Non-detect"),
    `Percent <LOD`=mean(detect=="Below LOD"),
    `CV`=sqrt(var(result_corrected,na.rm=T))/mean(result_corrected,na.rm=T)
  ) %>% write_csv("C:\\Users\\afossa\\Dropbox (Brown)\\Dust and EDC R21\\output\\blank_dust_desc.csv")

###Visualize blanks with house dust for reference----
####Import raw extract concentrations from house dust----
extracts<-read_xlsx(
  file.path(getwd(),"PFAS Data HOME Dust Sept 2023_9.27.2023.xlsx"),
  sheet="Extract Concentrations",
  range="A3:GB45",
  col_names = T,
  na=c("N/F","N/A")
)

####Melt and tidy raw extract concentrations----
extracts %>% 
  pivot_longer(
    cols=3:184,
    names_to="label",
    values_to = "result"
  ) %>% 
  rename(
    analyte=`Sample (Batch _PID)`,
    LOD=`LOD (ng/L)`
  ) %>% 
  mutate(
    batch=str_split_i(label,"_",1),
    Participant_ID=str_split_i(label,"_",2),
    batch_f=factor(
      case_when(
        batch=="B1"~1,
        batch=="B2"~2,
        batch=="B3"~3,
        batch=="B4"~4
      ),
      levels=c(1:4),
      labels=c("Batch 1","Batch 2","Batch 3","Batch 4")
    ),
    detect=factor(
      case_when(
        is.na(result)~1,
        !is.na(result) & (result<LOD)~2,
        !is.na(result) & (result>=LOD)~3
      ),
      levels=c(1:3),
      labels=c("Non-detect","Below LOD","Detect")
    ),
    result_corrected=case_when(
      detect %in% c("Detect","Below LOD")~result,
      detect=="Non-detect"~LOD/sqrt(2)
    ),
    source="House dust",
    class=case_when(
      analyte %in% c("PFBA","PFPeA","PFHxA","PFHpA","PFOA","PFNA","PFDA","PFUnA","PFDoA","PFTrDA","PFTeDA")~"PF carboxylic acids",
      analyte %in% c("PFBS","PFPeS","PFHxS","PFHpS","PFOS","PFNS","PFDS","PFDoS")~"PF sulfonic acids",
      analyte %in% c("4:2-FTS","6:2-FTS","8:2-FTS")~"FT sulfonic acids",
      analyte %in% c("PFOSA","NMeFOSA","NEtFOSA")~"PFO sulfonamides",
      analyte %in% c("NMeFOSAA","NEtFOSAA")~"PFO sulfonamidoacetic acids",
      analyte %in% c("NMeFOSE","NEtFOSE")~"PFO sulfonamide ethanols",
      analyte %in% c("HFPO-DA","ADONA","PFMPA","PFMBA","NFDHA")~"PFE carboxylic acids",
      analyte %in% c("9Cl-PF3ONS","11Cl-PF3OUdS","PFEESA")~"Ether sulfonic acids",
      analyte %in% c("3:3 FTCA","5:3 FTCA","7:3 FTCA")~"FT carboxylic acids",
      analyte %in% c("FBSA","FHxSA","NFDHA","HFPODA-GenX")~"Other"
    )
  )->extracts_long

####Add extracts to blanks data and plot----

blanks_long %>% 
  mutate(
    source="Blanks"
  ) %>% 
  bind_rows(
    .,
    extracts_long
  ) %>% 
  dplyr::filter(!(analyte %in% c("5:3 FTCA","7:3 FTCA","FBSA","PFUnA"))) %>% #dplyr::filter out 100% non-detects
  ggplot()+
  geom_jitter(
    aes(
      y=source,
      x=result_corrected,
      color=source
    ),
    alpha=0.4,
    width=0,
    height=0.1
  )+
  scale_x_continuous(trans='log',labels = round2,n.breaks = 5)+
  theme_classic()+
  theme(
    axis.text.x = element_text(angle=30,vjust=0.5,color='black'),
    axis.text.y = element_text(colour='black')
    )+
  labs(
    x="Analyte concentration (ng/L)",
    y="",
    #title="HOME House Dust and Blanks",
    #caption="Note: X-axis on log scale\nRed line is LOD"
  )+
  scale_color_manual(values=brewer.pal(4,"Paired")[c(2,4)],name="")+
  facet_wrap(~analyte,scales="free")+
  geom_vline(aes(xintercept=LOD),color="#E31A1C")

ggsave("C:\\Users\\afossa\\Dropbox (Brown)\\Dust and EDC R21\\output\\blanks_and_house_extracts.tiff",width = 1920,height=1080,units = 'px',scale=2)

ggsave("C:\\Users\\afossa\\Dropbox (Brown)\\Dust and EDC R21\\output\\blanks_and_house_extracts.tiff",dpi=300,scale=1)

