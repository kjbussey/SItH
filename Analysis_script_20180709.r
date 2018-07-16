####Needed Packages####
library(ggExtra)
library(ggplot2)
library(gridExtra)
library(Hmisc)
library(stats)
library(stats4)
library(MASS)
library(gdata)
library(plyr)
library(dplyr)
library(scales)
library(Biobase)
library(data.table)
library(kSamples)
library(survival)
library(survminer)
library(readr)
library(rlist)
library(stringr)
#################

###############load data###############
pcawg.descript.SID<-as.data.table(read_delim("~/PCAWG.SID.Descriptive.tsv","\t", escape_double = FALSE, trim_ws = TRUE))
clin.pcawg.sith.motif<-as.data.table(read_delim("~/PCAWG.SID.SITH.Motif.Clinical.tsv",sep="\t",escape_double = FALSE, trim_ws = TRUE))
cgi.descript.SID<-as.data.table(read_delim("~/CGI.SID.descriptive.tsv",sep="\t", escape_double = FALSE, trim_ws = TRUE))
CGI_SID_SITH_Motifs_15kb <- as.data.table(read_delim("~/CGI_SIDs_SITH_Motifs_Table_dsnv_max-15kb.tsv","\t", escape_double = FALSE, trim_ws = TRUE))

#################

##############Descriptive Statistics for motifs in clusters#####################
fivenum(pcawg.descript.SID$fraction_in_cltr_AID)
fivenum(pcawg.descript.SID$fraction_in_cltr_APO)
fivenum(pcawg.descript.SID$fraction_in_cltr)
fivenum(pcawg.descript.SID$fraction_in_cltr_TLS)
fivenum(cgi.descript.SID$fraction_in_cltr)
fivenum(cgi.descript.SID$fraction_in_cltr_TLS)
fivenum(cgi.descript.SID$fraction_APO_in_cltr)
fivenum(cgi.descript.SID$fraction_AID_in_cltr)

#################

##########Linear models of motif contribution to SITH###################
cgi.sith.lm<-lm(SITH_CM_All~SITH_CM_TLS+SITH_CM_APO+SITH_CM_AID, data=CGI_SID_SITH_Motifs_15kb)
summary(cgi.sith.lm)
pcawg.sith.lm<-lm(SITH~SITH_TLS+SITH_APO+SITH_AID+mult_tumor+is_maxSITH+strata(Organ), data=clin.pcawg.sith.motif)
summary(pcawg.sith.lm)

##############build table of organ based linear models#####################
organ.table<-unique(clinical$Organ)
organ.lms<-vector("list",length(organ.table))
for (i in 1:length(organ.table)) {
  x<-lm(SITH~SITH_TLS+SITH_APO+SITH_AID+mult_tumor+is_maxSITH, data=clin.pcawg.sith.motif[Organ==organ.table[i]])
  organ.lms[[i]]<-x
  names(organ.lms)[i]<-paste0(organ.table[i],".sith.lm")
}
organ.lm.summary<-lapply(organ.lms,summary)
for (i in 1:length(organ.table)) {
  if (i==1) {
    organ.coeff.table<-organ.lm.summary[[i]]$coefficients} else
      {organ.coeff.table<-rbind(organ.coeff.table,organ.lm.summary[[i]]$coefficients)}
}
organ.coeff.table<-data.frame(organ.coeff.table)
setDT(organ.coeff.table,keep.rownames = TRUE)
setnames(organ.coeff.table,"rn","Variable")
organ.coeff.table[grep("mult_tumor",organ.coeff.table$Variable),Estimate:=NA]
organ.coeff.table[grep("is_maxSITH",organ.coeff.table$Variable),Estimate:=NA]
organ.coeff.table<-organ.coeff.table[!is.na(Estimate),]
organ.names<-rep(organ.table,each=4)
organ.coeff.table<-data.table(cbind(organ.names,organ.coeff.table))
organ.coeff.table[grep("X.Intercept",organ.coeff.table$Variable),Variable:="X.Intercept"] 
organ.coeff.table[grep("SITH_TLS",organ.coeff.table$Variable),Variable:="SITH_TLS"] 
organ.coeff.table[grep("SITH_APO",organ.coeff.table$Variable),Variable:="SITH_APOBEC"] 
organ.coeff.table[grep("SITH_AID",organ.coeff.table$Variable),Variable:="SITH_AID"] 
organ.coeff.table[,adjusted_p:=p.adjust(organ.coeff.table$Pr...t..,method="BH")]
write.table(organ.coeff.table,"~/sith_linear_model_table.tsv",sep="\t",row.names = FALSE)

##########################################

############ANOVA for relevant clinical variables##################

variable.aov<-aov(SITH~Organ+mult_tumor+is_maxSITH, data = clin.pcawg.sith.motif)

############Cox proportional hazard analysis for survival#################
######for primary tumors:
primary.surv<-coxph(Surv(survival_time,survival_code)~SITH+mult_tumor+is_maxSITH+strata(Organ),data=clin.pcawg.sith.motif[sample_type=="primary",])
summary(primary.surv)
######for mets and recurrences
met_recur.surv<-coxph(Surv(survival_time,survival_code)~SITH+mult_tumor+is_maxSITH+strata(Organ),data=clin.pcawg.sith.motif[sample_type!="primary",])
summary(met_recur.surv)
######for IQR
SITH_IQR<-coxph(Surv(survival_time,survival_code)~INT_IQR+mult_tumor+is_maxINT_IQR+strata(Organ),data=clin.pcawg.sith.motif)
summary(SITH_IQR)
######for median IQR grouping
clin.pcawg.sith.motif[,IQR_group:=ifelse(INT_IQR<median(INT_IQR,na.rm=TRUE),1,2)]
SITH_IQRgrp_coxph<-coxph(Surv(survival_time,survival_code)~IQR_group + mult_tumor + is_maxINT_IQR + strata(Organ), data=clin.pcawg.sith.motif)
summary(SITH_IQRgrp_coxph)

################################################

##################plot survival curves for median IQR curves############################
theme1<-theme(text=element_text(family="sans",size=24,colour="black"),axis.title=element_text(face="bold",size=rel(1)),axis.text=element_text(size=rel(0.583),colour="black"),legend.title=element_text(face="bold",size=rel(0.75)),legend.text=element_text(size=rel(0.583)))

IQR_group_15kb_plot<-ggadjustedcurves(SITH_IQRgrp_coxph, data=clin.pcawg.sith.motif, variable = "IQR_group", method="conditional") + theme1 + scale_color_discrete(name= "",labels=c("Below the Median SItH IQR","Above the Median SItH IQR"))+labs(colour="IQR group", x="\nTime in Days", y="Survival Proportion\n")+annotate("text", x=5000, y=0.80, size=5, label = "HR = 1.26, p-value = 0.017") + geom_vline(xintercept = 1826.25, linetype="dashed")+annotate("text", x=2500, y=1.0, size = 4.5, label = "5 years")

IQR_group_15kb_plot
