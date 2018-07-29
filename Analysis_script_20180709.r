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
############# Note: because some of the data for individual calls is restricted in the public repositories, only summarized data is provided. Scripts for summarizing data will be provided upon request
pcawg.descript.SID<-as.data.table(read_delim("PCAWG.SID.Descriptive.tsv","\t", escape_double = FALSE, trim_ws = TRUE))
clin.pcawg.sith.motif<-as.data.table(read_delim("PCAWG.SID.SITH.Motif.Clinical.tsv","\t",escape_double = FALSE, trim_ws = TRUE))
cgi.descript.SID<-as.data.table(read_delim("CGI.SID.descriptive.tsv","\t", escape_double = FALSE, trim_ws = TRUE))
CGI_SID_SITH_Motifs_15kb <- as.data.table(read_delim("CGI_SIDs_SITH_Motifs_Table_dsnv_max-15kb.tsv","\t", escape_double = FALSE, trim_ws = TRUE))
cgi.summary<- as.data.table(read_delim("CGI.summarized.tsv","\t", escape_double = FALSE, trim_ws = TRUE))
pcawg.summary<-as.data.table(read_delim("PCAWG.summarized.tsv","\t", escape_double = FALSE, trim_ws = TRUE))
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

############# Generation of single SNV changes enrichment in clusters #####################
#fisher table building CGI#
# tables are row binding of c(in_cltr with motif, in_cltr without motif), c(not in cltr with motif, not in cltr without motif)
WT.MT<-c("A>C","A>G","A>T","C>A","C>G","C>T","G>A","G>C","G>T","T>A","T>C","T>G")
wt.mt<-str_replace(WT.MT,">",".")
z<-cgi.summary[,lapply(.SD[,2:37],sum, na.rm=TRUE)]

for (i in 1:12) {
k<-2*i+11
l<-2*i+12
quad1<-as.numeric(z[,..l])
quad2<-as.numeric(z[,2]-quad1)
quad3<-as.numeric(z[,..k]-quad1)
quad4<-as.numeric(z[,1]-(quad1+quad2+quad3))
mat1<-as.matrix(rbind(c(quad1,quad2),c(quad3,quad4)))
  assign(paste0(wt.mt[i],".cltr.2x2.OR"),fisher.test(mat1))
}
#####Make Supplemental table for paper#####
#WT.MT<-str_replace(wt.mt,".",">")
#WT.MT<-sort(WT.MT)
CGI.snv.fisher.table<-as.data.frame(rbind(c(A.C.cltr.2x2.OR$estimate,A.C.cltr.2x2.OR$conf.int),c(A.G.cltr.2x2.OR$estimate,A.G.cltr.2x2.OR$conf.int),c(A.T.cltr.2x2.OR$estimate,A.T.cltr.2x2.OR$conf.int),c(C.A.cltr.2x2.OR$estimate,C.A.cltr.2x2.OR$conf.int),c(C.G.cltr.2x2.OR$estimate,C.G.cltr.2x2.OR$conf.int),c(C.T.cltr.2x2.OR$estimate,C.T.cltr.2x2.OR$conf.int),c(G.A.cltr.2x2.OR$estimate,G.A.cltr.2x2.OR$conf.int),c(G.C.cltr.2x2.OR$estimate,G.C.cltr.2x2.OR$conf.int),c(G.T.cltr.2x2.OR$estimate,G.T.cltr.2x2.OR$conf.int),c(T.A.cltr.2x2.OR$estimate,T.A.cltr.2x2.OR$conf.int),c(T.C.cltr.2x2.OR$estimate,T.C.cltr.2x2.OR$conf.int),c(T.G.cltr.2x2.OR$estimate,T.G.cltr.2x2.OR$conf.int)))

CGI.snv.fisher.table<-cbind(WT.MT,CGI.snv.fisher.table)
names(CGI.snv.fisher.table)[3:4]<-c("CI_95_low","CI_95_high")

#fisher table building PCAWG#
# tables are row binding of c(in_cltr with motif, in_cltr without motif), c(not in cltr with motif, not in cltr without motif)
z1<-pcawg.summary[,lapply(.SD[,2:37],sum, na.rm=TRUE)]

for (i in 1:12) {
k<-2*i+11
l<-2*i+12
quad1<-as.numeric(z1[,..l])
quad2<-as.numeric(z1[,2]-quad1)
quad3<-as.numeric(z1[,..k]-quad1)
quad4<-as.numeric(z1[,1]-(quad1+quad2+quad3))
mat1<-as.matrix(rbind(c(quad1,quad2),c(quad3,quad4)))
  assign(paste0(wt.mt[i],".pcawg.cltr.2x2.OR"),fisher.test(mat1))
}
#####Make Supplemental table for paper#####
PCAWG.snv.fisher.table<-as.data.frame(rbind(c(A.C.pcawg.cltr.2x2.OR$estimate,A.C.pcawg.cltr.2x2.OR$conf.int),c(A.G.pcawg.cltr.2x2.OR$estimate,A.G.pcawg.cltr.2x2.OR$conf.int),c(A.T.pcawg.cltr.2x2.OR$estimate,A.T.pcawg.cltr.2x2.OR$conf.int),c(C.A.pcawg.cltr.2x2.OR$estimate,C.A.pcawg.cltr.2x2.OR$conf.int),c(C.G.pcawg.cltr.2x2.OR$estimate,C.G.pcawg.cltr.2x2.OR$conf.int),c(C.T.pcawg.cltr.2x2.OR$estimate,C.T.pcawg.cltr.2x2.OR$conf.int),c(G.A.pcawg.cltr.2x2.OR$estimate,G.A.pcawg.cltr.2x2.OR$conf.int),c(G.C.pcawg.cltr.2x2.OR$estimate,G.C.pcawg.cltr.2x2.OR$conf.int),c(G.T.pcawg.cltr.2x2.OR$estimate,G.T.pcawg.cltr.2x2.OR$conf.int),c(T.A.pcawg.cltr.2x2.OR$estimate,T.A.pcawg.cltr.2x2.OR$conf.int),c(T.C.pcawg.cltr.2x2.OR$estimate,T.C.pcawg.cltr.2x2.OR$conf.int),c(T.G.pcawg.cltr.2x2.OR$estimate,T.G.pcawg.cltr.2x2.OR$conf.int)))

PCAWG.snv.fisher.table<-cbind(WT.MT,PCAWG.snv.fisher.table)
names(PCAWG.snv.fisher.table)[3:4]<-c("CI_95_low","CI_95_high")

######## Fisher Test for Enrichment of clustered mutations in 3'-end of genes ##########
q1<-as.numeric(z[,10])
q2<-as.numeric(z[,9]-q1)
q3<-as.numeric(z[,12]-q1)
q4<-as.numeric(z[,11]-(q1+q2+q3))
cgi.in_3prime.mat<-as.matrix(rbind(c(q1,q2),c(q3,q4)))
cgi.in_3prime.OR<-fisher.test(cgi.in_3prime.mat)
cgi.in_3prime.OR

q1<-as.numeric(z1[,10])
q2<-as.numeric(z1[,9]-q1)
q3<-as.numeric(z1[,12]-q1)
q4<-as.numeric(z1[,11]-(q1+q2+q3))
pcawg.in_3prime.mat<-as.matrix(rbind(c(q1,q2),c(q3,q4)))
pcawg.in_3prime.OR<-fisher.test(pcawg.in_3prime.mat)
pcawg.in_3prime.OR

##########Linear models of motif contribution to SITH###################
cgi.sith.lm<-lm(SITH_CM_All~SITH_CM_TLS+SITH_CM_APO+SITH_CM_AID, data=CGI_SID_SITH_Motifs_15kb)
summary(cgi.sith.lm)
pcawg.sith.lm<-lm(SITH~SITH_TLS+SITH_APO+SITH_AID+mult_tumor+is_maxSITH+strata(Organ), data=clin.pcawg.sith.motif)
summary(pcawg.sith.lm)

##############build table of organ based linear models#####################
organ.table<-unique(clin.pcawg.sith.motif$Organ)
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
write.table(organ.coeff.table,"sith_linear_model_table.tsv",sep="\t",row.names = FALSE)

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
SITH_IQRgrp_coxph<-coxph(Surv(survival_time,survival_code)~IQR_group + mult_tumor + is_maxINT_IQR + strata(Organ), data=clin.pcawg.sith.motif)
summary(SITH_IQRgrp_coxph)

################################################

##################plot survival curves for median IQR curves############################
theme1<-theme(text=element_text(family="sans",size=24,colour="black"),axis.title=element_text(face="bold",size=rel(1)),axis.text=element_text(size=rel(0.583),colour="black"),legend.title=element_text(face="bold",size=rel(0.75)),legend.text=element_text(size=rel(0.583)))

fit1<-survfit(Surv(survival_time,survival_code)~IQR_group, data = clin.pcawg.sith.motif)

primary.HR.plot<-ggforest(primary.surv,data=clin.pcawg.sith.motif[sample_type=="primary"], fontsize = 1)

met.HR.plot<-ggforest(met_recur.surv,data=clin.pcawg.sith.motif[sample_type!="primary"], fontsize = 1)

SITH_IQR.HR.plot<-ggforest(SITH_IQR,data=clin.pcawg.sith.motif, fontsize = 1)

SITH_IQR.KM.plot<-ggsurvplot(fit1,data=clin.pcawg.sith.motif,risk.table=TRUE, break.time.by = 1500, xlim=c(0,9500),xlab="Time in Days", legend.title="",censor.shape=3, censor.size=4.5, palette="npg",legend.labs=c("Above Median IQR","Below Median IQR"), font.legend=12)

pdf("primary.HR.plot.pdf",width=8.5,height=11)
primary.HR.plot
dev.off()

pdf("met.HR.plot.pdf",width=8.5,height=11)
met.HR.plot
dev.off()

pdf("SITH_IQR.HR.plot.pdf",width=8.5,height=11)
SITH_IQR.HR.plot
dev.off()

pdf("primary.HR.plot.pdf",width=11,height=8.5)
SITH_IQR.KM.plot
dev.off()

##lables were edited in Adobe Acrobat for finished figures