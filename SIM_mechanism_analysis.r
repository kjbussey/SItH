rm(list = ls())

####Data access:
## Except for cell lines, only summarized data is provided because the data requires approval for restricted access from the appropriate study Data Access Committee. This is noted in the comments below. Aggregated individual level data files will be provided upon request and verification of Data Access Committee approval for access. 

####Data
library(tidyverse)
library(grid)
library(patchwork)

#####Normal
vcf.data.cltr  <- read_delim("normal_vcf.data.cltr.csv", 
                             delim = ",", 
                             escape_double = FALSE, trim_ws = TRUE)

###Normal and PCAWG
pcawg_cluster_enrich_calls <- read_delim("pcawg_cluster_enrich_calls.csv", 
                                         delim = ",", 
                                         escape_double = FALSE, trim_ws = TRUE)

#####Cell Lines
kcccg_snvs_motifs_15kb <- read_delim("kcccg_snvs_motifs_15kb.tsv", 
                                     delim = "\t", escape_double = FALSE, 
                                     trim_ws = TRUE)
kcccg_sids_sith_motifs_15kb <- read_delim("kcccg_sids_sith_motifs_15kb.tsv", 
                                          delim = "\t", escape_double = FALSE, 
                                          trim_ws = TRUE)
kcccg_sid_sith_15kb <- read_delim("kcccg_sid_sith_15kb.tsv", 
                                  delim = "\t", escape_double = FALSE, trim_ws = TRUE)
result_table <- read_delim("cell_line_fisher_results_table.tsv", 
                           delim = "\t", escape_double = FALSE, trim_ws = TRUE)


###plots



### Normal and PCAWG
norm_cancer_forest <- ggplot(data=pcawg_cluster_enrich_calls,aes(x = call,y = `odds ratio`, ymin = CI_95_high, ymax = CI_95_low, shape = type)) +
  geom_pointrange(aes(col=type)) +
  geom_hline(yintercept =1, linetype=2)+
  labs(title = "Mutational Enrichment in Clusters", x='Mutation, WT>MT', y="Odds Ratio (95% Confidence Interval)") +
  geom_errorbar(aes(ymin=CI_95_low, ymax=CI_95_high,col=type), width=0.5,cex=1)+ ylim(c(0,2.25)) +
  facet_wrap(~call,strip.position="left",nrow=9,scales = "free_y") +
  theme(plot.title=element_text(size=16,face="bold"),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.x=element_text(face="bold"),
        axis.title=element_text(size=12,face="bold"),
        strip.text.y = element_text(hjust=0,vjust = 1,angle=180,face="bold"))+
  coord_flip()


###Cell Lines

them <- theme(text         = element_text(size   = 24, 
                                          colour = "black", 
                                          face   = "bold"),
              axis.text    = element_text(size   = rel(0.75)),
              axis.title.x = element_text(size   = rel(1),
                                          margin = margin(1,0,0,0, unit = "lines")),
              axis.title.y = element_text(size   = rel(1),
                                          margin = margin(0,1,0,0, unit = "lines")),
              axis.text.y  = element_text(angle  = 90, 
                                          hjust  = 0.5),
              plot.margin  = unit(c(1, 1, 1, 1), "lines"), #(up, right, down, left)
              plot.title   = element_text(size   = rel(1), 
                                          hjust  = 0,
                                          margin = margin(0,0,1,0, unit = "lines")),
              legend.text           = element_text(size = 10)) 

#sith vs number of snvs
n_snvs_sith_plot <- ggplot(kcccg_sids_sith_motifs_15kb, 
                           aes(n_snv, sith_cm_all, colour = condition, shape = type)) + 
  geom_point(size   = 3) + labs(x = "Total Number of SNVs", y = "Overall SItH") + 
  ylim(c(0.7,0.9)) + them + ggtitle("(A)")
#tls sith
tls_sith_plot <- ggplot(kcccg_sids_sith_motifs_15kb, 
                        aes(sith_cm_tls, sith_cm_all,colour = condition, shape = type )) + 
  geom_point(size   = 3)  + labs(x = "TLS SItH", y = "Overall SItH") + 
  lims(x = c(0.7,0.9), y = c(0.7, 0.9)) + them  + ggtitle("(B)")
  #geom_abline(slope = 1, intercept = 0,linetype = "dashed") 
#apobec sith
apo_sith_plot <- ggplot(kcccg_sids_sith_motifs_15kb, 
                        aes(sith_cm_apo, sith_cm_all, colour = condition, shape = type)) + 
  geom_point(size   = 3) + labs(x = "APOBEC SItH", y = "Overall SItH") + 
  lims(x = c(0.7,0.9), y = c(0.7, 0.9)) + them + ggtitle("(C)")
  #geom_abline(slope = 1, intercept = 0,linetype = "dashed") 
#aid sith
aid_sith_plot <- ggplot(kcccg_sids_sith_motifs_15kb, 
                        aes(sith_cm_aid, sith_cm_all, colour = condition, shape = type)) + 
  geom_point(size   = 3) + labs(x = "AID SItH", y = "Overall SItH") + 
  lims(x = c(0.7,0.9), y = c(0.7, 0.9)) + them + ggtitle("(D)")
  #geom_abline(slope = 1, intercept = 0,linetype = "dashed") 
#iqr dispersion plot
iqr_plot <- ggplot(kcccg_sid_sith_15kb, aes(type,sith_iqr, fill = condition)) + 
  geom_boxplot() + labs(x = "", y = "SItH IQR") + 
  them + ylim(c(0.15,0.32)) + ggtitle("(E)")

layout <- "
ABE
CDE
"

fileout <- paste("kcccg_sith_plot.png", sep="")
png(file=fileout, width = 20 , height = 15, units = "in", res = 300, bg = "white", type = "cairo")
combined <- n_snvs_sith_plot + tls_sith_plot + apo_sith_plot + aid_sith_plot + iqr_plot + 
  plot_layout(design = layout, guides = "collect") & theme(legend.position = "bottom")
print(combined)
dev.off()

#OR forest plots
#subset data by condition
frap1kd_or_to_plot <- result_table %>% filter(condition == "FRAP1_kd") %>% select(type, call, odds_ratio, conf_int_min_95, conf_in_max_95)
tunicamycin_or_to_plot<- result_table %>% filter(condition == "tunicamycin") %>% select(type, call, odds_ratio, conf_int_min_95, conf_in_max_95)
vemurafenib_or_to_plot<- result_table %>% filter(condition == "vemurafenib") %>% select(type, call, odds_ratio, conf_int_min_95, conf_in_max_95)

### forrest plots
ggplot(data=tunicamycin_or_to_plot,aes(x = call,y = odds_ratio, ymin = conf_int_min_95, ymax = conf_in_max_95, shape = type)) +
  geom_pointrange(aes(col=type)) +
  geom_hline(yintercept =1, linetype=2)+
  labs(title = "Tunicamycin Selection: Enrichment in Clusters", x='Mutation, WT:MT', y="Odds Ratio (95% Confidence Interval)") +
  geom_errorbar(aes(ymin=conf_int_min_95, ymax=conf_in_max_95,col=type), width=0.5,cex=1)+ ylim(c(0,2)) +
  facet_wrap(~call,strip.position="left",nrow=9,scales = "free_y") +
  theme(plot.title=element_text(size=16,face="bold"),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.x=element_text(face="bold"),
        axis.title=element_text(size=12,face="bold"),
        strip.text.y = element_text(hjust=0,vjust = 1,angle=180,face="bold"))+
  coord_flip()

ggplot(data=vemurafenib_or_to_plot,aes(x = call,y = odds_ratio, ymin = conf_int_min_95, ymax = conf_in_max_95, shape = type)) +
  geom_pointrange(aes(col=type)) +
  geom_hline(yintercept =1, linetype=2)+
  labs(title = "Vemurafenib Selection: Enrichment in Clusters", x='Mutation, WT:MT', y="Odds Ratio (95% Confidence Interval)") +
  geom_errorbar(aes(ymin=conf_int_min_95, ymax=conf_in_max_95,col=type), width=0.5,cex=1)+ ylim(c(0,2)) +
  facet_wrap(~call,strip.position="left",nrow=9,scales = "free_y") +
  theme(plot.title=element_text(size=16,face="bold"),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.x=element_text(face="bold"),
        axis.title=element_text(size=12,face="bold"),
        strip.text.y = element_text(hjust=0,vjust = 1,angle=180,face="bold"))+
  coord_flip()

ggplot(data=frap1kd_or_to_plot,aes(x = call,y = odds_ratio, ymin = conf_int_min_95, ymax = conf_in_max_95, shape = type)) +
  geom_pointrange(aes(col=type)) +
  geom_hline(yintercept =1, linetype=2)+
  labs(title = "FRAP1 Knock-Down: Enrichment in Clusters", x='Mutation, WT:MT', y="Odds Ratio (95% Confidence Interval)") +
  geom_errorbar(aes(ymin=conf_int_min_95, ymax=conf_in_max_95,col=type), width=0.5,cex=1)+ ylim(c(0,2)) + 
  facet_wrap(~call,strip.position="left",nrow=9,scales = "free_y") +
  theme(plot.title=element_text(size=16,face="bold"),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.x=element_text(face="bold"),
        axis.title=element_text(size=12,face="bold"),
        strip.text.y = element_text(hjust=0,vjust = 1,angle=180,face="bold"))+
  coord_flip()


###Code for mutational signature analysis - Load the following libraries
#Based on the flow chart from Maura et al, Nat. Comm. 10:2969 (2019)
#Build SNV catalog
rm(list = ls())
source("http://bioconductor.org/biocLite.R")
library(MutationalPatterns)
library(data.table)
library(BSgenome.Hsapiens.UCSC.hg19) 
library(GenomicRanges)
library(IRanges)
library(ggplot2)
library(RColorBrewer)
library(GenomicFeatures)
library(BSgenome)
library("NMF")
library(gridExtra)
#### Because of size of data, recommend running each data set separately. Note that the code is written making this assumption
### PCAWG Data:
####NOTE: Requires individual level data. vcg.data.cltr provided upon request as indicated above
#####Analysis for clustered SNVs only######
#### Create GRanges object to serve as input into the rest of the analysis. Object includes genomic position, reference call, snv call, sample name, and study
aux <-GRanges(seqnames = vcf.data.cltr[,Chr], 
              ranges = IRanges(start = vcf.data.cltr[,Loc], 
                               end = vcf.data.cltr[,Loc]), 
              REF = vcf.data.cltr[,Wld], 
              ALT = vcf.data.cltr[,Mut], 
              sampleNames = vcf.data.cltr[,SID],
              study = vcf.data.cltr[,Organ])
genome(aux) <- "hg19"
types<-mut_type(aux)
context <- mut_context(aux,ref_genome)
type.context<- type_context(aux, ref_genome)
#### for each sample, create GRange file in a list
g_cltr <-list()
samples<-unique(vcf.data.cltr$SID)

for(i in (1:length(samples)))
{
  vcf.data.single <- vcf.data.cltr[vcf.data.cltr$SID==samples[i],]
  aux.single <- GRanges(seqnames = vcf.data.single[,Chr], 
                        ranges = IRanges(start = vcf.data.single[,Loc], 
                                         end = vcf.data.single[,Loc]), 
                        REF = vcf.data.single[,Wld], 
                        ALT = vcf.data.single[,Mut], 
                        study = vcf.data.single[,Organ])
  names(aux.single) <- vcf.data.single$SID
  genome(aux.single) <- "hg19"
  g_cltr[[i]] <- aux.single
}
names(g_cltr) <- samples
#####create the 96 mutational profile for all samples in order to do de novo extraction and then assignment to COSMIC signatures
mut_mat_cltr<- mut_matrix(vcf_list = g_cltr, ref_genome = ref_genome)
##extract 26 signatures for comparison to COSMIC and 82 reference signatures from SIGNAL
mut_mat_cltr<-mut_mat_cltr + 0.0001
nmf_res_cltr_26 <- extract_signatures(mut_mat_cltr, rank = 26, nrun = 50)
colnames(nmf_res_cltr_tls_26$signatures) <-sig.names.26 
rownames(nmf_res_cltr_tls_26$contribution) <- sig.names.26
nmf_res_cltr_26 <- rename_nmf_signatures(nmf_res_cltr_tls_26, cancer_signatures, cutoff = 0.85)
###compare de novo 26 signatures to COSMIC V2
#Import COSMIC v2 signatures
library(readr)
cancer_signatures <- read_delim("COSMIC_v2_SBS_GRCh37.txt", delim = "\t", escape_double = FALSE, trim_ws = TRUE)
View(cancer_signatures)
new_order = match(row.names(mut_mat_cltr), cancer_signatures$Type)
cancer_signatures = cancer_signatures[as.vector(new_order),]
cancer_signatures_rownames <- cancer_signatures$Type
cancer_signatures = as.matrix(cancer_signatures[,2:31])
row.names(cancer_signatures) = cancer_signatures_rownames

cos_sim_cltr26_cosmic <- cos_sim_matrix(nmf_res_cltr_26$signatures,cancer_signatures)
plot_cosine_heatmap(cos_sim_cltr_tls26_cosmic)
###compare de novo 26 to 82 ref sigs for clustered
hclust_ref82 = cluster_signatures(ref_sbs_sigs_82, method = "average")
ref82_order = colnames(ref_sbs_sigs_82)[hclust_ref82$order]
cos_sim_denovo26_cltr_ref82 <- cos_sim_matrix(nmf_res_cltr_26$signatures,ref_sbs_sigs_82)
plot_cosine_heatmap(cos_sim_denovo26_cltr_ref82, col_order = ref82_order)

###Stranded Analysis:
###strand bias analysis for denovo 26 from clustered SNVs#####
#make mutational count matrix with strand information
mut_mat_s_cltr <- mut_matrix_stranded(g_cltr, ref_genome, annote1)
#count the number of mutations on each strand
strand_counts_cltr <- strand_occurrences(mut_mat_s_cltr)
#perform poisson test for strand assymetry 
strand_bias_cltr <- strand_bias_test(strand_counts_cltr)
#plot mutation spectrum with strand distinction
ps1_cltr <- plot_strand(strand_counts_cltr, mode = "relative") + scale_alpha_discrete(range = c(1,0.4), name = "strand", breaks = c("transcribed","untranscribed"), labels = c("reference", "complement"))
# combine plots into one figure
ps2_cltr <- plot_strand_bias(strand_bias_cltr) + labs(y = "log2(reference/complement)")
grid.arrange(ps1_cltr, ps2_cltr)
#extract signatures with strand bias
nmf_res_strand_cltr_26 <- extract_signatures(mut_mat_s_cltr, rank = 26, nrun=50)
colnames(nmf_res_strand_cltr_26$signatures) <-sig.names.26
rownames(nmf_res_strand_cltr_26$contribution) <- sig.names.26
###match signatures
## collapse stranded back into 96
collapsed_cltr_stranded_26_sig <- read.csv("collapsed_cltr_stranded_26_sig.csv", row.names=1)
cos_sim_cltr_26_stranded_cltr_26 <-cos_sim_matrix(collapsed_cltr_stranded_26_sig,nmf_res_cltr_26$signatures)
cos_sim_stranded_cltr_26_cosmic <-cos_sim_matrix(collapsed_cltr_stranded_26_sig,cancer_signatures)
plot_cosine_heatmap(cos_sim_cltr_26_stranded_cltr_26)
plot_cosine_heatmap(cos_sim_stranded_cltr_26_cosmic)
sig26.names.stranded <- c("S13L","StrB","StrC","StrD","SBST","StrF","StrG","S26L","SBSC","StrJ","StrK","S28L","S17L","StrN","StrO","StrP","SBSD","StrR","SBSS","S1L","S7L","StrV","StrW","StrX","StrY","SBSG")
colnames(nmf_res_strand_cltr_26$signatures) <-sig26.names.stranded
rownames(nmf_res_strand_cltr_26$contribution) <- sig26.names.stranded
##plot for paper:
sigs_to_plot <- nmf_res_strand_cltr_26$signatures[,c(1,5,8,9,12,13,17,19:21,26)]
plot_192_profile(sigs_to_plot, condensed = TRUE) + scale_alpha_discrete(range = c(0.1,1), name = "strand", breaks = c("transcribed","untranscribed"), labels = c("reference", "complement"))
plot_signature_strand_bias(sigs_to_plot)+ labs(y = "log2(reference/complement)")

### Correlation of SITH and SITH IQR with Signature Contribution
###correlations of all 26 signatures from clustered mutations with SITH and IQR
sig_26_cltr_pivot <- sig_contrib_cltr_26 %>% pivot_longer(cols = FI660216:FI671314, names_to = "sid", values_to = "contribution") %>% glimpse()
sig_26_cltr_pivot <-sig_26_cltr_pivot %>% pivot_wider(names_from = "sig_names", values_from = "contribution") %>% glimpse()
all_sig_26_cltr_corr <- merge(sig_26_cltr_pivot, pcawg_for_sig_corr, by.x = "sid", by.y = "FileID")

#correlation table:
for (i in 2:27){
  sith_cor <- cor.test(all_sig_26_cltr_corr[,i],all_sig_26_cltr_corr$SITH,, alternative = "two.sided", method = "spearman", exact = NULL)
  iqr_cor <- cor.test(all_sig_26_cltr_corr[,i],all_sig_26_cltr_corr$INT_SITH_IQR, alternative = "two.sided", method = "spearman", exact = NULL)
  if(i == 2) 
  {spearman_results <- data.frame(Signature = names(all_sig_26_cltr_corr)[i], SITH_rho = sith_cor$estimate[[1]], SITH_pvalue = sith_cor$p.value, IQR_rho = iqr_cor$estimate[[1]], IQR_pvalue = iqr_cor$p.value)} else
  {next_row <- c(Signature = names(all_sig_26_cltr_corr)[i],SITH_rho = sith_cor$estimate[[1]], SITH_pvalue = sith_cor$p.value,IQR_rho = iqr_cor$estimate[[1]], IQR_pvalue = iqr_cor$p.value)
  spearman_results <- rbind(spearman_results,next_row)}
}
sith_p.adjust <- p.adjust(spearman_results$SITH_pvalue, method = "BY")
iqr_p.adjust <- p.adjust(spearman_results$IQR_pvalue, method = "BY")
spearman_results <- spearman_results %>% mutate(SITH_adjusted_p = p.adjust(SITH_pvalue, method = "BY"), IQR_adjusted_p = p.adjust(IQR_pvalue, method = "BY")) %>% relocate(SITH_adjusted_p, .after = SITH_pvalue) %>% glimpse()

### Normal Data:
####NOTE: Requires individual level data. vcg.data.cltr provided upon request as indicated above
###import data
read_delim(
  ##Set reference genome
  ref_genome<-"BSgenome.Hsapiens.UCSC.hg19"
  library(ref_genome, character.only = TRUE)
  ####Analysis of clustered mutations only
  #### Create GRanges object to serve as input into the rest of the analysis. Object includes genomic position, reference call, snv call, sample name, and study
  aux <-GRanges(seqnames = vcf.data.cltr[,Chrom], 
                ranges = IRanges(start = vcf.data.cltr[,Loc], 
                                 end = vcf.data.cltr[,Loc]), 
                REF = vcf.data.cltr[,Wild_sequence], 
                ALT = vcf.data.cltr[,Mutant_sequence], 
                sampleNames = vcf.data.cltr[,SID],
                study = "CGI")
  genome(aux) = "hg19"
  types<-mut_type(aux)
  context <- mut_context(aux,ref_genome)
  type.context<- type_context(aux, ref_genome)
  #### for each sample, create GRange file in a list
  g_cltr <-list()
  samples<-unique(vcf.data.cltr$SID)
  #studies<-unique(vcf.data$Organ)
  for(i in (1:length(samples)))
  {
    vcf.data.single <- vcf.data.cltr[vcf.data.cltr$SID==samples[i],]
    aux.single <- GRanges(seqnames = vcf.data.single[,Chrom], 
                          ranges = IRanges(start = vcf.data.single[,Loc], 
                                           end = vcf.data.single[,Loc]), 
                          REF = vcf.data.single[,Wild_sequence], 
                          ALT = vcf.data.single[,Mutant_sequence], 
                          study = "CGI")
    names(aux.single) <- vcf.data.single$SID
    genome(aux.single) <- "hg19"
    g_cltr[[i]] <- aux.single
  }
  names(g_cltr) <- samples
  #####compute and plot 6 classes and CpG prevalence for samples
  type_occurrences <- mut_type_occurrences(g_cltr,ref_genome) #computes spectrum for all samples in the GRanges list in a single graph
  
  #####create the 96 mutational profile for all samples in order to do de novo extraction and then assignment to COSMIC signatures
  mut_mat_cltr<- mut_matrix(vcf_list = g_cltr, ref_genome = ref_genome)
  ##Extract 26 signatures to compare to COSMIC version 2. Using 26 because it is what worked best in PCAWG data
  mut_mat_cltr<-mut_mat_cltr + 0.0001 #prevent dividing by 0
  cgi_nmf_res_cltr <- extract_signatures(mut_mat_cltr, rank = 26, nrun = 50)
  colnames(cgi_nmf_res_cltr$signatures) <-sig.names.26 
  rownames(cgi_nmf_res_cltr$contribution) <- sig.names.26
  cgi_nmf_res_cltr <- rename_nmf_signatures(cgi_nmf_res_cltr, cancer_signatures, cutoff = 0.80)
  cos_sim_denovo26_cosmic_cltr <- cos_sim_matrix(cgi_nmf_res_cltr$signatures, cancer_signatures)
  plot_cosine_heatmap(cos_sim_denovo26_cosmic_cltr, col_order = cosmic_order, cluster_rows = TRUE)
  
  ####Compare to PCAWG 26 de novo signatures
  pcawg_cltr_sig26 <- read_csv("cluster_signatures26.csv")
  names(pcawg_cltr_sig26)[1] <- "type"
  pcawg_cltr_sig26 <- pcawg_cltr_sig26 %>% column_to_rownames(var = "type") 
  cos_sim_denovo26_pcawg26_cltr <- cos_sim_matrix(cgi_nmf_res_cltr$signatures, pcawg_cltr_sig26)
  plot_cosine_heatmap(cos_sim_denovo26_pcawg26_cltr, cluster_rows = TRUE, cluster_cols = TRUE)
  
  ### KCCG Cell Line Data:
  read_delim(
    ###Prep initial data
    vcf.data <- data.table(kcccg_snvs_motifs_15kb) 
    vcf.data.cltr <- vcf.data[in_cltr == TRUE]
    vcf.data.noCID <- vcf.data[in_cltr == FALSE]
    ##Set reference genome
    ref_genome<-"BSgenome.Hsapiens.UCSC.hg19"
    library(ref_genome, character.only = TRUE)
    
    ####Analysis for clustered SNVs only######
    #### Create GRanges object to serve as input into the rest of the analysis. Object includes genomic position, reference call, snv call, sample name, and study
    aux <-GRanges(seqnames = vcf.data.cltr[,chrom], 
                  ranges = IRanges(start = vcf.data.cltr[,loc], 
                                   end = vcf.data.cltr[,loc]), 
                  REF = vcf.data.cltr[,wild_sequence], 
                  ALT = vcf.data.cltr[,mutant_sequence], 
                  sampleNames = vcf.data.cltr[,sid],
                  histology = vcf.data.cltr[,histology],
                  type = vcf.data.cltr[,type],
                  condition = vcf.data.cltr[,condition])
    types<-mut_type(aux)
    context <- mut_context(aux,ref_genome)
    type.context<- type_context(aux, ref_genome)
    #### for each sample, create GRange file in a list
    g <-list()
    samples<-unique(vcf.data.cltr$sid)
    
    for(i in (1:length(samples)))
    {
      vcf.data.single <- vcf.data.cltr[vcf.data.cltr$sid==samples[i],]
      aux.single <- GRanges(seqnames = vcf.data.single[,chrom], 
                            ranges = IRanges(start = vcf.data.single[,loc], 
                                             end = vcf.data.single[,loc]), 
                            REF = vcf.data.single[,wild_sequence], 
                            ALT = vcf.data.single[,mutant_sequence], 
                            study = vcf.data.single[,histology],
                            type = vcf.data.single[,type],
                            condition = vcf.data.single[,condition])
      names(aux.single) <- vcf.data.single$sid
      genome(aux.single) <- "hg19"
      g[[i]] <- aux.single
    }
    names(g) <- samples
    #####compute and plot 6 classes and CpG prevalence for samples and plot the spectrum of all 6 classes by organ/tissue separating out CpG from other C>T contexts
    type_occurrences <- mut_type_occurrences(g,ref_genome)
    p4.cltr <- plot_spectrum(type_occurrences, CT = TRUE, legend = TRUE)
    
    
    
