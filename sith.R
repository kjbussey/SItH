#########################################################################################################
#  Plots basic SItH Statistics as a function of mutational load ( Total [SNVs]): 
#    - Overal SItH, 
#    - Cluster SItH IQR, 
#
#########################################################################################################

rm(list = ls())

library(ggplot2)
library(data.table);
library(grid)
library(scales)
library(readr)
library(rstudioapi) 
setwd(dirname(getActiveDocumentContext()$path)) 

load("sith.RData") 

#########################################################################################################

S_N_scatter      <- ggplot()

if (ca1_flag) {
  S_N_scatter <- S_N_scatter + geom_point( data       = ca1.stats.tbl, 
              aes(x      = N_SNV, 
                  y      = SITH,
                  colour = "ca1"),
              shape  = 1,
              size   = 3)
}

if (cgi_flag) {
  S_N_scatter <- S_N_scatter + geom_point( data       = cgi.stats.tbl, 
              aes(x      = N_SNV, 
                  y      = SITH,
                  colour = "cgi"),
              shape  = 1,
              size   = 3)
}

if (sim_flag) {
  S_N_scatter <- S_N_scatter + geom_line( data       = sim.tbl, 
             aes(x      = N_SNV, 
                 y      = Av_SITH,
                 colour = "sim"),
             size       = 1,
             linetype   = 2) +
  geom_errorbar( data        = sim.tbl,
                 aes( x      = N_SNV,
                      ymin   = Av_SITH - Sd_SITH, 
                      ymax   = Av_SITH + Sd_SITH,
                      colour = "sim"),
                 size        = 1
  )
}

S_N_scatter <- S_N_scatter +
  coord_cartesian(xlim      = c(500, 500000),
                  ylim      = c(0.1, 1.05),
                  expand    = FALSE
  ) +
  scale_x_log10(breaks       = logbreaks,
                minor_breaks = logminor_breaks,
                labels       = comma
  ) + 
  scale_y_continuous(breaks       = seq(0,2,.2),
                     minor_breaks = seq(0,2,.05),
                     labels       = comma
  ) + 
  scale_colour_manual(name   = element_blank(),
                      values = cols,
                      labels = labs,
                      breaks = nams
  ) +
  theme(text         = element_text(size   = 24, 
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
        legend.title          = element_blank(),
        legend.position       = c(0.99, 0.99),
        legend.justification  = c(1, 1),
        legend.text           = element_text(size = 10),
        legend.box.just       = "right", 
        legend.box.background = element_rect()
  ) + 
  labs( x = "Total Number of SNVs",
        y = "Overal SItH"
  ) 

print(S_N_scatter)

#########################################################################################################

IQR_N_scatter  <- ggplot() 

if (ca1_flag) {
  IQR_N_scatter <-IQR_N_scatter + geom_point( data   = ca1.stats.tbl, 
              aes(x  = N_SNV, 
                  y  = INT_SITH_iqr,
                  colour = "ca1"),
              shape  = 1,
              size   = 3) 
}

if (cgi_flag) {
  IQR_N_scatter <-IQR_N_scatter + geom_point( data   = cgi.stats.tbl, 
              aes(x  = N_SNV, 
                  y  = INT_SITH_iqr,
                  colour = "cgi"),
              shape  = 1,
              size   = 3)
}

if (sim_flag) {
  IQR_N_scatter <-IQR_N_scatter + geom_line( data   = sim.tbl, 
             aes(x  = N_SNV, 
                 y  = Av_INT_SITH_iqr,
                 colour = "sim"),
             size   = 1,
             linetype = 2) + 
  geom_errorbar( data   = sim.tbl,
                 aes( x    = N_SNV,
                      ymin = Av_INT_SITH_iqr - Sd_INT_SITH_iqr, 
                      ymax = Av_INT_SITH_iqr + Sd_INT_SITH_iqr,
                      colour = "sim"),
                 size   = 1
  )
}

IQR_N_scatter <-IQR_N_scatter + 
  coord_cartesian(xlim      = c(500, 500000),
                  ylim      = c(-0.02, 0.5),
                  expand    = FALSE
  ) +
  scale_x_log10(breaks       = logbreaks,
                minor_breaks = logminor_breaks,
                labels       = comma
  ) + 
  scale_y_continuous(breaks       = seq(0,2,.2),
                     minor_breaks = seq(0,2,.05),
                     labels       = comma
  ) + 
  scale_colour_manual(name   = element_blank(),
                      values = cols,
                      labels = labs,
                      breaks = nams
  ) +
  theme(text         = element_text(size   = 24, 
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
        legend.title          = element_blank(),
        legend.position       = c(0.99, 0.01),
        legend.justification  = c(1, 0),
        legend.text           = element_text(size = 10),
        legend.box.just       = "right", 
        legend.box.background = element_rect()
  ) + 
  labs( x = "Total Number of SNVs",
        y = "Cluster SItH IQR"
  ) 

print(IQR_N_scatter)



#########################################################################################################

fileout <- paste("fig_sith_1A.tiff", sep="")
tiff(file=fileout, width = 7.5 , height = 7.5, units = "in", res = 300, bg = "white", type = "cairo")
print(S_N_scatter)
dev.off()

fileout <- paste("fig_sith_1B.tiff", sep="")
tiff(file=fileout, width = 7.5 , height = 7.5, units = "in", res = 300, bg = "white", type = "cairo")
print(IQR_N_scatter)
dev.off()



#########################################################################################################

