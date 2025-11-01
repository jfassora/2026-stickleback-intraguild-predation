#--------------------------------------------------------------------------------------------------------------------------------------------
# Analysis script for Fassora, Martin, Matthews. "Multivariate morphological divergence due to intraguild predation". 2026.
# 
# For further questions contact Jonathan Fassora: jonathanf@unis.no / https://github.com/jfassora
# 
# The script is organised in blocks delimited by #----------, using Rstudio, you can close each block to make reading easier
#                                                                           (click on the down arrow next to the line number)
#
# In this script, model fitting commands are commented out as we read the saved model outputs instead,
# to document and replicate the exact results discussed in the paper.
# These model outputs can be found in the github repo for the paper (where this script is also found) as a split zip archive.
#--------------------------------------------------------------------------------------------------------------------------------------------

setwd()

library(tidyverse)
library(ggpubr)
library(patchwork)
library(readxl)
library(lme4)
library(rstan)
library(cmdstanr)
library(loo)
library(shinystan)
library(brms)
library(tidybayes)
library(reshape2)
library(dplyr)
library(tidyr)
library(ggplot2)
library(cowplot)
library(scales)
library(igraph)
library(qgraph)
library(grid)

clrs = c("#02887f","#A30502") # color palette for figures

# Read data and do some cleaning
#--------------------------------------------------------------------------------------------------------------------------------------------
fish = read_csv("fish_traits_pub.csv")
lakes = read_csv("lake_info_pub.csv")
lakes$Area_m2 = lakes$Area_ha * 10000
lakes$area_perim_ratio = lakes$Area_m2/lakes$Perimeter_m

# for the paper we use more logical site names
# here's the mapping to the names we use internally and in the csv files
display_sites = c( 
  "AK2" = "A01",
  "AP3" = "A02",
  "AK11" = "A03",
  "AK10" = "A04",
  "AK9" = "A05",
  "AK4" = "A06",
  "AL4" = "A07",
  "AK8" = "A08",
  "AK7" = "A09",
  "L26" = "Q01",
  "EFL2" = "Q02",
  "Beach_lake" = "Q3",
  "ERL125" = "Q04",
  "ERL124" = "Q05",
  "ERL084" = "Q06",
  "ERL153" = "Q07",
  "ERL006" = "Q08",
  "ERL033" = "Q09",
  "ERL137" = "Q10",
  "ERL138" = "Q11",
  "ERL067" = "Q12",
  "ERL001" = "Q13",
  "ERL016" = "Q14",
  "ERL017" = "Q15",
  "Cedar_lake" = "Q16",
  "TULX10" = "T01",
  "TULX11" = "T02",
  "TULX8" = "T03",
  "TULX5" = "T04",
  "TULX6" = "T05",
  "TUL10" = "T06",
  "TUL8" = "T07",
  "TUL6" = "T08",
  "TUPX1" = "T09",
  "TUPX2" = "T10",
  "TUP19" = "T11"
)
#--------------------------------------------------------------------------------------------------------------------------------------------


# Calculate biomechanical traits
#--------------------------------------------------------------------------------------------------------------------------------------------
# suction index (pi*... part is the cross sectional area, approximated like this since we only have length measures from the side)
fish = mutate(fish, Suction_index=(pi*(Epax_height_open*EpaxWidth/2)/2)
                    *(Epax_height_open/Neuro_Outlever)/(Gape*Buccal_length))

# kinematic transmission code from https://academic.oup.com/evolut/article/71/11/2738/6725984#402569322
#              "CalculateKTCode.R" https://datadryad.org/stash/dataset/doi:10.5061/dryad.5b1k0
# slightly modified on the return step (as we only need the KT value)
calc_KT_operc = function(Fixed, In, Coupler, Out, Diagonal, rotation = 5){ # calculated at a given aperture (default 5°)
  Lf = 1
  Li = In/Fixed
  Lc = Coupler/Fixed
  Lo = Out/Fixed
  Ld = Diagonal/Fixed
  std.shape = data.frame(Lf, Li, Lc, Lo, Ld)
  # All angles in radians
  rotation_rad = rotation*pi/180
  
  # Test for a valid 4bar:
  if(Ld < Lc + Lo & Ld < Lf + Li) valid.start = T else valid.start = F
  
  if(valid.start){
    theta_in_start = acos((Li^2 + 1 - Ld^2)/(2*Li*1))
    theta_1_start =acos((1 + Ld^2 - Li^2)/(2*1*Ld))
    theta_2_start =acos((Lo^2 + Ld^2 - Lc^2)/(2*Lo*Ld))
    theta_out_start = theta_1_start + theta_2_start
  } else theta_out_start = NA
  
  # Calculations for ending output angle
  if(is.na(theta_out_start) == F ){
    theta_in_end = theta_in_start + rotation_rad
    Ld_end = sqrt(1 + Li^2 - 2*1*Li*cos(theta_in_end))
    if(Ld_end < Lc + Lo & Ld_end < Lf + Li) valid.end = T else valid.end = F
    if(valid.end){
      theta_1_end =acos((1 + Ld_end^2 - Li^2)/(2*1*Ld_end))
      theta_2_end =acos((Lo^2 + Ld_end^2 - Lc^2)/(2*Lo*Ld_end))
      theta_out_end = theta_1_end + theta_2_end
    } else theta_out_end = NA
    
    # calculate kt
    kt = abs(theta_out_end - theta_out_start)/rotation_rad
    
    return(kt)
  } else  return(NA)
}
fish = mutate(fish, KT=mapply(fish$Fixed, fish$Input, fish$Coupler, fish$Output, fish$Diagonal, FUN=calc_KT_operc))

# lever ratio
fish = mutate(fish, Lever_ratio=Jaw_Outlever/Jaw_Inlever)
#--------------------------------------------------------------------------------------------------------------------------------------------


# Trait correction + scaling
#--------------------------------------------------------------------------------------------------------------------------------------------
# calculate mean raker length and spacing
fish$Mean_raker_length = (fish$Raker_1_length + fish$Raker_2_length)/2
fish$Mean_spacing = (fish$Spacing_1 + fish$Spacing_2)/2

# trait correction following Paccard (2020), originally from Lleonart (2000)
Trait_list = c("Spine_length_dorsal_1", "Spine_length_dorsal_2", "Spine_length_pelvic",
                "Pelvic_height", "Pelvic_length",
                "Mean_raker_length", "Mean_spacing")
SL_mean = mean(fish$SLF)
O = lapply(Trait_list, function(trait) {
  X1 = fish[, c("SLF","Site","FishEc", trait)]
  colnames(X1) = c("SLF", "Site","FishEc", "Trait")
  X2 = X1 %>% 
    mutate(log_SLF = log10(SLF), 
           log_Trait = log10(Trait)) %>% 
    na.omit()
  mixed = lmer(log_Trait ~ log_SLF + (1 + log_SLF | Site ), data = X2,
                control=lmerControl(optimizer="bobyqa"))
  Slopes = as.data.frame(coef(mixed)$Site); Slopes$Site = rownames(Slopes); colnames(Slopes) = c("y","Slope","Site")
  X3 = X1 %>% right_join(Slopes, by = "Site")
  X3$Traits_cor =  X3$Trait*(SL_mean / X3$SLF)^ X3$Slope
  colnames(X3) = c("SLF","Site","FishEc", "Trait","y","Slope", paste(trait, "cor", sep="_"))
  X4 = X3 %>% dplyr::select(7)
  return(X4)
})
TraitsCor = as_tibble(do.call(cbind, O))
fish = cbind(fish, TraitsCor)

# trait correction of jaw protrusion and gape size on the head rather than body
# see methods section and Moosmann et al. (2025) for reasoning
Trait_list_head = c("Gape","Jaw_prot")
HL_mean = mean(fish$Head_length)
O = lapply(Trait_list_head, function(trait) {
  X1 = fish[, c("Head_length","Site","FishEc", trait)]
  colnames(X1) = c("Head_length", "Site","FishEc", "Trait")
  X2 = X1 %>%
    mutate(log_Head = log10(Head_length),
           log_Trait = log10(Trait)) %>%
    na.omit()
  mixed = lmer(log_Trait ~ log_Head + (1 + log_Head | Site ), data = X2,
                         control=lmerControl(optimizer="bobyqa"))
  Slopes = as.data.frame(coef(mixed)$Site); Slopes$Site = rownames(Slopes); colnames(Slopes) = c("y","Slope","Site")
  X3 = X1 %>% right_join(Slopes, by = "Site")
  X3$Traits_cor =  X3$Trait*(HL_mean / X3$Head_length)^ X3$Slope
  colnames(X3) = c("Head_length","Site","FishEc", "Trait","y","Slope", paste(trait, "cor_head", sep="_"))
  X4 = X3 %>% dplyr::select(7)
  return(X4)
})
TraitsCorHead = as_tibble(do.call(cbind, O))
fish = cbind(fish, TraitsCorHead)
TraitsCor = cbind(TraitsCor, TraitsCorHead)
#--------------------------------------------------------------------------------------------------------------------------------------------


# Final dataset
#--------------------------------------------------------------------------------------------------------------------------------------------
# prep data frame
fish = left_join(fish, lakes[,c("Lake","Char","area_perim_ratio","Drainage","Region")],
                  by=c("Site"="Lake"))
fish$FishEc = as.factor(fish$FishEc)
fish$Char = as.factor(fish$Char)
fish$Site = as.factor(fish$Site)
fish$Drainage = as.factor(fish$Drainage)
fish$Region = as.factor(fish$Region)

# remove redundant linear measures (final trait set)
trait_set = c("Spine_length_dorsal_1_cor", "Spine_length_dorsal_2_cor", "Spine_length_pelvic_cor",
              "Pelvic_height_cor", "Pelvic_length_cor",
              "Gape_cor_head", "Jaw_prot_cor_head",
              "Mean_raker_length_cor", "Mean_spacing_cor",
              "Suction_index", "KT", "Lever_ratio",
              "Plate_count", "Raker_count")

fish.scaled = fish[, c("FishEc", "Drainage", "Site", "area_perim_ratio", "Char", "SLF", trait_set)]
fish.scaled[,7:18] = scale(fish.scaled[,7:18])
fish.scaled$area_perim_ratio_scaled = scale(fish.scaled$area_perim_ratio)
fish.scaled$olre = 1:nrow(fish.scaled) #observation-level random effect to capture overdispersion (Poisson)
#--------------------------------------------------------------------------------------------------------------------------------------------


# Single trait comparisons
#--------------------------------------------------------------------------------------------------------------------------------------------
### Models
# this commented chunk (within the brackets) fits the models, de-comment if you wish to rerun them
# in the line following the chunk (list2env) we instead read the saved output
{

prior =
  c(prior("normal(0,1)", class = "Intercept"),
    prior("normal(0,1)", class = "b"),
    prior("exponential(2)", class = "sd"),
    prior("exponential(2)", class = "sigma"))

prior.pois =
  c(prior("normal(0,1)", class = "Intercept"),
    prior("normal(0,1)", class = "b"),
    prior("exponential(2)", class = "sd"))

### Defense
# mod.spine1 = brm(Spine_length_dorsal_1_cor ~ Char + (1|Drainage) + (1|Site), prior = prior,
#                   fish.scaled, iter=3000, control=list(adapt_delta=0.99), seed=1234, cores=4)
# mod.spine2 = brm(Spine_length_dorsal_2_cor ~ Char + (1|Drainage) + (1|Site), prior = prior,
#                   fish.scaled, iter=3000, control=list(adapt_delta=0.99), seed=1234, cores=4)
# mod.spinep = brm(Spine_length_pelvic_cor ~ Char + (1|Drainage) + (1|Site), prior = prior,
#                   fish.scaled, iter=3000, control=list(adapt_delta=0.99), seed=1234, cores=4)
# mod.pelvih = brm(Pelvic_height_cor ~ Char + (1|Drainage) + (1|Site), prior = prior,
#                   fish.scaled, iter=3000, control=list(adapt_delta=0.99), seed=1234, cores=4)
# mod.pelvil = brm(Pelvic_length_cor ~ Char + (1|Drainage) + (1|Site), prior = prior,
#                   fish.scaled, iter=3000, control=list(adapt_delta=0.99), seed=1234, cores=4)
# note that in the poisson model we add (1|olre), observation-level random effect, to soak up overdispersion
# and that in the gaussian model we scale the count
# mod.plates = brm(Plate_count ~ Char + (1|Drainage) + (1|Site) + (1|olre), prior = prior.pois,
#                   fish.scaled, family=poisson(link="log"),
#                   iter=3000, control=list(adapt_delta=0.99), seed=1234, cores=4)
# mod.platesg = brm(scale(Plate_count) ~ Char + (1|Drainage) + (1|Site), prior = prior,
#                   fish.scaled, iter=3000, control=list(adapt_delta=0.99), seed=1234, cores=4)

### Foraging (compared to models for defensive traits, here we add area:perimeter ratio as explanatory variable)
# mod.gape = brm(Gape_cor_head ~ Char*area_perim_ratio_scaled + (1|Drainage) + (1|Site), prior = prior,
#                 fish.scaled, iter=3000, control=list(adapt_delta=0.99), seed=1234, cores=4)
# mod.jp = brm(Jaw_prot_cor_head ~ Char*area_perim_ratio_scaled + (1|Drainage) + (1|Site), prior = prior,
#               fish.scaled, iter=3000, control=list(adapt_delta=0.99), seed=1234, cores=4)
# mod.grl = brm(Mean_raker_length_cor ~ Char*area_perim_ratio_scaled + (1|Drainage) + (1|Site), prior = prior,
#                 fish.scaled, iter=3000, control=list(adapt_delta=0.99), seed=1234, cores=4)
# mod.spacing = brm(Mean_spacing_cor ~ Char*area_perim_ratio_scaled + (1|Drainage) + (1|Site), prior = prior,
#                              fish.scaled, iter=3000, control=list(adapt_delta=0.99), seed=1234, cores=4)
# mod.rakcou = brm(Raker_count ~ Char*area_perim_ratio_scaled + (1|Drainage) + (1|Site) + (1|olre), prior = prior.pois,
#                    fish.scaled, family=poisson(link="log"),
#                    iter=3000, control=list(adapt_delta=0.99), seed=1234, cores=4)
# mod.rakcoug = brm(scale(Raker_count) ~ Char*area_perim_ratio_scaled + (1|Drainage) + (1|Site), prior = prior,
#                   fish.scaled, iter=3000, control=list(adapt_delta=0.99), seed=1234, cores=4)
# mod.si = brm(Suction_index ~ Char*area_perim_ratio_scaled + (1|Drainage) + (1|Site), prior = prior,
#                fish.scaled, iter=3000, control=list(adapt_delta=0.99), seed=1234, cores=4)
# mod.kt = brm(KT ~ Char*area_perim_ratio_scaled + (1|Drainage) + (1|Site), prior = prior,
#               fish.scaled, iter=3000, control=list(adapt_delta=0.99), seed=1234, cores=4)
# mod.lr = brm(Lever_ratio ~ Char*area_perim_ratio_scaled + (1|Drainage) + (1|Site), prior = prior,
#                fish.scaled, iter=3000, control=list(adapt_delta=0.99), seed=1234, cores=4)

# trait_models = list("mod.spine1" = mod.spine1,
#                       "mod.spine2" = mod.spine2,
#                       "mod.spinep" = mod.spinep,
#                       "mod.pelvih" = mod.pelvih,
#                       "mod.pelvil" = mod.pelvil,
#                       "mod.plates" = mod.plates,
#                       "mod.platesg" = mod.platesg,
#                       "mod.gape" = mod.gape,
#                       "mod.jp" = mod.jp,
#                       "mod.grl" = mod.grl,
#                       "mod.spacing" = mod.spacing,
#                       "mod.rakcou" = mod.rakcou,
#                       "mod.rakcoug" = mod.rakcoug,
#                       "mod.si" = mod.si,
#                       "mod.kt" = mod.kt,
#                       "mod.lr" = mod.lr)
# 
# saveRDS(trait_models, "trait_models.rds")
}
list2env(readRDS("trait_models.rds"), envir=globalenv())


ggplot(data=fish.scaled) + geom_histogram(aes(x=Raker_count, fill=Char), position="identity", alpha=.3)



# posterior prediction plots (diagnostic)
# Figure S1
ggsave(
   filename="figure_pp_plots.png", width=10, height=10,
   plot=ggarrange(
     ggarrange(pp_check(mod.spine1) + ggtitle("Dorsal spine 1"),
               pp_check(mod.spine2) + ggtitle("Dorsal spine 2"),
               pp_check(mod.spinep) + ggtitle("Pelvic spine"),
               nrow=1),
     ggarrange(pp_check(mod.pelvih) + ggtitle("Ascending process height"),
               pp_check(mod.pelvil) + ggtitle("Ascending process width"),
               pp_check(mod.gape) + ggtitle("Gape size"),
               nrow=1),
     ggarrange(pp_check(mod.jp) + ggtitle("Jaw protrusion"),
               pp_check(mod.grl) + ggtitle("Gill raker length"),
               pp_check(mod.spacing) + ggtitle("Gill raker spacing"),
               nrow=1),
     ggarrange(pp_check(mod.si) + ggtitle("Suction index"),
               pp_check(mod.kt) + ggtitle("Kinematic transmission"),
               pp_check(mod.lr) + ggtitle("Lever ratio"),
               nrow=1),
     ggarrange(pp_check(mod.plates) + ggtitle("Plate number (Poisson)"),
               pp_check(mod.platesg) + ggtitle("Plate number (Gaussian)"),
               pp_check(mod.rakcou) + ggtitle("Raker number (Poisson)"),
               pp_check(mod.rakcoug) + ggtitle("Raker number (Gaussian)"),
               nrow=1),
     nrow=5)
 )

# summarise results (reported in Table 2)
hypothesis(mod.spine1, "Char1 > 0", robust = TRUE)
hypothesis(mod.spine2, "Char1 > 0", robust = TRUE)
hypothesis(mod.spinep, "Char1 > 0", robust = TRUE)
hypothesis(mod.pelvih, "Char1 > 0", robust = TRUE)
hypothesis(mod.pelvil, "Char1 > 0", robust = TRUE)
hypothesis(mod.platesg, "Char1 < 0", robust = TRUE)

hypothesis(mod.gape, "Char1 < 0", robust = TRUE)
hypothesis(mod.gape, "area_perim_ratio_scaled > 0", robust = TRUE)
hypothesis(mod.gape, "Char1:area_perim_ratio_scaled < 0", robust = TRUE)

hypothesis(mod.jp, "Char1 < 0", robust = TRUE)
hypothesis(mod.jp, "area_perim_ratio_scaled > 0", robust = TRUE)
hypothesis(mod.jp, "Char1:area_perim_ratio_scaled < 0", robust = TRUE)

hypothesis(mod.grl, "Char1 < 0", robust = TRUE)
hypothesis(mod.grl, "area_perim_ratio_scaled > 0", robust = TRUE)
hypothesis(mod.grl, "Char1:area_perim_ratio_scaled < 0", robust = TRUE)

hypothesis(mod.spacing, "Char1 < 0", robust = TRUE)
hypothesis(mod.spacing, "area_perim_ratio_scaled > 0", robust = TRUE)
hypothesis(mod.spacing, "Char1:area_perim_ratio_scaled < 0", robust = TRUE)

hypothesis(mod.rakcoug, "Char1 < 0", robust = TRUE)
hypothesis(mod.rakcoug, "area_perim_ratio_scaled > 0", robust = TRUE)
hypothesis(mod.rakcoug, "Char1:area_perim_ratio_scaled < 0", robust = TRUE)

hypothesis(mod.si, "Char1 > 0", robust = TRUE)
hypothesis(mod.si, "area_perim_ratio_scaled < 0", robust = TRUE)
hypothesis(mod.si, "Char1:area_perim_ratio_scaled > 0", robust = TRUE)

hypothesis(mod.kt, "Char1 > 0", robust = TRUE)
hypothesis(mod.kt, "area_perim_ratio_scaled < 0", robust = TRUE)
hypothesis(mod.kt, "Char1:area_perim_ratio_scaled < 0", robust = TRUE)

hypothesis(mod.lr, "Char1 > 0", robust = TRUE)
hypothesis(mod.lr, "area_perim_ratio_scaled > 0", robust = TRUE)
hypothesis(mod.lr, "Char1:area_perim_ratio_scaled < 0", robust = TRUE)


### Figure 4
{
### Defense
#separating char and no char to avoid facet wrap
tmp = fish.scaled
tmp$Site = as.character(tmp$Site)
for(i in 1:length(display_sites)) tmp[tmp$Site==names(display_sites)[i],"Site"] = display_sites[i]
tmp = rbind(tmp, NA)
tmp[nrow(tmp),"Site"] = ""
tmp$Site = as.factor(tmp$Site)
tmp$Char = as.character(tmp$Char)
tmp[nrow(tmp),"Char"] = 0.5
tmp$Char = as.factor(tmp$Char)
tmp[nrow(tmp),"Spine_length_dorsal_1_cor"] = 1
tmp[nrow(tmp),"Pelvic_length_cor"] = 1
tmp = tmp %>%
  arrange(Char) %>%
  mutate(Site = factor(Site, levels = unique(Site)))
ordered_sites = unique(tmp$Site)
sep_line = which(levels(tmp$Site) == "")
levels(tmp$Char) = c("Nocharr", "x", "Charr")

spine.plot = ggplot(data=tmp, aes(x=Site, y=Spine_length_dorsal_1_cor, col=Char)) + geom_boxplot(outliers=F, key_glyph="point") +
  geom_jitter(aes(alpha=Char), width=0.1) +
  geom_vline(xintercept=sep_line, linetype="dotted") + scale_x_discrete(limits = ordered_sites) +
  scale_alpha_manual(values=c(Nocharr = 0.3, x = 0, Charr = 0.3), guide="none") +
  scale_color_manual(values=c(Nocharr = clrs[1], x = "#00000000", Charr = clrs[2]), breaks=c("Nocharr","Charr"), labels=c("Absent","Present")) +
  xlab("Lake") + ylab("First dorsal spine length") + ggtitle("Defence") + guides(col=guide_legend(title="Arctic charr")) +
  annotate(geom = "text", size = 5,
           label = "A",
           x = -Inf, hjust = 2.8,
           y = Inf, vjust = 1.5,
           fontface="bold") + coord_cartesian(clip="off") +
  theme_classic() +
  theme(legend.position = c(1.24,1),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        plot.title = element_text(hjust = 0.5, vjust=2, size=20))

spine.density = ggplot(data=tmp[tmp$Site!="",], aes(x=Spine_length_dorsal_1_cor, fill=Char, col=Char)) + 
  geom_density(alpha = 0.4) +
  geom_vline(xintercept=mean(fish.scaled[fish.scaled$Char==0,"Spine_length_dorsal_1_cor"], na.rm=T), col=clrs[1]) +
  geom_vline(xintercept=mean(fish.scaled[fish.scaled$Char==1,"Spine_length_dorsal_1_cor"], na.rm=T), col=clrs[2]) +
  scale_fill_manual(values=c(clrs[1], clrs[2])) +
  scale_color_manual(values=c(clrs[1], clrs[2])) +
  scale_y_continuous(expand=c(0,0), limits=c(0,0.6)) +
  theme_void() +
  theme(legend.position = "none") +
  coord_flip()

pelvic.plot = ggplot(data=tmp, aes(x=Site, y=Pelvic_length_cor, col=Char)) + geom_boxplot(outliers=F) +
  geom_jitter(aes(alpha=Char), width=0.1) +
  geom_vline(xintercept=sep_line, linetype="dotted") + scale_x_discrete(limits = ordered_sites) +
  scale_alpha_manual(values=c(0.3, 0, 0.3)) +
  scale_color_manual(values=c(clrs[1], "#00000000", clrs[2])) +
  xlab("Lake") + ylab("Ascending process width") +
  annotate(geom = "text", size = 5,
           label = "B",
           x = -Inf, hjust = 2.8,
           y = Inf, vjust = 1.5,
           fontface="bold") + coord_cartesian(clip="off") +
  theme_classic() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

pelvic.density = ggplot(data=tmp[tmp$Site!="",], aes(x=Pelvic_length_cor, fill=Char, col=Char)) + 
  geom_density(alpha = 0.4) +
  geom_vline(xintercept=mean(fish.scaled[fish.scaled$Char==0,"Pelvic_length_cor"], na.rm=T), col=clrs[1]) +
  geom_vline(xintercept=mean(fish.scaled[fish.scaled$Char==1,"Pelvic_length_cor"], na.rm=T), col=clrs[2]) +
  scale_fill_manual(values=c(clrs[1], clrs[2])) +
  scale_color_manual(values=c(clrs[1], clrs[2])) +
  scale_y_continuous(expand=c(0,0), limits=c(0,0.6)) +
  theme_void() +
  theme(legend.position = "none") +
  coord_flip()


### Foraging

## Gill raker length
post.grl = tibble(area_perim_ratio_scaled = c(
  seq(min(fish.scaled[fish.scaled$Char==0,]$area_perim_ratio_scaled), 
      max(fish.scaled[fish.scaled$Char==0,]$area_perim_ratio_scaled), length.out=100),
  seq(min(fish.scaled[fish.scaled$Char==1,]$area_perim_ratio_scaled), 
      max(fish.scaled[fish.scaled$Char==1,]$area_perim_ratio_scaled),length.out=100)),
  Char=as.factor(c(rep(0,100), rep(1,100)))) %>%
  add_epred_draws(object=mod.grl,
                   scale = "response", ndraws = 1e3,
                   re_formula =~1|Site + 1|Drainage,
                   allow_new_levels = TRUE)

grl.plot = ggplot() +
  geom_line(data=post.grl[post.grl$Char==0,],
            aes(x = area_perim_ratio_scaled, y = .epred, group = .draw), alpha = 0.015, col=clrs[1]) +
  geom_line(data=post.grl[post.grl$Char==1,],
            aes(x = area_perim_ratio_scaled, y = .epred, group = .draw), alpha = 0.015, col=clrs[2]) +
  geom_point(data=fish.scaled, aes(x=area_perim_ratio_scaled, y=Mean_raker_length_cor, col=Char), alpha=0.3) +
  scale_color_manual(values=clrs) + xlab("Lake Area:Perimeter") + ylab("Gill raker length") + ggtitle("Foraging") +
  annotate(geom = "text", size = 5,
           label = "C",
           x = -Inf, hjust = 2.8,
           y = Inf, vjust = 1.5,
           fontface="bold") + coord_cartesian(clip="off") +
  scale_y_continuous(breaks=c(-2, 0, 2, 4)) +
  theme_classic() +
  theme(legend.position = "none",
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        plot.title = element_text(hjust = 0.5, vjust=2, size=20))

grl.density = ggplot(data=fish.scaled, aes(x=Mean_raker_length_cor, fill=Char, col=Char)) + 
  geom_density(alpha = 0.4) +
  geom_vline(xintercept=mean(fish.scaled[fish.scaled$Char==0,"Mean_raker_length_cor"], na.rm=T), col=clrs[1]) +
  geom_vline(xintercept=mean(fish.scaled[fish.scaled$Char==1,"Mean_raker_length_cor"], na.rm=T), col=clrs[2]) +
  scale_fill_manual(values=c(clrs[1], clrs[2])) +
  scale_color_manual(values=c(clrs[1], clrs[2])) +
  scale_y_continuous(expand=c(0,0), limits=c(0,0.6)) +
  theme_void() +
  theme(legend.position = "none") +
  coord_flip()

## Suction index
post.si = tibble(area_perim_ratio_scaled = c(
  seq(min(fish.scaled[fish.scaled$Char==0,]$area_perim_ratio_scaled), 
      max(fish.scaled[fish.scaled$Char==0,]$area_perim_ratio_scaled), length.out=100),
  seq(min(fish.scaled[fish.scaled$Char==1,]$area_perim_ratio_scaled), 
      max(fish.scaled[fish.scaled$Char==1,]$area_perim_ratio_scaled),length.out=100)),
  Char=as.factor(c(rep(0,100), rep(1,100)))) %>%
  add_epred_draws(object=mod.si,
                   scale = "response", ndraws = 1e3,
                   re_formula =~1|Site + 1|Drainage,
                   allow_new_levels = TRUE)

si.plot = ggplot() +
  geom_line(data=post.si[post.si$Char==0,],
            aes(x = area_perim_ratio_scaled, y = .epred, group = .draw), alpha = 0.015, col=clrs[1]) +
  geom_line(data=post.si[post.si$Char==1,],
            aes(x = area_perim_ratio_scaled, y = .epred, group = .draw), alpha = 0.015, col=clrs[2]) +
  geom_point(data=fish.scaled, aes(x=area_perim_ratio_scaled, y=Suction_index, col=Char), alpha=0.3) +
  scale_color_manual(values=clrs) + xlab("Lake Area:Perimeter") + ylab("Suction index") +
  annotate(geom = "text", size = 5,
           label = "D",
           x = -Inf, hjust = 2.8,
           y = Inf, vjust = 1.5,
           fontface="bold") + coord_cartesian(clip="off") +
  theme_classic() +
  theme(legend.position = "none")

si.density = ggplot(data=fish.scaled, aes(x=Suction_index, fill=Char, col=Char)) + 
  geom_density(alpha = 0.4) +
  geom_vline(xintercept=mean(fish.scaled[fish.scaled$Char==0,"Suction_index"], na.rm=T), col=clrs[1]) +
  geom_vline(xintercept=mean(fish.scaled[fish.scaled$Char==1,"Suction_index"], na.rm=T), col=clrs[2]) +
  scale_fill_manual(values=c(clrs[1], clrs[2])) +
  scale_color_manual(values=c(clrs[1], clrs[2])) +
  scale_y_continuous(expand=c(0,0), limits=c(0,0.6)) +
  theme_void() +
  theme(legend.position = "none") +
  coord_flip()


### Full plot
fig.comparison = spine.plot + spine.density + plot_spacer() + grl.plot + grl.density +
  plot_spacer() + plot_spacer() + plot_spacer() + plot_spacer() + plot_spacer() +
  pelvic.plot + pelvic.density + plot_spacer() + si.plot + si.density +
  plot_layout(nrow=3, ncol=5, widths=c(4,1,0.3,4,1), heights=c(1,-0.101,1))
}

ggsave(filename="figure_means_comparison.png",
        plot=fig.comparison,
        device="png", dpi = 300,
        width=30, height=20, units="cm"
)

# Appendix figure for plate and raker number
# Figure S2
{
tmp = fish.scaled
tmp$Site = as.character(tmp$Site)
tmp = subset(tmp, select=-c(Plate_count, Raker_count))
tmp = left_join(tmp, fish[,c("FishEc","Plate_count","Raker_count")], by = "FishEc") # because we want them in true scale
for(i in 1:length(display_sites)) tmp[tmp$Site==names(display_sites)[i],"Site"] = display_sites[i]
tmp = rbind(tmp, NA)
tmp[nrow(tmp),"Site"] = ""
tmp$Site = as.factor(tmp$Site)
tmp$Char = as.character(tmp$Char)
tmp[nrow(tmp),"Char"] = 0.5
tmp$Char = as.factor(tmp$Char)
tmp[nrow(tmp),"Plate_count"] = 1
tmp[nrow(tmp),"Raker_count"] = 12
tmp = tmp %>%
  arrange(Char) %>%
  mutate(Site = factor(Site, levels = unique(Site)))
ordered_sites = unique(tmp$Site)
sep_line = which(levels(tmp$Site) == "")
levels(tmp$Char) = c("Nocharr", "x", "Charr")

plates.plot = ggplot(data=tmp, aes(x=Site, y=Plate_count, col=Char)) + geom_boxplot(outliers=F, key_glyph="point") +
  geom_jitter(aes(alpha=Char), width=0.1) +
  geom_vline(xintercept=sep_line, linetype="dotted") + scale_x_discrete(limits = ordered_sites) +
  scale_alpha_manual(values=c(Nocharr = 0.3, x = 0, Charr = 0.3), guide="none") +
  scale_color_manual(values=c(Nocharr = clrs[1], x = "#00000000", Charr = clrs[2]), breaks=c("Nocharr","Charr"), labels=c("Absent","Present")) +
  xlab("Lake") + ylab("Number of plates") + guides(col=guide_legend(title="Arctic charr")) +
  annotate(geom = "text", size = 5,
           label = "A",
           x = -Inf, hjust = 2.8,
           y = Inf, vjust = 1.5,
           fontface="bold") + coord_cartesian(clip="off") +
  theme_classic() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

grn.plot = ggplot(data=tmp, aes(x=Site, y=Raker_count, col=Char)) + geom_boxplot(outliers=F, key_glyph="point") +
  geom_jitter(aes(alpha=Char), width=0.1) +
  geom_vline(xintercept=sep_line, linetype="dotted") + scale_x_discrete(limits = ordered_sites) +
  scale_alpha_manual(values=c(Nocharr = 0.3, x = 0, Charr = 0.3), guide="none") +
  scale_color_manual(values=c(Nocharr = clrs[1], x = "#00000000", Charr = clrs[2]), breaks=c("Nocharr","Charr"), labels=c("Absent","Present")) +
  xlab("Lake") + ylab("Number of gill rakers") + guides(col=guide_legend(title="Arctic charr")) +
  annotate(geom = "text", size = 5,
           label = "B",
           x = -Inf, hjust = 2.8,
           y = Inf, vjust = 1.5,
           fontface="bold") + coord_cartesian(clip="off") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
}

ggsave(filename="figure_mean_comparisons_counts.png",
      plot=plates.plot + grn.plot + plot_layout(ncol=2),
      device="png", dpi = 300,
      width=30, height=15, units="cm"
)

#--------------------------------------------------------------------------------------------------------------------------------------------


# (Co)variance and correlation matrix estimation with test of char effect
#--------------------------------------------------------------------------------------------------------------------------------------------
#prep data
fish.scaled.nona = na.omit(fish.scaled) #no NAs in Stan
Xd = model.matrix(~ as.character(Char), data = fish.scaled.nona) #fe for defensive traits 
Xf = model.matrix(~ as.character(Char)*area_perim_ratio_scaled, data = fish.scaled.nona) #fe for foraging traits

#arrange traits (defensive then foraging)
traits = fish.scaled.nona[,c("Spine_length_dorsal_1_cor", "Spine_length_dorsal_2_cor", "Spine_length_pelvic_cor",
                   "Pelvic_height_cor", "Pelvic_length_cor", "Plate_count",
                   "Gape_cor_head", "Jaw_prot_cor_head", "Mean_raker_length_cor",
                   "Mean_spacing_cor", "Raker_count", "Suction_index", "KT", "Lever_ratio")]
#scale counts
traits$Plate_count = (traits$Plate_count - mean(traits$Plate_count)) / sd(traits$Plate_count)
traits$Raker_count = (traits$Raker_count - mean(traits$Raker_count)) / sd(traits$Raker_count)

#Make data list for Stan
fishdl = list(
  N = nrow(fish.scaled.nona),
  T = 14,
  Td = 6,
  Tf = 8,
  L = length(unique(fish.scaled.nona$Site)),
  D = length(unique(fish.scaled.nona$Drainage)),
  Fd = ncol(Xd), 
  Ff = ncol(Xf),
  charr_pres = Xd[,2],
  drain_id = as.integer(factor(fish.scaled.nona$Drainage)),
  lake_id = as.integer(factor(fish.scaled.nona$Site)),
  Xd = Xd,
  Xf = Xf,
  traits = traits
)

#custom function to extract posteriors
extract = function(fit_obj) {
  vars = fit_obj$metadata()$stan_variables
  draws = posterior::as_draws_rvars(fit_obj$draws())
  
  lapply(vars, \(var_name){  
    posterior::draws_of(draws[[var_name]], with_chains = FALSE)
  }) |> setNames(vars)
}

#cmdstan directory (change accordingly if necessary)
set_cmdstan_path(cmdstan_path())

#compile models
# vcovm1 = cmdstan_model(stan_file = "mvtraitmod1.stan") #single P matrix
# vcovm2 = cmdstan_model(stan_file = "mvtraitmod2.stan") #separate P matrices

#estimate models /!\ can take a few hours
# m1 = vcovm1$sample(
#   data = fishdl,
#   iter_sampling = 500,
#   iter_warmup = 1000,
#   init = 1e-4,
#   chains = 4,  
#   parallel_chains = 4,
#   adapt_delta = 0.90,
#   refresh = 10)
# 
# post1 = extract(m1)
# l1 = loo(m1$draws("log_lik"))
# saveRDS(m1,"m1_vcovm1.RDS")
# saveRDS(post1,"post1_vcovm1.RDS")
# saveRDS(l1,"l1_vcovm1.RDS")
# launch_shinystan(m1) #inspect
m1 = readRDS("m1_vcovm1.RDS")
l1 = readRDS("l1_vcovm1.RDS")
post1 = readRDS("post1_vcovm1.RDS")

# m2 = vcovm2$sample(
#   data = fishdl,
#   iter_sampling = 500,
#   iter_warmup = 1000,
#   init = 1e-4,
#   chains = 4,  
#   parallel_chains = 4,
#   adapt_delta = 0.90,
#   refresh = 10)
# 
# post2 = extract(m2)
# l2 = loo(m2$draws("log_lik"))
# saveRDS(m2,"m2_vcovm2.RDS")
# saveRDS(post2,"post_vcovm2.RDS")
# saveRDS(l2,"l2_vcovm2.RDS")
# launch_shinystan(m2) #inspect
m2 = readRDS("m2_vcovm2.RDS")
l2 = readRDS("l2_vcovm2.RDS")
post2 = readRDS("post2_vcovm2.RDS")

#compare LOO
sum(l1$pointwise[,"looic"] - l2$pointwise[,"looic"]) #diff
sd(l1$pointwise[,"looic"] - l2$pointwise[,"looic"]) * #se
  sqrt(nrow(l1$pointwise))
#--------------------------------------------------------------------------------------------------------------------------------------------


# Trait-specific and pairwise comparisons of (co)variances and correlations
#--------------------------------------------------------------------------------------------------------------------------------------------
#organize results
post2 = readRDS("post2_vcovm2.RDS")
Tname = c("1st dorsal spine length", "2nd dorsal spine length", "Pelvic spine length",
          "Ascending process height", "Ascending process width", "Plate number",
          "Gape size", "Jaw protrusion", "Gill raker length", "Gill raker gap",
          "Gill raker number", "Suction index", "Kinematic transmission", "Lever ratio")
RNC = post2$RNC #corr no char / stickleback only
RC = post2$RC #corr w char
sdNC = post2$sd_obsNC #standard dev
sdC = post2$sd_obsC

#compare SDs
lVR = log(sdC/sdNC)
colnames(lVR) = Tname
lVRp = 
  as_tibble(lVR) %>%
  pivot_longer(everything(), names_to = "trait", values_to = "lVR") %>%
  mutate(trait = factor(trait, levels = Tname)) %>% 
  group_by(trait) %>%
  summarise(
    median = median((exp(lVR)-1)*100),
    lCI = quantile((exp(lVR)-1)*100, 0.05),
    uCI = quantile((exp(lVR)-1)*100, 0.95),
    pd = mean(sign((exp(lVR)-1)*100) == 
                sign(median((exp(lVR)-1)*100))),
    .groups = "drop"
  )
lVRp[,2:5] = round(lVRp[,2:5], 2)
lVRp
write.csv(lVRp, "lVRp.csv")

#plot SDs
sdlNC = reshape2::melt(sdNC)
sdlC = reshape2::melt(sdC)
sdlNC$type= "Absent"
sdlC$type= "Present"
sds = rbind(sdlNC, sdlC)
sds$trait = Tname[sds$Var2]
sds$trait = factor(sds$trait, levels = Tname)
sds$type = factor(sds$type, level = c("Present","Absent"))

deltaSD = lVRp[lVRp$trait %in% Tname[c(3, 4, 5, 6, 7, 13)], "median", drop = T]
deltaSD[deltaSD>0] <- paste0("+",deltaSD[deltaSD>0])
labels_vec = setNames(
  paste0("'", Tname[c(3, 4, 5, 6, 7, 13)], "'~(", deltaSD, "*'%'['σ'])"),
  Tname[c(3, 4, 5, 6, 7, 13)])

#plot traits with larger effects
sd.plot = 
    ggplot(data=sds[sds$trait %in% Tname[c(3, 4, 5, 6, 7, 13)],], 
                  aes(x = value, group = type, fill = type, color = type)) + 
    geom_density(aes(y = after_stat(scaled)), alpha = 0.4)+
    scale_x_continuous(breaks = pretty_breaks(n = 3))+
    scale_color_manual(name = "Artic charr", guide = guide_legend(reverse = TRUE), values = clrs[c(2,1)])+
    scale_fill_manual(name = "Artic charr", guide = guide_legend(reverse = TRUE), values = clrs[c(2,1)])+
    facet_wrap(.~trait, ncol = 3, scales = "free",
               labeller = labeller(trait = as_labeller(labels_vec, label_parsed)))+
    theme_classic()+
    theme(legend.position = "top", strip.background = element_blank(),
          axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.line.y = element_blank(),
          plot.margin = margin(t = 10, r = 35, b = 10, l = 35))+
    labs(x = "Standard deviation", y = "Density")

sd.plot = 
    ggdraw(sd.plot) +
      draw_label("A", x = 0.02, y = 0.98, hjust = 0, vjust = 1, fontface = "bold", size = 16)  
  
#compare correlations
RNCl = 
    as.data.frame.table(RNC, responseName = "corr") %>%
    mutate(
      trait1 = Tname[as.integer(Var2)],
      trait2 = Tname[as.integer(Var3)],
      trait_pair = paste(trait1, trait2, sep = "___")) %>%
    filter(trait1 < trait2) %>%
    select(draw = Var1, trait_pair, corr) 
RNCl$type = "absent"

RCl = 
  as.data.frame.table(RC, responseName = "corr") %>%
  mutate(
    trait1 = Tname[as.integer(Var2)],
    trait2 = Tname[as.integer(Var3)],
    trait_pair = paste(trait1, trait2, sep = "___")) %>%
  filter(trait1 < trait2) %>%
  select(draw = Var1, trait_pair, corr) 
RCl$type = "present"

#combine
Rcomb = rbind(RNCl, RCl)

#summarize
aggregate(corr ~ type * trait_pair, data = Rcomb, FUN = median)
aggregate(corr ~ type * trait_pair, data = Rcomb, 
                          FUN = quantile, c(0.05, 0.95))
#compare (q)
Cohenq = 
  Rcomb %>%
  pivot_wider(
    names_from = type,
    values_from = corr
  ) %>%
  mutate(
    z_present = atanh(present),
    z_absent  = atanh(absent),
    q = z_present - z_absent
  ) %>%
  group_by(trait_pair) %>%
  summarise(
    median_q = median(q, na.rm = TRUE),
    lCI = quantile(q, 0.05, na.rm = TRUE),
    uCI = quantile(q, 0.95, na.rm = TRUE),
    pd = mean(sign(q) == sign(median(q))),
    .groups = "drop"
  )

#save q values
Cohenq[,2:5] = round(Cohenq[,2:5], 2)
Cohenq$q_size = ifelse(abs(Cohenq$median_q)<0.1, "negligible",
                       ifelse(abs(Cohenq$median_q)>=0.1 & abs(Cohenq$median_q)<0.3, "small",
                              ifelse(abs(Cohenq$median_q)>=0.3 & abs(Cohenq$median_q)<0.5, "medium",
                                     ifelse(abs(Cohenq$median_q)>=0.5, "large", NA))))
#write.csv(Cohenq, "cohenq.csv")

#non-negligible effects
nne = Cohenq[Cohenq$q_size!= "negligible",]
print(nne[order(nne$median_q),], n = Inf)

#plot correlation matrix
mRNC = round(apply(RNC, c(2, 3), median),2)
mRC = round(apply(RC, c(2, 3), median),2)
rownames(mRNC) = colnames(mRNC) = rownames(mRC) = colnames(mRC) = Tname

qmat = round(apply(atanh(RC) - atanh(RNC), c(2,3), median), 2)
q_labels = matrix("", nrow = 14, ncol = 14) # empty label matrix
sel = upper.tri(qmat) & abs(qmat) > 0.1
q_labels[sel] = sprintf("%.2f", qmat[sel]) # 2 decimal places
q_labels_rot = t(apply(q_labels, 2, rev))
label_df = reshape2::melt(q_labels_rot)
label_df[label_df$value>0,]$value <- paste0("+",label_df[label_df$value>0,]$value)

mRNC[upper.tri(mRNC)] = mRC[upper.tri(mRC)] #make upper triangle char present
mRNC = t(apply(mRNC, 2, rev)) # so diagonal is \ rather than / (rotated clockwise, upper tri is sbc)

corr.plot = 
 ggplot(data=melt(mRNC)) + geom_point(aes(x=Var1, y=Var2, col=value, size=abs(value)), shape=15) +
  scale_color_gradient2(limits=c(-1,1), name = "correlation\n ", low="#001177", mid="#FFFFFF", high="#10cd6b") +
  guides(color = guide_colorbar(ticks.colour = "NA")) +
  geom_text(data = label_df, aes(x = Var1, y = Var2, label = value),
            size = 3.1, col="#FF0000", nudge_x = 0.1, nudge_y = 0.2, fontface = "bold") +
  annotation_custom(
    grob = textGrob("Cohen's q\n(Charr effect)", gp = gpar(col = "#FF0000", fontface = "bold")),
    xmin = 13, xmax = 20, ymin = 12, ymax = 12
  ) +
  geom_abline(slope=-1, intercept=15, col="#333333", linewidth=0.7) +
  scale_size(limits=c(0,1), range=c(2,8), guide="none") +
  theme_bw() +
  labs(title = "") +
  theme(legend.position="right",
        legend.margin=margin(t = 0, r = 0, b = 0, l = 0), 
        legend.frame=element_rect(fill=NA, color="#333333"),
        legend.direction="vertical",
        plot.margin = margin(t = 0.5, r = 2, b = 2, l = 2), 
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        axis.title = element_blank()) +
  coord_fixed(clip = "off", ratio=1)
corr.plot =
  ggdraw(corr.plot) +
  draw_label("B", x = 0.02, y = 0.98, hjust = 0, vjust = 1, fontface = "bold", size = 16)

varcor.plot = sd.plot + corr.plot +
  plot_layout(ncol=1, heights = c(0.6,1))
# Figure 5
ggsave(filename="figure_varcor.png",
        plot=varcor.plot,
        device="png",
        width=20, height=25, units="cm")


#--------------------------------------------------------------------------------------------------------------------------------------------


# Whole matrix comparisons through various metrics
#--------------------------------------------------------------------------------------------------------------------------------------------
#organize results
post2 = readRDS("post2_vcovm2.RDS")
Tname = c("1st dorsal spine length", "2nd dorsal spine length", "Pelvic spine length",
          "Ascending process height", "Ascending process width", "Plate number",
          "Gape size", "Jaw protrusion", "Gill raker length", "Gill raker gap",
          "Gill raker number", "Suction index", "Kinematic transmission", "Lever ratio")
RNC = post2$RNC #corr no char / stickleback only
RC = post2$RC #corr w char
sdNC = post2$sd_obsNC #standard dev
sdC = post2$sd_obsC

#data
PNC = post2$PNC #cov no char / stickleback only
PC = post2$PC #cov w char

#median cor matrix
mRNC = round(apply(RNC, c(2, 3), median),2)
mRC = round(apply(RC, c(2, 3), median),2)

#median cov matrix
mPNC = round(apply(PNC, c(2, 3), median),2)
mPC = round(apply(PC, c(2, 3), median),2)

#functions
vect_angle = function(v1, v2){
  a = 180/pi * acos(v1 %*% v2 / (norm(matrix(v1), type="2")*norm(matrix(v2), type="2")))
  return(min(a, 180-a))
}

subspace_similarity = function(M, N, k){
  A = eigen(M)$vectors[,1:k]
  B = eigen(N)$vectors[,1:k]
  S = t(A) %*% B %*% t(B) %*% A
  return(sum(eigen(S)$values)/k)
}

ed = function(eig_vals) {
  return(1/prod((eig_vals/sum(eig_vals))**(eig_vals/sum(eig_vals))))
}

auton<-function(G){
  lambda<-function(G){
    #eigenvalues
    lambda<-eigen(G,only.values=TRUE)$values
    return(lambda)
  }
  H<-function(x) {
    #harmonic mean
    H=1/mean(1/x)
    return(H)
  }
  I<-function(x) {
    #mean.standardized variance
    var.x=(sum(x^2)/length(x))-(sum(x)/length(x))^2
    I=var.x/(mean(x)^2)
    return(I)
  }
  #below the generic autonomy function calling the functions above
  a<-function(G) {
    k= length(lambda(G))
    a<-( H(lambda(G)) / mean(lambda(G)) ) * (1+( (2*(I(lambda(G)) + I(1/lambda(G)) -1 + ( H(lambda(G)) / mean(lambda(G)) ) + 
                                                       ((2*I(lambda(G))*I(1/lambda(G)))/(k+2)) )) / (k+2)) )
    return(a)
  }
  gen.a=a(G)
  return(gen.a)
}

# amount of explained variance in p_max
eigen(mPC)$values[1] / sum(eigen(mPC)$values)
eigen(mPNC)$values[1] / sum(eigen(mPNC)$values)

# effective dimensionality
ed(eigen(mPC)$values)
ed(eigen(mPNC)$values)

# total phenotypic variance
sum(eigen(mPC)$values)
sum(eigen(mPNC)$values)

# angle between p_maxs
vect_angle(eigen(mPNC)$vectors[,1], eigen(mPC)$vectors[,1])

# compare pmax
eig1 = eigen(mPNC)
eig2 = eigen(mPC)
pmax1 = eig1$vectors[, 1]
pmax2 = eig2$vectors[, 1]
flip_sign = function(vec) { #make sure dom loadings +
  idx = which.max(abs(vec))
  if(vec[idx] < 0) vec = -vec
  return(vec)
}
pmax1 = flip_sign(pmax1)
pmax2 = flip_sign(pmax2)
df = data.frame(trait = rep(Tname, times = 2),
  type = rep(c("Charr absent", "Charr present"), each = 14),
  loading = c(as.vector(pmax1), as.vector(pmax2))
)
df$trait = factor(df$trait, levels = Tname)
df$type = factor(df$type, levels = c("Charr absent", "Charr present"))
pmax_plot =
  ggplot(df, aes(x = type, y = trait, fill = loading)) +
  geom_tile(color = "#ffffff") +
  scale_y_discrete(limits = rev)+
  geom_text(aes(label = round(loading, 2)), size = 3) +   
  scale_fill_gradient2(low = "#001177", mid = "#ffffff", high = "#10cd6b", 
                       midpoint = 0, limits = c(-1, 1), breaks = c(-1,0,1),
                       name = expression(p[max]~"loading   ")  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(face = "bold"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "top")
# Figure S4
ggsave(filename="figure_pmax.png",
       plot=pmax_plot,
       device="png",
       width=15, height=15, units="cm")

# how many PCs are expected beyond chance? (5)
psych::fa.parallel(mRC, n.obs = fishdl$N)

# similarity at k=5
subspace_similarity(mPC, mPNC, 5)

# Krzanowski similarity and angles
set.seed(510)
similarity.df = data.frame()
for(j in 1:nrow(PNC)){
  P_sb = PNC[j,,]
  P_sb_eig = eigen(P_sb)
  P_sbc = PC[j,,]
  P_sbc_eig = eigen(P_sbc)
  similarity.df = rbind(similarity.df, c(
    subspace_similarity(P_sb, P_sbc, 5),
    vect_angle(P_sb_eig$vectors[,1], P_sbc_eig$vectors[,1])
  ))
}
names(similarity.df) = c("similarity","angle")

# effective number of dimensions and total variance from all MCMC samples
eff_dim_sb = eff_dim_sbc = tot_var_sb = tot_var_sbc = c()
for(j in 1:nrow(PNC)){
  eig_vals = eigen(PNC[j,,])$values
  eff_dim_sb = c(eff_dim_sb, ed(eig_vals))
  tot_var_sb = c(tot_var_sb, sum(eig_vals))
  eig_vals = eigen(PC[j,,])$values
  eff_dim_sbc = c(eff_dim_sbc, ed(eig_vals))
  tot_var_sbc = c(tot_var_sbc, sum(eig_vals))
}
metrics.df = as.data.frame(cbind(eff_dim_sb, eff_dim_sbc, tot_var_sb, tot_var_sbc))

# trait autonomy
adf = data.frame(a_NC = rep(NA, nrow(PNC)), a_C = rep(NA, nrow(PNC)) )
for(j in 1:nrow(PNC)){
 adf$a_NC[[j]]  = auton(PNC[j,,])
 adf$a_C[[j]]  = auton(PC[j,,])
}

### Figure 6
simdist.plot = ggplot(data=similarity.df) +
  geom_density(aes(x=similarity, y = after_stat(scaled))) +
  xlab("Similarity coefficient (k=5)") + ylab("Density") +
  scale_x_continuous(limits=c(0,1)) + theme_classic() +
  scale_y_continuous(expand = c(0,0))+
  theme(  axis.title.y = element_blank(),
          axis.text.y  = element_blank(),
          axis.ticks.y = element_blank(),
          axis.line.y  = element_blank())

angle.plot = ggplot(data=similarity.df) +
  geom_density(aes(x=angle, y = after_stat(scaled))) +
  xlab(bquote("Angle between p"[max])) + ylab("Density") +
  scale_x_continuous(limits=c(0,91), breaks=c(0,45,90)) +
  scale_y_continuous(expand = c(0,0))+
  theme_classic()+
  theme(  axis.title.y = element_blank(),
          axis.text.y  = element_blank(),
          axis.ticks.y = element_blank(),
          axis.line.y  = element_blank())

dim.plot = ggplot(data=metrics.df) +
  geom_density(aes(x=eff_dim_sb, y = after_stat(scaled)), col=clrs[1], fill=clrs[1], alpha=0.4) +
  geom_density(aes(x=eff_dim_sbc, y = after_stat(scaled)), col=clrs[2], fill=clrs[2], alpha=0.4) +
  xlab("Effective number of dimensions") + ylab("Density") +
  scale_y_continuous(expand = c(0,0))+
  theme_classic()+  theme(  axis.title.y = element_blank(),
                            axis.text.y  = element_blank(),
                            axis.ticks.y = element_blank(),
                            axis.line.y  = element_blank())

var.plot = ggplot(data=metrics.df) +
  geom_density(aes(x=tot_var_sb, y = after_stat(scaled)), col=clrs[1], fill=clrs[1], alpha=0.4) +
  geom_density(aes(x=tot_var_sbc, y = after_stat(scaled)), col=clrs[2], fill=clrs[2], alpha=0.4) +
  xlab("Total phenotypic variance") + ylab("Density") +
  scale_y_continuous(expand = c(0,0))+
  theme_classic()+
  theme(  axis.title.y = element_blank(),
          axis.text.y  = element_blank(),
          axis.ticks.y = element_blank(),
          axis.line.y  = element_blank())

aut.plot = 
  ggplot(data= melt(adf), aes(x = value, color = variable, group = variable, fill = variable)) +
  geom_density(alpha = 0.4)+
  xlab("Average autonomy") + ylab("Density") +
  theme_classic()+
  scale_y_continuous(expand = c(0,0))+
  scale_color_manual(values = clrs[1:2])+
    scale_fill_manual(values = clrs[1:2])+
    theme(  axis.title.y = element_blank(),
            axis.text.y  = element_blank(),
            axis.ticks.y = element_blank(),
            axis.line.y  = element_blank())+
    guides(fill = "none", color = "none")

# plot just for legend
legend.plot = 
  ggplot(data.frame(x=1, y=1, grp=c("charr absent","charr present")),
                      aes(x, y, fill=grp)) +
  geom_col(alpha = 0.5) +
  scale_fill_manual(values = clrs[1:2],
                    name = "Arctic charr",
                    labels = c("Absent", "Present")) + theme_classic()
legend = cowplot::get_legend(legend.plot)
legend = ggdraw(legend) + theme(plot.tag = element_text(color = "NA"))

metrics.plot = 
  simdist.plot + angle.plot + var.plot + dim.plot + aut.plot + legend +
  plot_layout(ncol = 3, guides = "collect") +
  plot_annotation(tag_levels = 'A') & 
  theme(plot.tag = element_text(face = "bold", hjust = 0, vjust = 1))

ggsave(filename="figure_matrix_metrics.png",
        plot=metrics.plot,
        device="png",
        width=20, height=10, units="cm")


# Mdularity
# Figure S3
colnames(mRNC) = rownames(mRNC) = colnames(mRC) = rownames(mRC) = Tname
ggC = as.matrix(EBICglasso(mRC,n = fishdl$N, threshold = T))
ggmC = graph_from_adjacency_matrix(ggC, mode = "undirected", weighted = TRUE, diag = FALSE)
ggmaC = graph_from_adjacency_matrix(abs(ggC), mode = "undirected", weighted = TRUE, diag = FALSE)
cC = cluster_louvain(ggmaC)

ggNC = as.matrix(EBICglasso(mRNC, n = fishdl$N, threshold = T))
ggmNC = graph_from_adjacency_matrix(ggNC, mode = "undirected", weighted = TRUE, diag = FALSE)
ggmaNC = graph_from_adjacency_matrix(abs(ggNC), mode = "undirected", weighted = TRUE, diag = FALSE)
cNC = cluster_louvain(ggmaNC)

png(file = "figure_network.png", res=600, width=15, height=7.5, units="in")
layout(matrix(c(1,2),nrow=1))

colNC = c(rep("#B2DFEE",3), rep("#7FFFD4",2), "#FFFFFF", rep("#FFD700",3), 
          rep("#98FB98", 2), rep("white",2), "#FFD700")
colC = c(rep("#B2DFEE",3), rep("#7FFFD4",2), "#FFFFFF", rep("#FFD700",3), 
         rep("#98FB98", 2), "#FFD700", "#FFFFFF", "#FFD700")

qgraph::qgraph(mRNC, maximum=1, layout="circle", esize=10,vsize=10,colFactor=0.4,
               posCol="#35abd6", fade = T, groups = split(names(membership(cNC)), membership(cNC)),
               minimum=0.1, label.prop=0.5, color = colNC, 
               labels = abbreviate(Tname, minlength = 3), title = "Charr absent", title.cex = 2,
               legend = F, labels = TRUE, edge.labels=F, node.resolution=300, mar=c(2,2,4,2))

qgraph::qgraph(mRC, maximum=1, layout="circle", esize=10,vsize=10,colFactor=0.4,
               posCol="#35abd6", fade = T, groups = split(names(membership(cC)), membership(cC)),
               minimum=0.1, label.prop=0.5, color = colC, 
               labels = abbreviate(Tname, minlength = 3), title = "Charr present", title.cex = 2,
               legend = F, labels = TRUE, edge.labels=F, node.resolution=300, mar=c(2,2,4,2))
dev.off()
#--------------------------------------------------------------------------------------------------------------------------------------------
