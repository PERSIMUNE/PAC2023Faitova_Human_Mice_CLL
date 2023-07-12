#R script to perform statistical analysis and visualizations for the mice samples

##########################################################
#Load meta data for mice samples

meta_mouse = read.csv('meta_mouse_complete.csv', sep = ';')
rownames(meta_mouse) = meta_mouse$ID_metag

##########################################################
#Load taxonimical profiling data

split = function(string){unlist(str_split(string ,"_", n=2))[1]}

#species level
mouse_OTUS_rel = read.csv('motus_and_specI_relative_abundance.txt', sep = '\t', header = T)    
rownames(mouse_OTUS_rel) = mouse_OTUS_rel[,1]
mouse_OTUS_rel = mouse_OTUS_rel[ , -1]
rownames(mouse_OTUS_rel)[2] = 'Unclassified'
for (i in 1:length(colnames(mouse_OTUS_rel))){ 
  colnames(mouse_OTUS_rel)[i] = split(colnames(mouse_OTUS_rel)[i])
}

#genera level
mouse_genera_rel = read.csv('genus_relative_abundance.txt', sep = '\t', header = T)
rownames(mouse_genera_rel) = mouse_genera_rel[,1]
mouse_genera_rel = mouse_genera_rel[ , -1]
rownames(mouse_genera_rel)[1] = 'Unclassified'
for (i in 1:length(colnames(mouse_genera_rel))){ 
  colnames(mouse_genera_rel)[i] = split(colnames(mouse_genera_rel)[i])
}


##########################################################
#mouseID categories
mouse_low_hygiene = rownames(subset(meta_mouse, hygiene_level == 'Low Hygiene'))
mouse_high_hygiene = rownames(subset(meta_mouse, hygiene_level == 'High Hygiene'))

mouse_tumor = rownames(subset(meta_mouse, timepoint.x == 'after injection'))
mouse_notumor = rownames(subset(meta_mouse, timepoint.x == 'before injection'))

meta_mouse$ID_tumor = paste(meta_mouse$`Input.ID.x`,'.',meta_mouse$TM, sep = "")
meta_mouse$ID_tumor_num = as.numeric(meta_mouse$ID_tumor)

meta_mouse$hygiene_num = ifelse(meta_mouse$hygiene_level == 'Low Hygiene', 0, 1)

##########################################################
#Filter data & clr transform
filtered_mouse_names_no0 = names(which(rowSums(mouse_OTUS_rel) != 0))
mouse_filtered_no0 = mouse_OTUS_rel[filtered_mouse_names_no0, ]
mouse_filtered_no0_rel = mouse_filtered_no0[rownames(mouse_filtered_no0) != 'Unclassified', ]
mouse_clr_filtered_no0_rel = as.data.frame(compositions::clr(mouse_filtered_no0_rel))

filtered_mouse_names_no0 = names(which(rowSums(mouse_OTUS) != 0))
mouse_filtered_no0 = mouse_OTUS[filtered_mouse_names_no0, ]

filtered_mouse_names_no0 = names(which(rowSums(mouse_genera_rel) != 0))
mouse_genera_rel_no0 = mouse_genera_rel[filtered_mouse_names_no0, ]
mouse_genera_no0_clr_rel = as.data.frame(compositions::clr(mouse_genera_rel_no0))

##########################################################
#Alpha diversity low vs high hygiene
richness_LH <- specnumber(mouse_filtered_no0[,mouse_low_hygiene],MARGIN = 2,)
richness_HH <- specnumber(mouse_filtered_no0[,mouse_high_hygiene], MARGIN = 2,)

df_diversity_LH = as.data.frame(richness_LH)
df_diversity_HH = as.data.frame(richness_HH)

df_diversity_LH$Shannon <- vegan::diversity(mouse_filtered_no0[ ,mouse_low_hygiene], MARGIN = 2, index = "shannon")
df_diversity_LH$InvSimpson <- vegan::diversity(mouse_filtered_no0[ ,mouse_low_hygiene], MARGIN = 2, index = "invsimpson")
df_diversity_HH$Shannon <- vegan::diversity(mouse_filtered_no0[ ,mouse_high_hygiene], MARGIN = 2, index = "shannon")
df_diversity_HH$InvSimpson <- vegan::diversity(mouse_filtered_no0[ ,mouse_high_hygiene], MARGIN = 2, index = "invsimpson")

df_diversity_LH$hygiene = 'low hygiene'
df_diversity_HH$hygiene = 'high hygiene'

colnames(df_diversity_LH)[1] = 'Richness'
colnames(df_diversity_HH)[1] = 'Richness'
df = rbind(df_diversity_LH, df_diversity_HH)
df1 = merge(df, meta_mouse %>% dplyr::select(timepoint.x), by = 0)
rownames(df1) = df1$Row.names
df1$Row.names = NULL
dfmelt = reshape2::melt(df1)
colnames(dfmelt)[2] = 'tumor'

#Plot diversity
#stats 
stats_df = compare_means(value ~ hygiene, group.by = "variable", p.adjust.method = "BH", data = dfmelt)

#plot
ggplot(dfmelt, aes(x=hygiene, y=value)) + 
  geom_boxplot(aes(fill = tumor)) +
  labs(title="Alpha diversity measures - Species level", x ="Hygiene", y = "Diversity") +
  theme(strip.background = element_rect(colour="black", fill="white", size=0.7, linetype="solid"), 
        strip.text.x = element_text(size=10, color="black", face = 'bold')) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle=-40, hjust=.1, face = 'bold'))+
  ggsignif::geom_signif(data=stats_df, 
                        aes(xmin=group1, xmax=group2, annotations= paste(p.adj, p.signif), 
                            y_position=c(109, 3.5,  21)), 
                        manual=TRUE,  size = 0.9, tip_length = 0.005) +
  
  facet_wrap(~ variable, scales = "free_y", nrow = 1)+
  theme(strip.text.x = element_text(size = 11, color = "black", face ='bold'))


##########################################################
#Plot the composition of all mice microbiome samples 

#Include those genera being under 1% and merge them into one variable 'Others' to see how much they take up
#keep the unassigned genera
mouse_genera_rel_no0_t = as.data.frame(t(mouse_genera_rel_no0))

#mouse_genera_rel_no0_t$SampleID <- rownames(mouse_genera_rel_no0_t)
mouse_genera_rel_no0_t = merge(mouse_genera_rel_no0_t, meta_mouse %>% dplyr::select(ID_tumor), by = 0)
mouse_genera_rel_no0_t$Row.names = NULL

#add artificial fill for 12.0 missing sample
mouse_genera_rel_no0_t$fill = as.numeric(0)
df_fill = as.data.frame(t(c(rep(0,43), '12.0', 1)))
colnames(df_fill) = colnames(mouse_genera_rel_no0_t)
df_fill[ ,c(1:43, 45)] = as.numeric(df_fill[ ,c(1:43, 45)])
mouse_genera_rel_no0_t = rbind(mouse_genera_rel_no0_t, df_fill)

df <- reshape2::melt(mouse_genera_rel_no0_t, 'ID_tumor')


#assign own color to each genus - custom color scale created in R script 'data,packages.R'
colors43 = c("#999999", "#FF8C00", "#1E90FF", "#009E73", "#DAA520", "#0072B2", "#FF4500", "#CC79A7",
             "#B9D4DB", '#D8BFD8', "#E8175D", "#B8860B",'#FF1493',"darkblue", '#9ACD32',"#800000",'#F08080',
             '#98FB98','#8A2BE2','#F5DEB3','#2F4F4F','#FFFF80','#BA55D3',"lightcyan3",
             '#CD5C5C',"#00FECA", "#8B0000","#4B0082", "#FF85EA",'#00FFFF', '#FFD700','#808000','#00CED1','#FF6347',
             "#000000",'linen', 'red',"#027A9F","#BDBDFD","yellow2",'mediumorchid1','#0000FF',"olivedrab1",
             "deepskyblue1", "paleturquoise1")
show_col(colors43)
names(colors43) <- levels(factor(df$variable))

#add ranks - most abundant genera
sorted <- as.data.frame(sort(colMeans(mouse_genera_rel_no0_t[ ,1:43]), decreasing = T))
sorted$rank = c(1:43)

#fill becomes rank 44 automatically
df$rank = 44
for (i in 1:43){
  df$rank[which(df$variable %in% rownames(sorted)[i])] = sorted$rank[i]
}

#put Unclassified to the bottom of ranking
df[1:19,4] = max(df$rank)+1 


levels_xaxis = c("10.0","11.0", "12.0", "13.0", "10.1", "11.1", "12.1", "13.1", 
                 "18.0", "28.0", "29.0", "30.0", "34.0", "35.0", 
                 "18.1", "28.1", "29.1", "30.1", "34.1", "35.1")

#plot
p = ggplot(df, aes(x = factor(ID_tumor, level = levels_xaxis), y = value, fill = variable, group = plyr::desc(rank))) +
  geom_bar(aes(fill=variable), stat = "identity", width=0.9, colour = "#303030", size = 0.05) +
  scale_fill_manual(name = "variable", values = colors43) +
  guides(fill=guide_legend(title="Genus", ncol = 2)) +
  xlab("Mouse ID, before/after tumor injection (0/1)")+
  ylab("Relative Abundance Genus") +
  ggtitle("Composition of mice gut microbiome \n Biological classification: Genus \n") +
  theme_bw()+ 
  theme(axis.text.x = element_text(angle=-40, hjust=.1, color = 'black', face = 'bold'))+ #, legend.position='none') 
  theme(legend.title = element_text(size = 8), legend.text = element_text(size = 9))
p

ggplotly(p)

df_notumor = dcast(df[grep('\\.0', df$ID_tumor), ], ID_tumor ~ variable)
names = names(which(colSums(df_notumor[,-1])< 0.01))


##########################################################
#data merging
mouse_genera_no0_clr_rel_t = as.data.frame(t(mouse_genera_no0_clr_rel))
mouse_logistic = merge(mouse_genera_no0_clr_rel_t, meta_mouse %>% dplyr::select(hygiene_num, TM), by = 0)
rownames(mouse_logistic) = mouse_logistic$Row.names
mouse_logistic = mouse_logistic[,-1]
mouse_logistic$TM = as.numeric(as.character(mouse_logistic$TM))
mouse_logistic$cage = as.factor(mouse_logistic$cage)


##########################################################
#PERMANOVA - test tumor, hygiene and cage effect

mouse_genera_no0_rel_t = as.data.frame(t(mouse_genera_rel_no0))
mat = as.matrix(mouse_genera_no0_rel_t[intersect(mouse_low_hygiene, mouse_notumor), ])
dist.mat_LH_TM0 <- vegan::vegdist(mat, method="robust.aitchison")
dist.mat_LH <- vegan::vegdist(mat, method="robust.aitchison")
dist.mat_TM1 <- vegan::vegdist(mat, method="robust.aitchison")



#cage effect among low hygiene mice kept separately
meta_mouse$cage = as.factor(meta_mouse$cage)
set.seed(99)
perMANOVA <- vegan::adonis2(dist.mat_LH_TM0 ~ cage, 
                            data = meta_mouse[meta_mouse$hygiene_num == 0 & meta_mouse$TM == 0, ], na.action = 'na.omit')
perMANOVA
#          Df SumOfSqs      R2      F Pr(>F)
#cage      1   53.178 0.29887 1.7051    0.1
#no cage effect

#tumor effect among low hygiene mice
set.seed(100)
perMANOVA <- vegan::adonis2(dist.mat_LH ~ TM, 
                            data = meta_mouse[meta_mouse$hygiene_num == 0, ], na.action = 'na.omit', permutations = 999)
perMANOVA

#hygiene effect among tumor transplanted mice (at time point 2)
set.seed(102)
perMANOVA <- vegan::adonis2(dist.mat_TM1 ~ hygiene_num, 
                            data = meta_mouse[meta_mouse$TM == 1, ], na.action = 'na.omit', permutations = 999)

perMANOVA


##########################################################
#Run linear regression in loop
#We  use the variables TM and hygiene (binary) to predict the effect on abundance of a given bacteria
#positive coefficients of hygiene means that hygiene_level == 1 (high hygiene barrier) are positively affecting the abundance of a bacteria (increase) 
#EXAMPLE: Take the variable hygiene as an example. 
#         If the other variable (TM) remain constant, a change in hygiene (0 to 1/low to high) means that 
#         there will be an average change in the bacterial abundance of about (-)XX units. 
#negative coefficients of hygiene means that hygiene_level == 1 (high hygiene barrier) is negatively affecting the abundance of a bacteria (decrease in abundance) 
##positive coefficients of tumor means that TM == 1 (tumor injection) is positively affecting the abundance of a bacteria (increase) 
##negative coefficients of tumor means that TM == 1 (tumor injection) is negatively affecting the abundance of a bacteria (decrease in abundance) 

df_lm = as.data.frame(rownames(sorted))
df_lm$coef_hygiene = 1
df_lm$std_e_hygiene = 2
df_lm$p.val_hygiene = 3
df_lm$coef_TM = 4
df_lm$std_e_TM = 5
df_lm$p.val_TM = 6

for (i in 1:length(rownames(sorted))){
  lm = glm(get(rownames(sorted)[i]) ~ hygiene_num + TM, data = mouse_logistic)
  df_lm[i,2] = coef(summary(lm))[2,1]
  df_lm[i,3] = coef(summary(lm))[2,2]
  df_lm[i,4] = coef(summary(lm))[2,4]
  df_lm[i,5] = coef(summary(lm))[3,1]
  df_lm[i,6] = coef(summary(lm))[3,2]
  df_lm[i,7] = coef(summary(lm))[3,4]
}

plot(df_lm$coef_hygiene, df_lm$coef_TM)

#Absolute change in relative abundance no tumor --> tumor
mouse_genera_rel_no0_t = as.data.frame(t(mouse_genera_rel_no0))

mouse_before = mouse_genera_rel_no0_t[mouse_notumor, ]
mouse_before = merge(mouse_before, meta_mouse %>% dplyr::select('Input.ID.x'), by = 0)
rownames(mouse_before) =  mouse_before$`Input ID`
mouse_before = mouse_before[,-1]
mouse_before = mouse_before[ order(row.names(mouse_before)), ]

mouse_after = mouse_genera_rel_no0_t[mouse_tumor, ]
mouse_after = merge(mouse_after, meta_mouse %>% dplyr::select('Input.ID.x'), by = 0)
rownames(mouse_after) =  mouse_after$`Input ID`
mouse_after = mouse_after[-7,-1] #remove also unpaired 'after tumor' sample 
mouse_after = mouse_after[ order(row.names(mouse_after)), ]

mouse_logistic_delta = mouse_after[,1:43] - mouse_before[,1:43]
mouse_logistic_delta = abs(mouse_logistic_delta)
#delta_abs_tumor = as.data.frame(colMeans(mouse_logistic_delta))
delta_abs_tumor = as.data.frame(apply(mouse_logistic_delta, 2, median))


#PLOT
rownames(df_lm) = df_lm$`rownames(sorted)`
df_lm = merge(df_lm,delta_abs_tumor, by=0)
df_lm$`rownames(sorted)` = NULL
colnames(df_lm)[6] = 'delta_tumor'
rownames(df_lm) = df_lm$Row.names

ggplot(data = df_lm, aes(x = coef_hygiene, y = coef_TM)) +
  geom_point(color = "blue", alpha=0.7, size = 2)+
  scale_size(range = c(.1, 10), name="Absolute change in rel. abundance after tumor injection")+
  geom_label_repel(aes(label = ifelse(coef_TM <= -0.05 | coef_TM >= 0.05 | coef_hygiene <= -0.05 | coef_hygiene >= 0.05, rownames(df_lm), '')),
                   size = 3,
                   max.overlaps = 20,
                   box.padding   = 0.2, 
                   point.padding = 0.3,
                   segment.size  = 0.2,
                   nudge_x = -0.2,
                   segment.color = 'grey50') +
  coord_fixed(ratio = 1)+
  xlim(-3,3.1)+
  ylim(-3,3)



##########################################################
####SIAMCAT - not included in the manuscript

library("SIAMCAT")
tax.profiles = as.matrix(mouse_OTUS_rel)

#Only species that have a relative abundance of at least 1e-05 in at least one of the samples will be kept
species.max.value <- apply(tax.profiles, 1, max)
f.idx <- which(species.max.value > 1e-05)
tax.profiles_filtered <- tax.profiles[f.idx,]
tax.profiles_filtered <- tax.profiles_filtered[-1,]

#Wilcoxon - pairwise comparison 
p.vals <- rep_len(1, nrow(tax.profiles_filtered))
names(p.vals) <- rownames(tax.profiles_filtered)

#order metadata so that order of rownames matches order of colnames of rel data 
meta_mouse$order = match(rownames(meta_mouse), colnames(mouse_OTUS))
meta_mouse_ordered = meta_mouse[with(meta_mouse, order(order)), ]
meta_mouse_ordered$timepoint.x = as.factor(meta_mouse_ordered$timepoint.x)
meta_mouse_ordered$TM = as.factor(meta_mouse_ordered$TM)
meta_mouse_ordered$hygiene_level = as.factor(meta_mouse_ordered$hygiene_level)


stopifnot(all(rownames(meta_mouse_ordered) == colnames(tax.profiles)))
#DOES NOT WORK!!!
for (i in rownames(tax.profiles_filtered)){
  x <- tax.profiles_filtered[1,]
  y <- meta_mouse_ordered$TM
  t <- wilcox.test(x~y)
  p.vals[i] <- t$p.value
}
sort(p.vals)

sc.obj <- siamcat(feat=tax.profiles[,which(meta_mouse_ordered$TM == '0')], meta=meta_mouse_ordered, 
                  label='hygiene_level', case='High Hygiene')

sc.obj <- filter.features(sc.obj, filter.method = 'abundance', cutoff = 1e-05)

check.confounders(sc.obj, fn.plot = 'conf_check_TM0.pdf')

#check.associations(sc.obj, fn.plot = 'assoc.pdf')
sc.obj <- check.associations(siamcat = sc.obj, log.n0 = 1e-06, alpha = 0.1)

#association plot
association.plot(sc.obj, sort.by = 'fc', panels = NULL)

### Volcano plot associations
df.assoc <- associations(sc.obj)
write.csv2(df.assoc, 'L:/LovbeskyttetMapper/CLL-EPI/people/Tereza/Mouse_data Heidelberg/association_data_p0.1_TM0.csv')

options(ggrepel.max.overlaps = 50)
df.assoc %>% 
  ggplot(aes(x=fc, y=-log10(p.adj), label= rownames(df.assoc))) + 
  geom_point() + 
  geom_label_repel(aes(label=ifelse(-log10(p.adj)> 1, as.character(rownames(df.assoc)),''))) +
  xlab('Fold change')

#Check species one by one 
species <-'Helicobacter typhlonius [ref_mOTU_v25_04710]' #'Burkholderiales bacterium YL45 [meta_mOTU_v25_13139]'  'Firmicutes bacterium M10-2 [ref_mOTU_v25_07501]'  #'Lactobacillus intestinalis [ref_mOTU_v25_04390]' # Parabacteroides goldsteinii [ref_mOTU_v25_01679]
df.plot <- tibble(fuso=tax.profiles_filtered[species, ],
                  label=meta_mouse_ordered$hygiene_level)

ggplot(df.plot, (aes(x=label, y=fuso ))) + 
geom_boxplot() + 
xlab('') + 
ylab('Helicobacter typhlonius rel. ab.')

#################
# PCoA plot
library("vegan")
library("labdsv")

dist.mat <- vegdist(t(tax.profiles_filtered))
pco.results <- pco(dist.mat)

plot_data = as.data.frame(pco.results$points)
plot_data = merge(plot_data, meta_mouse_ordered, by = 0)
rownames(plot_data) = plot_data[,1]

percentage <- pco.results$eig[1:2]/sum(pco.results$eig)
percentage <- percentage*100
percentage <- paste0(sprintf(fmt='%.2f', percentage),'%')


plot_data %>%
ggplot(aes(x=V1, y=V2, col=TM, label = hygiene_level)) + 
  geom_point() + 
  geom_label() +
  xlab(paste0('PCo1 [', percentage[1], ']')) + 
  ylab(paste0('PCo2 [', percentage[2], ']'))

