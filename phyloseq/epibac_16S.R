#### Epibac bacteria community analysis ####

## Preparing the data

# read OTU data, calling this 'dat' 
# use the 'stringsAsFactors' argument to read OTU names as character
dat <- read.delim('ASV_counts_nosingle.txt', stringsAsFactors=F)

dat$timepoint <- factor(dat$timepoint, levels=c('T0','T1','T2'))
dat$temp <- factor(dat$temp, levels=c('A','H'))
dat$treatment <- factor(dat$treatment, levels=c('CTRL', 'BMC', 'PATH', 'BOTH'))
dat$group <- factor(dat$group, levels=c('T0_ACTRL', 'T0_ABMC', 'T0_APATH', 'T0_ABOTH','T0_HCTRL', 'T0_HBMC', 'T0_HPATH', 'T0_HBOTH',
                                        'T1_ACTRL', 'T1_ABMC','T1_APATH', 'T1_ABOTH','T1_HCTRL', 'T1_HBMC','T1_HPATH', 'T1_HBOTH',
                                        'T2_ACTRL', 'T2_ABMC','T2_APATH', 'T2_ABOTH','T2_HCTRL', 'T2_HBMC','T2_HPATH', 'T2_HBOTH'))
str(dat[, 1:5])

# read in the taxonomy, use 'stringsAsFactors=F' to keep text as 'character'
tax <- read.delim('tax_table_16S_nosingle.txt', stringsAsFactors=F)

## Creating phyloseq object

# We need to specify the different components to create a [phyloseq](https://github.com/joey711/phyloseq) object.

# load phyloseq, dplyr(to manipulate data)
library(phyloseq)
library(dplyr)

# extract the OTU data from 'dat' and convert to matrix format
# using 'starts_with' with 'select' allows us to chose all columns containing OTU counts
otus <- as.matrix(select(dat, starts_with('ASV')))

# rows need to be named with the sample names from 'dat', which we can do directly because they are in the same order
rownames(otus) <- dat$sample

# extract the sample data from from 'dat' and keep as dataframe format (mix of character and numeric variables)
samps <- select(dat, sample:group)

# rows need to be named with the sample names from 'dat', which we can do directly because they are in the same order
rownames(samps) <- dat$sample

# extract the taxonomy info for each OTU from 'tax' and convert to matrix format
taxonomy <- as.matrix(select(tax, kingdom:genus))

# rows need to be named with the OTU names from 'tax', which we can do directly because they are in the same order
rownames(taxonomy) <- tax$sample

# merge the three objects into a single phyloseq object using 'merge_phyloseq'
# each function nexted within the call to 'merge_phyloseq' creates a special object for that type of data
# because of how the OTU matrix is oriented, we have to specify 'taxa_are_rows=F' 
phy <- merge_phyloseq(otu_table(otus, taxa_are_rows=F), 
                sample_data(samps),
                tax_table(taxonomy))

# remove extra objects to keep workspace tidy
rm(otus, samps, taxonomy)

# remove archaea/eukarya
phy <- subset_taxa(phy, kingdom=="k_Bacteria")

# remove OTUs with abundance zero
phy <- prune_samples(sample_sums(phy)>0, phy)
                     
# remove OTUs that hit to Chloroplast or Mitochondria (because they are showing eukarya)
phy <- subset_taxa(phy, order != "o_Chloroplast")
phy <- subset_taxa(phy, family != "f_Mitochondria")

## Alpha diversity

library(microbiome)
library(ggpubr)

# Rarefaction using the phyloseq function rarefy_even_depth
set.seed(1)
phy.rarefied<- rarefy_even_depth(phy, sample.size = min(sample_sums(phy)), replace=FALSE, rngseed = TRUE)
phy
phy.rarefied

# we can look at the alpha diversity in each of the samples
set.seed(1)
rich<-estimate_richness(phy.rarefied)
even<-evenness(phy.rarefied)

# subset just the high temps to do 4 plots showing alpha diversity progression over time
phy.h <- subset_samples(phy.rarefied, temp=="H")
rich.h <- estimate_richness(phy.h)
even.h <- evenness(phy.h)

# we can plot at the alpha diversity in each of the samples
l<-plot_richness(phy.h, x="group", measures=c("Shannon")) +
  geom_boxplot(aes(fill=(sample_data(phy.h)$group)), outlier.shape = NA) +
  geom_point(size=1)+
  facet_wrap(~treatment, scales = "free_x", nrow = 1, strip.position = 'top',
             labeller = labeller(treatment = c('CTRL' = 'Control','BMC' = 'Cobetia',
                                           'PATH' = 'Vibrio','BOTH' = 'Both'))) +
  scale_x_discrete(labels=c('T0','T1','T2')) +
  theme_bw() +
  theme(axis.text=element_text(size=14, family = "sans", colour = "black"),axis.title=element_text(size=14, family = "sans", colour = "black"))+
  theme(legend.position="none") + 
  theme(axis.text.x=element_text(hjust=0.5,vjust=1)) +
  theme(axis.title.x = element_blank(), strip.text = element_blank()) +
  scale_y_continuous(name="Shannon diversity", limits=c(2,5.7), breaks = seq(2,6,1)) +
  theme(strip.placement = 'outside', strip.background = element_blank(), strip.text = element_text(size=14, family = "sans", colour = "black"))
l

ggsave("shannon_heats.pdf", plot = l, height = 5, width = 8)

# plot the evenness
k<-ggplot(even.h,aes(x=(sample_data(phy.h)$group), y=pielou,))+
  geom_boxplot(aes(fill=(sample_data(phy.h)$group)), outlier.shape = NA) +
  geom_point(size=1)+
  theme_bw()+
  facet_wrap((sample_data(phy.h)$treatment), scales = "free_x", nrow = 1, strip.position = 'top',
             labeller = labeller(treatment = c('CTRL' = 'Control','BMC' = 'Cobetia',
                                               'PATH' = 'Vibrio','BOTH' = 'Both'))) +
  scale_x_discrete(labels=c('T0','T1','T2')) +
  theme(axis.text=element_text(size=14, family = "sans", colour = "black"),axis.title=element_text(size=14, family = "sans", colour = "black"))+
  theme(legend.position="none")+ 
  theme(axis.text.x=element_text(hjust = 0.5, vjust = 1))+
  theme(axis.title.x = element_blank()) +
  scale_y_continuous(name="Evenness", limits=c(0.3,0.8), breaks = seq(0.3,0.8,0.1)) +
  theme(strip.placement = 'outside', strip.background = element_blank(), strip.text = element_text(size=14, family = "sans", colour = "black"))
k

ggsave("evenness_heats.pdf", plot = k, height = 5, width = 8)

# statistics for comparing across timepoints
hist(rich.h$Shannon, main="Shannon index", xlab="")
ggqqplot(rich.h$Shannon, main="Shannon index", xlab="")
shapiro.test(rich.h$Shannon) ##Not normal W = 0.88128, p-value = 0.000165

# subset for each treatment
phy.h.ctrl <- subset_samples(phy.h, treatment=="CTRL")
rich.h.ctrl <-estimate_richness(phy.h.ctrl)
even.h.ctrl <- evenness(phy.h.ctrl)
phy.h.bmc <- subset_samples(phy.h, treatment=="BMC")
rich.h.bmc <-estimate_richness(phy.h.bmc)
even.h.bmc <- evenness(phy.h.bmc)
phy.h.path <- subset_samples(phy.h, treatment=="PATH")
rich.h.path <-estimate_richness(phy.h.path)
even.h.path <- evenness(phy.h.path)
phy.h.both <- subset_samples(phy.h, treatment=="BOTH")
rich.h.both <-estimate_richness(phy.h.both)
even.h.both <- evenness(phy.h.both)

kruskal.test(even.h.ctrl$pielou ~ sample_data(phy.h.ctrl)$timepoint) # substitute richness/evenness
FSA::dunnTest(even.h.ctrl$pielou ~ sample_data(phy.h.ctrl)$timepoint, 
              data = sample_data(phy.h.ctrl)$timepoint, method = "bh") # substitute richness/evenness

## Abundance analysis barplots

library(ggplot2)

# create a dummy variable with both the categories to make a relative abundance plot
phy.dum<-phy
sample_data(phy.dum)$dummy <- paste0(sample_data(phy.dum)$timepoint, sample_data(phy.dum)$temp, sample_data(phy.dum)$treatment)
phym = merge_samples(phy.dum, "dummy")

# repair the variables that you just destroyed
sample_data(phym)$timepoint <- levels(sample_data(phy.dum)$timepoint)[get_variable(phym, "timepoint")]
sample_data(phym)$temp <- levels(sample_data(phy.dum)$temp)[get_variable(phym, "temp")]
sample_data(phym)$treatment <- levels(sample_data(phy.dum)$treatment)[get_variable(phym, "treatment")]

# to do relative abundance (out of 100%), we need to also transform the counts to percents
phy.prop <- transform_sample_counts(phym, function(otu) 100 * otu/sum(otu))

# convert phyloseq object to dataframe to make figures in ggplot2
phy.propdf<-psmelt(phy.prop)

# keep the same levels on the x axis
phy.propdf$genus <- factor(phy.propdf$genus,levels=rev(unique(phy.propdf$genus)))

# filter to just the top 20 taxa at the family level.

# agglomerate at family level
phyg <- tax_glom(phy.prop, "family", NArm = FALSE)

# get just the top 20 order in the samples
top20otus = names(sort(taxa_sums(phyg), TRUE)[1:20])
taxtab20 = cbind(tax_table(phyg), phy20 = NA)
taxtab20[top20otus, "fam20"] <- as(tax_table(phyg)[top20otus, "family"], 
    "character")
tax_table(phyg) <- tax_table(taxtab20)

# prune samples
phyg20 = prune_taxa(top20otus, phyg)

# convert phyloseq object to dataframe to make figures in ggplot2
phyg20df<-psmelt(phyg20)

# keep the same levels on the x axis
phyg20df$phy20 <- factor(phyg20df$phy20,levels=rev(unique(phyg20df$phy20)))
phyg20df$Sample <- factor(phyg20df$Sample, 
                          levels = c("T0ACTRL", "T0ABMC", "T0APATH", "T0ABOTH", "T0HCTRL", "T0HBMC", "T0HPATH", "T0HBOTH", 
                                     "T1ACTRL", "T1ABMC", "T1APATH", "T1ABOTH", "T1HCTRL", "T1HBMC", "T1HPATH", "T1HBOTH", 
                                     "T2ACTRL", "T2ABMC", "T2APATH", "T2ABOTH", "T2HCTRL", "T2HBMC", "T2HPATH", "T2HBOTH"))

# you can define the colors
colours <- c("#808080", "#0075DC","#FFCC99","#F0A3FF","#94FFB5","#8F7C00","#9DCC00","#C20088","#003380","#FFA405",
             "#FFA8BB","#426600","#FF0010","#5EF1F2","#00998F","#740AFF","#990000","#FF8000", "#993F00","#4C005C","#2BCE48")

# plot
n <- ggplot(phyg20df, aes(x=Sample, fill = phy20, y=Abundance)) +
  geom_bar(aes(color=phy20), stat = 'identity', position = 'stack', colour = NA) +
  facet_wrap(~timepoint, nrow = 1, scales = "free_x", strip.position = 'top') +
  scale_x_discrete(labels=c('A_Control','A_Cobetia','A_Vibrio','A_Both',
                            'H_Control','H_Cobetia','H_Vibrio','H_Both')) +
  theme_bw() +
  theme(axis.text=element_text(size=12, family = "sans", colour = "black"),axis.title=element_text(size=14, family = "sans", colour = "black"))+
  theme(legend.text = element_text(size=12, family = "sans", colour = "black"), legend.title = element_text(size=12, family = "sans", colour = "black")) +
  scale_fill_manual(name = "Top 20 families", values = colours) +
  scale_y_continuous(name= 'Abundance (%)', expand = c(0,0), limits=c(0,101)) +
  theme(axis.text.x=element_text(angle = 45, vjust = 1, hjust = 1))+
  theme(axis.title.x = element_blank()) +
  theme(strip.placement = 'outside', strip.text = element_text(size=14, family = "sans", colour = "black"))
n

ggsave("top20_families.pdf", plot = n, height = 6, width = 9)

## Beta diversity

# Check the beta diversity by making a PCoA comparing bacterial communities between experimental groups.

library(ggplot2)

# we can just use a subset of the samples so that we don't clog the plot
phy.t0 = subset_samples(phy, timepoint == "T0") # substitute for t0/t1/t2

# this is to normalize the abundances of the samples
phy.t0.norm = transform_sample_counts(phy.t0, function(x) 100 * x/sum(x)) # substitute for t0/t1/t2

# calculate the distance matrix
phy.dist.t0 = phyloseq::distance(phy.t0.norm, method="bray") # substitute for t0/t1/t2

# calculate the ordination (PCoA) using distance (from bray)
ord0 <- ordinate(phy.t0.norm, method="PCoA", distance=phy.dist.t0) # substitute for t0/t1/t2

# see the eigen value of the axes
plot_scree(ord0) # substitute for t0/t1/t2

# plot
p <- plot_ordination(phy.t0.norm, ord2, color = "temp", shape = "treatment")+
  geom_point(size=2, alpha=0.8)+
  theme_bw()+
  stat_ellipse(aes(color=phy.t0.norm@sam_data[["temp"]], group = phy.02.norm@sam_data[["temp"]]))
p

# calculate the beta diversity
library(vegan)
adonis2(phy.dist.t0 ~ sample_data(phy.t0.norm)$temp * sample_data(phy.t0.norm)$treatment, permutations = 1000) # substitute for t0/t1/t2

# run pairwise adonis on all the combinations to see which groups are different from each other
anova(betadisper(phy.dist.t0, sample_data(phy.t0.norm)$temp)) # substitute for t0/t1/t2
pairwiseAdonis::pairwise.adonis(phy.dist.t0, sample_data(phy.t0.norm)$group, perm = 1000, p.adjust.m = 'BH') # substitute for t0/t1/t2

## Multi boxplot showing the relative abundances of the incoulated bacteria

# to do relative abundance (out of 100%), we need to also transform the counts to percents
phy.rel <- transform_sample_counts(phy, function(otu) 100 * otu/sum(otu))

Taxa = c("ASV176", "ASV319")
# "ASV36", "ASV37", "ASV57" # V. coralliilyticus ASVs with 100% identity to inoculated bacteria
# "ASV176", "ASV319" # Cobetia. sp ASVs with 100% identity to inoculated bacteria
allTaxa = taxa_names(phy.rel)
myTaxa <- allTaxa[(allTaxa %in% Taxa)]
phy.rel.inoc = prune_taxa(myTaxa, phy.rel)

phy.rel.inoc<-tax_glom(phy.rel.inoc, "genus", NArm = FALSE)

phy.reldf.inoc<-psmelt(phy.rel.inoc)

phy.reldf.inoc<-phy.reldf.inoc

# plot the boxplot
q<-ggplot(phy.reldf.inoc,aes(x=group, fill = genus, y=Abundance))+
  geom_boxplot(aes(fill=genus), outlier.shape = NA) +
  geom_point(size=0.8,position=position_dodge(width=0.75)) +
  facet_wrap(~timepoint, scales = "free_x") +
  theme_bw()+
  theme(axis.text.x=element_text(angle = 45, vjust = 1, hjust = 1))
q

# we can just use a subset of the samples so we can do stats within timepoints
phy.rel.inoc.1 = subset_samples(phy.rel.inoc, timepoint == "T0") # substitute for t0/t1/t2
phy.rel.inoc.1 = subset_taxa(phy.rel.inoc.1, genus == "g_Vibrio") # substitute for Vibrio/Cobetia
phy.rel.inoc.1a = subset_samples(phy.rel.inoc.1, temp == "A") # substitute for ambient/heat-stressed

# convert to data frame
phy.reldf.inoc.1<-psmelt(phy.rel.inoc.1a)

# testing for normality
hist(phy.reldf.inoc.1<b$Abundance)
shapiro.test(phy.reldf.inoc.1<$Abundance)

# significance testing
kruskal.test(phy.reldf.inoc.1$Abundance,phy.reldf.inoc.1$treatment)
FSA::dunnTest(phy.reldf.inoc.1$Abundance~phy.reldf.inoc.1$treatment, data = phy.reldf.inoc.1, method='bh')
