source('D:/MAQC/20200508/Analysis/code_library/00_initialization.R')
# P5
p5.deep.expnames = meta_data$Exp_code[which(
  meta_data$Person == 'P5' & meta_data$Purpose == 'DeepCoverage'
)]
p5.deep.df = get_fot_df_by_expnames(fot_data, p5.deep.expnames)


p5.single.expnames = meta_data$Exp_code[which(
  meta_data$Person == 'P5' & meta_data$Purpose != 'DeepCoverage'
)]
p5.single.df = get_fot_df_by_expnames(fot_data, p5.single.expnames)


p5.deep.plot_df <- get_mean_and_rank(p5.deep.df)
p5.single.plot_df <- get_mean_and_rank(p5.single.df)


# P6
p6.deep.expnames = meta_data$Exp_code[which(
  meta_data$Person == 'P6' & meta_data$Purpose == 'DeepCoverage'
)]
p6.deep.df = get_fot_df_by_expnames(fot_data, p6.deep.expnames)


p6.single.expnames = meta_data$Exp_code[which(
  meta_data$Person == 'P6' & meta_data$Purpose != 'DeepCoverage'
)]
p6.single.df = get_fot_df_by_expnames(fot_data, p6.single.expnames)


p6.deep.plot_df <- get_mean_and_rank(p6.deep.df)
p6.single.plot_df <- get_mean_and_rank(p6.single.df)


# P7
p7.deep.expnames = meta_data$Exp_code[which(
  meta_data$Person == 'P7' & meta_data$Purpose == 'DeepCoverage'
)]
p7.deep.df = get_fot_df_by_expnames(fot_data, p7.deep.expnames)


p7.single.expnames = meta_data$Exp_code[which(
  meta_data$Person == 'P7' & meta_data$Purpose != 'DeepCoverage'
)]
p7.single.df = get_fot_df_by_expnames(fot_data, p7.single.expnames)


p7.deep.plot_df <- get_mean_and_rank(p7.deep.df)
p7.single.plot_df <- get_mean_and_rank(p7.single.df)


# P8
p8.deep.expnames = meta_data$Exp_code[which(
  meta_data$Person == 'P8' & meta_data$Purpose == 'DeepCoverage'
)]
p8.deep.df = get_fot_df_by_expnames(fot_data, p8.deep.expnames)


p8.single.expnames = meta_data$Exp_code[which(
  meta_data$Person == 'P8' & meta_data$Purpose != 'DeepCoverage'
)]
p8.single.df = get_fot_df_by_expnames(fot_data, p8.single.expnames)


p8.deep.plot_df <- get_mean_and_rank(p8.deep.df)
p8.single.plot_df <- get_mean_and_rank(p8.single.df)



# plot 
par(
  mar = c(5,5,3,4)
)
plot(
  NA, 
  # log = 'y',
  xlab = 'Rank', ylab = expression(log10(FOT[ave])),
  xlim = c(0, 10000), ylim = c(-5, 5), 
  yaxt = 'n', xaxt = 'n'
)
axis(
  1, at = c(1, 2000, 4000, 6000, 8000, 10000), labels = c(1, 2000, 4000, 6000, 8000, 10000)
)
axis(
  2, 
  at = seq(-4, 4, 2), labels = seq(-4, 4, 2)
)
alpha_per = 0.5
points(
  p5.single.plot_df$Rank,
  log10(p5.single.plot_df$FOT), 
  pch = 16, col = alpha(p_colors[1],alpha_per)
)
points(
  p6.single.plot_df$Rank,
  log10(p6.single.plot_df$FOT), 
  pch = 16, col = alpha(p_colors[2],alpha_per)
)
points(
  p7.single.plot_df$Rank,
  log10(p7.single.plot_df$FOT), 
  pch = 16, col = alpha(p_colors[3],alpha_per)
)
points(
  p8.single.plot_df$Rank,
  log10(p8.single.plot_df$FOT), 
  pch = 16, col = alpha(p_colors[4],alpha_per)
)



points(
  p5.deep.plot_df$Rank,
  log10(p5.deep.plot_df$FOT), 
  pch = 17, col = alpha(p_colors[1],alpha_per)
)
points(
  p6.deep.plot_df$Rank,
  log10(p6.deep.plot_df$FOT), 
  pch = 17, col = alpha(p_colors[2],alpha_per)
)
points(
  p7.deep.plot_df$Rank,
  log10(p7.deep.plot_df$FOT), 
  pch = 17, col = alpha(p_colors[3],alpha_per)
)
points(
  p8.deep.plot_df$Rank,
  log10(p8.deep.plot_df$FOT), 
  pch = 17, col = alpha(p_colors[4],alpha_per)
)

as.vector(p5.single.plot_df$GeneSymbol[1:10])
as.vector(p6.single.plot_df$GeneSymbol[1:10])
as.vector(p7.single.plot_df$GeneSymbol[1:10])
as.vector(p8.single.plot_df$GeneSymbol[1:10])

as.vector(p5.deep.plot_df$GeneSymbol[1:10])
as.vector(p6.deep.plot_df$GeneSymbol[1:10])
as.vector(p7.deep.plot_df$GeneSymbol[1:10])
as.vector(p8.deep.plot_df$GeneSymbol[1:10])

# Epstein-Barr virus infection
ebv_infection = 'BCL2|CD19|CD38|CDK1|CR2|CSNK2B|HLA-DPA1|HLA-DQA1|HLA-DQB1|HLA-DRB1|HLA-DRB3|HLA-DRB4|HLA-DRB5|HLA-E|HSPA1L|HSPA6|HSPB1|ICAM1|JAK3|LYN|PLCG2|POLR2C|POLR2E|POLR2H|PRKACA|PTMA|RB1|RELA|SPN|SYK|TRAF1|VIM|POLR1C|DDX58'
ebv_infection = unlist(strsplit(ebv_infection, split = '|', fixed  = T))
# Protein processing in endoplasmic reticulum
pp_in_er = 'BCL2|CALR|CAPN2|DDOST|DNAJB2|HSPA1L|HSPA5|HSPA6|DNAJB1|STT3A|LMAN1|P4HB|DNAJC3|RPN1|RPN2|RRBP1|SEL1L|SSR1|SSR4|HSP90B1|PDIA4|SEC24C|SEC24D|PREB|PDIA6|HYOU1|SEC24A|HSPH1|CKAP4|TRAM1|SEC61G|ERLEC1|SEC61A1|DNAJB11|DNAJC10|NPLOC4|ERO1B|TXNDC5|SYVN1'
pp_in_er = unlist(strsplit(pp_in_er, split = '|', fixed  = T))
# Adaptive Immune System
ais = 'B2M|BTK|CALR|CBLB|CD19|CD22|CD74|CTSC|CLTA|CTSH|CYBB|DYNC1LI2|HLA-DMB|HLA-DPA1|HLA-DQA1|HLA-DQB1|HLA-DRB1|HLA-DRB3|HLA-DRB4|HLA-DRB5|HLA-E|HSPA5|ICAM1|INPP5D|ITGB1|ITGB7|ITPR2|LCK|LMO7|LYN|CD99|NCF4|NFATC2|PLCG2|PPP2R5C|PRKACA|PRKCB|PSMB6|PSMD9|PTPN11|RAC1|REL|RELA|STIM1|SYK|CUL2|IFITM1|AP1S2|BCL10|SEC24C|SEC24D|PSME3|ATG7|CD226|SEC24A|MALT1|KIF2C|CD300A|SEC61G|CD274|BLNK|SEC61A1|ERAP1|TUBA8|PAG1|RNF213|SLAMF7|ERAP2|KLC2|UBA5|TUBB6|SLAMF6|PIK3AP1|DTX3L|TUBB2B'
# lymphocyte activation
la = 'ANXA1|B2M|BCL2|BTK|CAMK4|CD19|CD22|CD38|CD47|CD74|CDK6|CR1|CR2|DDOST|FLOT2|GSN|MSH6|HLA-DMB|HLA-DPA1|HLA-E|HMGB3|ICAM1|INPP5D|ITGB1|JAK3|LCK|LGALS1|LYN|MNDA|NFATC2|PNP|PAWR|PLCG2|POU2F2|PRKCB|PRKDC|PTPN11|RAC1|RAC2|SPN|VAMP2|VAMP7|SYK|LAT2|BCL10|SART1|DNAJA3|EBI3|HSPH1|MALT1|CD300A|IKZF3|PYCARD|CD274|BLNK|MZB1|CYRIB|DOCK10|PCID2|PAG1|PREX1|SLAMF7|SAMSN1|SEMA4A|SLAMF6|IGLL5'
la = unlist(strsplit(la, split = '|', fixed  = T))

df.ebv_infection = p8.deep.plot_df[na.omit(match(ebv_infection,p8.deep.plot_df$GeneSymbol)),]
df.pp_in_er = p8.deep.plot_df[na.omit(match(pp_in_er,p8.deep.plot_df$GeneSymbol)),]
df.la = p8.deep.plot_df[na.omit(match(la,p8.deep.plot_df$GeneSymbol)),]


ebv_infection = c('PTMA', 'HLA-DQA1', 'HLA-DRB1', 'HLA-DRB5', 'HLA-DQB1', 'VIM', 'HLA-E', 'HLA-DPA1', 'HSPA1L', 'HLA-DRB3')
la = c('CD38', 'IGLL5', 'CDK6', 'ANXA1', 'JAK3', 'CD274', 'CD47', 'BTK', 'LCK')
pp_in_er = c('HSPA5', 'HSP90B', 'CALR', 'PDIA6', 'PDIA4', 'P4HB', 'TXNDC5', 'RPN1','HSPA1L')


df.ebv_infection = p8.deep.plot_df[na.omit(match(ebv_infection,p8.deep.plot_df$GeneSymbol)),]
df.pp_in_er = p8.deep.plot_df[na.omit(match(pp_in_er,p8.deep.plot_df$GeneSymbol)),]
df.la = p8.deep.plot_df[na.omit(match(la,p8.deep.plot_df$GeneSymbol)),]

points(
  df.ebv_infection$Rank, log10(df.ebv_infection$FOT), pch = 16, col = '#FF6600'
)


points(
  df.pp_in_er$Rank, log10(df.pp_in_er$FOT), pch = 16, col = 'red'
)

points(
  df.la$Rank, log10(df.la$FOT), pch = 16, col = 'firebrick'
)




points(
  df.ebv_infection$Rank, log10(df.ebv_infection$FOT), pch = 16, col = '#009966'
)


points(
  df.pp_in_er$Rank, log10(df.pp_in_er$FOT), pch = 16, col = '#336666'
)

points(
  df.la$Rank, log10(df.la$FOT), pch = 16, col = '#336699'
)




################################################################################

plot(
  NA, 
  # log = 'y',
  xlab = 'Rank', ylab = expression(log10(FOT[ave])),
  xlim = c(0, 8000), ylim = c(-5, 5), 
  yaxt = 'n', xaxt = 'n'
)
axis(
  1, at = c(1, 2000, 4000, 6000, 8000, 10000), labels = c(1, 2000, 4000, 6000, 8000, 10000)
)
axis(
  2, 
  at = seq(-4, 4, 2), labels = seq(-4, 4, 2)
)
alpha_per = 0.5
points(
  p5.single.plot_df$Rank,
  log10(p5.single.plot_df$FOT), 
  pch = 16, col = alpha(p_colors[1],alpha_per)
)
points(
  p6.single.plot_df$Rank,
  log10(p6.single.plot_df$FOT), 
  pch = 16, col = alpha(p_colors[2],alpha_per)
)
points(
  p7.single.plot_df$Rank,
  log10(p7.single.plot_df$FOT), 
  pch = 16, col = alpha(p_colors[3],alpha_per)
)
points(
  p8.single.plot_df$Rank,
  log10(p8.single.plot_df$FOT), 
  pch = 16, col = alpha(p_colors[4],alpha_per)
)

ebv_infection = c('PTMA', 'HLA-DQA1', 'HLA-DRB1', 'HLA-DRB5', 'HLA-DQB1', 'VIM', 'HLA-E', 'HLA-DPA1', 'HSPA1L', 'HLA-DRB3')
la = c('CD38', 'IGLL5', 'CDK6', 'ANXA1', 'JAK3', 'CD274', 'CD47', 'BTK', 'LCK')
pp_in_er = c('HSPA5', 'HSP90B', 'CALR', 'PDIA6', 'PDIA4', 'P4HB', 'TXNDC5', 'RPN1','HSPA1L')


df.ebv_infection = p8.deep.plot_df[na.omit(match(ebv_infection,p8.deep.plot_df$GeneSymbol)),]
df.pp_in_er = p8.deep.plot_df[na.omit(match(pp_in_er,p8.deep.plot_df$GeneSymbol)),]
df.la = p8.deep.plot_df[na.omit(match(la,p8.deep.plot_df$GeneSymbol)),]

points(
  df.ebv_infection$Rank, log10(df.ebv_infection$FOT), pch = 16, col = '#FF6600'
)


points(
  df.pp_in_er$Rank, log10(df.pp_in_er$FOT), pch = 16, col = 'red'
)

points(
  df.la$Rank, log10(df.la$FOT), pch = 16, col = 'firebrick'
)




points(
  df.ebv_infection$Rank, log10(df.ebv_infection$FOT), pch = 16, col = '#009966'
)


points(
  df.pp_in_er$Rank, log10(df.pp_in_er$FOT), pch = 16, col = '#336666'
)

points(
  df.la$Rank, log10(df.la$FOT), pch = 16, col = '#336699'
)
























