# meta_data
meta_data.lm = meta_data[which(
  meta_data$Purpose == 'LongMonitoring' #& meta_data$Preservation == 'protein'
),]

batch = meta_data.lm$Batch
batch_str = ifelse(batch<10, paste('0', batch, sep = ''), batch)
meta_data.lm.uids = paste(
  meta_data.lm$Preservation, batch_str, sep = '_'
)
uids = sort(unique(meta_data.lm.uids))
uids_count = length(uids)


# peptide
px = 'P5'
px.df1 = meta_data.lm[meta_data.lm$Preservation=='Peptide' & meta_data.lm$Person==px,]
px.df1$Batch = factor(px.df1$Batch, levels = seq(1,15))
px.lm.df1 = get_fot_df_by_expnames(fot_data, px.df1$Exp_code)
df = px.lm.df1
ref_index = which(df[,2]>0)
df1 = df[ref_index,]
df1.row_id = apply(df1[,-1], 1, function(x){
  length(which(x>0))
})
df2 = df1[df1.row_id==15,]
df2.row = nrow(df2)
lm_time = seq(1,15)
pval_lst = c()
cor_lst = c()
fc_lst = c()
fc_min_lst = c()
for (i in 1:df2.row) {
  x = unlist(df2[i,-1])
  x_test = cor.test(x, lm_time, method = c("pearson", "kendall", "spearman")[3])
  x_pval = x_test$p.value
  x_cor = x_test$estimate
  pval_lst = c(pval_lst, x_pval)
  cor_lst = c(cor_lst, x_cor)
  fc = x/x[1]
  fc = fc[-1]
  fc = fc[fc>0]
  fc_min = min(fc)
  fc = fc[length(fc)]
  fc_lst = c(fc_lst, fc)
  fc_min_lst = c(fc_min_lst, fc_min) 
}
res_df = data.frame(
  Symbol = df2$GeneSymbol,
  ID = apply(df2[,-1], 1, function(x){
    length(which(x>0))
  }),
  Pvalue = pval_lst,
  Cor = cor_lst,
  FC = fc_lst,
  FOT = apply(df2[,-1], 1, function(x){
    mean(x)
  })
)

degrade_index = which(
  res_df$FC<0.5 & res_df$Cor< 0 & res_df$Pvalue<0.1
)
res_df1 = res_df[degrade_index,]

constant_index = which(
  res_df$FC>0.75 & res_df$FC<1.5 & 
  res_df$Cor>0.6 & res_df$Pvalue<0.05
)
res_df2 = res_df[constant_index,]


i = which(df2$GeneSymbol=='TAGLN2')
x = unlist(df2[i,-1])
           
lm_proteins1 = c('HLA-DRA', 'SNRPC', 'KRT1', 'KRT9')
lm_proteins2 = c('EEF2', 'CCT8', 'PSMC6', 'DDB1')
lm_proteins = c(lm_proteins2, lm_proteins1)


if(T){
  lm_proteins = lm_proteins1
  px_list = list()
  p_flags = sort(unique(meta_data.lm$UPerson))
  melt_df = NULL
  for (i in 1:4) {
    px = p_flags[i]
    px.df = meta_data.lm[meta_data.lm$Preservation=='Peptide' & meta_data.lm$UPerson==px,]
    px.df$Batch = factor(px.df$Batch, levels = seq(1,15))
    px.fot = get_fot_df_by_expnames(fot_data, px.df$Exp_code)
    px.fot.subset = px.fot[match(lm_proteins, px.fot$GeneSymbol),]
    px.fot.subset_melt = melt(px.fot.subset)
    colnames(px.fot.subset_melt) = c('GeneSymbol', 'Month', 'FOT')
    px.fot.subset_melt$Month = rep(seq(1,15), each = 4)
    px.fot.subset_melt$Sample = px
    px_list[[i]] = px.fot.subset_melt
    melt_df = rbind.data.frame(melt_df, px.fot.subset_melt)
  }
  melt_df$Month = factor(melt_df$Month, levels = seq(1,15))
  melt_df$GeneSymbol = factor(melt_df$GeneSymbol, levels = lm_proteins)
  
  
  
  df1 = melt_df
  p <- ggplot(df1, aes(x = Month, y = (FOT)))
  p <- p + geom_point(aes(shape=GeneSymbol, color = Sample),  size = 3)
  # p <- p + geom_line()
  p <- p + scale_shape_manual(values=c(15, 16, 17, 18)) 
  p <- p #+ scale_y_log10() 
  p <- p + ylab('FOT') + theme_classic2() +  ggtitle('Degradative')
  p <- p + facet_grid(GeneSymbol ~ .)
  print(p)
  p1 <- p
}



if(T){
  lm_proteins = lm_proteins2
  px_list = list()
  p_flags = sort(unique(meta_data.lm$UPerson))
  melt_df = NULL
  for (i in 1:4) {
    px = p_flags[i]
    px.df = meta_data.lm[meta_data.lm$Preservation=='Peptide' & meta_data.lm$UPerson==px,]
    px.df$Batch = factor(px.df$Batch, levels = seq(1,15))
    px.fot = get_fot_df_by_expnames(fot_data, px.df$Exp_code)
    px.fot.subset = px.fot[match(lm_proteins, px.fot$GeneSymbol),]
    px.fot.subset_melt = melt(px.fot.subset)
    colnames(px.fot.subset_melt) = c('GeneSymbol', 'Month', 'FOT')
    px.fot.subset_melt$Month = rep(seq(1,15), each = 4)
    px.fot.subset_melt$Sample = px
    px_list[[i]] = px.fot.subset_melt
    melt_df = rbind.data.frame(melt_df, px.fot.subset_melt)
  }
  melt_df$Month = factor(melt_df$Month, levels = seq(1,15))
  melt_df$GeneSymbol = factor(melt_df$GeneSymbol, levels = lm_proteins)
  
  
  
  df1 = melt_df
  p <- ggplot(df1, aes(x = Month, y = (FOT)))
  p <- p + geom_point(aes(shape=GeneSymbol, color = Sample),  size = 3)
  # p <- p + geom_line()
  p <- p + scale_shape_manual(values=c(15, 16, 17, 18)) 
  p <- p + scale_y_log10()
  p <- p + ylab('FOT') + theme_classic2() +  ggtitle('Constant')
  p <- p + facet_grid(GeneSymbol ~ .)
  print(p)
  p2 <- p
  
  ggarrange(
    p2, p1, nrow = 2
  )
}

















for (i in 1:8) {
  protein = lm_proteins[i]
  df1 = melt_df[melt_df$GeneSymbol==protein,]
  p <- ggplot(df1, aes(x = Month, y = FOT))
  p <- p + geom_point(aes(shape=GeneSymbol, color = Sample),  size = 3)
  # p <- p + geom_line()
  p <- p + scale_shape_manual(values=c(16)) 
  p <- p + scale_y_log10()
  p <- p + ylab('FOT') + theme_classic2() +  ggtitle(protein)
  p <- p + facet_grid(GeneSymbol ~ .)
  print(p)
}



degrade_df_list = list()
constant_df_list = list()
for (j in 1:4) {
  # peptide
  px = p_flags[j]
  px.df1 = meta_data.lm[meta_data.lm$Preservation=='Peptide' & meta_data.lm$UPerson==px,]
  px.df1$Batch = factor(px.df1$Batch, levels = seq(1,15))
  px.lm.df1 = get_fot_df_by_expnames(fot_data, px.df1$Exp_code)
  df = px.lm.df1
  ref_index = which(df[,2]>0)
  df1 = df[ref_index,]
  df1.row_id = apply(df1[,-1], 1, function(x){
    length(which(x>0))
  })
  df2 = df1[df1.row_id==15,]
  df2.row = nrow(df2)
  lm_time = seq(1,15)
  pval_lst = c()
  cor_lst = c()
  fc_lst = c()
  fc_min_lst = c()
  for (i in 1:df2.row) {
    x = unlist(df2[i,-1])
    x_test = cor.test(x, lm_time, method = c("pearson", "kendall", "spearman")[3])
    x_pval = x_test$p.value
    x_cor = x_test$estimate
    pval_lst = c(pval_lst, x_pval)
    cor_lst = c(cor_lst, x_cor)
    fc = x/x[1]
    fc = fc[-1]
    fc = fc[fc>0]
    fc_min = min(fc)
    fc = fc[length(fc)]
    fc_lst = c(fc_lst, fc)
    fc_min_lst = c(fc_min_lst, fc_min) 
  }
  res_df = data.frame(
    Symbol = df2$GeneSymbol,
    ID = apply(df2[,-1], 1, function(x){
      length(which(x>0))
    }),
    Pvalue = pval_lst,
    Cor = cor_lst,
    FC = fc_lst,
    FOT = apply(df2[,-1], 1, function(x){
      mean(x)
    })
  )
  degrade_index = which(
    res_df$FC<0.5 & res_df$Cor< (-0.5) & res_df$Pvalue<0.05
  )
  degrade_df = res_df[degrade_index,]
  degrade_df_list[[j]] = degrade_df
  
  # constant_index = which(
  #   res_df$Cor>(0) & res_df$Cor<0.1 #& res_df$FC>1/1.2 & res_df$FC<1.2 #& res_df$Pvalue<0.05 #
  # )
  # constant_df = res_df[constant_index,]
  # constant_df_list[[j]] = constant_df
}
degrade_df_merge = degrade_df_list[[1]]
for (i in 2:4) {
  degrade_df_merge = merge(degrade_df_merge, degrade_df_list[[i]], by = 'Symbol', all = T)
}
degrade_df_merge = na.omit(degrade_df_merge)



if(T){
  
  px.df2 = meta_data.lm[meta_data.lm$Preservation=='Peptide',]
  px.df2.fot = get_fot_df_by_expnames(fot_data, px.df2$Exp_code)
  px.df2.fot.rowid = apply(px.df2.fot[,-1], 1, function(x){length(which(x>0))})
  px.df2.fot.subset = px.df2.fot[px.df2.fot.rowid==60,]
  
  constant_df_list = list()
  for (j in 1:4) {
    # peptide
    px = p_flags[j]
    px.df1 = meta_data.lm[meta_data.lm$Preservation=='Peptide' & meta_data.lm$UPerson==px,]
    px.df1$Batch = factor(px.df1$Batch, levels = seq(1,15))
    px.lm.df1 = get_fot_df_by_expnames(px.df2.fot.subset, px.df1$Exp_code)
    df2 = px.lm.df1
    res_df = data.frame(
      Symbol = df2$GeneSymbol,
      ID = apply(df2[,-1], 1, function(x){
        length(which(x>0))
      }),
      FOT = apply(df2[,-1], 1, function(x){
        mean(x)
      }),
      CV = apply(df2[,-1], 1, function(x){
        sd(x)/mean(x)
      })
    )
    
    
    constant_index = which(
      res_df$CV<0.25
    )
    constant_df = res_df[constant_index,]
    constant_df_list[[j]] = constant_df

  }
  
  constant_df_merge = constant_df_list[[1]]
  for (i in 2:4) {
    constant_df_merge = merge(constant_df_merge, constant_df_list[[i]], by = 'Symbol', all = T)
  }
  constant_df_merge = na.omit(constant_df_merge)
}






if(T){
  lm_proteins = lm_proteins1
  px_list = list()
  p_flags = sort(unique(meta_data.lm$UPerson))
  melt_df = NULL
  for (i in 1:4) {
    px = p_flags[i]
    px.df = meta_data.lm[meta_data.lm$Preservation=='Peptide' & meta_data.lm$UPerson==px,]
    px.df$Batch = factor(px.df$Batch, levels = seq(1,15))
    px.fot = get_fot_df_by_expnames(fot_data, px.df$Exp_code)
    px.fot.subset = px.fot[match(lm_proteins, px.fot$GeneSymbol),]
    px.fot.subset_melt = melt(px.fot.subset)
    colnames(px.fot.subset_melt) = c('GeneSymbol', 'Month', 'FOT')
    px.fot.subset_melt$Month = rep(seq(1,15), each = 4)
    px.fot.subset_melt$Sample = px
    px_list[[i]] = px.fot.subset_melt
    melt_df = rbind.data.frame(melt_df, px.fot.subset_melt)
  }
  melt_df$Month = factor(melt_df$Month, levels = seq(1,15))
  melt_df$GeneSymbol = factor(melt_df$GeneSymbol, levels = lm_proteins)
  
  
  
  df1 = melt_df#[melt_df$GeneSymbol=='SNRPC',]
  p <- ggplot(df1, aes(x = Month, y = (FOT)))
  p <- p + geom_point(aes(shape=GeneSymbol, color = Sample),  size = 3)
  # p <- p + geom_line()
  p <- p + scale_shape_manual(values=c(15, 16, 17, 18)) 
  p <- p + scale_y_log10() 
  p <- p + ylab('FOT') + theme_classic2() +  ggtitle('Degradative')
  p <- p + facet_grid(GeneSymbol ~ ., scales = "free")
  print(p)
  p1 <- p
}



