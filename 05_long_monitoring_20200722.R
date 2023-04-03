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

#### Identification #### 
if(T){
  # peptide
  df1 = meta_data.lm[meta_data.lm$Preservation=='Peptide',]
  df1$Batch = factor(df1$Batch, levels = seq(1,15))
  df1$Sample = df1$Person
  p1 <- ggplot(data=df1, aes(x=Batch, y=GPs, fill=Sample)) 
  p1 <- p1 + geom_bar(stat="identity", color="black", position=position_dodge(), width = 0.8)
  p1 <- p1 + scale_fill_manual(values=p_colors) + theme_classic()
  p1 <- p1 + ggtitle('Peptide standards')
  
  df2 = meta_data.lm[meta_data.lm$Preservation=='Protein',]
  df2$Batch = factor(df2$Batch, levels = seq(1,15))
  df2$Sample = df2$Person
  p2 <- ggplot(data=df2, aes(x=Batch, y=GPs, fill=Sample)) 
  p2 <- p2 + geom_bar(stat="identity", color="black", position=position_dodge(0.8), width = 0.8)
  p2 <- p2 + scale_fill_manual(values=p_colors) + theme_classic()
  p2 <- p2 + ggtitle('Protein standards')
  
  ggarrange(p1, p2, nrow = 2)
}

#### Reproducibility #### 
if(T){
  if(T){
    meta_data.lm.subset = meta_data.lm[meta_data.lm$Preservation=='Peptide',]
    
    px_rep_list = NULL
    for (i in 1:4) {
      # px = 'P5'
      px = p_flags[i]
      px.lm.expnames = meta_data.lm.subset$Exp_code[meta_data.lm.subset$Person==px]
      px.lm.df = get_fot_df_by_expnames(fot_data, px.lm.expnames)
      px.lm.df.reproducibility = get_lm_reproducibility(px.lm.df)
      px_rep_list[[i]] = px.lm.df.reproducibility
    }
    names(px_rep_list) = p_flags
    
    px.g.df = data.frame(
      Batch = rep(seq(2,15),4),
      Reproducibility = c(
        px_rep_list$P5[,1], px_rep_list$P6[,1], px_rep_list$P7[,1], px_rep_list$P8[,1]
      ),
      Sample = rep(p_flags, each = 14),
      Flag = 'Global'
    )
    
    px.l.df = data.frame(
      Batch = rep(seq(2,15),4),
      Reproducibility = c(
        px_rep_list$P5[,2], px_rep_list$P6[,2], px_rep_list$P7[,2], px_rep_list$P8[,2]
      ),
      Sample = rep(p_flags, each = 14),
      Flag = 'Low'
    )
    
    px.m.df = data.frame(
      Batch = rep(seq(2,15),4),
      Reproducibility = c(
        px_rep_list$P5[,3], px_rep_list$P6[,3], px_rep_list$P7[,3], px_rep_list$P8[,3]
      ),
      Sample = rep(p_flags, each = 14),
      Flag = 'Medium'
    )
    
    px.h.df = data.frame(
      Batch = rep(seq(2,15),4),
      Reproducibility = c(
        px_rep_list$P5[,4], px_rep_list$P6[,4], px_rep_list$P7[,4], px_rep_list$P8[,4]
      ),
      Sample = rep(p_flags, each = 14),
      Flag = 'High'
    )
    
    px.df = rbind.data.frame(
      px.l.df, px.m.df, px.h.df, px.g.df
    )
    px.df$Flag = factor(px.df$Flag, levels = c('Low', 'Medium', 'High', 'Global'))
    px.df$Batch = factor(px.df$Batch, levels = seq(2, 15))
    
    df1 = px.df#[px.df$Flag=='Low',]
    p1 <- ggplot(data=df1, aes(x=Batch, y=Reproducibility, group = Sample, color = Sample)) #  fill=Sample,
    p1 <- p1 + geom_point(show.legend = T)
    p1 <- p1+ scale_fill_manual(values=p_colors) + theme_classic() + facet_grid(Flag ~ .) 
    p1 <- p1 + ggtitle('Peptide standards (Reference = Batch 1)')
    p1 
  }
  
  if(T){
    meta_data.lm.subset = meta_data.lm[meta_data.lm$Preservation=='Protein',]
    
    px_rep_list = NULL
    for (i in 1:4) {
      # px = 'P5'
      px = p_flags[i]
      px.lm.expnames = meta_data.lm.subset$Exp_code[meta_data.lm.subset$Person==px]
      px.lm.df = get_fot_df_by_expnames(fot_data, px.lm.expnames)
      px.lm.df.reproducibility = get_lm_reproducibility(px.lm.df)
      px_rep_list[[i]] = px.lm.df.reproducibility
    }
    names(px_rep_list) = p_flags
    
    px.g.df = data.frame(
      Batch = rep(seq(2,15),4),
      Reproducibility = c(
        px_rep_list$P5[,1], px_rep_list$P6[,1], px_rep_list$P7[,1], px_rep_list$P8[,1]
      ),
      Sample = rep(p_flags, each = 14),
      Flag = 'Global'
    )
    
    px.l.df = data.frame(
      Batch = rep(seq(2,15),4),
      Reproducibility = c(
        px_rep_list$P5[,2], px_rep_list$P6[,2], px_rep_list$P7[,2], px_rep_list$P8[,2]
      ),
      Sample = rep(p_flags, each = 14),
      Flag = 'Low'
    )
    
    px.m.df = data.frame(
      Batch = rep(seq(2,15),4),
      Reproducibility = c(
        px_rep_list$P5[,3], px_rep_list$P6[,3], px_rep_list$P7[,3], px_rep_list$P8[,3]
      ),
      Sample = rep(p_flags, each = 14),
      Flag = 'Medium'
    )
    
    px.h.df = data.frame(
      Batch = rep(seq(2,15),4),
      Reproducibility = c(
        px_rep_list$P5[,4], px_rep_list$P6[,4], px_rep_list$P7[,4], px_rep_list$P8[,4]
      ),
      Sample = rep(p_flags, each = 14),
      Flag = 'High'
    )
    
    px.df = rbind.data.frame(
      px.l.df, px.m.df, px.h.df, px.g.df
    )
    px.df$Flag = factor(px.df$Flag, levels = c('Low', 'Medium', 'High', 'Global'))
    px.df$Batch = factor(px.df$Batch, levels = seq(2, 15))
    
    if(T){
      df1 = px.df#[px.df$Flag=='Low',]
      p2 <- ggplot(data=df1, aes(x=Batch, y=Reproducibility, group = Sample, color = Sample)) #  
      p2 <- p2 + geom_point(shape=17, show.legend = T)
      p2 <- p2+ scale_fill_manual(values=p_colors) + theme_classic() + facet_grid(Flag ~ .) 
      p2 <- p2 + ggtitle('Protein standards (Reference = Batch 1)')
      p2
    }
  }
  
  ggarrange(p1, p2, ncol = 2)
}

#### Cluster #### 
if(T){# 120
  # meta_data
  meta_data.lm = meta_data[which(
    meta_data$Purpose == 'LongMonitoring' #& meta_data$Preservation == 'protein'
  ),]
  # meta_data.lm = meta_data.lm[meta_data.lm$Preservation=='Peptide',]
  batch = meta_data.lm$Batch
  batch_str = ifelse(batch<10, paste('0', batch, sep = ''), batch)
  meta_data.lm.uids = paste(
    meta_data.lm$Preservation, batch_str, sep = '_'
  )
  uids = sort(unique(meta_data.lm.uids))
  uids_count = length(uids)
  
  
  uid.df.scale.list = list()
  for (i in 1:uids_count) {
    uid = uids[i] 
    uid.metadata = meta_data.lm[meta_data.lm.uids==uid,]
    uid.metadata = uid.metadata[order(uid.metadata$Person),]
    uid.expnames = uid.metadata$Exp_code
    uid.df = get_fot_df_by_expnames(fot_data, uid.expnames)
    uid.df.genesymbol = as.vector(uid.df$GeneSymbol)
    uid.df.mat = as.vector(uid.df[,-1])
    uid.df.mat.scale = t(scale(t(uid.df.mat)))
    
    uid.df.scale = data.frame(
      GeneSymbol = uid.df.genesymbol,
      uid.df.mat.scale
    )
    uid.df.scale.list[[i]] = uid.df.scale
  }
  
  # merge
  lm.scale_df = uid.df.scale.list[[1]]
  for (i in 2:uids_count) {
    tmp_df = uid.df.scale.list[[i]]
    lm.scale_df = merge(lm.scale_df, tmp_df, by = 'GeneSymbol', all = T)
  }
  lm.scale_df = na.omit(lm.scale_df)
  lm.scale_df.mat = lm.scale_df[,-1]
  
  # ph <- get_fot_df_by_expnames(fot_data, meta_data.lm$Exp_code)
  # lm.scale_df.mat <- as.matrix(ph[,-1])
  library(pheatmap)
  meta_data.lm.new = meta_data.lm[match(colnames(lm.scale_df.mat), meta_data.lm$Exp_code),]
  annotation_col = data.frame(
    Sample = meta_data.lm.new$UPerson,
    Standards = meta_data.lm.new$Preservation,
    Batch = factor(meta_data.lm.new$Batch, levels = seq(1,15))
  )
  rownames(annotation_col) = colnames(lm.scale_df.mat)
  batch_col = viridis(15)
  names(batch_col) = seq(1,15)
  annotation_colors = list(
    Sample = c('D5'=p_colors[1], 'D6'=p_colors[2], 'F7'=p_colors[3], 'M8'=p_colors[4]),
    Standards = c(Peptide="#A6CEE3", Protein="#1F78B4"),
    Batch = batch_col
  )
  
  pdf(normalizePath(
    file.path(BASE_DIR, 'temp_figures', '05', '04_pheatmap_pep.pdf'), mustWork = F
  ), onefile = F)
  pheatmap(
    lm.scale_df.mat,
    scale = 'row',
    annotation_col = annotation_col,
    annotation_colors = annotation_colors,
    show_rownames = F,
    show_colnames = F
  )
  dev.off()
  
  
  
  
  
}




if(T){ # 60 Protein
  # meta_data
  meta_data.lm = meta_data[which(
    meta_data$Purpose == 'LongMonitoring' #& meta_data$Preservation == 'protein'
  ),]
  meta_data.lm = meta_data.lm[meta_data.lm$Preservation=='Protein',]
  batch = meta_data.lm$Batch
  batch_str = ifelse(batch<10, paste('0', batch, sep = ''), batch)
  meta_data.lm.uids = paste(
    meta_data.lm$Preservation, batch_str, sep = '_'
  )
  uids = sort(unique(meta_data.lm.uids))
  uids_count = length(uids)
  
  
  uid.df.scale.list = list()
  for (i in 1:uids_count) {
    uid = uids[i] 
    uid.metadata = meta_data.lm[meta_data.lm.uids==uid,]
    uid.metadata = uid.metadata[order(uid.metadata$Person),]
    uid.expnames = uid.metadata$Exp_code
    uid.df = get_fot_df_by_expnames(fot_data, uid.expnames)
    uid.df.genesymbol = as.vector(uid.df$GeneSymbol)
    uid.df.mat = as.vector(uid.df[,-1])
    uid.df.mat.scale = t(scale(t(uid.df.mat)))
    
    uid.df.scale = data.frame(
      GeneSymbol = uid.df.genesymbol,
      uid.df.mat.scale
    )
    uid.df.scale.list[[i]] = uid.df.scale
  }
  
  # merge
  lm.scale_df = uid.df.scale.list[[1]]
  for (i in 2:uids_count) {
    tmp_df = uid.df.scale.list[[i]]
    lm.scale_df = merge(lm.scale_df, tmp_df, by = 'GeneSymbol', all = T)
  }
  lm.scale_df = na.omit(lm.scale_df)
  lm.scale_df.mat = lm.scale_df[,-1]
  
  # ph <- get_fot_df_by_expnames(fot_data, meta_data.lm$Exp_code)
  # lm.scale_df.mat <- as.matrix(ph[,-1])
  library(pheatmap)
  meta_data.lm.new = meta_data.lm[match(colnames(lm.scale_df.mat), meta_data.lm$Exp_code),]
  annotation_col = data.frame(
    Sample = meta_data.lm.new$UPerson,
    # Standards = meta_data.lm.new$Preservation,
    Batch = factor(meta_data.lm.new$Batch, levels = seq(1,15))
  )
  rownames(annotation_col) = colnames(lm.scale_df.mat)
  batch_col = viridis(15)
  names(batch_col) = seq(1,15)
  annotation_colors = list(
    Sample = c('D5'=p_colors[1], 'D6'=p_colors[2], 'F7'=p_colors[3], 'M8'=p_colors[4]),
    # Standards = c(Peptide="#A6CEE3", Protein="#1F78B4"),
    Batch = batch_col
  )
  
  pdf(normalizePath(
    file.path(BASE_DIR, 'temp_figures', '05', '04_pheatmap_pro.pdf'), mustWork = F
  ), onefile = T)
  pheatmap(
    lm.scale_df.mat,
    scale = 'row',
    annotation_col = annotation_col,
    annotation_colors = annotation_colors,
    show_rownames = F,
    show_colnames = F
  )
  
  
  lm.colors = rep(p_colors, 15)
  pca <- prcomp(((t(lm.scale_df.mat))), center = F, scale = F)
  importance <- summary(pca)$importance
  pc1 = paste('PC1 (', round(importance[2,1]*100, 2), '%)', sep = '')
  pc2 = paste('PC2 (', round(importance[2,2]*100, 2), '%)', sep = '')
  pca_predict <- predict(pca)
  pca_predict.2d = pca_predict[,c(1,2)]
  main = paste('PCA; Z-score within batch') 
  plot(
    NA, 
    # t = 'n', 
    # col = lm.colors,
    main = main, 
    xlab = pc1, ylab = pc2, 
    xlim = get_xlim(pca_predict.2d), 
    # xlim = c(-20,30), 
    ylim = get_ylim(pca_predict.2d),
  )
  points(pca_predict.2d, col = lm.colors, pch = 16)
  
  dev.off()
  
  
}



if(T){ # 60 Peptide
  # meta_data
  meta_data.lm = meta_data[which(
    meta_data$Purpose == 'LongMonitoring' #& meta_data$Preservation == 'protein'
  ),]
  meta_data.lm = meta_data.lm[meta_data.lm$Preservation=='Peptide',]
  batch = meta_data.lm$Batch
  batch_str = ifelse(batch<10, paste('0', batch, sep = ''), batch)
  meta_data.lm.uids = paste(
    meta_data.lm$Preservation, batch_str, sep = '_'
  )
  uids = sort(unique(meta_data.lm.uids))
  uids_count = length(uids)
  
  
  uid.df.scale.list = list()
  for (i in 1:uids_count) {
    uid = uids[i] 
    uid.metadata = meta_data.lm[meta_data.lm.uids==uid,]
    uid.metadata = uid.metadata[order(uid.metadata$Person),]
    uid.expnames = uid.metadata$Exp_code
    uid.df = get_fot_df_by_expnames(fot_data, uid.expnames)
    uid.df.genesymbol = as.vector(uid.df$GeneSymbol)
    uid.df.mat = as.vector(uid.df[,-1])
    uid.df.mat.scale = t(scale(t(uid.df.mat)))
    
    uid.df.scale = data.frame(
      GeneSymbol = uid.df.genesymbol,
      uid.df.mat.scale
    )
    uid.df.scale.list[[i]] = uid.df.scale
  }
  
  # merge
  lm.scale_df = uid.df.scale.list[[1]]
  for (i in 2:uids_count) {
    tmp_df = uid.df.scale.list[[i]]
    lm.scale_df = merge(lm.scale_df, tmp_df, by = 'GeneSymbol', all = T)
  }
  lm.scale_df = na.omit(lm.scale_df)
  lm.scale_df.mat = lm.scale_df[,-1]
  
  # ph <- get_fot_df_by_expnames(fot_data, meta_data.lm$Exp_code)
  # lm.scale_df.mat <- as.matrix(ph[,-1])
  library(pheatmap)
  meta_data.lm.new = meta_data.lm[match(colnames(lm.scale_df.mat), meta_data.lm$Exp_code),]
  annotation_col = data.frame(
    Sample = meta_data.lm.new$UPerson,
    # Standards = meta_data.lm.new$Preservation,
    Batch = factor(meta_data.lm.new$Batch, levels = seq(1,15))
  )
  rownames(annotation_col) = colnames(lm.scale_df.mat)
  batch_col = viridis(15)
  names(batch_col) = seq(1,15)
  annotation_colors = list(
    Sample = c('D5'=p_colors[1], 'D6'=p_colors[2], 'F7'=p_colors[3], 'M8'=p_colors[4]),
    # Standards = c(Peptide="#A6CEE3", Protein="#1F78B4"),
    Batch = batch_col
  )
  
  pdf(normalizePath(
    file.path(BASE_DIR, 'temp_figures', '05', '04_pheatmap_pep.pdf'), mustWork = F
  ), onefile = T)
  pheatmap(
    lm.scale_df.mat,
    scale = 'row',
    annotation_col = annotation_col,
    annotation_colors = annotation_colors,
    show_rownames = F,
    show_colnames = F
  )
  
  
  
  
  
  lm.colors = rep(p_colors, 15)
  pca <- prcomp(((t(lm.scale_df.mat))), center = F, scale = F)
  importance <- summary(pca)$importance
  pc1 = paste('PC1 (', round(importance[2,1]*100, 2), '%)', sep = '')
  pc2 = paste('PC2 (', round(importance[2,2]*100, 2), '%)', sep = '')
  pca_predict <- predict(pca)
  pca_predict.2d = pca_predict[,c(1,2)]
  main = paste('PCA; Z-score within batch') 
  plot(
    NA, 
    # t = 'n', 
    # col = lm.colors,
    main = main, 
    xlab = pc1, ylab = pc2, 
    xlim = get_xlim(pca_predict.2d), 
    # xlim = c(-20,30), 
    ylim = get_ylim(pca_predict.2d),
  )
  points(pca_predict.2d, col = lm.colors, pch = 16)
  dev.off()
  
}























if(F){
  
  # PCA
  
  uid.df.scale.list = list()
  for (i in 1:uids_count) {
    uid = uids[i] 
    uid.metadata = meta_data.lm[meta_data.lm.uids==uid,]
    uid.metadata = uid.metadata[order(uid.metadata$Person),]
    uid.expnames = uid.metadata$Exp_code
    uid.df = get_fot_df_by_expnames(fot_data, uid.expnames)
    uid.df.genesymbol = as.vector(uid.df$GeneSymbol)
    uid.df.mat = as.vector(uid.df[,-1])
    uid.df.mat.scale = t(scale(t(uid.df.mat)))
    
    uid.df.scale = data.frame(
      GeneSymbol = uid.df.genesymbol,
      uid.df.mat.scale
    )
    uid.df.scale.list[[i]] = uid.df.scale
  }
  
  # merge
  lm.scale_df = uid.df.scale.list[[1]]
  for (i in 2:uids_count) {
    tmp_df = uid.df.scale.list[[i]]
    lm.scale_df = merge(lm.scale_df, tmp_df, by = 'GeneSymbol', all = T)
  }
  lm.scale_df = na.omit(lm.scale_df)
  
  # PCA
  df_mat_subset = as.matrix(lm.scale_df[,-1])
  df_mat_subset.metadata = meta_data.lm[match(colnames(df_mat_subset),meta_data.lm$Exp_code),]
  
  lm.colors = rep(p_colors, 30)
  pca <- prcomp(((t(df_mat_subset))), center = F, scale = F)
  importance <- summary(pca)$importance
  pc1 = paste('PC1 (', round(importance[2,1]*100, 2), '%)', sep = '')
  pc2 = paste('PC2 (', round(importance[2,2]*100, 2), '%)', sep = '')
  pca_predict <- predict(pca)
  pca_predict.2d = pca_predict[,c(1,2)]
  main = paste('PCA; Z-score within batch') 
  plot(
    NA, 
    # t = 'n', 
    # col = lm.colors,
    main = main, 
    xlab = pc1, ylab = pc2, 
    xlim = get_xlim(pca_predict.2d), 
    # xlim = c(-20,30), 
    ylim = get_ylim(pca_predict.2d),
  )
  points(pca_predict.2d[1:60,], col = lm.colors[1:60], pch = 16)
  points(pca_predict.2d[61:120,], col = lm.colors[61:120], pch = 17)
  a = pca_predict.2d[1:60,]
  
}