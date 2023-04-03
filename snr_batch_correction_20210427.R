library(factoextra)
library(dendextend)

meta_data.cp = meta_data[which(meta_data$Purpose == 'Replications'),]
meta_data.cp.uids = paste(
  meta_data.cp$Labcode, 
  meta_data.cp$SInstrument, 
  meta_data.cp$Batch, sep = '_'
)
uids = sort(unique(meta_data.cp.uids))
uids_count = length(uids)
uid.df.scale.list = list()
snr_list = NULL
snr_list_new = NULL
for (i in 1:uids_count) {
  uid = uids[i] 
  uid.metadata = meta_data.cp[meta_data.cp.uids==uid,]
  # uid.metadata = uid.metadata[uid.metadata$Person!='P8',]
  # uid.metadata = uid.metadata[uid.metadata$Person=='P5' | uid.metadata$Person=='P6',]
  uid.metadata = uid.metadata[order(uid.metadata$Person),]
  uid.expnames = uid.metadata$Exp_code
  uid.df = get_fot_df_by_expnames(fot_data, uid.expnames)
  uid.df.genesymbol = as.vector(uid.df$GeneSymbol)
  uid.df.mat = as.vector(uid.df[,-1])
  uid.df.mat.scale = t(scale(t(uid.df.mat)))
  
  mat = uid.df.mat.scale
  fot_group = uid.metadata$Person
  snr = round(snr_function(mat, fot_group),2)
  colors = rep(p_colors, each=3)
  pca <- prcomp(((t(mat))), center = F, scale = F)
  importance <- summary(pca)$importance
  pc1 = paste('PC1 (', round(importance[2,1]*100, 2), '%)', sep = '')
  pc2 = paste('PC2 (', round(importance[2,2]*100, 2), '%)', sep = '')
  pca_predict <- predict(pca)
  pca_predict.2d = pca_predict[,c(1,2)]
  main = paste(uid, '; All', ' Proteins (', nrow(mat), ')\nSNR = ', snr, sep = '') 
  plot(
    pca_predict.2d, 
    t = 'n', 
    main = main, 
    xlab = pc1, ylab = pc2, 
    xlim = get_xlim(pca_predict.2d), 
    ylim = get_ylim(pca_predict.2d),
  )
  points(pca_predict.2d, pch = 16, col = colors)
  print(i)
  snr_list = c(snr_list, snr)
  
  snr_new = round(snrdb_function(mat, fot_group),2)
  snr_list_new = c(snr_list_new, snr_new)
}
names(snr_list) = uids

df = data.frame(
  UID = uids,
  SNR_old = snr_list,
  SNR = snr_list_new
)
df = df[order(df$SNR),]
df$SNR = as.numeric(df$SNR)





Company =  apply(df, 1, function(x){
  strsplit(x, '_')[[1]][1]
})
Instrument = apply(df, 1, function(x){
  strsplit(x, '_')[[1]][2]
})
Batch = apply(df, 1, function(x){
  strsplit(x, '_')[[1]][3]
})


df1 = data.frame(
  df, Company, Instrument, Batch, Label = paste('T', seq(1,26), sep = '')
)
write.csv(
  df1, '26_SNR_Figures_20210427.csv', row.names = F
)

















if(T){
  df = read_excel('26_SNR_Figures_20210427.xlsx')
  df = df[!is.na(df$Code),]
  df = df[order(df$Code),]
  df$SNR = as.numeric(df$SNR)
  
  jco_col = pal_jco("default", alpha = 0.6)(10)#[1:5]
  jco_col.list = c(
    rep(jco_col[1], 2),
    rep(jco_col[2], 5),
    rep(jco_col[4], 2)
  )
  jco_col.list = rep(jco_col[1], nrow(df))
  
  ggdotchart(
    df, x = "UID", y = "SNR",
    color = jco_col.list,                                # Color by groups
    palette = jco_col.list, #c("#0073C2FF", "#EFC000FF"), # Custom color palette
    sorting = 'none',                        # Sort value in descending order
    add = "segments",                             # Add segments from y = 0 to dots
    ggtheme = theme_pubr(),                        # ggplot2 theme
    # group = "Level",                                # Order by groups
    dot.size = 5,                                 # Large dot size
    label = round(df$SNR,1),
    add.params = list(color = jco_col.list, size = 2), # Change segment color and size,
    order = df$UID,
    xlab = '',
    ylab = 'SNR',
    ylim = c(0, 30),
    font.label = list(color = "black", size = 10, vjust = -1)
  )
  
  
}










# 24

df = read_excel('26_SNR_Figures_20210427.xlsx')
df = df[is.na(df$Complementary),]
df = df[order(df$SNR),]
df$SNR = as.numeric(df$SNR)


alpha_cols = rgb(
  raw_cols_rgb[1L, ], raw_cols_rgb[2L, ], 
  raw_cols_rgb[3L, ], alpha = 0.7 * 255L, names = 'zdd_firebrick', 
  maxColorValue = 255L
)

jco_col.list = c(
  rep(jco_col[3], 4),
  rep(jco_col[5], 5),
  rep(jco_col[1], 3),
  rep(jco_col[4], 9),
  rep(alpha_cols, 3)
)

ggdotchart(
  df, x = "UID", y = "SNR",
  color = jco_col.list,                                # Color by groups
  palette = jco_col.list, #c("#0073C2FF", "#EFC000FF"), # Custom color palette
  # sorting = "ascending",                        # Sort value in descending order
  add = "segments",                             # Add segments from y = 0 to dots
  ggtheme = theme_pubr(),                        # ggplot2 theme
  # group = "Level",                                # Order by groups
  dot.size = 6,                                 # Large dot size
  label = round(df$SNR,1),
  add.params = list(color = jco_col.list, size = 1.5), # Change segment color and size,
  xlab = '',
  ylab = 'SNR',
  ylim = c(0, 30),
  font.label = list(color = "black", size = 10, vjust = -1)
)



snrdb_function<-function(exprMat, group){
  
  library(data.table)
  
  exprMat = uid.df.mat.scale
  group = uid.metadata$Person
  
  IDs<- colnames(exprMat)
  IDs.group.mat<-data.table(
    IDs=IDs,
    group=group) 
  
  pca_prcomp <- prcomp(t(exprMat),retx=T)
  pcs <- as.data.frame(predict(pca_prcomp))
  pcs$Sample_id <- rownames(pcs)
  
  dt.perc.pcs <- data.table(PCX=1:nrow(pcs),
                            Percent=summary(pca_prcomp)$importance[2,],
                            AccumPercent=summary(pca_prcomp)$importance[3,])
  
  dt.dist <- data.table(ID.A = rep(IDs,each=length(IDs)),
                        ID.B = rep(IDs,time=length(IDs)))
  
  dt.dist$group.A <- IDs.group.mat[match(dt.dist$ID.A,IDs.group.mat$IDs)]$group
  dt.dist$group.B <- IDs.group.mat[match(dt.dist$ID.B,IDs.group.mat$IDs)]$group
  
  dt.dist[,Type:=ifelse(ID.A==ID.B,'Same',
                        ifelse(group.A==group.B,'Intra','Inter'))]
  
  dt.dist[,Dist:=(dt.perc.pcs[1]$Percent*(pcs[ID.A,1]-pcs[ID.B,1])^2+dt.perc.pcs[2]$Percent*(pcs[ID.A,2]-pcs[ID.B,2])^2)]
  
  dt.dist.stats <- dt.dist[,.(Avg.Dist=mean(Dist)),by=.(Type)]
  setkey(dt.dist.stats,Type)
  signoise <- dt.dist.stats['Inter']$Avg.Dist/dt.dist.stats['Intra']$Avg.Dist  
  
  signoise_db <- 10*log10(signoise)
  return(signoise_db)
  
}

