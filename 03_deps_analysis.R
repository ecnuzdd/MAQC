# fot_data
# meta_data
meta_data.cp = meta_data[which(meta_data$Purpose == 'Replications'),]

meta_data.cp$uids = paste(
  meta_data.cp$Labcode, meta_data.cp$SInstrument, meta_data.cp$Batch, sep = '_'
)
uids_list = sort(unique(meta_data.cp$uids))

com_list = list(
  'D5_D6' = c('D5', 'D6'),
  'D5_F7' = c('D5', 'F7'),
  'P5_M8' = c('D5', 'M8'),
  'D6_F7' = c('D6', 'F7'),
  'D6_P8' = c('D6', 'M8'),
  'F7_M8' = c('F7', 'M8')
)


for(j in 1:length(com_list)){
  # j = 1
  com = com_list[[j]]
  com_names = names(com_list)[j]
  com.deps_df.list = list()
  for (i in 1:length(uids_list)) {
    if(T){
      # i = 1
      uid = uids_list[i]
      uid.index = which(meta_data.cp$uids==uid)
      uid.meta_data = meta_data.cp[uid.index,]
      com1.index = which(uid.meta_data$UPerson==com[1])
      com2.index = which(uid.meta_data$UPerson==com[2])
      com.meta_data = uid.meta_data[c(com1.index, com2.index),]
      com.grps = com.meta_data$UPerson
      com.fot = get_fot_df_by_expnames(fot_data, com.meta_data$Exp_code)
      com.deps_df = get_deps(com.fot, com.grps)
      colnames(com.deps_df) = c('GeneSymbol', paste(uid, '_FC', sep = ''))
      com.deps_df.list[[i]] = com.deps_df
      print(paste(i, ': ', uid, sep = ''))
    }
  }
  merge_deps = merge_deps_list_1(com.deps_df.list)
  merge_deps_path = normalizePath(
    file.path(
      BASE_DIR, 'deps_20200708',
      paste(com_names, '_DEPs.xlsx', sep = '')
    ),
    mustWork = F
  )
  write.xlsx(merge_deps, merge_deps_path, row.names = F)
  print('*******************************')
  print(paste(j, ': ', com_names, sep = ''))
  print('*******************************')
}

stable.genesymbol_list = NULL
for(j in 1:length(com_list)){
  if(T){
    # j = 1
    com_names = names(com_list)[j]
    merge_deps_path = normalizePath(
      file.path(
        BASE_DIR, 'deps_20200708',
        paste(com_names, '_DEPs.xlsx', sep = '')
      ),
      mustWork = F
    )
    merge_deps = read_excel(merge_deps_path)
    
    merge_deps.row_id = apply(merge_deps[,-1], 1, function(x){
      x = x[!is.na(x)]
      length(x)
    })
    merge_deps.row_id.up = apply(merge_deps[,-1], 1, function(x){
      x = x[!is.na(x)]
      x = x[x>2]
      length(x)
    })
    merge_deps.row_id.dn = apply(merge_deps[,-1], 1, function(x){
      x = x[!is.na(x)]
      x = x[x<1/2]
      length(x)
    })
    
    stable.up_index = which(merge_deps.row_id.up == 24)
    stable.dn_index = which(merge_deps.row_id.dn == 24)
    stable.genesymbol = as.vector(merge_deps$GeneSymbol[c(stable.up_index, stable.dn_index)])
    stable.genesymbol_list = c(stable.genesymbol_list, stable.genesymbol)
  }
  print('*******************************')
  print(paste(j, ': ', com_names, sep = ''))
  print('*******************************')
}

st_proteins = unique(stable.genesymbol_list)
snr_vector = NULL
pdf(
  normalizePath(
    file.path(
      BASE_DIR, 'deps_20200708',
      'stable_deps_pca.pdf'
    ),
    mustWork = F
  )
)
for (i in 1:24) {
  if(T){
    # i = 3
    uid = uids_list[i]
    uid.index = which(meta_data.cp$uids==uid)
    uid.meta_data = meta_data.cp[uid.index,]
    uid.fot = get_fot_df_by_expnames(fot_data, uid.meta_data$Exp_code)
    uid.fot.subset = uid.fot[match(st_proteins, uid.fot$GeneSymbol),]
    if(T){
      uid.df.mat = as.vector(uid.fot.subset[,-1])
      uid.df.mat.scale = t(scale(t(uid.df.mat)))
      mat = uid.df.mat.scale
      fot_group = uid.meta_data$UPerson
      uid.snr = round(snr_function(mat, fot_group),2)
      colors = rep(p_colors, each=3)
      pca <- prcomp(((t(mat))), center = F, scale = F)
      importance <- summary(pca)$importance
      pc1 = paste('PC1 (', round(importance[2,1]*100, 2), '%)', sep = '')
      pc2 = paste('PC2 (', round(importance[2,2]*100, 2), '%)', sep = '')
      pca_predict <- predict(pca)
      pca_predict.2d = pca_predict[,c(1,2)]
      main = paste(uid, '; DEPs', ' (', nrow(mat), ')\nSNR = ', uid.snr, sep = '') 
      plot(
        pca_predict.2d, 
        t = 'n', 
        main = main, 
        xlab = pc1, ylab = pc2, 
        xlim = get_xlim(pca_predict.2d), 
        ylim = get_ylim(pca_predict.2d),
      )
      points(pca_predict.2d, pch = 16, col = colors)
    }
  }
  snr_vector = c(snr_vector, uid.snr)
}
dev.off()


if(T){
  df = data.frame(
    UID = uids_list,
    SNR = snr_vector
  )
  df = df[order(df$SNR),]
  df$SNR = as.numeric(df$SNR)
  
  jco_col = pal_jco("default", alpha = 0.6)(10)#[1:5]
  jco_col.list = c(
    rep(jco_col[3], 4),
    rep(jco_col[5], 7),
    rep(jco_col[1], 9),
    rep(jco_col[4], 3),
    rep(jco_col[9], 1)
  )
  
  ggdotchart(
    df, x = "UID", y = "SNR",
    color = jco_col.list,                                # Color by groups
    palette = jco_col.list, #c("#0073C2FF", "#EFC000FF"), # Custom color palette
    # sorting = "ascending",                        # Sort value in descending order
    add = "segments",                             # Add segments from y = 0 to dots
    ggtheme = theme_pubr(),                        # ggplot2 theme
    # group = "Level",                                # Order by groups
    dot.size = 5,                                 # Large dot size
    label = round(df$SNR,1),
    add.params = list(color = jco_col.list, size = 2), # Change segment color and size,
    xlab = '',
    ylab = 'SNR',
    ylim = c(0, 30),
    font.label = list(color = "black", size = 10, vjust = -1)
  )
  
}








merge_deps_list_1 <- function(deps_list){
  deps_list_count = length(deps_list)
  merge_df = deps_list[[1]]
  for (i in 2:deps_list_count) {
    tmp_df = deps_list[[i]]
    merge_df = merge(merge_df, tmp_df, by = 'GeneSymbol', all = T)
  }
  return(merge_df)
}


get_deps <- function(df, grps){
  # df = com.fot
  # grps = com.grps
  
  fc = apply(df[,-1], 1, function(x, y){
    x_mean = tapply(x, y, mean)
    x_mean[1]/x_mean[2]
  }, y = grps)
  
  id1 = apply(df[,-1], 1, function(x, y){
    x_id = tapply(x, y, function(z){
      length(which(z>0))
    })
    x_id[1]
  }, y = grps)
  
  id2 = apply(df[,-1], 1, function(x, y){
    x_id = tapply(x, y, function(z){
      length(which(z>0))
    })
    x_id[2]
  }, y = grps)
  
  deps_index.up = which(fc>2 & id1>=2)
  deps_index.dn = which(fc<1/2 & id2>=2)
  
  deps.df = data.frame(
    GeneSymbol = df$GeneSymbol,
    FC = fc
  )[c(deps_index.up, deps_index.dn),]
  
  return(deps.df)
}







