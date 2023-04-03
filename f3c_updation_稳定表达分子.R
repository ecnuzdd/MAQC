# meta_data
if(T){
  px = 'P5'
  meta_data.px = meta_data[which(meta_data$Person == px & meta_data$Purpose == 'Replications'),]
  meta_data.px.uid = paste(meta_data.px$Labcode, meta_data.px$SInstrument, meta_data.px$Loading, meta_data.px$Batch, sep = '_')
  uids = sort(unique(meta_data.px.uid))
  
  
  px.fot = get_fot_df_by_expnames(fot_data, meta_data.px$Exp_code)
  px.fot.row_id = apply(px.fot[,-1], 1, function(x){
    length(which(x>0))
  })
  px.fot.subset = px.fot[px.fot.row_id==72,]
  
  
  uid.mean.df = NULL
  for (i in 1:24) {
    uid = uids[i]
    uid.expcode = meta_data.px$Exp_code[meta_data.px.uid==uid]
    uid.fot_df = get_fot_df_by_expnames(px.fot.subset, uid.expcode)
    uid.mean = apply(uid.fot_df[,-1], 1, function(x){mean(x)})
    uid.mean.df = cbind(uid.mean.df, uid.mean)
  }
  colnames(uid.mean.df) = uids
  uid.mean.df = data.frame(
    GeneSymbol = px.fot.subset$GeneSymbol,
    uid.mean.df
  )
  
  px.df = data.frame(
    gene = uid.mean.df$GeneSymbol,
    cv = apply(uid.mean.df[,-1], 1, function(x){sd(x)/mean(x)}),
    mean = apply(uid.mean.df[,-1], 1, function(x){mean(x)})
  )
  
  px.df.subset = px.df[px.df$cv<0.25,]
  px.df.subset1 = px.df.subset[px.df.subset$mean>1000 & px.df.subset$mean<10000,]
  
  
  gene_subset = c(
    'HUWE1',
    'HSPH1', 'ACTN4', 'PARP1',
    'MYH9', 
    'HSPA5',
    'LDHB',
    'PKM'
  )
  gene_subset.df1 = uid.mean.df[match(gene_subset, uid.mean.df$GeneSymbol),]
  colnames(gene_subset.df1) = c('GeneSymbol', paste('S', seq(1, 24), sep = ''))
  
  library(reshape2)
  df0 = gene_subset.df1
  df = melt(df0)
  colnames(df) = c('GeneSymbol', 'Site', 'Value')
  df$Site = factor(df$Site, levels =  paste('S', seq(1, 24), sep = ''))
  df$GeneSymbol = factor(df$GeneSymbol, levels = gene_subset)
  # df$Value = log10(df$Value)
  p <- ggplot(df, aes(x = Site, y = Value))
  p <- p + geom_point(aes(shape=GeneSymbol), color = p_colors[1], fill = p_colors[3], size = 3)
  p <- p + scale_shape_manual(values=c(15, 0, 16, 1, 17, 2, 18, 5)) + scale_y_log10()
  p <- p + ylab('FOT') + theme_classic2() +  ggtitle('P5')
  print(p)
  df5 = df
}








# meta_data
if(T){
  px = 'P6'
  meta_data.px = meta_data[which(meta_data$Person == px & meta_data$Purpose == 'Replications'),]
  meta_data.px.uid = paste(meta_data.px$Labcode, meta_data.px$SInstrument, meta_data.px$Loading, meta_data.px$Batch, sep = '_')
  uids = sort(unique(meta_data.px.uid))
  
  
  px.fot = get_fot_df_by_expnames(fot_data, meta_data.px$Exp_code)
  px.fot.row_id = apply(px.fot[,-1], 1, function(x){
    length(which(x>0))
  })
  px.fot.subset = px.fot[px.fot.row_id==72,]
  
  
  uid.mean.df = NULL
  for (i in 1:24) {
    uid = uids[i]
    uid.expcode = meta_data.px$Exp_code[meta_data.px.uid==uid]
    uid.fot_df = get_fot_df_by_expnames(px.fot.subset, uid.expcode)
    uid.mean = apply(uid.fot_df[,-1], 1, function(x){mean(x)})
    uid.mean.df = cbind(uid.mean.df, uid.mean)
  }
  colnames(uid.mean.df) = uids
  uid.mean.df = data.frame(
    GeneSymbol = px.fot.subset$GeneSymbol,
    uid.mean.df
  )
  
  px.df = data.frame(
    gene = uid.mean.df$GeneSymbol,
    cv = apply(uid.mean.df[,-1], 1, function(x){sd(x)/mean(x)}),
    mean = apply(uid.mean.df[,-1], 1, function(x){mean(x)})
  )
  
  px.df.subset = px.df[px.df$cv<0.25,]
  px.df.subset1 = px.df.subset[px.df.subset$mean>1000 & px.df.subset$mean<10000,]
  
  
  gene_subset = c(
    'HUWE1',
    'HSPH1', 'ACTN4', 'PARP1',
    'MYH9', 
    'HSPA5',
    'LDHB',
    'PKM'
  )
  gene_subset.df1 = uid.mean.df[match(gene_subset, uid.mean.df$GeneSymbol),]
  colnames(gene_subset.df1) = c('GeneSymbol', paste('S', seq(1, 24), sep = ''))
  
  library(reshape2)
  df0 = gene_subset.df1
  df = melt(df0)
  colnames(df) = c('GeneSymbol', 'Site', 'Value')
  df$Site = factor(df$Site, levels =  paste('S', seq(1, 24), sep = ''))
  df$GeneSymbol = factor(df$GeneSymbol, levels = gene_subset)
  # df$Value = log10(df$Value)
  p <- ggplot(df, aes(x = Site, y = Value))
  p <- p + geom_point(aes(shape=GeneSymbol), color = p_colors[2], fill = p_colors[3], size = 3)
  p <- p + scale_shape_manual(values=c(15, 0, 16, 1, 17, 2, 18, 5)) + scale_y_log10()
  p <- p + ylab('FOT') + theme_classic2() +  ggtitle('P6')
  print(p)
  df6 = df
}


if(T){
  px = 'P7'
  meta_data.px = meta_data[which(meta_data$Person == px & meta_data$Purpose == 'Replications'),]
  meta_data.px.uid = paste(meta_data.px$Labcode, meta_data.px$SInstrument, meta_data.px$Loading, meta_data.px$Batch, sep = '_')
  uids = sort(unique(meta_data.px.uid))
  
  
  px.fot = get_fot_df_by_expnames(fot_data, meta_data.px$Exp_code)
  px.fot.row_id = apply(px.fot[,-1], 1, function(x){
    length(which(x>0))
  })
  px.fot.subset = px.fot[px.fot.row_id==72,]
  
  
  uid.mean.df = NULL
  for (i in 1:24) {
    uid = uids[i]
    uid.expcode = meta_data.px$Exp_code[meta_data.px.uid==uid]
    uid.fot_df = get_fot_df_by_expnames(px.fot.subset, uid.expcode)
    uid.mean = apply(uid.fot_df[,-1], 1, function(x){mean(x)})
    uid.mean.df = cbind(uid.mean.df, uid.mean)
  }
  colnames(uid.mean.df) = uids
  uid.mean.df = data.frame(
    GeneSymbol = px.fot.subset$GeneSymbol,
    uid.mean.df
  )
  
  px.df = data.frame(
    gene = uid.mean.df$GeneSymbol,
    cv = apply(uid.mean.df[,-1], 1, function(x){sd(x)/mean(x)}),
    mean = apply(uid.mean.df[,-1], 1, function(x){mean(x)})
  )
  
  px.df.subset = px.df[px.df$cv<0.25,]
  px.df.subset1 = px.df.subset[px.df.subset$mean>1000 & px.df.subset$mean<10000,]
  
  
  gene_subset = c(
    'HUWE1',
    'HSPH1', 'ACTN4', 'PARP1',
    'MYH9', 
    'HSPA5',
    'LDHB',
    'PKM'
  )
  gene_subset.df1 = uid.mean.df[match(gene_subset, uid.mean.df$GeneSymbol),]
  colnames(gene_subset.df1) = c('GeneSymbol', paste('S', seq(1, 24), sep = ''))
  
  library(reshape2)
  df0 = gene_subset.df1
  df = melt(df0)
  colnames(df) = c('GeneSymbol', 'Site', 'Value')
  df$Site = factor(df$Site, levels =  paste('S', seq(1, 24), sep = ''))
  df$GeneSymbol = factor(df$GeneSymbol, levels = gene_subset)
  # df$Value = log10(df$Value)
  p <- ggplot(df, aes(x = Site, y = Value))
  p <- p + geom_point(aes(shape=GeneSymbol), color = p_colors[3], fill = p_colors[3], size = 3)
  p <- p + scale_shape_manual(values=c(15, 0, 16, 1, 17, 2, 18, 5)) + scale_y_log10()
  p <- p + ylab('FOT') + theme_classic2() +  ggtitle('P7')
  print(p)
  df7 = df
}


if(T){
  px = 'P8'
  meta_data.px = meta_data[which(meta_data$Person == px & meta_data$Purpose == 'Replications'),]
  meta_data.px.uid = paste(meta_data.px$Labcode, meta_data.px$SInstrument, meta_data.px$Loading, meta_data.px$Batch, sep = '_')
  uids = sort(unique(meta_data.px.uid))
  
  
  px.fot = get_fot_df_by_expnames(fot_data, meta_data.px$Exp_code)
  px.fot.row_id = apply(px.fot[,-1], 1, function(x){
    length(which(x>0))
  })
  px.fot.subset = px.fot[px.fot.row_id==72,]
  
  
  uid.mean.df = NULL
  for (i in 1:24) {
    uid = uids[i]
    uid.expcode = meta_data.px$Exp_code[meta_data.px.uid==uid]
    uid.fot_df = get_fot_df_by_expnames(px.fot.subset, uid.expcode)
    uid.mean = apply(uid.fot_df[,-1], 1, function(x){mean(x)})
    uid.mean.df = cbind(uid.mean.df, uid.mean)
  }
  colnames(uid.mean.df) = uids
  uid.mean.df = data.frame(
    GeneSymbol = px.fot.subset$GeneSymbol,
    uid.mean.df
  )
  
  px.df = data.frame(
    gene = uid.mean.df$GeneSymbol,
    cv = apply(uid.mean.df[,-1], 1, function(x){sd(x)/mean(x)}),
    mean = apply(uid.mean.df[,-1], 1, function(x){mean(x)})
  )
  
  px.df.subset = px.df[px.df$cv<0.25,]
  px.df.subset1 = px.df.subset[px.df.subset$mean>1000 & px.df.subset$mean<10000,]
  
  
  gene_subset = c(
    'HUWE1',
    'HSPH1', 'ACTN4', 'PARP1',
    'MYH9', 
    'HSPA5',
    'LDHB',
    'PKM'
  )
  gene_subset.df1 = uid.mean.df[match(gene_subset, uid.mean.df$GeneSymbol),]
  colnames(gene_subset.df1) = c('GeneSymbol', paste('S', seq(1, 24), sep = ''))
  
  library(reshape2)
  df0 = gene_subset.df1
  df = melt(df0)
  colnames(df) = c('GeneSymbol', 'Site', 'Value')
  df$Site = factor(df$Site, levels =  paste('S', seq(1, 24), sep = ''))
  df$GeneSymbol = factor(df$GeneSymbol, levels = gene_subset)
  # df$Value = log10(df$Value)
  p <- ggplot(df, aes(x = Site, y = Value))
  p <- p + geom_point(aes(shape=GeneSymbol), color = p_colors[4], fill = p_colors[3], size = 3)
  p <- p + scale_shape_manual(values=c(15, 0, 16, 1, 17, 2, 18, 5)) + scale_y_log10()
  p <- p + ylab('FOT') + theme_classic2() +  ggtitle('P8')
  print(p)
  df8 = df
}




if(T){
  df2 = rbind.data.frame(
    df5, df6, df7, df8
  )
  df2$Sample = rep(unique(meta_data$UPerson), each = 192)
  df = df2
  p <- ggplot(df, aes(x = Site, y = Value))
  p <- p + geom_point(aes(shape=GeneSymbol, color = Sample),  size = 3)
  p <- p + scale_shape_manual(values=c(15, 0, 16, 1, 17, 2, 18, 5)) + scale_y_log10()
  p <- p + ylab('FOT') + theme_classic2() #+  ggtitle('P8')
  p <- p + facet_grid(GeneSymbol ~ .)
  print(p)

}



