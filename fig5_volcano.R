meta_data.NPB.Lumos = meta_data[which(
  meta_data$Labcode == 'NPB' 
  & meta_data$SInstrument == 'Lumos'
  & meta_data$Purpose == 'Replications'
  & meta_data$Batch == 1
),]


meta_data.BGI.Lumos = meta_data[which(
  meta_data$Labcode == 'BGI' 
  & meta_data$SInstrument == 'Lumos'
  & meta_data$Purpose == 'Replications'
  & meta_data$Batch == 1
),]

if(T){
  # single sample
  # D5
  px_list = c('D5', 'D6', 'F7', 'M8')
  px = px_list[1]
  # npb
  npb.px.index = which(meta_data.NPB.Lumos$UPerson == px)
  npb.px.fot = get_fot_df_by_expnames(fot_data, meta_data.NPB.Lumos$Exp_code[npb.px.index])
  npb.px.fot.row_id = apply(npb.px.fot[,-1], 1, function(x){length(which(x>0))})
  npb.px.fot.subset = npb.px.fot[npb.px.fot.row_id>=2,]
  
  # bgi
  bgi.px.index = which(meta_data.BGI.Lumos$UPerson == px)
  bgi.px.fot = get_fot_df_by_expnames(fot_data, meta_data.BGI.Lumos$Exp_code[bgi.px.index])
  bgi.px.fot.row_id = apply(bgi.px.fot[,-1], 1, function(x){length(which(x>0))})
  bgi.px.fot.subset = bgi.px.fot[bgi.px.fot.row_id>=2,]
  
  px_merge_df = merge(npb.px.fot.subset, bgi.px.fot.subset, by = 'GeneSymbol', all = T)
  px_merge_df.mat = as.matrix(px_merge_df[,-1])
  na_idnex = which(is.na(px_merge_df.mat))
  px_merge_df.mat[na_idnex] = 0
  # min_v = min(px_merge_df.mat[px_merge_df.mat>0])
  # zero_index = which(px_merge_df.mat==0)
  # px_merge_df.mat[zero_index] = min_v/10
  px_merge_df1 = data.frame(
    GeneSymbol = px_merge_df$GeneSymbol,
    px_merge_df.mat
  )
  px_group = rep(c('NPB', 'BGI'), each = 3)
  px_group = factor(px_group, levels = c('NPB', 'BGI'))
  px_merge_df1.row = nrow(px_merge_df1)
  pva_lst = c()
  fc_lst = c()
  for (i in 1:px_merge_df1.row) {
    x = unlist(px_merge_df1[i,-1])
    x = x + 1
    x_test = t.test(log2(x)~px_group)
    x.pval = x_test$p.value
    x.fc = mean(x[1:3])/mean(x[4:6])
    # if(x.fc == 0){break}
    pva_lst = c(pva_lst, x.pval)
    fc_lst = c(fc_lst, x.fc)
  }
  deps_df = data.frame(
    GeneSymbol = px_merge_df1$GeneSymbol,
    Pvalue = pva_lst,
    FC = fc_lst,
    FDR = p.adjust(pva_lst, method = 'BH')
  )
  index = which(deps_df$FC>0 & deps_df$FC<Inf)
  fc.min = min(deps_df$FC[index])
  fc.max = max(deps_df$FC[index])
  deps_df$FC[deps_df$FC==0] = fc.min
  deps_df$FC[is.infinite(deps_df$FC)] = fc.max
  h_index = which(deps_df$FC>2 & deps_df$Pvalue<0.05)
  l_index = which(deps_df$FC<1/2 & deps_df$Pvalue<0.05)
  deps_df$Flag = 'Constant'
  deps_df$Flag[h_index] = 'UP'
  deps_df$Flag[l_index] = 'Down'
  deps_df$Flag = factor(deps_df$Flag, levels = c('UP', 'Constant', 'Down'))
  quantile(log2(deps_df$FC))
  
  if(T){
    # volcano
    df = deps_df
    p <- ggplot(data = df, aes(x = log2(FC), y = -log10(Pvalue)))
    p <- p + geom_point(aes(color = Flag), size = 2)
    p <- p + scale_color_manual(values = c("firebrick", "grey", "steelblue"))
    p <- p + scale_x_continuous(
      limits = c(-11, 11), breaks = seq(-10, 10, 5)
    )
    p <- p + geom_hline(yintercept = -log10(0.05), color = 'black', linetype = "dashed")
    p <- p + geom_vline(xintercept = c(-1, 1), color = 'black', linetype = "dashed")
    p <- p + labs(title = px) + theme_classic2()
    print(p)
    p5_1 <- p
  }
  
  if(T){
    # histgram
    df = deps_df[deps_df$Flag!='Constant',]
    df$Flag = factor(df$Flag, levels = c('UP', 'Down'))
    df = as.data.frame(table(df$Flag))
    colnames(df) = c('Flag', 'Count')
    
    p <- ggplot(data = df, aes(x = Flag, y = Count, fill = Flag))
    p <- p + geom_bar(stat="identity", width = 0.5,  position=position_dodge(0.7))
    p <- p + scale_fill_manual(values = c("firebrick", "steelblue"))
    p <- p + geom_text(
      aes(label=Count), vjust=1.6, 
      color="yellow", size=3.5
    )
    p <- p + scale_y_continuous(
      limits = c(0, 300), breaks = seq(0, 300, 100)
    )
    p <- p + labs(title = px) + theme_classic2()
    print(p)
    p5_2 <- p
  }
  
  if(T){
    # ggarrange(
    #   p5_1, p5_2, widths = c(10,5)
    # )
    
    p5_1 + annotation_custom(
      ggplotGrob(p5_2),
      xmin = 3, ymin = 5,
      xmax = 11, ymax = 8
    ) 
    # 10 * 8
  }
  
}


if(T){
  # single sample
  # D5
  px_list = c('D5', 'D6', 'F7', 'M8')
  px = px_list[2]
  # npb
  npb.px.index = which(meta_data.NPB.Lumos$UPerson == px)
  npb.px.fot = get_fot_df_by_expnames(fot_data, meta_data.NPB.Lumos$Exp_code[npb.px.index])
  npb.px.fot.row_id = apply(npb.px.fot[,-1], 1, function(x){length(which(x>0))})
  npb.px.fot.subset = npb.px.fot[npb.px.fot.row_id>=2,]
  
  # bgi
  bgi.px.index = which(meta_data.BGI.Lumos$UPerson == px)
  bgi.px.fot = get_fot_df_by_expnames(fot_data, meta_data.BGI.Lumos$Exp_code[bgi.px.index])
  bgi.px.fot.row_id = apply(bgi.px.fot[,-1], 1, function(x){length(which(x>0))})
  bgi.px.fot.subset = bgi.px.fot[bgi.px.fot.row_id>=2,]
  
  px_merge_df = merge(npb.px.fot.subset, bgi.px.fot.subset, by = 'GeneSymbol', all = T)
  px_merge_df.mat = as.matrix(px_merge_df[,-1])
  na_idnex = which(is.na(px_merge_df.mat))
  px_merge_df.mat[na_idnex] = 0
  # min_v = min(px_merge_df.mat[px_merge_df.mat>0])
  # zero_index = which(px_merge_df.mat==0)
  # px_merge_df.mat[zero_index] = min_v/10
  px_merge_df1 = data.frame(
    GeneSymbol = px_merge_df$GeneSymbol,
    px_merge_df.mat
  )
  px_group = rep(c('NPB', 'BGI'), each = 3)
  px_group = factor(px_group, levels = c('NPB', 'BGI'))
  px_merge_df1.row = nrow(px_merge_df1)
  pva_lst = c()
  fc_lst = c()
  for (i in 1:px_merge_df1.row) {
    x = unlist(px_merge_df1[i,-1])
    x = x + 1
    x_test = t.test(log2(x)~px_group)
    x.pval = x_test$p.value
    x.fc = mean(x[1:3])/mean(x[4:6])
    # if(x.fc == 0){break}
    pva_lst = c(pva_lst, x.pval)
    fc_lst = c(fc_lst, x.fc)
  }
  deps_df = data.frame(
    GeneSymbol = px_merge_df1$GeneSymbol,
    Pvalue = pva_lst,
    FC = fc_lst,
    FDR = p.adjust(pva_lst, method = 'BH')
  )
  index = which(deps_df$FC>0 & deps_df$FC<Inf)
  fc.min = min(deps_df$FC[index])
  fc.max = max(deps_df$FC[index])
  deps_df$FC[deps_df$FC==0] = fc.min
  deps_df$FC[is.infinite(deps_df$FC)] = fc.max
  h_index = which(deps_df$FC>2 & deps_df$Pvalue<0.05)
  l_index = which(deps_df$FC<1/2 & deps_df$Pvalue<0.05)
  deps_df$Flag = 'Constant'
  deps_df$Flag[h_index] = 'UP'
  deps_df$Flag[l_index] = 'Down'
  deps_df$Flag = factor(deps_df$Flag, levels = c('UP', 'Constant', 'Down'))
  quantile(log2(deps_df$FC))
  
  if(T){
    # volcano
    df = deps_df
    p <- ggplot(data = df, aes(x = log2(FC), y = -log10(Pvalue)))
    p <- p + geom_point(aes(color = Flag), size = 2)
    p <- p + scale_color_manual(values = c("firebrick", "grey", "steelblue"))
    p <- p + scale_x_continuous(
      limits = c(-11, 11), breaks = seq(-10, 10, 5)
    ) + scale_y_continuous(
      limits = c(0, 8), breaks = seq(0, 8, 2)
    )
    p <- p + geom_hline(yintercept = -log10(0.05), color = 'black', linetype = "dashed")
    p <- p + geom_vline(xintercept = c(-1, 1), color = 'black', linetype = "dashed")
    p <- p + labs(title = px) + theme_classic2()
    print(p)
    p5_1 <- p
  }
  
  if(T){
    # histgram
    df = deps_df[deps_df$Flag!='Constant',]
    df$Flag = factor(df$Flag, levels = c('UP', 'Down'))
    df = as.data.frame(table(df$Flag))
    colnames(df) = c('Flag', 'Count')
    
    p <- ggplot(data = df, aes(x = Flag, y = Count, fill = Flag))
    p <- p + geom_bar(stat="identity", width = 0.5,  position=position_dodge(0.7))
    p <- p + scale_fill_manual(values = c("firebrick", "steelblue"))
    p <- p + geom_text(
      aes(label=Count), vjust=1.6, 
      color="yellow", size=3.5
    )
    p <- p + scale_y_continuous(
      limits = c(0, 400), breaks = seq(0, 400, 100)
    )
    p <- p + labs(title = px) + theme_classic2()
    print(p)
    p5_2 <- p
  }
  
  if(T){
    # ggarrange(
    #   p5_1, p5_2, widths = c(10,5)
    # )
    
    p5_1 + annotation_custom(
      ggplotGrob(p5_2),
      xmin = 3, ymin = 5,
      xmax = 11, ymax = 8
    ) 
    # 10 * 8
  }
  
}



if(T){
  # single sample
  # D5
  px_list = c('D5', 'D6', 'F7', 'M8')
  px = px_list[3]
  # npb
  npb.px.index = which(meta_data.NPB.Lumos$UPerson == px)
  npb.px.fot = get_fot_df_by_expnames(fot_data, meta_data.NPB.Lumos$Exp_code[npb.px.index])
  npb.px.fot.row_id = apply(npb.px.fot[,-1], 1, function(x){length(which(x>0))})
  npb.px.fot.subset = npb.px.fot[npb.px.fot.row_id>=2,]
  
  # bgi
  bgi.px.index = which(meta_data.BGI.Lumos$UPerson == px)
  bgi.px.fot = get_fot_df_by_expnames(fot_data, meta_data.BGI.Lumos$Exp_code[bgi.px.index])
  bgi.px.fot.row_id = apply(bgi.px.fot[,-1], 1, function(x){length(which(x>0))})
  bgi.px.fot.subset = bgi.px.fot[bgi.px.fot.row_id>=2,]
  
  px_merge_df = merge(npb.px.fot.subset, bgi.px.fot.subset, by = 'GeneSymbol', all = T)
  px_merge_df.mat = as.matrix(px_merge_df[,-1])
  na_idnex = which(is.na(px_merge_df.mat))
  px_merge_df.mat[na_idnex] = 0
  # min_v = min(px_merge_df.mat[px_merge_df.mat>0])
  # zero_index = which(px_merge_df.mat==0)
  # px_merge_df.mat[zero_index] = min_v/10
  px_merge_df1 = data.frame(
    GeneSymbol = px_merge_df$GeneSymbol,
    px_merge_df.mat
  )
  px_group = rep(c('NPB', 'BGI'), each = 3)
  px_group = factor(px_group, levels = c('NPB', 'BGI'))
  px_merge_df1.row = nrow(px_merge_df1)
  pva_lst = c()
  fc_lst = c()
  for (i in 1:px_merge_df1.row) {
    x = unlist(px_merge_df1[i,-1])
    x = x + 1
    x_test = t.test(log2(x)~px_group)
    x.pval = x_test$p.value
    x.fc = mean(x[1:3])/mean(x[4:6])
    # if(x.fc == 0){break}
    pva_lst = c(pva_lst, x.pval)
    fc_lst = c(fc_lst, x.fc)
  }
  deps_df = data.frame(
    GeneSymbol = px_merge_df1$GeneSymbol,
    Pvalue = pva_lst,
    FC = fc_lst,
    FDR = p.adjust(pva_lst, method = 'BH')
  )
  index = which(deps_df$FC>0 & deps_df$FC<Inf)
  fc.min = min(deps_df$FC[index])
  fc.max = max(deps_df$FC[index])
  deps_df$FC[deps_df$FC==0] = fc.min
  deps_df$FC[is.infinite(deps_df$FC)] = fc.max
  h_index = which(deps_df$FC>2 & deps_df$Pvalue<0.05)
  l_index = which(deps_df$FC<1/2 & deps_df$Pvalue<0.05)
  deps_df$Flag = 'Constant'
  deps_df$Flag[h_index] = 'UP'
  deps_df$Flag[l_index] = 'Down'
  deps_df$Flag = factor(deps_df$Flag, levels = c('UP', 'Constant', 'Down'))
  quantile(log2(deps_df$FC))
  
  if(T){
    # volcano
    df = deps_df
    p <- ggplot(data = df, aes(x = log2(FC), y = -log10(Pvalue)))
    p <- p + geom_point(aes(color = Flag), size = 2)
    p <- p + scale_color_manual(values = c("firebrick", "grey", "steelblue"))
    p <- p + scale_x_continuous(
      limits = c(-11, 11), breaks = seq(-10, 10, 5)
    ) + scale_y_continuous(
      limits = c(0, 8), breaks = seq(0, 8, 2)
    )
    p <- p + geom_hline(yintercept = -log10(0.05), color = 'black', linetype = "dashed")
    p <- p + geom_vline(xintercept = c(-1, 1), color = 'black', linetype = "dashed")
    p <- p + labs(title = px) + theme_classic2()
    print(p)
    p5_1 <- p
  }
  
  if(T){
    # histgram
    df = deps_df[deps_df$Flag!='Constant',]
    df$Flag = factor(df$Flag, levels = c('UP', 'Down'))
    df = as.data.frame(table(df$Flag))
    colnames(df) = c('Flag', 'Count')
    
    p <- ggplot(data = df, aes(x = Flag, y = Count, fill = Flag))
    p <- p + geom_bar(stat="identity", width = 0.5,  position=position_dodge(0.7))
    p <- p + scale_fill_manual(values = c("firebrick", "steelblue"))
    p <- p + geom_text(
      aes(label=Count), vjust=1.6, 
      color="yellow", size=3.5
    )
    p <- p + scale_y_continuous(
      limits = c(0, 400), breaks = seq(0, 400, 100)
    )
    p <- p + labs(title = px) + theme_classic2()
    print(p)
    p5_2 <- p
  }
  
  if(T){
    # ggarrange(
    #   p5_1, p5_2, widths = c(10,5)
    # )
    
    p5_1 + annotation_custom(
      ggplotGrob(p5_2),
      xmin = 3, ymin = 5,
      xmax = 11, ymax = 8
    ) 
    # 10 * 8
  }
  
}



if(T){
  # single sample
  # D5
  px_list = c('D5', 'D6', 'F7', 'M8')
  px = px_list[4]
  # npb
  npb.px.index = which(meta_data.NPB.Lumos$UPerson == px)
  npb.px.fot = get_fot_df_by_expnames(fot_data, meta_data.NPB.Lumos$Exp_code[npb.px.index])
  npb.px.fot.row_id = apply(npb.px.fot[,-1], 1, function(x){length(which(x>0))})
  npb.px.fot.subset = npb.px.fot[npb.px.fot.row_id>=2,]
  
  # bgi
  bgi.px.index = which(meta_data.BGI.Lumos$UPerson == px)
  bgi.px.fot = get_fot_df_by_expnames(fot_data, meta_data.BGI.Lumos$Exp_code[bgi.px.index])
  bgi.px.fot.row_id = apply(bgi.px.fot[,-1], 1, function(x){length(which(x>0))})
  bgi.px.fot.subset = bgi.px.fot[bgi.px.fot.row_id>=2,]
  
  px_merge_df = merge(npb.px.fot.subset, bgi.px.fot.subset, by = 'GeneSymbol', all = T)
  px_merge_df.mat = as.matrix(px_merge_df[,-1])
  na_idnex = which(is.na(px_merge_df.mat))
  px_merge_df.mat[na_idnex] = 0
  # min_v = min(px_merge_df.mat[px_merge_df.mat>0])
  # zero_index = which(px_merge_df.mat==0)
  # px_merge_df.mat[zero_index] = min_v/10
  px_merge_df1 = data.frame(
    GeneSymbol = px_merge_df$GeneSymbol,
    px_merge_df.mat
  )
  px_group = rep(c('NPB', 'BGI'), each = 3)
  px_group = factor(px_group, levels = c('NPB', 'BGI'))
  px_merge_df1.row = nrow(px_merge_df1)
  pva_lst = c()
  fc_lst = c()
  for (i in 1:px_merge_df1.row) {
    x = unlist(px_merge_df1[i,-1])
    x = x + 1
    x_test = t.test(log2(x)~px_group)
    x.pval = x_test$p.value
    x.fc = mean(x[1:3])/mean(x[4:6])
    # if(x.fc == 0){break}
    pva_lst = c(pva_lst, x.pval)
    fc_lst = c(fc_lst, x.fc)
  }
  deps_df = data.frame(
    GeneSymbol = px_merge_df1$GeneSymbol,
    Pvalue = pva_lst,
    FC = fc_lst,
    FDR = p.adjust(pva_lst, method = 'BH')
  )
  index = which(deps_df$FC>0 & deps_df$FC<Inf)
  fc.min = min(deps_df$FC[index])
  fc.max = max(deps_df$FC[index])
  deps_df$FC[deps_df$FC==0] = fc.min
  deps_df$FC[is.infinite(deps_df$FC)] = fc.max
  h_index = which(deps_df$FC>2 & deps_df$Pvalue<0.05)
  l_index = which(deps_df$FC<1/2 & deps_df$Pvalue<0.05)
  deps_df$Flag = 'Constant'
  deps_df$Flag[h_index] = 'UP'
  deps_df$Flag[l_index] = 'Down'
  deps_df$Flag = factor(deps_df$Flag, levels = c('UP', 'Constant', 'Down'))
  quantile(log2(deps_df$FC))
  
  if(T){
    # volcano
    df = deps_df
    p <- ggplot(data = df, aes(x = log2(FC), y = -log10(Pvalue)))
    p <- p + geom_point(aes(color = Flag), size = 2)
    p <- p + scale_color_manual(values = c("firebrick", "grey", "steelblue"))
    p <- p + scale_x_continuous(
      limits = c(-11, 11), breaks = seq(-10, 10, 5)
    ) + scale_y_continuous(
      limits = c(0, 8), breaks = seq(0, 8, 2)
    )
    p <- p + geom_hline(yintercept = -log10(0.05), color = 'black', linetype = "dashed")
    p <- p + geom_vline(xintercept = c(-1, 1), color = 'black', linetype = "dashed")
    p <- p + labs(title = px) + theme_classic2()
    print(p)
    p5_1 <- p
  }
  
  if(T){
    # histgram
    df = deps_df[deps_df$Flag!='Constant',]
    df$Flag = factor(df$Flag, levels = c('UP', 'Down'))
    df = as.data.frame(table(df$Flag))
    colnames(df) = c('Flag', 'Count')
    
    p <- ggplot(data = df, aes(x = Flag, y = Count, fill = Flag))
    p <- p + geom_bar(stat="identity", width = 0.5,  position=position_dodge(0.7))
    p <- p + scale_fill_manual(values = c("firebrick", "steelblue"))
    p <- p + geom_text(
      aes(label=Count), vjust=1.6, 
      color="yellow", size=3.5
    )
    p <- p + scale_y_continuous(
      limits = c(0, 400), breaks = seq(0, 400, 100)
    )
    p <- p + labs(title = px) + theme_classic2()
    print(p)
    p5_2 <- p
  }
  
  if(T){
    # ggarrange(
    #   p5_1, p5_2, widths = c(10,5)
    # )
    
    p5_1 + annotation_custom(
      ggplotGrob(p5_2),
      xmin = 3, ymin = 5,
      xmax = 11, ymax = 8
    ) 
    # 10 * 8
  }
  
}




if(T){
  # inter
  # D5
  px_list = c('D5', 'D6', 'F7', 'M8')
  px1 = px_list[1]
  px2 = px_list[2]
  px = paste(px1, 'vs', px2)
  # npb
  npb.px.index = which(meta_data.NPB.Lumos$UPerson == px1)
  npb.px.fot = get_fot_df_by_expnames(fot_data, meta_data.NPB.Lumos$Exp_code[npb.px.index])
  npb.px.fot.row_id = apply(npb.px.fot[,-1], 1, function(x){length(which(x>0))})
  npb.px.fot.subset1 = npb.px.fot[npb.px.fot.row_id>=2,]
  
  
  # npb
  npb.px.index = which(meta_data.NPB.Lumos$UPerson == px2)
  npb.px.fot = get_fot_df_by_expnames(fot_data, meta_data.NPB.Lumos$Exp_code[npb.px.index])
  npb.px.fot.row_id = apply(npb.px.fot[,-1], 1, function(x){length(which(x>0))})
  npb.px.fot.subset2 = npb.px.fot[npb.px.fot.row_id>=2,]
  
  
  px_merge_df = merge(npb.px.fot.subset1, npb.px.fot.subset2, by = 'GeneSymbol', all = T)
  px_merge_df.mat = as.matrix(px_merge_df[,-1])
  na_idnex = which(is.na(px_merge_df.mat))
  px_merge_df.mat[na_idnex] = 0
  # min_v = min(px_merge_df.mat[px_merge_df.mat>0])
  # zero_index = which(px_merge_df.mat==0)
  # px_merge_df.mat[zero_index] = min_v/10
  px_merge_df1 = data.frame(
    GeneSymbol = px_merge_df$GeneSymbol,
    px_merge_df.mat
  )
  px_group = rep(c('PX1', 'PX2'), each = 3)
  px_group = factor(px_group, levels = c('PX1', 'PX2'))
  px_merge_df1.row = nrow(px_merge_df1)
  pva_lst = c()
  fc_lst = c()
  for (i in 1:px_merge_df1.row) {
    x = unlist(px_merge_df1[i,-1])
    x = x+1
    x_test = t.test(log2(x)~px_group)
    x.pval = x_test$p.value
    x.fc = mean(x[1:3])/mean(x[4:6])
    # if(x.fc == 0){break}
    pva_lst = c(pva_lst, x.pval)
    fc_lst = c(fc_lst, x.fc)
  }
  deps_df = data.frame(
    GeneSymbol = px_merge_df1$GeneSymbol,
    Pvalue = pva_lst,
    FC = fc_lst,
    FDR = p.adjust(pva_lst, method = 'BH')
  )
  index = which(deps_df$FC>0 & deps_df$FC<Inf)
  fc.min = min(deps_df$FC[index])
  fc.max = max(deps_df$FC[index])
  deps_df$FC[deps_df$FC==0] = fc.min
  deps_df$FC[is.infinite(deps_df$FC)] = fc.max
  h_index = which(deps_df$FC>2 & deps_df$Pvalue<0.05)
  l_index = which(deps_df$FC<1/2 & deps_df$Pvalue<0.05)
  deps_df$Flag = 'Constant'
  deps_df$Flag[h_index] = 'UP'
  deps_df$Flag[l_index] = 'Down'
  deps_df$Flag = factor(deps_df$Flag, levels = c('UP', 'Constant', 'Down'))
  quantile(log2(deps_df$FC))
  
  if(T){
    # volcano
    df = deps_df
    p <- ggplot(data = df, aes(x = log2(FC), y = -log10(Pvalue)))
    p <- p + geom_point(aes(color = Flag), size = 2)
    p <- p + scale_color_manual(values = c("firebrick", "grey", "steelblue"))
    p <- p + scale_x_continuous(
      limits = c(-11, 11), breaks = seq(-10, 10, 5)
    ) + scale_y_continuous(
      limits = c(0, 8), breaks = seq(0, 8, 2)
    )
    p <- p + geom_hline(yintercept = -log10(0.05), color = 'black', linetype = "dashed")
    p <- p + geom_vline(xintercept = c(-1, 1), color = 'black', linetype = "dashed")
    p <- p + labs(title = px) + theme_classic2()
    print(p)
    p5_1 <- p
  }
  
  if(T){
    # histgram
    df = deps_df[deps_df$Flag!='Constant',]
    df$Flag = factor(df$Flag, levels = c('UP', 'Down'))
    df = as.data.frame(table(df$Flag))
    colnames(df) = c('Flag', 'Count')
    
    p <- ggplot(data = df, aes(x = Flag, y = Count, fill = Flag))
    p <- p + geom_bar(stat="identity", width = 0.5,  position=position_dodge(0.7))
    p <- p + scale_fill_manual(values = c("firebrick", "steelblue"))
    p <- p + geom_text(
      aes(label=Count), vjust=1.6, 
      color="yellow", size=3.5
    )
    p <- p + scale_y_continuous(
      limits = c(0, 100), breaks = seq(0, 100, 50)
    )
    p <- p + labs(title = px) + theme_classic2()
    print(p)
    p5_2 <- p
  }
  
  if(T){
    # ggarrange(
    #   p5_1, p5_2, widths = c(10,5)
    # )
    
    p5_1 + annotation_custom(
      ggplotGrob(p5_2),
      xmin = 3, ymin = 5,
      xmax = 11, ymax = 8
    ) 
    # 10 * 8
  }
  
}










if(T){
  # inter
  # D5
  px_list = c('D5', 'D6', 'F7', 'M8')
  px1 = px_list[3]
  px2 = px_list[4]
  px = paste(px1, 'vs', px2)
  # npb
  npb.px.index = which(meta_data.NPB.Lumos$UPerson == px1)
  npb.px.fot = get_fot_df_by_expnames(fot_data, meta_data.NPB.Lumos$Exp_code[npb.px.index])
  npb.px.fot.row_id = apply(npb.px.fot[,-1], 1, function(x){length(which(x>0))})
  npb.px.fot.subset1 = npb.px.fot[npb.px.fot.row_id>=2,]
  
  
  # npb
  npb.px.index = which(meta_data.NPB.Lumos$UPerson == px2)
  npb.px.fot = get_fot_df_by_expnames(fot_data, meta_data.NPB.Lumos$Exp_code[npb.px.index])
  npb.px.fot.row_id = apply(npb.px.fot[,-1], 1, function(x){length(which(x>0))})
  npb.px.fot.subset2 = npb.px.fot[npb.px.fot.row_id>=2,]
  
  
  px_merge_df = merge(npb.px.fot.subset1, npb.px.fot.subset2, by = 'GeneSymbol', all = T)
  px_merge_df.mat = as.matrix(px_merge_df[,-1])
  na_idnex = which(is.na(px_merge_df.mat))
  px_merge_df.mat[na_idnex] = 0
  # min_v = min(px_merge_df.mat[px_merge_df.mat>0])
  # zero_index = which(px_merge_df.mat==0)
  # px_merge_df.mat[zero_index] = min_v/10
  px_merge_df1 = data.frame(
    GeneSymbol = px_merge_df$GeneSymbol,
    px_merge_df.mat
  )
  px_group = rep(c('PX1', 'PX2'), each = 3)
  px_group = factor(px_group, levels = c('PX1', 'PX2'))
  px_merge_df1.row = nrow(px_merge_df1)
  pva_lst = c()
  fc_lst = c()
  for (i in 1:px_merge_df1.row) {
    x = unlist(px_merge_df1[i,-1])
    x = x+1
    x_test = t.test(log2(x)~px_group)
    x.pval = x_test$p.value
    x.fc = mean(x[1:3])/mean(x[4:6])
    # if(x.fc == 0){break}
    pva_lst = c(pva_lst, x.pval)
    fc_lst = c(fc_lst, x.fc)
  }
  deps_df = data.frame(
    GeneSymbol = px_merge_df1$GeneSymbol,
    Pvalue = pva_lst,
    FC = fc_lst,
    FDR = p.adjust(pva_lst, method = 'BH')
  )
  index = which(deps_df$FC>0 & deps_df$FC<Inf)
  fc.min = min(deps_df$FC[index])
  fc.max = max(deps_df$FC[index])
  deps_df$FC[deps_df$FC==0] = fc.min
  deps_df$FC[is.infinite(deps_df$FC)] = fc.max
  h_index = which(deps_df$FC>2 & deps_df$Pvalue<0.05)
  l_index = which(deps_df$FC<1/2 & deps_df$Pvalue<0.05)
  deps_df$Flag = 'Constant'
  deps_df$Flag[h_index] = 'UP'
  deps_df$Flag[l_index] = 'Down'
  deps_df$Flag = factor(deps_df$Flag, levels = c('UP', 'Constant', 'Down'))
  quantile(log2(deps_df$FC))
  
  if(T){
    # volcano
    df = deps_df
    p <- ggplot(data = df, aes(x = log2(FC), y = -log10(Pvalue)))
    p <- p + geom_point(aes(color = Flag), size = 2)
    p <- p + scale_color_manual(values = c("firebrick", "grey", "steelblue"))
    p <- p + scale_x_continuous(
      limits = c(-11, 11), breaks = seq(-10, 10, 5)
    ) + scale_y_continuous(
      limits = c(0, 8), breaks = seq(0, 8, 2)
    )
    p <- p + geom_hline(yintercept = -log10(0.05), color = 'black', linetype = "dashed")
    p <- p + geom_vline(xintercept = c(-1, 1), color = 'black', linetype = "dashed")
    p <- p + labs(title = px) + theme_classic2()
    print(p)
    p5_1 <- p
  }
  
  if(T){
    # histgram
    df = deps_df[deps_df$Flag!='Constant',]
    df$Flag = factor(df$Flag, levels = c('UP', 'Down'))
    df = as.data.frame(table(df$Flag))
    colnames(df) = c('Flag', 'Count')
    
    p <- ggplot(data = df, aes(x = Flag, y = Count, fill = Flag))
    p <- p + geom_bar(stat="identity", width = 0.5,  position=position_dodge(0.7))
    p <- p + scale_fill_manual(values = c("firebrick", "steelblue"))
    p <- p + geom_text(
      aes(label=Count), vjust=1.6, 
      color="yellow", size=3.5
    )
    p <- p + scale_y_continuous(
      limits = c(0, 120), breaks = seq(0, 100, 50)
    )
    p <- p + labs(title = px) + theme_classic2()
    print(p)
    p5_2 <- p
  }
  
  if(T){
    # ggarrange(
    #   p5_1, p5_2, widths = c(10,5)
    # )
    
    p5_1 + annotation_custom(
      ggplotGrob(p5_2),
      xmin = 3, ymin = 5,
      xmax = 11, ymax = 8
    ) 
    # 10 * 8
  }
  
}







if(T){
  # inter
  # D5
  px_list = c('D5', 'D6', 'F7', 'M8')
  px1 = px_list[2]
  px2 = px_list[3]
  px = paste(px1, 'vs', px2)
  
  # npb
  npb.px.index = which(meta_data.BGI.Lumos$UPerson == px1)
  npb.px.fot = get_fot_df_by_expnames(fot_data, meta_data.BGI.Lumos$Exp_code[npb.px.index])
  npb.px.fot.row_id = apply(npb.px.fot[,-1], 1, function(x){length(which(x>0))})
  npb.px.fot.subset1 = npb.px.fot[npb.px.fot.row_id>=2,]
  
  
  # npb
  npb.px.index = which(meta_data.BGI.Lumos$UPerson == px2)
  npb.px.fot = get_fot_df_by_expnames(fot_data, meta_data.BGI.Lumos$Exp_code[npb.px.index])
  npb.px.fot.row_id = apply(npb.px.fot[,-1], 1, function(x){length(which(x>0))})
  npb.px.fot.subset2 = npb.px.fot[npb.px.fot.row_id>=2,]
  
  
  px_merge_df = merge(npb.px.fot.subset1, npb.px.fot.subset2, by = 'GeneSymbol', all = T)
  px_merge_df.mat = as.matrix(px_merge_df[,-1])
  na_idnex = which(is.na(px_merge_df.mat))
  px_merge_df.mat[na_idnex] = 0
  # min_v = min(px_merge_df.mat[px_merge_df.mat>0])
  # zero_index = which(px_merge_df.mat==0)
  # px_merge_df.mat[zero_index] = min_v/10
  px_merge_df1 = data.frame(
    GeneSymbol = px_merge_df$GeneSymbol,
    px_merge_df.mat
  )
  px_group = rep(c('PX1', 'PX2'), each = 3)
  px_group = factor(px_group, levels = c('PX1', 'PX2'))
  px_merge_df1.row = nrow(px_merge_df1)
  pva_lst = c()
  fc_lst = c()
  for (i in 1:px_merge_df1.row) {
    x = unlist(px_merge_df1[i,-1])
    x = x+1
    x_test = t.test(log2(x)~px_group)
    x.pval = x_test$p.value
    x.fc = mean(x[1:3])/mean(x[4:6])
    # if(x.fc == 0){break}
    pva_lst = c(pva_lst, x.pval)
    fc_lst = c(fc_lst, x.fc)
  }
  deps_df = data.frame(
    GeneSymbol = px_merge_df1$GeneSymbol,
    Pvalue = pva_lst,
    FC = fc_lst,
    FDR = p.adjust(pva_lst, method = 'BH')
  )
  index = which(deps_df$FC>0 & deps_df$FC<Inf)
  fc.min = min(deps_df$FC[index])
  fc.max = max(deps_df$FC[index])
  deps_df$FC[deps_df$FC==0] = fc.min
  deps_df$FC[is.infinite(deps_df$FC)] = fc.max
  h_index = which(deps_df$FC>2 & deps_df$Pvalue<0.05)
  l_index = which(deps_df$FC<1/2 & deps_df$Pvalue<0.05)
  deps_df$Flag = 'Constant'
  deps_df$Flag[h_index] = 'UP'
  deps_df$Flag[l_index] = 'Down'
  deps_df$Flag = factor(deps_df$Flag, levels = c('UP', 'Constant', 'Down'))
  quantile(log2(deps_df$FC))
  
  if(T){
    # volcano
    df = deps_df
    p <- ggplot(data = df, aes(x = log2(FC), y = -log10(Pvalue)))
    p <- p + geom_point(aes(color = Flag), size = 2)
    p <- p + scale_color_manual(values = c("firebrick", "grey", "steelblue"))
    p <- p + scale_x_continuous(
      limits = c(-11, 11), breaks = seq(-10, 10, 5)
    ) + scale_y_continuous(
      limits = c(0, 8), breaks = seq(0, 8, 2)
    )
    p <- p + geom_hline(yintercept = -log10(0.05), color = 'black', linetype = "dashed")
    p <- p + geom_vline(xintercept = c(-1, 1), color = 'black', linetype = "dashed")
    p <- p + labs(title = px) + theme_classic2()
    print(p)
    p5_1 <- p
  }
  
  if(T){
    # histgram
    df = deps_df[deps_df$Flag!='Constant',]
    df$Flag = factor(df$Flag, levels = c('UP', 'Down'))
    df = as.data.frame(table(df$Flag))
    colnames(df) = c('Flag', 'Count')
    
    p <- ggplot(data = df, aes(x = Flag, y = Count, fill = Flag))
    p <- p + geom_bar(stat="identity", width = 0.5,  position=position_dodge(0.7))
    p <- p + scale_fill_manual(values = c("firebrick", "steelblue"))
    p <- p + geom_text(
      aes(label=Count), vjust=1.6, 
      color="yellow", size=3.5
    )
    p <- p + scale_y_continuous(
      limits = c(0, 130), breaks = seq(0, 100, 50)
    )
    p <- p + labs(title = px) + theme_classic2()
    print(p)
    p5_2 <- p
  }
  
  if(T){
    # ggarrange(
    #   p5_1, p5_2, widths = c(10,5)
    # )
    
    p5_1 + annotation_custom(
      ggplotGrob(p5_2),
      xmin = 3, ymin = 5,
      xmax = 11, ymax = 8
    ) 
    # 10 * 8
  }
  
}




