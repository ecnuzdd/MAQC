meta_data.cp = meta_data[which(meta_data$Purpose == 'Replications'),]
meta_data.cp.uids = paste(
  meta_data.cp$Labcode, 
  meta_data.cp$SInstrument, 
  meta_data.cp$Batch, 
  sep = '_'
)
uids = sort(unique(meta_data.cp.uids))

meta_data_1 = meta_data.cp[which(meta_data.cp.uids=="FDU_QE-HFX_1"),]
#######################
meta_data_4 = meta_data.cp[which(meta_data.cp.uids=="FDU_QE-HFX_3"),]
meta_data_3 = meta_data.cp[which(meta_data.cp.uids=="FDU_QE-HFX_4"),]


f1.fot = get_fot_df_by_expnames(fot_data, meta_data_1$Exp_code)
f3.fot = get_fot_df_by_expnames(fot_data, meta_data_3$Exp_code)
f4.fot = get_fot_df_by_expnames(fot_data, meta_data_4$Exp_code)

#### Identification #### 
if(T){
  
  # peptide
  df0 = data.frame(
    Sample = meta_data_1$UPerson,
    GPs = apply(f1.fot[,-1], 2, function(x){length(which(x>0))})
  )
  df0.new = data.frame(
    Sample = unique(df0$Sample),
    Mean = tapply(df0$GPs, df0$Sample, mean),
    SD = tapply(df0$GPs, df0$Sample, sd)
  )
  p0 <- ggplot(data=df0.new, aes(x=Sample, y=Mean, fill=Sample)) 
  p0 <- p0 + geom_bar(stat="identity", color="black", position=position_dodge(), width=0.8)
  p0 <- p0 + scale_fill_manual(values=p_colors) + theme_classic()
  p0 <- p0 + ggtitle('Load order: 5678-5678-5678 (random)') + ylim(0, 4000)
  p0 <- p0 + geom_errorbar(
    aes(ymin=Mean-SD, ymax=Mean+SD), position=position_dodge(0.9), width=0.3
  )
  p0
  
  
  # peptide
  df1 = data.frame(
    Sample = meta_data_3$UPerson,
    GPs = apply(f3.fot[,-1], 2, function(x){length(which(x>0))})
  )
  df1.new = data.frame(
    Sample = unique(df1$Sample),
    Mean = tapply(df1$GPs, df1$Sample, mean),
    SD = tapply(df1$GPs, df1$Sample, sd)
  )
  p1 <- ggplot(data=df1.new, aes(x=Sample, y=Mean, fill=Sample)) 
  p1 <- p1 + geom_bar(stat="identity", color="black", position=position_dodge(), width=0.8)
  p1 <- p1 + scale_fill_manual(values=p_colors) + theme_classic()
  p1 <- p1 + ggtitle('Load order: 555-666-777-888 (continuous)') + ylim(0, 4000)
  p1 <- p1 + geom_errorbar(
    aes(ymin=Mean-SD, ymax=Mean+SD), position=position_dodge(0.9), width=0.3
  )
  p1
  
  df2 = data.frame(
    Sample = meta_data_4$UPerson,
    GPs = apply(f4.fot[,-1], 2, function(x){length(which(x>0))})
  )
  df2.new = data.frame(
    Sample = unique(df2$Sample),
    Mean = tapply(df2$GPs, df2$Sample, mean),
    SD = tapply(df2$GPs, df2$Sample, sd)
  )
  p2 <- ggplot(data=df2.new, aes(x=Sample, y=Mean, fill=Sample)) 
  p2 <- p2 + geom_bar(stat="identity", color="black", position=position_dodge(), width=0.8)
  p2 <- p2 + scale_fill_manual(values=p_colors) + theme_classic()
  p2 <- p2 + ggtitle('Load order: 5678-5678-5678 (continuous)') + ylim(0, 4000)
  p2 <- p2 + geom_errorbar(
    aes(ymin=Mean-SD, ymax=Mean+SD), position=position_dodge(0.9), width=0.3
  )
  p2
  
  ggarrange(p0, p1, p2, ncol = 3)
}


if(F){
  library(VennDiagram)
  library(tidyverse)
  
  f1.fot.5 = get_fot_df_by_expnames(fot_data, meta_data_1$Exp_code[meta_data_1$UPerson=='D5'])
  f1.fot.6 = get_fot_df_by_expnames(fot_data, meta_data_1$Exp_code[meta_data_1$UPerson=='D6'])
  f1.fot.7 = get_fot_df_by_expnames(fot_data, meta_data_1$Exp_code[meta_data_1$UPerson=='F7'])
  f1.fot.8 = get_fot_df_by_expnames(fot_data, meta_data_1$Exp_code[meta_data_1$UPerson=='M8'])
  
  f3.fot.5 = get_fot_df_by_expnames(fot_data, meta_data_3$Exp_code[meta_data_3$UPerson=='D5'])
  f3.fot.6 = get_fot_df_by_expnames(fot_data, meta_data_3$Exp_code[meta_data_3$UPerson=='D6'])
  f3.fot.7 = get_fot_df_by_expnames(fot_data, meta_data_3$Exp_code[meta_data_3$UPerson=='F7'])
  f3.fot.8 = get_fot_df_by_expnames(fot_data, meta_data_3$Exp_code[meta_data_3$UPerson=='M8'])
  
  f4.fot.5 = get_fot_df_by_expnames(fot_data, meta_data_4$Exp_code[meta_data_3$UPerson=='D5'])
  f4.fot.6 = get_fot_df_by_expnames(fot_data, meta_data_4$Exp_code[meta_data_3$UPerson=='D6'])
  f4.fot.7 = get_fot_df_by_expnames(fot_data, meta_data_4$Exp_code[meta_data_3$UPerson=='F7'])
  f4.fot.8 = get_fot_df_by_expnames(fot_data, meta_data_4$Exp_code[meta_data_3$UPerson=='M8'])
  
  jco_col = pal_jco("default", alpha = 0.5)(10)[c(4,2,1)]
  temp.d5 <- venn.diagram(
    x = list(
      T1 = f1.fot.5$GeneSymbol,
      T2 = f3.fot.5$GeneSymbol,
      T3 = f4.fot.5$GeneSymbol
    ),
    filename = NULL,
    main = 'D5',
    # main.just = c(0, 5),
    main.pos = c(0.5,1.05),
    main.cex = 0.5,
    output = TRUE ,
    imagetype="png" ,
    height = 480 , 
    width = 480 , 
    resolution = 300,
    compression = "lzw",
    lwd = 1,
    col = jco_col,
    fill = jco_col,
    cex = 0.5,
    fontfamily = "sans",
    cat.cex = 0.5,
    cat.default.pos = "outer",
    # cat.pos = c(-2, 27, -135, 135),
    # cat.dist = c(0.055, 0.055, 0.085),
    cat.fontfamily = "sans"
    # cat.col = p_colors
    # rotation = 1,
    # print.mode = 'percent'
  )
  
  
  
  temp.d6 <- venn.diagram(
    x = list(
      T1 = f1.fot.6$GeneSymbol,
      T2 = f3.fot.6$GeneSymbol,
      T3 = f4.fot.6$GeneSymbol
    ),
    filename = NULL,
    main = 'D6',
    # main.just = c(0, 5),
    main.pos = c(0.5,1.05),
    main.cex = 0.5,
    output = TRUE ,
    imagetype="png" ,
    height = 480 , 
    width = 480 , 
    resolution = 300,
    compression = "lzw",
    lwd = 1,
    col = jco_col,
    fill = jco_col,
    cex = 0.5,
    fontfamily = "sans",
    cat.cex = 0.5,
    cat.default.pos = "outer",
    # cat.pos = c(-2, 27, -135, 135),
    # cat.dist = c(0.055, 0.055, 0.085),
    cat.fontfamily = "sans"
    # cat.col = p_colors
    # rotation = 1,
    # print.mode = 'percent'
  )
  
  
  
  
  
  temp.f7 <- venn.diagram(
    x = list(
      T1 = f1.fot.7$GeneSymbol,
      T2 = f3.fot.7$GeneSymbol,
      T3 = f4.fot.7$GeneSymbol
    ),
    filename = NULL,
    main = 'F7',
    # main.just = c(0, 5),
    main.pos = c(0.5,1.05),
    main.cex = 0.5,
    output = TRUE ,
    imagetype="png" ,
    height = 480 , 
    width = 480 , 
    resolution = 300,
    compression = "lzw",
    lwd = 1,
    col = jco_col,
    fill = jco_col,
    cex = 0.5,
    fontfamily = "sans",
    cat.cex = 0.5,
    cat.default.pos = "outer",
    # cat.pos = c(-2, 27, -135, 135),
    # cat.dist = c(0.055, 0.055, 0.085),
    cat.fontfamily = "sans"
    # cat.col = p_colors
    # rotation = 1,
    # print.mode = 'percent'
  )
  
  
  temp.m8 <- venn.diagram(
    x = list(
      T1 = f1.fot.8$GeneSymbol,
      T2 = f3.fot.8$GeneSymbol,
      T3 = f4.fot.8$GeneSymbol
    ),
    filename = NULL,
    main = 'M8',
    # main.just = c(0, 5),
    main.pos = c(0.5,1.05),
    main.cex = 0.5,
    output = TRUE ,
    imagetype="png" ,
    height = 480 , 
    width = 480 , 
    resolution = 300,
    compression = "lzw",
    lwd = 1,
    col = jco_col,
    fill = jco_col,
    cex = 0.5,
    fontfamily = "sans",
    cat.cex = 0.5,
    cat.default.pos = "outer",
    # cat.pos = c(-2, 27, -135, 135),
    # cat.dist = c(0.055, 0.055, 0.085),
    cat.fontfamily = "sans"
    # cat.col = p_colors
    # rotation = 1,
    # print.mode = 'percent'
  )
  
  ggarrange(
    temp.d5, temp.d6, 
    temp.f7, temp.m8, 
    ncol = 2, nrow = 2
  )
  
}

#### correlation ####
if(T){
  f1.fot.5 = get_fot_df_by_expnames(fot_data, meta_data_1$Exp_code[meta_data_1$UPerson=='D5'])
  f1.fot.6 = get_fot_df_by_expnames(fot_data, meta_data_1$Exp_code[meta_data_1$UPerson=='D6'])
  f1.fot.7 = get_fot_df_by_expnames(fot_data, meta_data_1$Exp_code[meta_data_1$UPerson=='F7'])
  f1.fot.8 = get_fot_df_by_expnames(fot_data, meta_data_1$Exp_code[meta_data_1$UPerson=='M8'])
  
  f3.fot.5 = get_fot_df_by_expnames(fot_data, meta_data_3$Exp_code[meta_data_3$UPerson=='D5'])
  f3.fot.6 = get_fot_df_by_expnames(fot_data, meta_data_3$Exp_code[meta_data_3$UPerson=='D6'])
  f3.fot.7 = get_fot_df_by_expnames(fot_data, meta_data_3$Exp_code[meta_data_3$UPerson=='F7'])
  f3.fot.8 = get_fot_df_by_expnames(fot_data, meta_data_3$Exp_code[meta_data_3$UPerson=='M8'])
  
  f4.fot.5 = get_fot_df_by_expnames(fot_data, meta_data_4$Exp_code[meta_data_3$UPerson=='D5'])
  f4.fot.6 = get_fot_df_by_expnames(fot_data, meta_data_4$Exp_code[meta_data_3$UPerson=='D6'])
  f4.fot.7 = get_fot_df_by_expnames(fot_data, meta_data_4$Exp_code[meta_data_3$UPerson=='F7'])
  f4.fot.8 = get_fot_df_by_expnames(fot_data, meta_data_4$Exp_code[meta_data_3$UPerson=='M8'])
  
  get_ggplot_cor_density <- function(dat, label){
    # dat = log10(dat)
    df <- data.frame(
      x = dat[,1],
      y = dat[,2]
    )
    df.cor = cor(df$x, df$y)
    df.cor = round(df.cor, 2)
    df.cor.text = paste(label, '(N=', nrow(df), ')', ': Cor=', df.cor, sep = '')
    
    df = log10(df)
    p <- ggplot(df, aes(x=x, y=y))
    p <- p + stat_density2d(aes(fill=..density..), geom="tile", contour = F)
    p <- p + scale_fill_gradientn(colours=magma(100)[10:100])
    p <- p + xlab(paste("log10(FOT)/", 'F3', sep = '')) + ylab(paste("log10(FOT)/", 'F4', sep = '')) 
    p <- p + scale_x_continuous(breaks=seq(-5, 5, 1), limits = c(-5, 5)) 
    p <- p + scale_y_continuous(breaks=seq(-5, 5, 1), limits = c(-5, 5))
    p <- p + theme_classic()
    p <- p + annotate("text", x = 0, y = -3.5, label = df.cor.text, color='yellow', size=4)
    p
    return(p)
  }
  
  if(T){
    df1 = f1.fot.5
    f1.df = data.frame(
      GeneSymbol = df1$GeneSymbol,
      FOT1 = apply(df1[,-1], 1, mean)
    )
    df3 = f3.fot.5
    f3.df = data.frame(
      GeneSymbol = df3$GeneSymbol,
      FOT3 = apply(df3[,-1], 1, mean)
    )
    df4 = f4.fot.5
    f4.df = data.frame(
      GeneSymbol = df4$GeneSymbol,
      FOT4 = apply(df4[,-1], 1, mean)
    )
    
    f3_f4.df = merge(f3.df, f4.df, by = 'GeneSymbol', all = T)
    f3_f4.df.mat = as.matrix(f3_f4.df[,-1])
    f3_f4.df.mat[which(is.na(f3_f4.df.mat))] = 0.00001
    
    p5 = get_ggplot_cor_density(
      f3_f4.df.mat, label = 'D5'
    )
  }
  
  if(T){
    df3 = f3.fot.6
    f3.df = data.frame(
      GeneSymbol = df3$GeneSymbol,
      FOT3 = apply(df3[,-1], 1, mean)
    )
    df4 = f4.fot.6
    f4.df = data.frame(
      GeneSymbol = df4$GeneSymbol,
      FOT4 = apply(df4[,-1], 1, mean)
    )
    f3_f4.df = merge(f3.df, f4.df, by = 'GeneSymbol', all = T)
    f3_f4.df.mat = as.matrix(f3_f4.df[,-1])
    f3_f4.df.mat[which(is.na(f3_f4.df.mat))] = 0.00001
    
    p6 = get_ggplot_cor_density(
      f3_f4.df.mat, label = 'D6'
    )
  }
  
  
  if(T){
    df3 = f3.fot.7
    f3.df = data.frame(
      GeneSymbol = df3$GeneSymbol,
      FOT3 = apply(df3[,-1], 1, mean)
    )
    df4 = f4.fot.7
    f4.df = data.frame(
      GeneSymbol = df4$GeneSymbol,
      FOT4 = apply(df4[,-1], 1, mean)
    )
    f3_f4.df = merge(f3.df, f4.df, by = 'GeneSymbol', all = T)
    f3_f4.df.mat = as.matrix(f3_f4.df[,-1])
    f3_f4.df.mat[which(is.na(f3_f4.df.mat))] = 0.00001
    
    p7 = get_ggplot_cor_density(
      f3_f4.df.mat, label = 'F7'
    )
  }
  
  if(T){
    df3 = f3.fot.8
    f3.df = data.frame(
      GeneSymbol = df3$GeneSymbol,
      FOT3 = apply(df3[,-1], 1, mean)
    )
    df4 = f4.fot.8
    f4.df = data.frame(
      GeneSymbol = df4$GeneSymbol,
      FOT4 = apply(df4[,-1], 1, mean)
    )
    f3_f4.df = merge(f3.df, f4.df, by = 'GeneSymbol', all = T)
    f3_f4.df.mat = as.matrix(f3_f4.df[,-1])
    f3_f4.df.mat[which(is.na(f3_f4.df.mat))] = 0.00001
    
    p8 = get_ggplot_cor_density(
      f3_f4.df.mat, label = 'M8'
    )
  }
  
  ggarrange(#8*6
    p5, p6, p7, p8, 
    nrow = 2, ncol = 2,
    labels = c('A', 'B', 'C', 'D')
  )
  
}



#### correlation ####
if(T){
  f1.fot.5 = get_fot_df_by_expnames(fot_data, meta_data_1$Exp_code[meta_data_1$UPerson=='D5'])
  f1.fot.6 = get_fot_df_by_expnames(fot_data, meta_data_1$Exp_code[meta_data_1$UPerson=='D6'])
  f1.fot.7 = get_fot_df_by_expnames(fot_data, meta_data_1$Exp_code[meta_data_1$UPerson=='F7'])
  f1.fot.8 = get_fot_df_by_expnames(fot_data, meta_data_1$Exp_code[meta_data_1$UPerson=='M8'])
  
  f3.fot.5 = get_fot_df_by_expnames(fot_data, meta_data_3$Exp_code[meta_data_3$UPerson=='D5'])
  f3.fot.6 = get_fot_df_by_expnames(fot_data, meta_data_3$Exp_code[meta_data_3$UPerson=='D6'])
  f3.fot.7 = get_fot_df_by_expnames(fot_data, meta_data_3$Exp_code[meta_data_3$UPerson=='F7'])
  f3.fot.8 = get_fot_df_by_expnames(fot_data, meta_data_3$Exp_code[meta_data_3$UPerson=='M8'])
  
  f4.fot.5 = get_fot_df_by_expnames(fot_data, meta_data_4$Exp_code[meta_data_3$UPerson=='D5'])
  f4.fot.6 = get_fot_df_by_expnames(fot_data, meta_data_4$Exp_code[meta_data_3$UPerson=='D6'])
  f4.fot.7 = get_fot_df_by_expnames(fot_data, meta_data_4$Exp_code[meta_data_3$UPerson=='F7'])
  f4.fot.8 = get_fot_df_by_expnames(fot_data, meta_data_4$Exp_code[meta_data_3$UPerson=='M8'])
  
  get_ggplot_cor_density <- function(dat, label, t1, t2){
    # dat = log10(dat)
    df <- data.frame(
      x = dat[,1],
      y = dat[,2]
    )
    df.cor = cor(df$x, df$y)
    df.cor = round(df.cor, 2)
    df.cor.text = paste(label, '(N=', nrow(df), ')', ': Cor=', df.cor, sep = '')
    
    df = log10(df)
    p <- ggplot(df, aes(x=x, y=y))
    p <- p + stat_density2d(aes(fill=..density..), geom="tile", contour = F)
    p <- p + scale_fill_gradientn(colours=magma(100)[10:100])
    p <- p + xlab(paste("log10(FOT)/", t1, sep = '')) + ylab(paste("log10(FOT)/", t2, sep = '')) 
    p <- p + scale_x_continuous(breaks=seq(-5, 5, 1), limits = c(-5, 5)) 
    p <- p + scale_y_continuous(breaks=seq(-5, 5, 1), limits = c(-5, 5))
    p <- p + theme_classic()
    p <- p + annotate("text", x = 0, y = -3.5, label = df.cor.text, color='yellow', size=4)
    p
    return(p)
  }
  
  if(T){
    df1 = f1.fot.5
    f1.df = data.frame(
      GeneSymbol = df1$GeneSymbol,
      FOT1 = apply(df1[,-1], 1, mean)
    )
    df3 = f3.fot.5
    f3.df = data.frame(
      GeneSymbol = df3$GeneSymbol,
      FOT3 = apply(df3[,-1], 1, mean)
    )
    df4 = f4.fot.5
    f4.df = data.frame(
      GeneSymbol = df4$GeneSymbol,
      FOT4 = apply(df4[,-1], 1, mean)
    )
    
    f1_f3.df = merge(f1.df, f3.df, by = 'GeneSymbol', all = T)
    f1_f3.df.mat = as.matrix(f1_f3.df[,-1])
    f1_f3.df.mat[which(is.na(f1_f3.df.mat))] = 0.00001
    p5_13 = get_ggplot_cor_density(
      f1_f3.df.mat, label = 'D5', 'T1', 'T2'
    )
    
    f1_f4.df = merge(f1.df, f4.df, by = 'GeneSymbol', all = T)
    f1_f4.df.mat = as.matrix(f1_f4.df[,-1])
    f1_f4.df.mat[which(is.na(f1_f4.df.mat))] = 0.00001
    p5_14 = get_ggplot_cor_density(
      f1_f4.df.mat, label = 'D5', 'T1', 'T3'
    )
    
    
    f3_f4.df = merge(f3.df, f4.df, by = 'GeneSymbol', all = T)
    f3_f4.df.mat = as.matrix(f3_f4.df[,-1])
    f3_f4.df.mat[which(is.na(f3_f4.df.mat))] = 0.00001
    p5_34 = get_ggplot_cor_density(
      f3_f4.df.mat, label = 'D5', 'T2', 'T3'
    )
    
    ggarrange(#8*6
      p5_13, p5_14, p5_34, NULL,
      nrow = 2, ncol = 2,
      labels = c('A', 'B', 'C', 'D')
    )

  }
  
  
  
  if(T){
    df1 = f1.fot.6
    f1.df = data.frame(
      GeneSymbol = df1$GeneSymbol,
      FOT1 = apply(df1[,-1], 1, mean)
    )
    df3 = f3.fot.6
    f3.df = data.frame(
      GeneSymbol = df3$GeneSymbol,
      FOT3 = apply(df3[,-1], 1, mean)
    )
    df4 = f4.fot.6
    f4.df = data.frame(
      GeneSymbol = df4$GeneSymbol,
      FOT4 = apply(df4[,-1], 1, mean)
    )
    
    f1_f3.df = merge(f1.df, f3.df, by = 'GeneSymbol', all = T)
    f1_f3.df.mat = as.matrix(f1_f3.df[,-1])
    f1_f3.df.mat[which(is.na(f1_f3.df.mat))] = 0.00001
    p6_13 = get_ggplot_cor_density(
      f1_f3.df.mat, label = 'D6', 'T1', 'T2'
    )
    
    f1_f4.df = merge(f1.df, f4.df, by = 'GeneSymbol', all = T)
    f1_f4.df.mat = as.matrix(f1_f4.df[,-1])
    f1_f4.df.mat[which(is.na(f1_f4.df.mat))] = 0.00001
    p6_14 = get_ggplot_cor_density(
      f1_f4.df.mat, label = 'D6', 'T1', 'T3'
    )
    
    
    f3_f4.df = merge(f3.df, f4.df, by = 'GeneSymbol', all = T)
    f3_f4.df.mat = as.matrix(f3_f4.df[,-1])
    f3_f4.df.mat[which(is.na(f3_f4.df.mat))] = 0.00001
    p6_34 = get_ggplot_cor_density(
      f3_f4.df.mat, label = 'D6', 'T2', 'T3'
    )
    
    ggarrange(#8*6
      p6_13, p6_14, p6_34, NULL,
      nrow = 2, ncol = 2,
      labels = c('A', 'B', 'C', 'D')
    )
    
  }
 
  if(T){
    df1 = f1.fot.7
    f1.df = data.frame(
      GeneSymbol = df1$GeneSymbol,
      FOT1 = apply(df1[,-1], 1, mean)
    )
    df3 = f3.fot.7
    f3.df = data.frame(
      GeneSymbol = df3$GeneSymbol,
      FOT3 = apply(df3[,-1], 1, mean)
    )
    df4 = f4.fot.7
    f4.df = data.frame(
      GeneSymbol = df4$GeneSymbol,
      FOT4 = apply(df4[,-1], 1, mean)
    )
    
    f1_f3.df = merge(f1.df, f3.df, by = 'GeneSymbol', all = T)
    f1_f3.df.mat = as.matrix(f1_f3.df[,-1])
    f1_f3.df.mat[which(is.na(f1_f3.df.mat))] = 0.00001
    p7_13 = get_ggplot_cor_density(
      f1_f3.df.mat, label = 'F7', 'T1', 'T2'
    )
    
    f1_f4.df = merge(f1.df, f4.df, by = 'GeneSymbol', all = T)
    f1_f4.df.mat = as.matrix(f1_f4.df[,-1])
    f1_f4.df.mat[which(is.na(f1_f4.df.mat))] = 0.00001
    p7_14 = get_ggplot_cor_density(
      f1_f4.df.mat, label = 'F7', 'T1', 'T3'
    )
    
    
    f3_f4.df = merge(f3.df, f4.df, by = 'GeneSymbol', all = T)
    f3_f4.df.mat = as.matrix(f3_f4.df[,-1])
    f3_f4.df.mat[which(is.na(f3_f4.df.mat))] = 0.00001
    p7_34 = get_ggplot_cor_density(
      f3_f4.df.mat, label = 'F7', 'T2', 'T3'
    )
    
    ggarrange(#8*6
      p7_13, p7_14, p7_34, NULL,
      nrow = 2, ncol = 2,
      labels = c('A', 'B', 'C', 'D')
    )
    
  }
 
  
  if(T){
    df1 = f1.fot.8
    f1.df = data.frame(
      GeneSymbol = df1$GeneSymbol,
      FOT1 = apply(df1[,-1], 1, mean)
    )
    df3 = f3.fot.8
    f3.df = data.frame(
      GeneSymbol = df3$GeneSymbol,
      FOT3 = apply(df3[,-1], 1, mean)
    )
    df4 = f4.fot.8
    f4.df = data.frame(
      GeneSymbol = df4$GeneSymbol,
      FOT4 = apply(df4[,-1], 1, mean)
    )
    
    f1_f3.df = merge(f1.df, f3.df, by = 'GeneSymbol', all = T)
    f1_f3.df.mat = as.matrix(f1_f3.df[,-1])
    f1_f3.df.mat[which(is.na(f1_f3.df.mat))] = 0.00001
    p8_13 = get_ggplot_cor_density(
      f1_f3.df.mat, label = 'M8', 'T1', 'T2'
    )
    
    f1_f4.df = merge(f1.df, f4.df, by = 'GeneSymbol', all = T)
    f1_f4.df.mat = as.matrix(f1_f4.df[,-1])
    f1_f4.df.mat[which(is.na(f1_f4.df.mat))] = 0.00001
    p8_14 = get_ggplot_cor_density(
      f1_f4.df.mat, label = 'M8', 'T1', 'T3'
    )
    
    
    f3_f4.df = merge(f3.df, f4.df, by = 'GeneSymbol', all = T)
    f3_f4.df.mat = as.matrix(f3_f4.df[,-1])
    f3_f4.df.mat[which(is.na(f3_f4.df.mat))] = 0.00001
    p8_34 = get_ggplot_cor_density(
      f3_f4.df.mat, label = 'M8', 'T2', 'T3'
    )
    
    ggarrange(#8*6
      p8_13, p8_14, p8_34, NULL,
      nrow = 2, ncol = 2,
      labels = c('A', 'B', 'C', 'D')
    )
    
  }
  
  
  ggarrange(#8*6
    p5_13, p5_14, p5_34,
    p6_13, p6_14, p6_34,
    p7_13, p7_14, p7_34,
    p8_13, p8_14, p8_34,
    nrow = 4, ncol = 3
    # labels = c('A', 'B', 'C', 'D')
  )
}



#### correlation ####
if(T){
  f1.fot.5 = get_fot_df_by_expnames(fot_data, meta_data_1$Exp_code[meta_data_1$UPerson=='D5'])
  f1.fot.6 = get_fot_df_by_expnames(fot_data, meta_data_1$Exp_code[meta_data_1$UPerson=='D6'])
  f1.fot.7 = get_fot_df_by_expnames(fot_data, meta_data_1$Exp_code[meta_data_1$UPerson=='F7'])
  f1.fot.8 = get_fot_df_by_expnames(fot_data, meta_data_1$Exp_code[meta_data_1$UPerson=='M8'])
  
  f3.fot.5 = get_fot_df_by_expnames(fot_data, meta_data_3$Exp_code[meta_data_3$UPerson=='D5'])
  f3.fot.6 = get_fot_df_by_expnames(fot_data, meta_data_3$Exp_code[meta_data_3$UPerson=='D6'])
  f3.fot.7 = get_fot_df_by_expnames(fot_data, meta_data_3$Exp_code[meta_data_3$UPerson=='F7'])
  f3.fot.8 = get_fot_df_by_expnames(fot_data, meta_data_3$Exp_code[meta_data_3$UPerson=='M8'])
  
  f4.fot.5 = get_fot_df_by_expnames(fot_data, meta_data_4$Exp_code[meta_data_3$UPerson=='D5'])
  f4.fot.6 = get_fot_df_by_expnames(fot_data, meta_data_4$Exp_code[meta_data_3$UPerson=='D6'])
  f4.fot.7 = get_fot_df_by_expnames(fot_data, meta_data_4$Exp_code[meta_data_3$UPerson=='F7'])
  f4.fot.8 = get_fot_df_by_expnames(fot_data, meta_data_4$Exp_code[meta_data_3$UPerson=='M8'])
  
  f.metadata = rbind.data.frame(
    meta_data_1, meta_data_3, meta_data_4
  )
  px = 'D6'
  f.metadata.px = f.metadata[f.metadata$UPerson==px,]
  f.fot.px = get_fot_df_by_expnames(fot_data, f.metadata.px$Exp_code)
  
  if(T){
    px.fot = f.fot.px
    px.fot.mat = as.matrix(px.fot[,-1])
    # colnames(px.fot.mat) = c(
    #   paste('NPB', c('R1', 'R2', 'R3'), sep = '_'),
    #   paste('BGI', c('R1', 'R2', 'R3'), sep = '_')
    # )
    colnames(px.fot.mat) = c(
      paste('RI', c('R1', 'R2', 'R3'), sep = '_'),
      paste('CI1', c('R1', 'R2', 'R3'), sep = '_'),
      paste('CI2', c('R1', 'R2', 'R3'), sep = '_')
    )
    na.index = which(is.na(px.fot.mat))
    px.fot.mat[na.index] = 0
    cor_method = 'pearson' 
    
    cormat <- round(cor(px.fot.mat, method = cor_method),2)
    
    cormat <- matrix(0, 9, 9)
    rownames(cormat) = c(
      paste('RI', c('R1', 'R2', 'R3'), sep = '_'),
      paste('CI1', c('R1', 'R2', 'R3'), sep = '_'),
      paste('CI2', c('R1', 'R2', 'R3'), sep = '_')
    )
    colnames(cormat) = c(
      paste('RI', c('R1', 'R2', 'R3'), sep = '_'),
      paste('CI1', c('R1', 'R2', 'R3'), sep = '_'),
      paste('CI2', c('R1', 'R2', 'R3'), sep = '_')
    )
    for (i in 1:9) {
      x = px.fot.mat[,i]
      for (j in 1:9) {
        y = px.fot.mat[,j]
        index = which(x>0 & y>0)
        xy.cor = round(cor(x[index], y[index], method = cor_method), 2)
        cormat[i,j] = xy.cor
      }
      
    }
    
    
    library(reshape2)
    melted_cormat <- melt(cormat)
    library(ggplot2)
    p1 <- ggplot(data = melted_cormat, aes(x=Var1, y=Var2, fill=value)) 
    p1 <- p1 +  geom_tile() 
    p1 <- p1 + xlab('Site-Replicates') + ylab('Site-Replicates')
    p1 <- p1 + theme(
      panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
      panel.background = element_blank(), axis.line = element_line(colour = "black")
    )
    # p1 <- p1 + scale_fill_gradient2(
    #   low = "#132B43", high ="#56B1F7", mid = "white",
    #   limits = c(0.5, 1), breaks = seq(0.5,1,0.1),
    #   guide = guide_colourbar(nbin=100, draw.ulim = FALSE, draw.llim = FALSE)
    # )
    p1 <- p1 + scale_fill_continuous(
      limits = c(0.7, 1), breaks = seq(0.7, 1, 0.1),
      guide = guide_colourbar(nbin=100, draw.ulim = FALSE, draw.llim = FALSE)
    ) + labs(fill = "PCorr")
    p1 <- p1 + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
    p1 + geom_text(aes(label = round(value, 2))) + ggtitle(px)
    # 6*5
    # col<- colorRampPalette(c("blue", "white", "red"))(20)
    # heatmap(x = cormat, col = col, symm = TRUE)
  }
  
}



#### PCA ####
if(T){
  if(T){
    uid.metadata = meta_data_1
    uid.metadata = uid.metadata[order(uid.metadata$Person),]
    uid.expnames = uid.metadata$Exp_code
    uid.df = get_fot_df_by_expnames(fot_data, uid.expnames)
    uid.df.genesymbol = as.vector(uid.df$GeneSymbol)
    uid.df.mat = as.vector(uid.df[,-1])
    uid.df.mat.col_id = apply(uid.df.mat, 2, function(x){length(which(x>0))})
    uid.df.mat.scale = t(scale(t(uid.df.mat)))
    mat = uid.df.mat.scale
    uid.df1 = data.frame(
      GeneSymbol = uid.df.genesymbol,
      mat
    )
    fot_group = uid.metadata$UPerson
    fdu_new.snr = round(snrdb_function(mat, fot_group),2)
    colors = rep(p_colors, each=3)
    pca <- prcomp(((t(mat))), center = F, scale = F)
    importance <- summary(pca)$importance
    pc1 = paste('PC1 (', round(importance[2,1]*100, 2), '%)', sep = '')
    pc2 = paste('PC2 (', round(importance[2,2]*100, 2), '%)', sep = '')
    pca_predict <- predict(pca)
    pca_predict.2d = pca_predict[,c(1,2)]
    main = paste('T1', '; All', ' Proteins (', nrow(mat), ')\nSNR = ', fdu_new.snr, sep = '') 
    
    df = data.frame(
      x = pca_predict.2d[,1],
      y = pca_predict.2d[,2],
      label = rownames(pca_predict.2d),
      Sample = factor(uid.metadata$UPerson),
      Order = uid.metadata$Order[match(rownames(pca_predict.2d), uid.metadata$Exp_code)]
    )
    df$Order = as.factor(df$Order)
    
    library(ggrepel)
    main = paste('T1', '; All', ' Proteins (', nrow(mat), ')\nSNR = ', fdu_new.snr, sep = '') 
    p <- ggplot(df, aes(x=x, y=y))
    p <- p + geom_point(aes(color = Sample, shape = Order), size=3, position = "identity")
    p <- p + scale_color_manual(values=p_colors)
    
    p <- p + theme_classic(base_size = 16)
    p <- p + xlim(-80, 60) + ylim(-60, 60) + xlab(pc1) + ylab(pc2) 
    
    p1 <- p + ggtitle(main) +  theme(plot.title = element_text(size = 15, face = "bold"))
    p1
    
  }
  
  
  if(T){
    uid.metadata = meta_data_3
    uid.metadata = uid.metadata[order(uid.metadata$Person),]
    uid.expnames = uid.metadata$Exp_code
    uid.df = get_fot_df_by_expnames(fot_data, uid.expnames)
    uid.df.genesymbol = as.vector(uid.df$GeneSymbol)
    uid.df.mat = as.vector(uid.df[,-1])
    uid.df.mat.col_id = apply(uid.df.mat, 2, function(x){length(which(x>0))})
    uid.df.mat.scale = t(scale(t(uid.df.mat)))
    mat = uid.df.mat.scale
    uid.df3 = data.frame(
      GeneSymbol = uid.df.genesymbol,
      mat
    )
    fot_group = uid.metadata$UPerson
    fdu_new.snr = round(snrdb_function(mat, fot_group),2)
    colors = rep(p_colors, each=3)
    pca <- prcomp(((t(mat))), center = F, scale = F)
    importance <- summary(pca)$importance
    pc1 = paste('PC1 (', round(importance[2,1]*100, 2), '%)', sep = '')
    pc2 = paste('PC2 (', round(importance[2,2]*100, 2), '%)', sep = '')
    pca_predict <- predict(pca)
    pca_predict.2d = pca_predict[,c(1,2)]
    main = paste('F3', '; All', ' Proteins (', nrow(mat), ')\nSNR = ', fdu_new.snr, sep = '') 
    
    df = data.frame(
      x = pca_predict.2d[,1],
      y = pca_predict.2d[,2],
      label = rownames(pca_predict.2d),
      Sample = factor(uid.metadata$UPerson),
      Order = uid.metadata$Order[match(rownames(pca_predict.2d), uid.metadata$Exp_code)]
    )
    df$Order = as.factor(df$Order)
    
    library(ggrepel)
    main = paste('T2', '; All', ' Proteins (', nrow(mat), ')\nSNR = ', fdu_new.snr, sep = '') 
    p <- ggplot(df, aes(x=x, y=y))
    p <- p + geom_point(aes(color = Sample, shape = Order), size=3, position = "identity")
    p <- p + scale_color_manual(values=p_colors)
    
    p <- p + theme_classic(base_size = 16)
    p <- p + xlim(-80, 60) + ylim(-60, 60) + xlab(pc1) + ylab(pc2) 
    
    p3 <- p + ggtitle(main) +  theme(plot.title = element_text(size = 15, face = "bold"))
    p3
    
  }
  
  
  if(T){
    uid.metadata = meta_data_4
    uid.metadata = uid.metadata[order(uid.metadata$Person),]
    uid.expnames = uid.metadata$Exp_code
    uid.df = get_fot_df_by_expnames(fot_data, uid.expnames)
    uid.df.genesymbol = as.vector(uid.df$GeneSymbol)
    uid.df.mat = as.vector(uid.df[,-1])
    uid.df.mat.col_id = apply(uid.df.mat, 2, function(x){length(which(x>0))})
    uid.df.mat.scale = t(scale(t(uid.df.mat)))
    mat = uid.df.mat.scale
    uid.df4 = data.frame(
      GeneSymbol = uid.df.genesymbol,
      mat
    )
    fot_group = uid.metadata$UPerson
    fdu_new.snr = round(snrdb_function(mat, fot_group),2)
    colors = rep(p_colors, each=3)
    pca <- prcomp(((t(mat))), center = F, scale = F)
    importance <- summary(pca)$importance
    pc1 = paste('PC1 (', round(importance[2,1]*100, 2), '%)', sep = '')
    pc2 = paste('PC2 (', round(importance[2,2]*100, 2), '%)', sep = '')
    pca_predict <- predict(pca)
    pca_predict.2d = pca_predict[,c(1,2)]
    main = paste('F4', '; All', ' Proteins (', nrow(mat), ')\nSNR = ', fdu_new.snr, sep = '') 
    
    df = data.frame(
      x = pca_predict.2d[,1],
      y = pca_predict.2d[,2],
      label = rownames(pca_predict.2d),
      Sample = factor(uid.metadata$UPerson),
      Order = uid.metadata$Order[match(rownames(pca_predict.2d), uid.metadata$Exp_code)]
    )
    df$Order = as.factor(df$Order)
    
    library(ggrepel)
    main = paste('T3', '; All', ' Proteins (', nrow(mat), ')\nSNR = ', fdu_new.snr, sep = '') 
    p <- ggplot(df, aes(x=x, y=y))
    p <- p + geom_point(aes(color = Sample, shape = Order), size=3, position = "identity")
    p <- p + scale_color_manual(values=p_colors)
    
    p <- p + theme_classic(base_size = 16)
    p <- p + xlim(-80, 60) + ylim(-60, 60) + xlab(pc1) + ylab(pc2) 
    
    p4 <- p + ggtitle(main) +  theme(plot.title = element_text(size = 15, face = "bold"))
    p4
    
  }
  
  if(T){
    uid.metadata = rbind.data.frame(meta_data_1, meta_data_3, meta_data_4)
    uid.metadata$Flag = rep(c('T1', 'T2', 'T3'), each = 12)
    uid.expnames = uid.metadata$Exp_code
    uid.df = get_fot_df_by_expnames(fot_data, uid.expnames)
    uid.df.genesymbol = as.vector(uid.df$GeneSymbol)
    uid.df.mat = as.vector(uid.df[,-1])
    uid.df.mat.col_id = apply(uid.df.mat, 2, function(x){length(which(x>0))})
    uid.df.mat.scale = t(scale(t(uid.df.mat)))
    mat = uid.df.mat.scale
    
    fot_group = paste(uid.metadata$Flag, uid.metadata$UPerson, sep = '_')
    fot_group = uid.metadata$UPerson
    fdu_new.snr = round(snrdb_function(mat, fot_group),2)
    colors = rep(p_colors, each=3)
    pca <- prcomp(((t(mat))), center = F, scale = F)
    importance <- summary(pca)$importance
    pc1 = paste('PC1 (', round(importance[2,1]*100, 2), '%)', sep = '')
    pc2 = paste('PC2 (', round(importance[2,2]*100, 2), '%)', sep = '')
    pca_predict <- predict(pca)
    pca_predict.2d = pca_predict[,c(1,2)]
    main = paste('F3_4', '; All', ' Proteins (', nrow(mat), ')\nSNR = ', fdu_new.snr, ' (Before correction)', sep = '') 
    
    df = data.frame(
      x = pca_predict.2d[,1],
      y = pca_predict.2d[,2],
      label = rownames(pca_predict.2d),
      Sample = factor(uid.metadata$UPerson),
      Batch = uid.metadata$Flag
    )
    df$Batch = as.factor(df$Batch)
    
    library(ggrepel)
    p <- ggplot(df, aes(x=x, y=y))
    p <- p + geom_point(aes(color = Sample, shape = Batch), size=3, position = "identity")
    p <- p + scale_color_manual(values=p_colors)
    
    p <- p + theme_classic(base_size = 16)
    p <- p + xlim(-80, 60) + ylim(-60, 60) + xlab(pc1) + ylab(pc2) 
    
    pb <- p + ggtitle(main) +  theme(plot.title = element_text(size = 15, face = "bold"))
    pb
    
  }
  
  
  
  
  
  
  
  
  
  
  
  
  
 
  
  
  if(T){
    uid.metadata = rbind.data.frame(meta_data_1, meta_data_3, meta_data_4)
    uid.metadata$Flag = rep(c('T1', 'T2', 'T3'), each = 12)
    
    uid.df = merge(
      uid.df1, uid.df3, by='GeneSymbol', all = T
    ) 
    
    uid.df = merge(
      uid.df, uid.df4, by='GeneSymbol', all = T
    ) 
    mat = na.omit(as.matrix(uid.df[,-1]))
    
    
    fot_group = uid.metadata$UPerson
    fdu_new.snr = round(snrdb_function(mat, fot_group),2)
    colors = rep(p_colors, each=3)
    pca <- prcomp(((t(mat))), center = F, scale = F)
    importance <- summary(pca)$importance
    pc1 = paste('PC1 (', round(importance[2,1]*100, 2), '%)', sep = '')
    pc2 = paste('PC2 (', round(importance[2,2]*100, 2), '%)', sep = '')
    pca_predict <- predict(pca)
    pca_predict.2d = pca_predict[,c(1,2)]
    main = paste('T1/2/3', '; All', ' Proteins (', nrow(mat), ')\nSNR = ', fdu_new.snr, sep = '') 
    
    df = data.frame(
      x = pca_predict.2d[,1],
      y = pca_predict.2d[,2],
      label = rownames(pca_predict.2d),
      Sample = factor(uid.metadata$UPerson),
      Batch = uid.metadata$Flag
    )
    df$Batch = as.factor(df$Batch)
    
    library(ggrepel)
    p <- ggplot(df, aes(x=x, y=y))
    p <- p + geom_point(aes(color = Sample, shape = Batch), size=3, position = "identity")
    p <- p + scale_color_manual(values=p_colors)
    
    p <- p + theme_classic(base_size = 16)
    p <- p + xlim(-80, 60) + ylim(-60, 60) + xlab(pc1) + ylab(pc2) 
    
    pa <- p + ggtitle(main) +  theme(plot.title = element_text(size = 15, face = "bold"))
    pa
    
  }
  
  ggarrange(#8*6
    p1, p3, p4, pa, 
    nrow = 2, ncol = 2,
    labels = c('A', 'B', 'C', 'D')
  )
  
  
}
