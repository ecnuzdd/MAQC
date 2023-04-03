# meta_data
if(T){
  px = 'P5'
  meta_data.px = meta_data[which(meta_data$Person == 'P5'),]
  px.fot = get_fot_df_by_expnames(fot_data, meta_data.px$Exp_code)
  px.fot.row_id = apply(px.fot[,-1], 1, function(x){
    length(which(x>0))
  })
  
  id_cutoff = seq(0, 1, 0.1)
  i_end = length(id_cutoff)-1
  id_v = NULL
  down_v = NULL
  up_v = NULL
  for (i in 1:i_end) {
    down_i = id_cutoff[i]
    up_i = id_cutoff[i+1]
    
    down_i_count = down_i*110
    up_i_count = up_i*110
    
    if(i!=i_end){
      index = which(px.fot.row_id>down_i_count & px.fot.row_id<=up_i_count)
    }else{
      index = which(px.fot.row_id>down_i_count & px.fot.row_id<up_i_count)
    }
    
    id_v = c(id_v, length(index))
    down_v = c(down_v, down_i)
    up_v = c(up_v, up_i)
  }
  df = data.frame(
    down = c(down_v, 1),
    up = c(up_v,1),
    count = c(id_v, length(which(px.fot.row_id==110)))
  )
  df$cumsum = cumsum(df$count)
  df$labels = c(
    '(0%, 10%]', '(10%, 20%]', '(20%, 30%]', '(30%, 40%]', '(40%, 50%]', '(50%, 60%]',
    '(60%, 70%]', '(70%, 80%]', '(80%, 90%]', '(90%, 100%)','[100%, 100%]'
  )
  
  p <- ggbarplot(
    df, x = "labels", y = "count",
    fill = "blue",               # change fill color by cyl
    color = "white",            # Set bar border colors to white
    # palette = "jco",            # jco journal color palett. see ?ggpar
    # sort.val = NULL,           # Sort the value in dscending order
    # sort.by.groups = TRUE,      # Sort inside each group
    x.text.angle = 90,           # Rotate vertically x axis texts,
    rotate = TRUE
    # ggtheme = theme_minimal()
  )
  
  df$labels = factor(df$labels, levels = df$labels)
  
  p <- ggplot(df, aes(x=labels, y=count)) 
  p <- p + geom_bar(aes(fill = (down)), stat="identity", width=0.7) 
  # p <- p + coord_flip() + theme_minimal()
  p <- p +  scale_fill_gradient(low = "yellow", high = "red", na.value = NA)
  p <- p + ylab('Protein count') + xlab('Confidence interval')
  p <- p + scale_y_continuous(breaks=seq(0, 5000, 1000), limits=c(0, 5000))
  p <- p + geom_text(data=df, aes(x=labels,y=count,label=count,), vjust= -0.5, hjust = 0.5, size=4)
  p <- p + guides(fill=FALSE) + theme_classic() 
  p5 <- p + ggtitle('P5')
  print(p5)
  
}


if(T){
  px = 'P6'
  meta_data.px = meta_data[which(meta_data$Person == px),]
  px.fot = get_fot_df_by_expnames(fot_data, meta_data.px$Exp_code)
  px.fot.row_id = apply(px.fot[,-1], 1, function(x){
    length(which(x>0))
  })
  
  id_cutoff = seq(0, 1, 0.1)
  i_end = length(id_cutoff)-1
  id_v = NULL
  down_v = NULL
  up_v = NULL
  for (i in 1:i_end) {
    down_i = id_cutoff[i]
    up_i = id_cutoff[i+1]
    
    down_i_count = down_i*110
    up_i_count = up_i*110
    
    if(i!=i_end){
      index = which(px.fot.row_id>down_i_count & px.fot.row_id<=up_i_count)
    }else{
      index = which(px.fot.row_id>down_i_count & px.fot.row_id<up_i_count)
    }
    
    id_v = c(id_v, length(index))
    down_v = c(down_v, down_i)
    up_v = c(up_v, up_i)
  }
  df = data.frame(
    down = c(down_v, 1),
    up = c(up_v,1),
    count = c(id_v, length(which(px.fot.row_id==110)))
  )
  df$cumsum = cumsum(df$count)
  df$labels = c(
    '(0%, 10%]', '(10%, 20%]', '(20%, 30%]', '(30%, 40%]', '(40%, 50%]', '(50%, 60%]',
    '(60%, 70%]', '(70%, 80%]', '(80%, 90%]', '(90%, 100%)','[100%, 100%]'
  )
  
  p <- ggbarplot(
    df, x = "labels", y = "count",
    fill = "blue",               # change fill color by cyl
    color = "white",            # Set bar border colors to white
    # palette = "jco",            # jco journal color palett. see ?ggpar
    # sort.val = NULL,           # Sort the value in dscending order
    # sort.by.groups = TRUE,      # Sort inside each group
    x.text.angle = 90,           # Rotate vertically x axis texts,
    rotate = TRUE
    # ggtheme = theme_minimal()
  )
  
  df$labels = factor(df$labels, levels = df$labels)
  
  p <- ggplot(df, aes(x=labels, y=count)) 
  p <- p + geom_bar(aes(fill = (down)), stat="identity", width=0.7) 
  # p <- p + coord_flip() + theme_minimal()
  p <- p +  scale_fill_gradient(low = "yellow", high = "red", na.value = NA)
  p <- p + ylab('Protein count') + xlab('Confidence interval')
  p <- p + scale_y_continuous(breaks=seq(0, 5000, 1000), limits=c(0, 5000))
  p <- p + geom_text(data=df, aes(x=labels,y=count,label=count,), vjust= -0.5, hjust = 0.5, size=4)
  p <- p + guides(fill=FALSE) + theme_classic() 
  p6 <- p + ggtitle('P6')
  print(p6)
  
}


if(T){
  px = 'P7'
  meta_data.px = meta_data[which(meta_data$Person == px),]
  px.fot = get_fot_df_by_expnames(fot_data, meta_data.px$Exp_code)
  px.fot.row_id = apply(px.fot[,-1], 1, function(x){
    length(which(x>0))
  })
  
  id_cutoff = seq(0, 1, 0.1)
  i_end = length(id_cutoff)-1
  id_v = NULL
  down_v = NULL
  up_v = NULL
  for (i in 1:i_end) {
    down_i = id_cutoff[i]
    up_i = id_cutoff[i+1]
    
    down_i_count = down_i*110
    up_i_count = up_i*110
    
    if(i!=i_end){
      index = which(px.fot.row_id>down_i_count & px.fot.row_id<=up_i_count)
    }else{
      index = which(px.fot.row_id>down_i_count & px.fot.row_id<up_i_count)
    }
    
    id_v = c(id_v, length(index))
    down_v = c(down_v, down_i)
    up_v = c(up_v, up_i)
  }
  df = data.frame(
    down = c(down_v, 1),
    up = c(up_v,1),
    count = c(id_v, length(which(px.fot.row_id==110)))
  )
  df$cumsum = cumsum(df$count)
  df$labels = c(
    '(0%, 10%]', '(10%, 20%]', '(20%, 30%]', '(30%, 40%]', '(40%, 50%]', '(50%, 60%]',
    '(60%, 70%]', '(70%, 80%]', '(80%, 90%]', '(90%, 100%)','[100%, 100%]'
  )
  
  p <- ggbarplot(
    df, x = "labels", y = "count",
    fill = "blue",               # change fill color by cyl
    color = "white",            # Set bar border colors to white
    # palette = "jco",            # jco journal color palett. see ?ggpar
    # sort.val = NULL,           # Sort the value in dscending order
    # sort.by.groups = TRUE,      # Sort inside each group
    x.text.angle = 90,           # Rotate vertically x axis texts,
    rotate = TRUE
    # ggtheme = theme_minimal()
  )
  
  df$labels = factor(df$labels, levels = df$labels)
  
  p <- ggplot(df, aes(x=labels, y=count)) 
  p <- p + geom_bar(aes(fill = (down)), stat="identity", width=0.7) 
  # p <- p + coord_flip() + theme_minimal()
  p <- p +  scale_fill_gradient(low = "yellow", high = "red", na.value = NA)
  p <- p + ylab('Protein count') + xlab('Confidence interval')
  p <- p + scale_y_continuous(breaks=seq(0, 5000, 1000), limits=c(0, 5000))
  p <- p + geom_text(data=df, aes(x=labels,y=count,label=count,), vjust= -0.5, hjust = 0.5, size=4)
  p <- p + guides(fill=FALSE) + theme_classic() 
  p7 <- p + ggtitle('P7')
  print(p7)
  
}

if(T){
  px = 'P8'
  meta_data.px = meta_data[which(meta_data$Person == px),]
  px.fot = get_fot_df_by_expnames(fot_data, meta_data.px$Exp_code)
  px.fot.row_id = apply(px.fot[,-1], 1, function(x){
    length(which(x>0))
  })
  
  id_cutoff = seq(0, 1, 0.1)
  i_end = length(id_cutoff)-1
  id_v = NULL
  down_v = NULL
  up_v = NULL
  for (i in 1:i_end) {
    down_i = id_cutoff[i]
    up_i = id_cutoff[i+1]
    
    down_i_count = down_i*110
    up_i_count = up_i*110
    
    if(i!=i_end){
      index = which(px.fot.row_id>down_i_count & px.fot.row_id<=up_i_count)
    }else{
      index = which(px.fot.row_id>down_i_count & px.fot.row_id<up_i_count)
    }
    
    id_v = c(id_v, length(index))
    down_v = c(down_v, down_i)
    up_v = c(up_v, up_i)
  }
  df = data.frame(
    down = c(down_v, 1),
    up = c(up_v,1),
    count = c(id_v, length(which(px.fot.row_id==110)))
  )
  df$cumsum = cumsum(df$count)
  df$labels = c(
    '(0%, 10%]', '(10%, 20%]', '(20%, 30%]', '(30%, 40%]', '(40%, 50%]', '(50%, 60%]',
    '(60%, 70%]', '(70%, 80%]', '(80%, 90%]', '(90%, 100%)','[100%, 100%]'
  )
  
  p <- ggbarplot(
    df, x = "labels", y = "count",
    fill = "blue",               # change fill color by cyl
    color = "white",            # Set bar border colors to white
    # palette = "jco",            # jco journal color palett. see ?ggpar
    # sort.val = NULL,           # Sort the value in dscending order
    # sort.by.groups = TRUE,      # Sort inside each group
    x.text.angle = 90,           # Rotate vertically x axis texts,
    rotate = TRUE
    # ggtheme = theme_minimal()
  )
  
  df$labels = factor(df$labels, levels = df$labels)
  
  p <- ggplot(df, aes(x=labels, y=count)) 
  p <- p + geom_bar(aes(fill = (down)), stat="identity", width=0.7) 
  # p <- p + coord_flip() + theme_minimal()
  p <- p +  scale_fill_gradient(low = "yellow", high = "red", na.value = NA)
  p <- p + ylab('Protein count') + xlab('Confidence interval')
  p <- p + scale_y_continuous(breaks=seq(0, 5000, 1000), limits=c(0, 5000))
  p <- p + geom_text(data=df, aes(x=labels,y=count,label=count,), vjust= -0.5, hjust = 0.5, size=4)
  p <- p + guides(fill=FALSE) + theme_classic() 
  p8 <- p + ggtitle('P8')
  print(p8)
  
}


if(T){
  meta_data.px = meta_data#[which(meta_data$Person == px),]
  px.fot = get_fot_df_by_expnames(fot_data, meta_data.px$Exp_code)
  px.fot.row_id = apply(px.fot[,-1], 1, function(x){
    length(which(x>0))
  })
  
  id_cutoff = seq(0, 1, 0.1)
  i_end = length(id_cutoff)-1
  id_v = NULL
  down_v = NULL
  up_v = NULL
  for (i in 1:i_end) {
    down_i = id_cutoff[i]
    up_i = id_cutoff[i+1]
    
    down_i_count = down_i*110*4
    up_i_count = up_i*110*4
    
    if(i!=i_end){
      index = which(px.fot.row_id>down_i_count & px.fot.row_id<=up_i_count)
    }else{
      index = which(px.fot.row_id>down_i_count & px.fot.row_id<up_i_count)
    }
    
    id_v = c(id_v, length(index))
    down_v = c(down_v, down_i)
    up_v = c(up_v, up_i)
  }
  df = data.frame(
    down = c(down_v, 1),
    up = c(up_v,1),
    count = c(id_v, length(which(px.fot.row_id==440)))
  )
  df$cumsum = cumsum(df$count)
  df$labels = c(
    '(0%, 10%]', '(10%, 20%]', '(20%, 30%]', '(30%, 40%]', '(40%, 50%]', '(50%, 60%]',
    '(60%, 70%]', '(70%, 80%]', '(80%, 90%]', '(90%, 100%)','[100%, 100%]'
  )
  
  p <- ggbarplot(
    df, x = "labels", y = "count",
    fill = "blue",               # change fill color by cyl
    color = "white",            # Set bar border colors to white
    # palette = "jco",            # jco journal color palett. see ?ggpar
    # sort.val = NULL,           # Sort the value in dscending order
    # sort.by.groups = TRUE,      # Sort inside each group
    x.text.angle = 90,           # Rotate vertically x axis texts,
    rotate = TRUE
    # ggtheme = theme_minimal()
  )
  
  df$labels = factor(df$labels, levels = df$labels)
  
  p <- ggplot(df, aes(x=labels, y=count)) 
  p <- p + geom_bar(aes(fill = (down)), stat="identity", width=0.7) 
  # p <- p + coord_flip() + theme_minimal()
  p <- p +  scale_fill_gradient(low = "yellow", high = "red", na.value = NA)
  p <- p + ylab('Protein count') + xlab('Confidence interval')
  p <- p + scale_y_continuous(breaks=seq(0, 6500, 1000), limits=c(0, 6500))
  p <- p + geom_text(data=df, aes(x=labels,y=count,label=count,), vjust= -0.5, hjust = 0.5, size=4)
  p <- p + guides(fill=FALSE) + theme_classic() 
  p <- p + ggtitle('All_440')
  print(p)
  
  
}


# ggarrange(
#   p,                                                 # First row with scatter plot
#   ggarrange(p5, p6, p7, p8, ncol = 2, nrow = 2), # Second row with box and dot plots
#   nrow = 2
#   # labels = "A"                                        # Labels of the scatter plot
# ) 

ggarrange(p, p5, p6, p7, p8, ncol = 1, nrow = 5) # Second row with box and dot plots

