source('D:/MAQC/20200508/Analysis/code_library/00_initialization.R')

meta_data.replicates = meta_data[meta_data$Purpose=='Replications',]

meta_data.replicates.uids = paste(
  meta_data.replicates$Affiliation, 
  meta_data.replicates$SInstrument,
  meta_data.replicates$Batch, sep = '_'
)

meta_data.replicates.uids.unique = sort(unique(meta_data.replicates.uids))

uid.deps_list.list = list()

for (i in 1:length(meta_data.replicates.uids.unique)) {
  uid = meta_data.replicates.uids.unique[i]
  uid.meta_data = meta_data.replicates[meta_data.replicates.uids==uid,]
  uid.fot_data = get_fot_df_by_expnames(fot_data, uid.meta_data$Exp_code)
  uid.fot_data.imputed = impute_df_with_value(uid.fot_data, 0)
  uid.deps_list = get_deps_df(uid.fot_data.imputed, uid.meta_data$Person, minFC=2, minPvalue = 0.05)
  uid.deps_list.list[[i]] = uid.deps_list
  print('##########################')
  print(i)
  print('##########################')
}

# saveRDS(uid.deps_list.list, "uid.deps_list.list.rds")
# uid.deps_list.list <- readRDS("uid.deps_list.list.rds")
uid.deps_list.list_count = length(uid.deps_list.list)

if(T){
  # P5_P6
  p5_p6.deps_list = list()
  for (i in 1:uid.deps_list.list_count) {
    tmp_list = uid.deps_list.list[[i]][[1]]
    p5_p6.deps_list[[i]] = tmp_list
  }
  p5_p6.mdf = merge_deps_list(p5_p6.deps_list, meta_data.replicates.uids.unique)
  p5_p6.deps_df = check_real_deps_df(p5_p6.mdf, min_id=2)
  p5_p6.col_id = apply(p5_p6.mdf[,-1], 2, function(x){
    x = as.vector(unlist(x))
    length(which(!is.na(x)))
  })
  # lapply(p5_p6.deps_list, colnames)
  
  # P5_P7
  p5_p7.deps_list = list()
  for (i in 1:uid.deps_list.list_count) {
    tmp_list = uid.deps_list.list[[i]][[2]]
    p5_p7.deps_list[[i]] = tmp_list
  }
  p5_p7.mdf = merge_deps_list(p5_p7.deps_list, meta_data.replicates.uids.unique)
  p5_p7.deps_df = check_real_deps_df(p5_p7.mdf, min_id=2)
  p5_p7.col_id = apply(p5_p7.mdf[,-1], 2, function(x){
    x = as.vector(unlist(x))
    length(which(!is.na(x)))
  })
  
  
  # P5_P8
  p5_p8.deps_list = list()
  for (i in 1:uid.deps_list.list_count) {
    tmp_list = uid.deps_list.list[[i]][[3]]
    p5_p8.deps_list[[i]] = tmp_list
  }
  p5_p8.mdf = merge_deps_list(p5_p8.deps_list, meta_data.replicates.uids.unique)
  p5_p8.deps_df = check_real_deps_df(p5_p8.mdf, min_id=2)
  p5_p8.col_id = apply(p5_p8.mdf[,-1], 2, function(x){
    x = as.vector(unlist(x))
    length(which(!is.na(x)))
  })
  
  # P6_P7
  p6_p7.deps_list = list()
  for (i in 1:uid.deps_list.list_count) {
    tmp_list = uid.deps_list.list[[i]][[4]]
    p6_p7.deps_list[[i]] = tmp_list
  }
  p6_p7.mdf = merge_deps_list(p6_p7.deps_list, meta_data.replicates.uids.unique)
  p6_p7.deps_df = check_real_deps_df(p6_p7.mdf, min_id=2)
  p6_p7.col_id = apply(p6_p7.mdf[,-1], 2, function(x){
    x = as.vector(unlist(x))
    length(which(!is.na(x)))
  })
  
  # P6_P8
  p6_p8.deps_list = list()
  for (i in 1:uid.deps_list.list_count) {
    tmp_list = uid.deps_list.list[[i]][[5]]
    p6_p8.deps_list[[i]] = tmp_list
  }
  p6_p8.mdf = merge_deps_list(p6_p8.deps_list, meta_data.replicates.uids.unique)
  p6_p8.deps_df = check_real_deps_df(p6_p8.mdf, min_id=2)
  p6_p8.col_id = apply(p6_p8.mdf[,-1], 2, function(x){
    x = as.vector(unlist(x))
    length(which(!is.na(x)))
  })
  
  # P7_P8
  p7_p8.deps_list = list()
  for (i in 1:uid.deps_list.list_count) {
    tmp_list = uid.deps_list.list[[i]][[6]]
    p7_p8.deps_list[[i]] = tmp_list
  }
  p7_p8.mdf = merge_deps_list(p7_p8.deps_list, meta_data.replicates.uids.unique)
  p7_p8.deps_df = check_real_deps_df(p7_p8.mdf, min_id=2)
  p7_p8.col_id = apply(p7_p8.mdf[,-1], 2, function(x){
    x = as.vector(unlist(x))
    length(which(!is.na(x)))
  })
}

# plot

deps_groups = c(
  'P5 VS P6', 'P5 VS P7', 'P5 VS P8', 
  'P6 VS P7', 'P6 VS P8', 'P7 VS P8'
)


#### 高置信差异蛋白 ####
write.csv(p5_p6.deps_df, '01_p5_p6.deps_df.csv', row.names = F)
write.csv(p5_p7.deps_df, '01_p5_p7.deps_df.csv', row.names = F)
write.csv(p5_p8.deps_df, '01_p5_p8.deps_df.csv', row.names = F)
write.csv(p6_p7.deps_df, '01_p6_p7.deps_df.csv', row.names = F)
write.csv(p6_p8.deps_df, '01_p6_p8.deps_df.csv', row.names = F)
write.csv(p7_p8.deps_df, '01_p7_p8.deps_df.csv', row.names = F)

union_deps = union(
  p5_p6.deps_df$GeneSymbol, union(
    p5_p7.deps_df$GeneSymbol, union(
      p5_p8.deps_df$GeneSymbol, union(
        p6_p7.deps_df$GeneSymbol, union(
          p6_p8.deps_df$GeneSymbol, p7_p8.deps_df$GeneSymbol
        )
      )
    )
  )
)
write.csv(union_deps, 'union_deps.csv', row.names = F)
intersect_deps = intersect(
  p5_p6.deps_df$GeneSymbol, intersect(
    p5_p7.deps_df$GeneSymbol, intersect(
      p5_p8.deps_df$GeneSymbol, intersect(
        p6_p7.deps_df$GeneSymbol, intersect(
          p6_p8.deps_df$GeneSymbol, p7_p8.deps_df$GeneSymbol
        )
      )
    )
  )
)


hc_deps_count = c(
  length(p5_p6.deps_df$GeneSymbol),
  length(p5_p7.deps_df$GeneSymbol),
  length(p5_p8.deps_df$GeneSymbol),
  length(p6_p7.deps_df$GeneSymbol),
  length(p6_p8.deps_df$GeneSymbol),
  length(p7_p8.deps_df$GeneSymbol),
  length(union_deps)
)
names(hc_deps_count) = c(deps_groups, 'Total')
par(
  mar = c(4,10,4,10)
)
bp = barplot(
  hc_deps_count,
  space = 0.5,
  ylab = 'DEPs',
  main = 'DEPs from pariwise comparison\n(Identification >= 2)',
  xaxt = 'n',
  ylim = c(0,870)
)
axis(1, bp, labels = c(deps_groups, 'Total'), cex.axis = 0.8)
text(bp, as.vector(hc_deps_count), labels = as.vector(hc_deps_count), pos = 3, cex = 0.8)

#### 高置信差异蛋白被检出的数目 ####
if(T){
  # P5 vs P6
  p5p6.hc_deps = table(p5_p6.deps_df$id)
  
  # par(mfrow=c(2,1))
  opar <- par(no.readonly = TRUE)
  # par(mai = c(2,5,1,2))
  par(
    fig = c(0,1,0.6,1), 
    mar = c(4,5,2,2)
  )
  bp = barplot(
    p5p6.hc_deps, 
    ylab = 'DEPs',
    xlab = 'Support evidence',
    col = '#EDEDED',
    yaxt = 'n',
    ylim = c(0,150),
    main = 'P5 VS P6'
  )
  axis(2,at=seq(0, 155, 25), labels = seq(0, 150, 25), cex.axis=0.8)
  text(
    bp, as.vector(p5p6.hc_deps),
    labels = as.vector(p5p6.hc_deps), 
    pos = 3, cex = 0.8, offset = 0.4
  )
  
  par(
    fig = c(0,1,0,0.65), 
    mar = c(5,5,2,2),
    new = TRUE
  )
  bp = barplot(
    p5_p6.col_id, 
    horiz = T, 
    xlim = c(0, 350),
    yaxt = 'n',
    col = '#EDEDED',
    xlab = 'DEPs',
    ylab = 'Location_MS_TestNo.'
    
  )
  axis(2, at = bp, labels = LETTERS[1:24], cex.axis=0.8)
  text(
    as.vector(p5_p6.col_id), bp, 
    labels = as.vector(p5_p6.col_id), 
    pos = 4, cex = 0.8, offset = 0.4
  )
}



if(T){
  # P5 vs P7
  p5p7.hc_deps = table(p5_p7.deps_df$id)
  
  # par(mfrow=c(2,1))
  opar <- par(no.readonly = TRUE)
  # par(mai = c(2,5,1,2))
  par(
    fig = c(0,1,0.6,1), 
    mar = c(4,5,2,2)
  )
  bp_ymax = 180
  bp = barplot(
    p5p7.hc_deps, 
    ylab = 'DEPs',
    xlab = 'Support evidence',
    col = '#EDEDED',
    yaxt = 'n',
    ylim = c(0,bp_ymax),
    main = 'P5 VS P7'
  )
  axis(2,at=seq(0, bp_ymax, 25), labels = seq(0, bp_ymax, 25), cex.axis=0.8)
  text(
    bp, as.vector(p5p7.hc_deps),
    labels = as.vector(p5p7.hc_deps), 
    pos = 3, cex = 0.8, offset = 0.4
  )
  
  par(
    fig = c(0,1,0,0.65), 
    mar = c(5,5,2,2),
    new = TRUE
  )
  bp = barplot(
    p5_p7.col_id, 
    horiz = T, 
    xlim = c(0, 350),
    yaxt = 'n',
    col = '#EDEDED',
    xlab = 'DEPs',
    ylab = 'Location_MS_TestNo.'
    
  )
  axis(2, at = bp, labels = LETTERS[1:24], cex.axis=0.8)
  text(
    as.vector(p5_p7.col_id), bp, 
    labels = as.vector(p5_p7.col_id), 
    pos = 4, cex = 0.8, offset = 0.4
  )
}





if(T){
  # P5 vs P8
  p5p8.hc_deps = table(p5_p8.deps_df$id)
  
  # par(mfrow=c(2,1))
  opar <- par(no.readonly = TRUE)
  # par(mai = c(2,5,1,2))
  par(
    fig = c(0,1,0.6,1), 
    mar = c(4,5,2,2)
  )
  bp_ymax = 220
  bp = barplot(
    p5p8.hc_deps, 
    ylab = 'DEPs',
    xlab = 'Support evidence',
    col = '#EDEDED',
    yaxt = 'n',
    ylim = c(0,bp_ymax),
    main = 'P5 VS P8'
  )
  axis(2,at=seq(0, bp_ymax, 25), labels = seq(0, bp_ymax, 25), cex.axis=0.8)
  text(
    bp, as.vector(p5p8.hc_deps),
    labels = as.vector(p5p8.hc_deps), 
    pos = 3, cex = 0.8, offset = 0.4
  )
  
  par(
    fig = c(0,1,0,0.65), 
    mar = c(5,5,2,2),
    new = TRUE
  )
  bp = barplot(
    p5_p8.col_id, 
    horiz = T, 
    xlim = c(0, 375),
    yaxt = 'n',
    col = '#EDEDED',
    xlab = 'DEPs',
    ylab = 'Location_MS_TestNo.'
    
  )
  axis(2, at = bp, labels = LETTERS[1:24], cex.axis=0.8)
  text(
    as.vector(p5_p8.col_id), bp, 
    labels = as.vector(p5_p8.col_id), 
    pos = 4, cex = 0.8, offset = 0.4
  )
}


if(T){
  # P6 vs P7
  p6p7.hc_deps = table(p6_p7.deps_df$id)
  
  # par(mfrow=c(2,1))
  opar <- par(no.readonly = TRUE)
  # par(mai = c(2,5,1,2))
  par(
    fig = c(0,1,0.6,1), 
    mar = c(4,5,2,2)
  )
  bp_ymax = 65
  bp = barplot(
    p6p7.hc_deps, 
    ylab = 'DEPs',
    xlab = 'Support evidence',
    col = '#EDEDED',
    yaxt = 'n',
    ylim = c(0,bp_ymax),
    main = 'P6 VS P7'
  )
  axis(2,at=seq(0, bp_ymax, 25), labels = seq(0, bp_ymax, 25), cex.axis=0.8)
  text(
    bp, as.vector(p6p7.hc_deps),
    labels = as.vector(p6p7.hc_deps), 
    pos = 3, cex = 0.8, offset = 0.4
  )
  
  par(
    fig = c(0,1,0,0.65), 
    mar = c(5,5,2,2),
    new = TRUE
  )
  bp = barplot(
    p6_p7.col_id, 
    horiz = T, 
    xlim = c(0, 350),
    yaxt = 'n',
    col = '#EDEDED',
    xlab = 'DEPs',
    ylab = 'Location_MS_TestNo.'
    
  )
  axis(2, at = bp, labels = LETTERS[1:24], cex.axis=0.8)
  text(
    as.vector(p6_p7.col_id), bp, 
    labels = as.vector(p6_p7.col_id), 
    pos = 4, cex = 0.8, offset = 0.4
  )
}




if(T){
  # P6 vs P8
  p6p8.hc_deps = table(p6_p8.deps_df$id)
  
  # par(mfrow=c(2,1))
  opar <- par(no.readonly = TRUE)
  # par(mai = c(2,5,1,2))
  par(
    fig = c(0,1,0.6,1), 
    mar = c(4,5,2,2)
  )
  bp_ymax = 125
  bp = barplot(
    p6p8.hc_deps, 
    ylab = 'DEPs',
    xlab = 'Support evidence',
    col = '#EDEDED',
    yaxt = 'n',
    ylim = c(0,bp_ymax),
    main = 'P6 VS P8'
  )
  axis(2,at=seq(0, bp_ymax, 25), labels = seq(0, bp_ymax, 25), cex.axis=0.8)
  text(
    bp, as.vector(p6p8.hc_deps),
    labels = as.vector(p6p8.hc_deps), 
    pos = 3, cex = 0.8, offset = 0.4
  )
  
  par(
    fig = c(0,1,0,0.65), 
    mar = c(5,5,2,2),
    new = TRUE
  )
  bp = barplot(
    p6_p8.col_id, 
    horiz = T, 
    xlim = c(0, 350),
    yaxt = 'n',
    col = '#EDEDED',
    xlab = 'DEPs',
    ylab = 'Location_MS_TestNo.'
    
  )
  axis(2, at = bp, labels = LETTERS[1:24], cex.axis=0.8)
  text(
    as.vector(p6_p8.col_id), bp, 
    labels = as.vector(p6_p8.col_id), 
    pos = 4, cex = 0.8, offset = 0.4
  )
}



if(T){
  # P7 vs P8
  p7p8.hc_deps = table(p7_p8.deps_df$id)
  
  # par(mfrow=c(2,1))
  opar <- par(no.readonly = TRUE)
  # par(mai = c(2,5,1,2))
  par(
    fig = c(0,1,0.6,1), 
    mar = c(4,5,2,2)
  )
  bp_ymax = 155
  bp = barplot(
    p7p8.hc_deps, 
    ylab = 'DEPs',
    xlab = 'Support evidence',
    col = '#EDEDED',
    yaxt = 'n',
    ylim = c(0,bp_ymax),
    main = 'P7 VS P8'
  )
  axis(2,at=seq(0, bp_ymax, 25), labels = seq(0, bp_ymax, 25), cex.axis=0.8)
  text(
    bp, as.vector(p7p8.hc_deps),
    labels = as.vector(p7p8.hc_deps), 
    pos = 3, cex = 0.8, offset = 0.4
  )
  
  par(
    fig = c(0,1,0,0.65), 
    mar = c(5,5,2,2),
    new = TRUE
  )
  bp = barplot(
    p7_p8.col_id, 
    horiz = T, 
    xlim = c(0, 350),
    yaxt = 'n',
    col = '#EDEDED',
    xlab = 'DEPs',
    ylab = 'Location_MS_TestNo.'
    
  )
  axis(2, at = bp, labels = LETTERS[1:24], cex.axis=0.8)
  text(
    as.vector(p7_p8.col_id), bp, 
    labels = as.vector(p7_p8.col_id), 
    pos = 4, cex = 0.8, offset = 0.4
  )
}





#### 不同中心不同仪器可以被检测出来的差异蛋白 ####
if(T){
  par(mfrow = c(3,2))
  p5p6.h = hist(
    p5_p6.col_id,
    breaks = 25,
    xlim = c(0, 350),
    col = '#CCCCCC',
    xlab = 'DEPs', 
    main = "P5 VS P6"
  )
  p5p7.h = hist(
    p5_p7.col_id,
    breaks = 25,
    xlim = c(0, 350),
    col = '#CCCCCC',
    xlab = 'DEPs', 
    main = "P5 VS P7"
  )
  p5p8.h = hist(
    p5_p8.col_id,
    breaks = 25,
    xlim = c(0, 350),
    col = '#CCCCCC',
    xlab = 'DEPs', 
    main = "P5 VS P8"
  )
  p6p7.h = hist(
    p6_p7.col_id,
    breaks = 25,
    xlim = c(0, 350),
    col = '#CCCCCC',
    xlab = 'DEPs', 
    main = "P6 VS P7"
  )
  p6p8.h = hist(
    p6_p8.col_id,
    breaks = 25,
    xlim = c(0, 350),
    col = '#CCCCCC',
    xlab = 'DEPs', 
    main = "P6 VS P8"
  )
  p7p8.h = hist(
    p7_p8.col_id,
    breaks = 25,
    xlim = c(0, 350),
    col = '#CCCCCC',
    xlab = 'DEPs', 
    main = "P7 VS P8"
  )
  
}
















