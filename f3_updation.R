get_deps_set_by_cx <- function(df, min_id = 8){
  
  # df = deps_up_df
  df.row_id = apply(df[,-1], 1, function(x){
    length(which(x>0))
  })
  df.index = which(df.row_id>=min_id)
  df.genes = as.vector(df$GeneSymbol[df.index])
  return(df.genes)
}

merge_deps_set <- function(df, step=3){
  df.row = nrow(df)
  df.grp_count = df.row/3
  df.grps = rep(seq(1,df.grp_count), each = step)
  df.grps = factor(df.grps, levels = seq(1,df.grp_count))
  df.sum_id = tapply(df$Count, df.grps, sum)
  df.names = seq(1, df.row, step)
  df.names = apply(data.frame(df.names), 1, function(x){
    x_str = paste('[', x, ', ', x+2, ']', sep = '')
    return(x_str)
  })
  df = data.frame(
    ID = df.names,
    Count = df.sum_id
  )
  return(df)
  
}

#############################################
meta_data.cp = meta_data[which(meta_data$Purpose == 'Replications'),]
meta_data.cp.uids = paste(
  meta_data.cp$Labcode, 
  meta_data.cp$SInstrument, 
  meta_data.cp$Batch, 
  sep = '_'
)
uids = unique(meta_data.cp.uids)[-c(25, 26)]
uids = sort(uids)
px_list = c('D5', 'D6', 'F7', 'M8')
px = px_list[1]

load("./deps_0815/deps_list_up.RData")
deps_up_list = cx.df.list

load("./deps_0815/deps_list_down.RData")
deps_dn_list = cx.df.list

#### D5 vs D6 ####
deps_list = list()
min_id = 1
for (i in 1:6) {
  deps_up_df = deps_up_list[[i]]
  deps_dn_df = deps_dn_list[[i]]
  deps_up = get_deps_set_by_cx(deps_up_df, min_id = min_id)
  deps_down = get_deps_set_by_cx(deps_dn_df, min_id = min_id)
  deps = setdiff(union(deps_up, deps_down), intersect(deps_up, deps_down))
  deps_list[[i]] = deps
}
deps_list_vector = unique(unlist(deps_list))



scale_df_list = list()
metadata_list = list()
metadata = NULL
for (i in 1:length(uids)) {
  uid = uids[i] 
  uid.metadata = meta_data.cp[meta_data.cp.uids==uid,]
  uid.metadata = uid.metadata[order(uid.metadata$Person),]
  metadata = rbind.data.frame(metadata, uid.metadata)
  uid.expnames = uid.metadata$Exp_code
  uid.df = get_fot_df_by_expnames(fot_data, uid.expnames)
  uid.df_index = na.omit(match(deps_list_vector, uid.df$GeneSymbol))
  uid.df = uid.df[uid.df_index,]
  uid.df.genesymbol = as.vector(uid.df$GeneSymbol)
  uid.df.mat = as.vector(uid.df[,-1])
  uid.df.mat.scale = t(scale(t(uid.df.mat)))
  uid.df.df.scale = data.frame(GeneSymbol = uid.df.genesymbol, uid.df.mat.scale)
  scale_df_list[[i]] = uid.df.df.scale
  metadata_list[[i]] = uid.metadata
  
}

merge_df = scale_df_list[[1]]
for (i in 2:length(uids)) {
  merge_df = merge(merge_df, scale_df_list[[i]], by = 'GeneSymbol', all = T)
}

row_na_index = apply(merge_df[,-1], 1, function(x){
  if(length(which(is.na(x)))>0){
    return(T)
  }else{
    return(F)
  }
})

merge_df.non_na = merge_df[-which(row_na_index),]
merge_df.na = merge_df[which(row_na_index),]

ph_metadata = metadata
# order_index = order(ph_metadata$UPerson)
# ph_metadata = ph_metadata[order_index,]

# 蛋白至少在24个中心都被鉴定到一次1/12
#### merge_df.non_na ####
library(pheatmap)
ph_mat = as.matrix(merge_df.non_na[,-1])
# ph_mat = ph_mat[,order_index]
ph_mat = na.omit(ph_mat)
annotation_col = data.frame(
  Sample = ph_metadata$UPerson
)
rownames(annotation_col) = colnames(ph_mat)
annotation_colors = list(
  Sample = c(D5=p_colors[1], D6=p_colors[2], F7=p_colors[3], M8=p_colors[4])
)
main = paste('merge_df.non_na', nrow(ph_mat))
ph_non_na = pheatmap(
  ph_mat, 
  scale = 'none',
  main = main,
  annotation_col = annotation_col,
  annotation_colors = annotation_colors,
  clustering_method = 'ward.D2',
  cluster_cols = T,
  cluster_rows = T,
  show_rownames = F,
  show_colnames = F
)
ph_non_na.row_order = ph_non_na$tree_row$order
ph_non_na_1 = pheatmap(
  ph_mat[ph_non_na.row_order,], 
  scale = 'none',
  main = main,
  annotation_col = annotation_col,
  annotation_colors = annotation_colors,
  clustering_method = 'ward.D2',
  cluster_cols = F,
  cluster_rows = F,
  show_rownames = F,
  show_colnames = F
) 

ph_mat1 <- ph_mat[ph_non_na.row_order,]


#### merge_df.na ####
library(pheatmap)
ph_mat = as.matrix(merge_df.na[,-1])
ph_mat = ph_mat[,order_index]
ph_mat2 = ph_mat

na_index = which(is.na(ph_mat))
ph_mat[na_index] = 0
annotation_col = data.frame(
  Sample = ph_metadata$UPerson
)
rownames(annotation_col) = colnames(ph_mat)
annotation_colors = list(
  Sample = c(D5=p_colors[1], D6=p_colors[2], F7=p_colors[3], M8=p_colors[4])
)
main = 'merge_df.na'
ph_na_1 = pheatmap(
  ph_mat, 
  scale = 'none',
  main = main,
  annotation_col = annotation_col,
  annotation_colors = annotation_colors,
  clustering_method = 'ward.D2',
  cluster_cols = F,
  cluster_rows = T,
  show_rownames = F,
  show_colnames = F
) 

#### merge_df ####
library(pheatmap)

ph_mat3 = ph_mat2[ph_na_1$tree_row$order,]

ph_mat = rbind(
  ph_mat1, ph_mat3
)
annotation_col = data.frame(
  Sample = ph_metadata$UPerson
)
rownames(annotation_col) = colnames(ph_mat)
annotation_colors = list(
  Sample = c(D5=p_colors[1], D6=p_colors[2], F7=p_colors[3], M8=p_colors[4])
)
main = 'merge_df'
ph_na_1 = pheatmap(
  ph_mat, 
  scale = 'none',
  main = main,
  annotation_col = annotation_col,
  annotation_colors = annotation_colors,
  clustering_method = 'ward.D2',
  cluster_cols = F,
  cluster_rows = F,
  show_rownames = F,
  show_colnames = F
) 



####################################################################################################
#### Reproducibility ####
#### up ####
deps_up_list_names = names(deps_up_list)
df_up_set = NULL
for (i in 1:6) {
  deps_up_df = deps_up_list[[i]]
  deps_up_df.row_id = apply(deps_up_df[,-1], 1, function(x){length(which(x>0))})
  df = data.frame(table(deps_up_df.row_id))
  colnames(df) = c('ID', 'Count')
  if(deps_up_list_names[i]=='D6 VS F7'){
    df$ID = as.vector(df$ID)
    df = rbind(df, c(21, 0))
    df = df[order(as.numeric(df$ID)),]
  }
  df = merge_deps_set(df, step=3)
  df$Condition = deps_up_list_names[i]
  # df$ID = factor(df$ID, levels = seq(1,24))
  lev = unique(as.vector(df$ID))
  lev = lev[length(lev):1]
  
  df$ID = factor(df$ID, levels = lev)
  df$Count = as.numeric(df$Count)
  df_up_set = rbind(df_up_set, df)
}
tapply(df_up_set$Count, (df_up_set$Condition), sum)

# Stacked + percent
df1 = df_up_set
df1$Code = rep(seq(1,8), 6)
if(T){
  p <- ggplot(df1, aes(fill = Code, y = Count, x = Condition))
  # p <- p + coord_flip() 
  # p <- p + scale_color_gradient(low="blue", high="red")
  p <- p + geom_bar(position = "stack", stat = "identity", width = 0.6) 
  p <- p + scale_y_continuous(limits = c(0, 6000), breaks = seq(0,6000,1000))
  print(p + theme_classic2())
}

#### down ####
deps_dn_list_names = names(deps_dn_list)
df_dn_set = NULL
for (i in 1:6) {
  deps_up_df = deps_dn_list[[i]]
  deps_up_df.row_id = apply(deps_up_df[,-1], 1, function(x){length(which(x>0))})
  df = data.frame(table(deps_up_df.row_id))
  colnames(df) = c('ID', 'Count')
  df = merge_deps_set(df, step=3)
  df$Condition = deps_dn_list_names[i]
  lev = unique(as.vector(df$ID))
  lev = lev[length(lev):1]
  df$ID = factor(df$ID, levels = lev)
  df$Count = as.numeric(df$Count)
  df_dn_set = rbind(df_dn_set, df)
}
tapply(df_dn_set$Count, (df_dn_set$Condition), sum)

# Stacked + percent
df2 = df_dn_set
df2$Code = rep(seq(1,8), 6)
if(T){
  p <- ggplot(df2, aes(fill = Code, y = Count, x = Condition))
  # p <- p + coord_flip() 
  # p <- p + scale_color_gradient(low="blue", high="red")
  p <- p + geom_bar(position = "stack", stat = "identity", width = 0.6) 
  p <- p + scale_y_continuous(limits = c(0, 6000), breaks = seq(0,6000,1000))
  print(p + theme_classic2())
}

df2_1 = df2
df3 = rbind.data.frame(
  df1, df2
)
df3$Reg = rep(c('UP', 'DN'), each = nrow(df3)/2)
df3$Count[df3$Reg=='DN'] = df3$Count[df3$Reg=='DN']*(-1)
df3$Code[df3$Reg=='DN'] = df3$Code[df3$Reg=='DN']*(-1)

if(T){
  label_str =  c('[1, 3]', '[4, 6]', '[7, 9]', '[10, 12]', '[13, 15]', '[16, 18]', '[19, 21]', '[22, 24]')
  p <- ggplot(df3, aes(fill = Code, y = Count, x = Condition))
  p <- p + coord_flip() 
  p <- p + scale_fill_gradient2(
    low = "blue", mid = 'white', high = "red", midpoint = 0,
    # color bar setting
    limits = c(-8, 8), 
    breaks = setdiff(seq(-8, 8, 1), 0), 
    labels = c(label_str[8:1], label_str), 
    # color bar setting
    
    name = 'Reproducibility\nFrequency',
    guide = F
  )
  p <- p + geom_bar(position = "stack", stat = "identity", width = 0.6) 
  p <- p + scale_y_continuous(
    limits = c(-6000, 6000), 
    breaks = seq(-6000,6000,1000), 
    labels = c(seq(6000,0,-1000), seq(1000,6000,1000))
  )
  
  p <- p + geom_hline(yintercept=0, color = "white", size=1)
  
  p <- p + guides(fill = guide_colorbar(
    title = "Reproducibility\nFrequency",
    barwidth = 1,
    barheight = 15
  )) 
  p <- p + theme_classic2()
  print(p)
}

write.xlsx(
  df3,
  normalizePath(file.path(
    BASE_DIR, 'deps_0815', 'Rplot_repro_frequency.xlsx'
  ), mustWork = F),
  row.names = F
)

df_dn_set.mat = matrix(
  df_dn_set$Count, 8, 6, byrow = F, 
  dimnames = list(unique(df_dn_set$ID), unique(df_dn_set$Condition))
)
write.xlsx(
  df_dn_set.mat,
  normalizePath(file.path(
    BASE_DIR, 'deps_0815', 'Rplot_repro_frequency.xlsx'
  ), mustWork = F),
  row.names = T,
  sheetName = 'DN', append = T
)

df_up_set.mat = matrix(
  df_up_set$Count, 8, 6, byrow = F, 
  dimnames = list(unique(df_up_set$ID), unique(df_up_set$Condition))
)
write.xlsx(
  df_up_set.mat,
  normalizePath(file.path(
    BASE_DIR, 'deps_0815', 'Rplot_repro_frequency.xlsx'
  ), mustWork = F),
  row.names = T,
  sheetName = 'UP', append = T
)

stable.up <- data.frame(
  IDF = rownames(df_up_set.mat),
  df_up_set.mat
)
stable.up.p <- ggtexttable(
  stable.up, rows = NULL, 
  theme = ttheme("mRed")
)
print(stable.up.p)

stable.dn <- data.frame(
  IDF = rownames(df_dn_set.mat),
  df_dn_set.mat
)
stable.dn.p <- ggtexttable(
  stable.dn, rows = NULL, 
  theme = ttheme("mBlue")
)
print(stable.dn.p)

library("gridExtra")
grid.arrange(
  p,                                                      # bar plot spaning two columns
  stable.up.p, stable.dn.p,                               # box plot and scatter plot
  ncol = 2, nrow = 2, 
  layout_matrix = rbind(c(1,1), c(2,3))
)


ggarrange(
  p,                                                 # First row with scatter plot
  ggarrange(stable.dn.p, stable.up.p, ncol = 2, labels = c("B", "C")), # Second row with box and dot plots
  nrow = 2, 
  labels = "A"                                        # Labels of the scatter plot
) 
####################################################################################################




####################################################################################################



























