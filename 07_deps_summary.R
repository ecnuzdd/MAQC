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
for (i in 1:6) {
  deps_up_df = deps_up_list[[i]]
  deps_dn_df = deps_dn_list[[i]]
  deps_up = get_deps_set_by_cx(deps_up_df, min_id = 8)
  deps_down = get_deps_set_by_cx(deps_dn_df, min_id = 8)
  deps = setdiff(union(deps_up, deps_down), intersect(deps_up, deps_down))
  deps_list[[i]] = deps
}
deps_list_vector = unique(unlist(deps_list))


scale_df_list = list()
metadata_list = list()
metadata = NULL
for (i in 1:24) {
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
for (i in 2:24) {
  merge_df = merge(merge_df, scale_df_list[[i]], by = 'GeneSymbol', all = T)
}
merge_df.non_na = na.omit(merge_df)


library(pheatmap)
ph_mat = as.matrix(merge_df.non_na[,-1])
ph_metadata = metadata

order_index = order(ph_metadata$UPerson)
ph_metadata = ph_metadata[order_index,]
ph_mat = ph_mat[,order_index]

ph_metadata.uids = paste(
  ph_metadata$Labcode, 
  ph_metadata$SInstrument, 
  ph_metadata$Batch, 
  sep = '_'
)
ph_uids = sort(unique(ph_metadata.uids))

# i = 2
# ph_mat = as.matrix(scale_df_list[[i]][,-1])
# ph_metadata = metadata_list[[i]]
annotation_col = data.frame(
  Company = ph_uids,
  Sample = ph_metadata$UPerson
)
rownames(annotation_col) = colnames(ph_mat)
annotation_colors = list(
  Sample = c(D5=p_colors[1], D6=p_colors[2], F7=p_colors[3], M8=p_colors[4])
)
pheatmap(
  ph_mat, 
  scale = 'none',
  annotation_col = annotation_col,
  annotation_colors = annotation_colors,
  clustering_method = 'ward.D2',
  cluster_cols = F,
  show_rownames = F,
  show_colnames = F
)




library(pheatmap)
ph_metadata = metadata
ph_df = get_fot_df_by_expnames(fot_data, metadata$Exp_code)
ph_mat = as.matrix(ph_df[,-1])
annotation_col = data.frame(
  Sample = ph_metadata$UPerson
)
rownames(annotation_col) = colnames(ph_mat)
annotation_colors = list(
  Sample = c(D5=p_colors[1], D6=p_colors[2], F7=p_colors[3], M8=p_colors[4])
)
ph_mat.row_id = apply(ph_mat, 1, function(x){length(which(x>0))})
ph_mat1 = ph_mat[ph_mat.row_id>216,]

pheatmap(
  ph_mat1, 
  scale = 'row',
  annotation_col = annotation_col,
  annotation_colors = annotation_colors,
  show_rownames = F,
  show_colnames = F
)


get_deps_set_by_cx <- function(df, min_id = 8){
  
  # df = deps_up_df
  df.row_id = apply(df[,-1], 1, function(x){
    length(which(x>0))
  })
  df.index = which(df.row_id>=min_id)
  df.genes = as.vector(df$GeneSymbol[df.index])
  return(df.genes)
}


#### global ####
if(T){
  scale_df_list = list()
  metadata_list = list()
  metadata = NULL
  for (i in 1:24) {
    uid = uids[i] 
    uid.metadata = meta_data.cp[meta_data.cp.uids==uid,]
    uid.metadata = uid.metadata[order(uid.metadata$Person),]
    metadata = rbind.data.frame(metadata, uid.metadata)
    uid.expnames = uid.metadata$Exp_code
    uid.df = get_fot_df_by_expnames(fot_data, uid.expnames)
    # uid.df_index = na.omit(match(deps_list_vector, uid.df$GeneSymbol))
    # uid.df = uid.df[uid.df_index,]
    uid.df.genesymbol = as.vector(uid.df$GeneSymbol)
    uid.df.mat = as.vector(uid.df[,-1])
    uid.df.mat.scale = t(scale(t(uid.df.mat)))
    uid.df.df.scale = data.frame(GeneSymbol = uid.df.genesymbol, uid.df.mat.scale)
    scale_df_list[[i]] = uid.df.df.scale
    metadata_list[[i]] = uid.metadata
    
  }
  
  merge_df = scale_df_list[[1]]
  for (i in 2:24) {
    merge_df = merge(merge_df, scale_df_list[[i]], by = 'GeneSymbol', all = T)
  }
  merge_df.non_na = na.omit(merge_df)
  
  
  library(pheatmap)
  ph_mat = as.matrix(merge_df.non_na[,-1])
  ph_metadata = metadata
  
  order_index = order(ph_metadata$UPerson)
  ph_metadata = ph_metadata[order_index,]
  ph_mat = ph_mat[,order_index]
  
  ph_metadata.uids = paste(
    ph_metadata$Labcode, 
    ph_metadata$SInstrument, 
    ph_metadata$Batch, 
    sep = '_'
  )
  ph_uids = sort(unique(ph_metadata.uids))
  
  annotation_col = data.frame(
    Company = ph_uids,
    Sample = ph_metadata$UPerson
  )
  rownames(annotation_col) = colnames(ph_mat)
  annotation_colors = list(
    Sample = c(D5=p_colors[1], D6=p_colors[2], F7=p_colors[3], M8=p_colors[4])
  )
  pheatmap(
    ph_mat, 
    scale = 'none',
    annotation_col = annotation_col,
    annotation_colors = annotation_colors,
    clustering_method = 'ward.D2',
    cluster_cols = F,
    show_rownames = F,
    show_colnames = F
  )
  
}


