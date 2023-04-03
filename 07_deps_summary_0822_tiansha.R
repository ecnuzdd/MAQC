meta_data.cp = meta_data[1:440,]
meta_data.cp.uids = paste(
  meta_data.cp$Labcode, 
  meta_data.cp$SInstrument, 
  meta_data.cp$Batch, 
  sep = '_'
)
uids = unique(meta_data.cp.uids)
uids = sort(uids)
px_list = c('D5', 'D6', 'F7', 'M8')
px = px_list[1]


expnames = meta_data.cp$Exp_code
grps = meta_data.cp$UPerson
fot = get_fot_df_by_expnames(fot_data, expnames)
fot.row = nrow(fot)
deps_info = NULL
for (i in 1:fot.row) {
  g = as.vector(fot$GeneSymbol[i])
  x = as.vector(unlist(fot[i,-1]))
  x.mean = tapply(x, grps, mean)
  x.mean.max = max(x.mean)
  x.mean.fc = x.mean.max/x.mean
  fc_count = length(which(x.mean.fc>2))
  # if(fc_count==3){
  #   break
  # }
  
  
  if(fc_count==3){
    max_index = which(x.mean==x.mean.max)
    case_flag = names(max_index)
    pvalue = NULL
    case_index = which(grps==case_flag)
    for (j in 1:4) {
      ctl_flag = px_list[j]
      if(ctl_flag != case_flag){
        ctrl_index = which(grps==ctl_flag)
        wt = wilcox.test(x[case_index], x[ctrl_index])
        pvalue = c(pvalue, wt$p.value)
      }
    }
    pvalue_count = length(which(pvalue<0.05))
    # if(pvalue_count==3){
    if(T){
      max_v = c(g, case_flag)
      deps_info = rbind(deps_info, max_v)
    }
    
    
  }
  print(i)
}

deps_info = as.data.frame(deps_info)
deps_df = fot[match(deps_info$V1, fot$GeneSymbol),]


order_index = order(meta_data.cp$UPerson)
ph_metadata = meta_data.cp[order_index,]

library(pheatmap)
ph_mat = as.matrix(deps_df[,-1])
ph_mat = ph_mat[,order_index]
annotation_col = data.frame(
  Sample = ph_metadata$UPerson
)
rownames(annotation_col) = colnames(ph_mat)
annotation_colors = list(
  Sample = c(D5=p_colors[1], D6=p_colors[2], F7=p_colors[3], M8=p_colors[4])
)
main = 'DEPs'
ph_non_na = pheatmap(
  ph_mat, 
  scale = 'row',
  main = main,
  annotation_col = annotation_col,
  annotation_colors = annotation_colors,
  clustering_method = 'ward.D2',
  cluster_cols = F,
  cluster_rows = T,
  show_rownames = F,
  show_colnames = F
)












