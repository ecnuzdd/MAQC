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


i = 1
deps.list = list()
for (i in 1:length(uids)) {
  uid = uids[i]
  uid.metadata = meta_data.cp[meta_data.cp.uids==uid,]
  uid.fot = get_fot_df_by_expnames(fot_data, uid.metadata$Exp_code)
  flags = NULL
  flags.count = 0
  uid.deps.list = list()
  for (j in 1:3) {
    ctrl.flag = px_list[j]
    ctrl.index = which(uid.metadata$UPerson==ctrl.flag)
    k.start = j + 1
    for (k in k.start:4) {
      case.flag = px_list[k]
      case.index = which(uid.metadata$UPerson==case.flag)
      flag.jk = paste(ctrl.flag, ' VS ', case.flag, sep = '')
      print(flag.jk)
      flags = c(flags, flag.jk)
      uid.metadata.subset = uid.metadata[c(ctrl.index, case.index),]
      uid.fot.subset = get_fot_df_by_expnames(fot_data, uid.metadata.subset$Exp_code)
      uid.group.subset = as.vector(uid.metadata.subset$UPerson)
      flag.deps = get_deps_results(uid.fot.subset, min_fc = 2)
      flags.count = flags.count + 1
      uid.deps.list[[flags.count]] = flag.deps
    }
  }
  names(uid.deps.list) = flags
  deps.list[[i]] = uid.deps.list
}



#### D5 VS D6 ####

cx.df.list = list()
for (i in 1:6) {
  cx.list = list()
  for (j in 1:length(uids)){
    cx.list[[j]] = deps.list[[j]][[i]]
  }
  cx.df = merge_list2df(cx.list)
  colnames(cx.df) = c('GeneSymbol', uids)
  cx.df.list[[i]] = cx.df
}
names(cx.df.list) = flags
save(cx.df.list, file = "./deps_0815/deps_list_up.RData")




merge_list2df <- function(lst){
  lst = cx.list
  merge_df = lst[[1]]
  for (i in 2:length(lst)) {
    merge_df = merge(merge_df, lst[[i]], by = 'GeneSymbol', all = T)
  }
  merge_df = na2zero(merge_df)
  return(merge_df)
}




get_deps_results <- function(df, min_fc = 2){
  # df = uid.fot.subset
  # groups = uid.group.subset
  genes = as.vector(df$GeneSymbol)
  mat = as.matrix(df[,-1])
  mat.row = nrow(mat)
  deps.df = NULL
  for (i in 1:mat.row) {
    x1 = mat[i,1:3]
    x2 = mat[i,4:6]
    x1.id = length(x1)
    if(x1.id>=2){
      fc = mean(x1)/mean(x2)
      if(fc>min_fc){
        gene = genes[i]
        deps.i = c(gene, fc)
        deps.df = rbind(deps.df, deps.i)
      }
    }
  }
  colnames(deps.df) = c('GeneSymbol', 'FC')
  return(deps.df)
}















