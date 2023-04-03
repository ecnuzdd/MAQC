library(readr)
library(xlsx)
library(readxl)
library(stringr)
library(scales)
library(ggplot2)

#### get_fot_df_by_expnames ####
# expnames_index = which(
#   meta_data$Affiliation == 'NCPSB'
#   & meta_data$Purpose == 'Replications'
#   & meta_data$SInstrument == 'Lumos'
#   & meta_data$Batch == 1
# )
# expnames = meta_data$Exp_code[expnames_index]
get_fot_df_by_expnames <- function(fot_df, expnames){
  GeneSymbol = as.vector(fot_df$GeneSymbol)
  fot_mat = as.matrix(fot_df[,-1])
  fot_mat_colnames = colnames(fot_mat)
  expnames_match_index = match(expnames, fot_mat_colnames)
  if(length(which(is.na(expnames_match_index)))>0){
    print('NAs')
    return(0)
  }else{
    fot_mat_subset = fot_mat[,expnames_match_index]
    fot_mat_subset_rowsums = rowSums(fot_mat_subset)
    kept_index = which(fot_mat_subset_rowsums>0)
    return_df = data.frame(
      GeneSymbol, fot_mat_subset
    )[kept_index,]
    return(return_df)
  }
}



get_fot_df_by_expnames_with_na <- function(fot_df, expnames){
  GeneSymbol = as.vector(fot_df$GeneSymbol)
  fot_mat = as.matrix(fot_df[,-1])
  fot_mat_colnames = colnames(fot_mat)
  expnames_match_index = match(expnames, fot_mat_colnames)
  if(length(which(is.na(expnames_match_index)))>0){
    print('NAs')
    return(0)
  }else{
    fot_mat_subset = fot_mat[,expnames_match_index]
    fot_mat_subset_rowid = apply(fot_mat_subset, 1, function(x){
      # length(which(!is.na(x)))
      length(which(x != 'NA'))
    })
    kept_index = which(fot_mat_subset_rowid>0)
    return_df = data.frame(
      GeneSymbol, fot_mat_subset
    )[kept_index,]
    return(return_df)
  }
}

#### impute_df_with_value ####
# a_df = get_fot_df_by_expnames(fot_data, expnames)
impute_df_with_value <- function(fot_df, value = 0){
  GeneSymbol = as.vector(fot_df$GeneSymbol)
  fot_mat = as.matrix(fot_df[,-1])
  fot_mat.zero_index = which(fot_mat==0)
  if(value == 0){
    value = min(fot_mat[which(fot_mat>0)])/10
    fot_mat[fot_mat.zero_index] = value
  }
  return_df = data.frame(
    GeneSymbol, fot_mat
  )
  return(return_df)
}



#### get_deps_df ####
get_deps_df <- function(df, groups, minFC=2, minPvalue = 0.05){
  # df = uid.fot_data.imputed
  df.genesymbol = as.vector(df$GeneSymbol)
  df.mat = as.matrix(df[,-1])
  df.mat.min = min(df.mat)
  rownames(df.mat) = df.genesymbol
  # groups = uid.meta_data$Person
  groups.unique = sort(unique(groups))
  groups.unique_count = length(groups.unique)
  
  
  df.mat.row = nrow(df.mat)
  deps_list = list()
  deps_list_index = 0
  j_end = groups.unique_count-1
  for (j in 1:j_end) {
    x1.index = which(groups==groups.unique[j])
    k_start = j+1
    for (k in k_start:groups.unique_count) {
      flag = paste(groups.unique[j], groups.unique[k], sep = ' VS ')
      x2.index = which(groups==groups.unique[k])
      jk.deps.df = NULL
      jk.deps = NULL
      for (i in 1:df.mat.row) {
        gene = df.genesymbol[i]
        x = df.mat[i,]
        
        x1 = as.vector(x[x1.index])
        x1.id = length(which(x1 > df.mat.min))
        x2 = as.vector(x[x2.index])
        x2.id = length(which(x2 > df.mat.min))
        if(x1.id>=2 | x2.id>=2){
          log2FC = log2(mean(x1)) - log2(mean(x2))
          ttpvalue = t.test(log2(x1), log2(x2))$p.value
          jk.deps = c(jk.deps, gene)
          jk.deps.df = rbind.data.frame(jk.deps.df, c(log2FC, ttpvalue))
          
        }
        
      }
      deps.df = data.frame(
        jk.deps, jk.deps.df
      )
      deps.df$FDR = p.adjust(as.vector(unlist(deps.df[,3])), method = 'BH')
      colnames(deps.df) = c(
        'GeneSymbol',
        'log2FC',
        'Pvalue',
        'FDR'
      )
      deps_index = which(
        abs(deps.df$log2FC)>log2(minFC) & deps.df$FDR<0.1
      )
      deps.df.real = deps.df[deps_index,]
      colnames(deps.df.real) = c(
        'GeneSymbol',
        paste(groups.unique[j], groups.unique[k], 'log2FC', sep = '_'),
        paste(groups.unique[j], groups.unique[k], 'Pvalue', sep = '_'),
        paste(groups.unique[j], groups.unique[k], 'FDR', sep = '_')
      )
      deps_list_index = deps_list_index + 1
      # print(deps_list_index)
      deps_list[[deps_list_index]] = deps.df.real
      
    }
  }
  return(deps_list)
}



#### merge_deps_list ####
merge_deps_list <- function(deps_list, list_names){
  # deps_list = p5_p6.deps_list
  # list_names = meta_data.replicates.uids.unique
  deps_list_count = length(deps_list)
  merge_df = deps_list[[1]][,c(1,2)]
  for (i in 2:deps_list_count) {
    tmp_df = deps_list[[i]][,c(1,2)]
    merge_df = merge(merge_df, tmp_df, by = 'GeneSymbol', all = T)
  }
  colnames(merge_df) = c('GeneSymbol', list_names)
  return(merge_df)
}

#### check_real_deps ####
check_real_deps_df <- function(deps_df, min_id=2){
  # deps_df = p5_p6.mdf
  deps_df.mat = as.matrix(deps_df[,-1])
  deps_df.row.gt0.count = apply(deps_df.mat, 1, function(x){
    x = as.vector(unlist(x))
    length(which(x>0))
  })
  deps_df.row.lt0.count = apply(deps_df.mat, 1, function(x){
    x = as.vector(unlist(x))
    length(which(x<0))
  })
  deps_df.row.id = deps_df.row.gt0.count + deps_df.row.lt0.count
  df = data.frame(
    GeneSymbol = deps_df$GeneSymbol,
    id = deps_df.row.id,
    id_gt0 = deps_df.row.gt0.count,
    id_lt0 = deps_df.row.lt0.count
  )
  fake_deps_index = which(
    df$id_gt0!=0 & df$id_lt0!=0
  )
  if(length(fake_deps_index)>0){
    df = df[-fake_deps_index,]
  }
  ret_df = df[,c(1,2)]
  ret_df = ret_df[ret_df$id>=min_id,]
  return(ret_df)
}


#### na2zero ####
na2zero <- function(df){
  # df = fot_data
  GeneSymbol = df$GeneSymbol
  mat = as.matrix(df[,-1])
  na_index = which(is.na(mat) | mat == 'NA')
  if(length(na_index)>0){
    mat[na_index] = 0
  }
  return_df = data.frame(
    GeneSymbol, mat
  )
  return(return_df)
}


#### get_cumulative_curve_vector ####
get_cumulative_curve_vector <- function(df){
  # df = p5.df
  mat = df[,-1]
  mat.col_id = apply(mat, 2, function(x){
    length(which(x>0))
  })
  mat.col_id.asec.order = order(mat.col_id)
  mat = mat[,mat.col_id.asec.order]
  mat.col = ncol(mat)
  id_vector = NULL
  for (i in 1:mat.col) {
    if(F){
      x.id.all = NULL
      for (j in 1:50) {
        set.seed(2)
        sample_index = sample(mat.col, i)
        x = as.matrix(mat[,sample_index])
        x.rowsums = rowSums(x)
        x.id = length(which(x.rowsums>0))
        x.id.all = c(x.id.all, x.id)
      }
      id_vector = c(id_vector, floor(mean(x.id.all)))
      print(i)
    }
    
    if(T){
      x = as.matrix(mat[,c(1:i)])
      x.rowsums = rowSums(x)
      x.id = length(which(x.rowsums>0))
      id_vector = c(id_vector, x.id)
      print(i)
    }
  }
  
  # print(quantile(id_vector))
  # plot(seq(1,mat.col), id_vector)
  # points(seq(1,102), id_vector[1:102], pch = 16, col = 'blue')
  # points(seq(103,110), id_vector[103:110], pch = 16, col = 'red')
  names(id_vector) = colnames(mat)
  return(id_vector)
}


#### get_mean_and_rank ####
get_mean_and_rank <- function(df){
  mean_fot = apply(df[,-1], 1, function(x){
    x = as.vector(unlist(x))
    x = x[x>0]
    mean(x)
  })
  merge_df = data.frame(
    GeneSymbol = df$GeneSymbol,
    FOT = mean_fot
  )
  merge_df = merge_df[order(merge_df$FOT, decreasing = T),]
  merge_df$Rank = seq(1, nrow(merge_df))
  return(merge_df)
}


#### get_xlim #### 
get_xlim <- function(mat_2d){
  xlim <- c(floor(min(mat_2d[,1]))-5, ceiling(max(mat_2d[,1]))+5)
  return(xlim)
}

#### get_ylim #### 
get_ylim <- function(mat_2d){
  ylim <- c(floor(min(mat_2d[,2]))-5, ceiling(max(mat_2d[,2]))+5)
  return(ylim)
}


#### get_lm_reproducibility ####
get_lm_reproducibility <- function(lm.df, min_id=2){
  
  # lm.df = p5.lm.df
  lm.df.row_id = apply(lm.df[,-1], 1, function(x){length(which(x>0))})
  lm.df.row_kept_index = which(lm.df.row_id>=min_id)
  lm.df = lm.df[lm.df.row_kept_index,]
  
  lm.df.gps = as.vector(unlist(lm.df$GeneSymbol))
  lm.df.mat = as.matrix(lm.df[,-1])
  
  t1.proteins = lm.df.gps[lm.df.mat[,1]>0]
  t1.proteins.fot = lm.df.mat[,1][lm.df.mat[,1]>0]
  t1.proteins.fot_q1 = quantile(t1.proteins.fot, probs = 1/3)
  t1.proteins.fot_q2 = quantile(t1.proteins.fot, probs = 0.50)
  t1.proteins.fot_q3 = quantile(t1.proteins.fot, probs = 2/3)
  t1.proteins.l = t1.proteins[t1.proteins.fot<=t1.proteins.fot_q1]
  t1.proteins.m = t1.proteins[t1.proteins.fot>t1.proteins.fot_q1 & t1.proteins.fot<=t1.proteins.fot_q3]
  t1.proteins.h = t1.proteins[t1.proteins.fot>t1.proteins.fot_q3]
  
  
  t_seq = seq(2,ncol(lm.df.mat))
  t_reproducibility = NULL
  for (i in t_seq) {
    ti_proteins = lm.df.gps[lm.df.mat[,i]>0]
    
    t1_ti_intersect_proteins = intersect(t1.proteins, ti_proteins)
    t1_ti_intersect_proteins.l = intersect(t1.proteins.l, ti_proteins)
    t1_ti_intersect_proteins.m = intersect(t1.proteins.m, ti_proteins)
    t1_ti_intersect_proteins.h = intersect(t1.proteins.h, ti_proteins)
    
    ti_reproducibility = length(t1_ti_intersect_proteins)/length(t1.proteins)
    ti_reproducibility.l = length(t1_ti_intersect_proteins.l)/length(t1.proteins.l)
    ti_reproducibility.m = length(t1_ti_intersect_proteins.m)/length(t1.proteins.m)
    ti_reproducibility.h = length(t1_ti_intersect_proteins.h)/length(t1.proteins.h)
    
    ti_v = c(ti_reproducibility, ti_reproducibility.l, ti_reproducibility.m, ti_reproducibility.h)
    t_reproducibility = rbind(t_reproducibility, ti_v)
  }
  rownames(t_reproducibility) = t_seq
  colnames(t_reproducibility) = c(
    'Global', 'Low-intensity', 'Medium-intensity', 'High-intensity'
  )
  return(t_reproducibility)
}


# snr_function
snr_function<-function(exprMat, group){
  
  library(data.table)
  
  IDs<- colnames(exprMat)
  IDs.group.mat<-data.table(
    IDs=IDs,
    group=group) 
  
  pca_prcomp <- prcomp(t(exprMat),retx=T)
  pcs <- as.data.frame(predict(pca_prcomp))
  pcs$Sample_id <- rownames(pcs)
  
  dt.perc.pcs <- data.table(PCX=1:nrow(pcs),
                            Percent=summary(pca_prcomp)$importance[2,],
                            AccumPercent=summary(pca_prcomp)$importance[3,])
  
  dt.dist <- data.table(ID.A = rep(IDs,each=length(IDs)),
                        ID.B = rep(IDs,time=length(IDs)))
  
  dt.dist$group.A <- IDs.group.mat[match(dt.dist$ID.A,IDs.group.mat$IDs)]$group
  dt.dist$group.B <- IDs.group.mat[match(dt.dist$ID.B,IDs.group.mat$IDs)]$group
  
  dt.dist[,Type:=ifelse(ID.A==ID.B,'Same',
                        ifelse(group.A==group.B,'Intra','Inter'))]
  
  dt.dist[,Dist:=sqrt(dt.perc.pcs[1]$Percent*(pcs[ID.A,1]-pcs[ID.B,1])^2+dt.perc.pcs[2]$Percent*(pcs[ID.A,2]-pcs[ID.B,2])^2)]
  
  dt.dist.stats <- dt.dist[,.(Avg.Dist=mean(Dist)),by=.(Type)]
  setkey(dt.dist.stats,Type)
  signoise <- dt.dist.stats['Inter']$Avg.Dist/dt.dist.stats['Intra']$Avg.Dist  
  return(signoise)
}



################function
snrdb_function<-function(exprMat, group){
  
  library(data.table)
  
  IDs<- colnames(exprMat)
  IDs.group.mat<-data.table(
    IDs=IDs,
    group=group) 
  
  pca_prcomp <- prcomp(t(exprMat),retx=T)
  pcs <- as.data.frame(predict(pca_prcomp))
  pcs$Sample_id <- rownames(pcs)
  
  dt.perc.pcs <- data.table(PCX=1:nrow(pcs),
                            Percent=summary(pca_prcomp)$importance[2,],
                            AccumPercent=summary(pca_prcomp)$importance[3,])
  
  dt.dist <- data.table(ID.A = rep(IDs,each=length(IDs)),
                        ID.B = rep(IDs,time=length(IDs)))
  
  dt.dist$group.A <- IDs.group.mat[match(dt.dist$ID.A,IDs.group.mat$IDs)]$group
  dt.dist$group.B <- IDs.group.mat[match(dt.dist$ID.B,IDs.group.mat$IDs)]$group
  
  dt.dist[,Type:=ifelse(ID.A==ID.B,'Same',
                        ifelse(group.A==group.B,'Intra','Inter'))]
  
  dt.dist[,Dist:=(dt.perc.pcs[1]$Percent*(pcs[ID.A,1]-pcs[ID.B,1])^2+dt.perc.pcs[2]$Percent*(pcs[ID.A,2]-pcs[ID.B,2])^2)]
  
  dt.dist.stats <- dt.dist[,.(Avg.Dist=mean(Dist)),by=.(Type)]
  setkey(dt.dist.stats,Type)
  signoise <- dt.dist.stats['Inter']$Avg.Dist/dt.dist.stats['Intra']$Avg.Dist  
  
  signoise_db <- 10*log10(signoise)
  return(signoise_db)
  
}
