if(T){
  
  get_quantile_normalization_UID_df <- function(df){
    library(preprocessCore)
    GeneSymbol = as.vector(df$GeneSymbol)
    Mat = as.matrix(df[,-1])
    Mat_colnames = colnames(Mat)
    Mat_Q =  normalize.quantiles(Mat, copy=TRUE)
    colnames(Mat_Q) = Mat_colnames
    df_new = data.frame(GeneSymbol, Mat_Q)
    return(df_new)
  }
  
  get_features_from_replicate_df <- function(df, group){
    GeneSymbol = as.vector(df$GeneSymbol)
    Mat = as.matrix(df[,-1])
    replicate_flag = apply(Mat, 1, function(x, y){
      x_y_id = tapply(x, y, function(x){length(which(x>0))})
      if(length(which(x_y_id>=2))){
        return(TRUE)
      }else{
        return(FALSE)
      }
    }, y = group)
    return(df[replicate_flag,])
  }
  
  get_features_from_df <- function(df, row_max = 1, row_id = 8, row_cv = 1){
    GeneSymbol = as.vector(df$GeneSymbol)
    Mat = as.matrix(df[,-1])
    Mat_minV = min(Mat[Mat>0])
    Mat_maxV = max(Mat[Mat>0])
    Mat_row_max = apply(Mat, 1, max)
    Mat_row_id = apply(Mat, 1, function(x){length(which(x>0))})
    Mat_row_cv = apply(Mat, 1, function(x){sd(x)/mean(x)})
    Mat_row_iqr = apply(Mat, 1, IQR)
    index = which(
      Mat_row_max >= row_max
      & Mat_row_id >= row_id
      & Mat_row_cv >= row_cv
    )
    df_subset = df[index,]
    return(df_subset)
  }
  
  # 9
  get_xlim <- function(mat_2d){
    xlim <- c(floor(min(mat_2d[,1]))-5, ceiling(max(mat_2d[,1]))+5)
    return(xlim)
  }
  
  # 10
  get_ylim <- function(mat_2d){
    ylim <- c(floor(min(mat_2d[,2]))-5, ceiling(max(mat_2d[,2]))+5)
    return(ylim)
  }
  
  # 6
  get_topx_df <- function(df, topx){
    mat_value = df[,-1]
    mat_value_col = ncol(mat_value)
    union_index = NULL
    topx_seq = seq(1, topx)
    for(i in 1:mat_value_col){
      x = as.vector(unlist(mat_value[,i]))
      x_order = order(x, decreasing = T)
      # x_order_index = match(topx_seq, x_order)
      x_order_index = x_order[topx_seq]
      if(2 %in% x_order_index){
        print(i)
      }
      union_index = union(union_index, x_order_index)
    }
    topx_df = df[union_index,]
    return(topx_df)
  }
  
  # 6
  get_topx_df_1 <- function(df, topx){
    mat_value = df[,-1]
    mat_value_col = ncol(mat_value)
    topx_seq = seq(1, topx)
    for(i in 1:mat_value_col){
      x = as.vector(unlist(mat_value[,i]))
      x_order = order(x, decreasing = T)
      x_order_index = x_order[topx_seq]
      y = x
      y[-x_order_index] = 0
      mat_value[,i] = y
    }
    kept_index = which(rowSums(mat_value)>0)
    topx_df = data.frame(
      Symbol = as.vector(unlist(df[,1])),
      mat_value
    )[kept_index,]
    return(topx_df)
  }
  
  
  get_gg_plot <- function(data_df, geom_text_x, geom_text_y1, lim_v, title){
    theme1 <- theme(panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    panel.background = element_blank(),
                    axis.line.x = element_line(color="#919191", size = 0.1),
                    axis.line.y = element_line(color="#919191", size = 0.1))
    
    data = data.frame(
      x = log10(data_df[,1]),
      y = log10(data_df[,2])
    )
    m <- lm(y ~ x, data)
    r2 <- format(summary(m)$r.squared, digits = 3)
    pvalue <- format(summary(m)$coefficients[2,4], digits = 3)
    gg <- ggplot(data, aes(x = x, y = y))
    gg <- gg + xlim(lim_v[1], lim_v[2]) + ylim(lim_v[1], lim_v[2])
    gg <- gg + geom_hex(bins=80)
    gg <- gg + scale_fill_gradient(low='blue', high='orange')
    gg <- gg + xlab("log10(iFOT)") + ylab("log10(iFOT)") + theme1
    gg <- gg + geom_smooth(method="lm", colour = 'black')
    gg <- gg + geom_text(x = geom_text_x, y = geom_text_y1, label = paste('R^2 = ', r2, sep = ''), colour = 'black') + 
      geom_text(x = geom_text_x, y = geom_text_y1-0.4, label = paste('Pvalue = ', pvalue, sep = ''), colour = 'black')
    gg <- gg + guides(fill=FALSE) + ggtitle(label=title)
    return(gg)
  }
  
  
  plot_correaltion_matrix <- function(df_mat_subset, geom_text_x, geom_text_y1, lim_v){
    require(ggplot2)
    theme2 <- theme(axis.line=element_blank(),axis.text.x=element_blank(),
                    axis.text.y=element_blank(),axis.ticks=element_blank(),
                    axis.title.x=element_blank(),
                    axis.title.y=element_blank(),legend.position="none",
                    panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
                    panel.grid.minor=element_blank(),plot.background=element_blank())
    
    data_df = data.frame(
      x = seq(lim_v[1], lim_v[2]),
      y = seq(lim_v[1], lim_v[2])
    )
    plot_null <- ggplot(data_df, aes(x = x, y = y)) + geom_blank()  + theme2
    
    
    require(gridExtra)
    plot1 <- get_gg_plot(df_mat_subset[,c(1,2)], geom_text_x, geom_text_y1, lim_v, 'R1 VS R2')
    plot4 <- plot_null
    plot2 <- get_gg_plot(df_mat_subset[,c(1,3)], geom_text_x, geom_text_y1, lim_v, 'R1 VS R3')
    plot3 <- get_gg_plot(df_mat_subset[,c(2,3)], geom_text_x, geom_text_y1, lim_v, 'R2 VS R3')
    print(grid.arrange(plot1, plot2, plot3, plot4, ncol=2))
  }
  
  get_cor_lst <- function(df){
    # df = UID_df_0_P5
    df_mat = df[,-1]
    df_mat_rowid = apply(df_mat, 1, function(x){length(which(x>0))})
    df_mat_row_kept_index = which(df_mat_rowid >= 2)
    df_mat_subset = as.matrix(df_mat[df_mat_row_kept_index,])
    zero_index = which(df_mat_subset==0)
    min_v = floor(min(log10(df_mat_subset[-zero_index])))
    max_v = ceiling(max(log10(df_mat_subset[-zero_index])))
    geom_text_x = min_v + 1
    geom_text_y1 = max_v - 1
    lim_v = c(min_v, max_v) 
    df_mat_subset[zero_index] = NA
    cor_matrix <- cor(
      df_mat_subset, 
      method = 'spearman', 
      use = c("everything", "all.obs", "complete.obs", "na.or.complete", "pairwise.complete.obs")[5]
    )
    cor_v = unique(as.vector(cor_matrix)[as.vector(cor_matrix)<1])
    cor_v_m = format(mean(cor_v), digits = 3)
    cor_v_sd = format(sd(cor_v), digits = 3)
    plot_correaltion_matrix(df_mat_subset, geom_text_x, geom_text_y1, lim_v)
    print(paste('Average(cor) = ', cor_v_m, sep = ''))
    print(paste('Std(cor) = ', cor_v_sd, sep = ''))
    lst = list(
      cor_v_m = cor_v_m,
      cor_v_sd =cor_v_sd
    )
    return(lst)
    
  }
  
  
  format_nairen_data <- function(df){
    # df = dim_ifot_data
    Symbol = as.vector(unlist(df[,1]))
    Mat = as.matrix(df[,-1])
    na_index = which(is.na(Mat))
    if(length(na_index)>0){
      Mat[na_index] = 0
    }
    df_new = data.frame(Symbol, Mat)
    return(df_new)
  }
}


#########################################################################
BASE_DIR = normalizePath('D:/my projects/MAQC/MAQC_firmiana_data/T_analysis/MAQC_Analysis_crossplatform/cross_platform_data')
setwd(BASE_DIR)
source('00_initialization.R')
library('stringr')



protein_fot_path = normalizePath(
  file.path(BASE_DIR, 'merge_data_cross_platform_20190821', 'merge_fot_0_df.csv'), mustWork = F
)
protein_fot_df = format_nairen_data(read.csv(protein_fot_path))


summary_table_path = normalizePath(
  file.path(BASE_DIR, 'cross_platform_local_grouper_20190821.xlsx'), mustWork = F
)

cp_summary_table_data = read_excel(summary_table_path, sheet = 'ALL')
cp_exps_id = paste(
  cp_summary_table_data$location,
  cp_summary_table_data$instrument,
  # cp_summary_table_data$person,
  cp_summary_table_data$loading,
  cp_summary_table_data$code,
  sep = '_'
)
cp_summary_table_data$UID = cp_exps_id



UID_set = unique(cp_summary_table_data$UID)



#########################
P5_cor_df <- NULL
P6_cor_df <- NULL
P7_cor_df <- NULL
P8_cor_df <- NULL



#########################
correlation_dir = normalizePath(file.path(BASE_DIR, '01_correlation'), mustWork = F)

i = 1
for(i in 1:length(UID_set)){
  UID_i = UID_set[i]
  cp_summary_table_data_subset = cp_summary_table_data[cp_summary_table_data$UID==UID_i,]
  loading = unique(cp_summary_table_data_subset$loading)
  UID_df_0 = get_ifot_df_by_expnames(protein_fot_df, cp_summary_table_data_subset$name)
  UID_df_0_colnames = colnames(UID_df_0)
  UID_df_0_colnames = str_replace_all(UID_df_0_colnames, '_FOT', '')
  colnames(UID_df_0) = UID_df_0_colnames
  
  pdf_path = normalizePath(file.path(correlation_dir, paste(i, '_', UID_i, '.pdf', sep = '')), mustWork = F)
  pdf(pdf_path, height = 10, width = 10)
  UID_df_0_P5_expnames = cp_summary_table_data_subset$name[cp_summary_table_data_subset$person=='P5']
  UID_df_0_P5 = get_ifot_df_by_expnames(protein_fot_df, UID_df_0_P5_expnames)
  UID_df_0_P5_cor_lst = get_cor_lst(UID_df_0_P5)
  P5_cor_df <- rbind(P5_cor_df, as.numeric(unlist(UID_df_0_P5_cor_lst)))
  
  UID_df_0_P6_expnames = cp_summary_table_data_subset$name[cp_summary_table_data_subset$person=='P6']
  UID_df_0_P6 = get_ifot_df_by_expnames(protein_fot_df, UID_df_0_P6_expnames)
  UID_df_0_P6_cor_lst = get_cor_lst(UID_df_0_P6)
  P6_cor_df <- rbind(P6_cor_df, as.numeric(unlist(UID_df_0_P6_cor_lst)))
  
  UID_df_0_P7_expnames = cp_summary_table_data_subset$name[cp_summary_table_data_subset$person=='P7']
  UID_df_0_P7 = get_ifot_df_by_expnames(protein_fot_df, UID_df_0_P7_expnames)
  UID_df_0_P7_cor_lst = get_cor_lst(UID_df_0_P7)
  P7_cor_df <- rbind(P7_cor_df, as.numeric(unlist(UID_df_0_P7_cor_lst)))
  
  UID_df_0_P8_expnames = cp_summary_table_data_subset$name[cp_summary_table_data_subset$person=='P8']
  UID_df_0_P8 = get_ifot_df_by_expnames(protein_fot_df, UID_df_0_P8_expnames)
  UID_df_0_P8_cor_lst = get_cor_lst(UID_df_0_P8)
  P8_cor_df <- rbind(P8_cor_df, as.numeric(unlist(UID_df_0_P8_cor_lst)))
  dev.off()
  print(i)
}
rownames(P5_cor_df) = UID_set
colnames(P5_cor_df) = c('cor_v_m', 'cor_v_sd')
P5_cor_df = data.frame(P5_cor_df)
P5_cor_df$R2 = P5_cor_df$cor_v_m*P5_cor_df$cor_v_m

rownames(P6_cor_df) = UID_set
colnames(P6_cor_df) = c('cor_v_m', 'cor_v_sd')
P6_cor_df = data.frame(P6_cor_df)
P6_cor_df$R2 = P6_cor_df$cor_v_m*P6_cor_df$cor_v_m

rownames(P7_cor_df) = UID_set
colnames(P7_cor_df) = c('cor_v_m', 'cor_v_sd')
P7_cor_df = data.frame(P7_cor_df)
P7_cor_df$R2 = P7_cor_df$cor_v_m*P7_cor_df$cor_v_m

rownames(P8_cor_df) = UID_set
colnames(P8_cor_df) = c('cor_v_m', 'cor_v_sd')
P8_cor_df = data.frame(P8_cor_df)
P8_cor_df$R2 = P8_cor_df$cor_v_m*P8_cor_df$cor_v_m

write.csv(P5_cor_df, normalizePath(file.path(correlation_dir, 'P5_cor_df.csv'), mustWork = F))
write.csv(P6_cor_df, normalizePath(file.path(correlation_dir, 'P6_cor_df.csv'), mustWork = F))
write.csv(P7_cor_df, normalizePath(file.path(correlation_dir, 'P7_cor_df.csv'), mustWork = F))
write.csv(P8_cor_df, normalizePath(file.path(correlation_dir, 'P8_cor_df.csv'), mustWork = F))

if(T){
  library(pracma)
  d = density(P5_cor_df$R2)
  plot(d, main = 'Spearman correlation (R^2)')
  rug(P5_cor_df$R2)
  dx = d$x
  dy = d$y
  for(ii in 2:250){
    a = trapz(dx[1:ii], dy[1:ii])
    if(a>0.025){
      break
    }
  }
  ii = ii - 1
  rug(dx[ii], col = 'blue', lwd = 2)
  
  for(ii in 251:512){
    a = trapz(dx[ii:512], dy[ii:512])
    if(a<0.025){
      break
    }
  }
  ii = ii - 1
  rug(dx[ii], col = 'red', lwd = 2)
}







# UID_df_0
QconCAT_proteins = c(
 'ACACA	ACAD10	ACAT2	ACLY	ACO1	ACSS2	ALDOA	BPGM	CPT1A	CS	DGAT1	DLAT	DLST	ECHS1	ENO1	FASN	FH	GAPDH	GPAM	GPI	HADH	HADHA	IDH3A	MDH2	ME1	OGDH	PCK1	PDHB	PFKFB2	PFKP	PGK1	PKM	SCD'
)
QconCAT_proteins = strsplit(QconCAT_proteins, split = '\t')[[1]]

df = data.frame(
  Symbol = UID_df_0$Symbol,
  FOT_ave = apply(UID_df_0[,-1], 1, mean)
)
df$Rank = 3343 - rank(df$FOT_ave) + 1
df$Flag = rep(NA, 3343)
df$Flag[na.omit(match(QconCAT_proteins, df$Symbol))] = TRUE

df_subset = df[which(df$Flag),]

plot(
  df$Rank, log10(df$FOT_ave), 
  main = 'The dynamic ranges of Quartet proteomes\nand proteins as internal standards',
  xlab = 'Rank', ylab = 'log10(FOT)', 
  type = 'p',
  ylim = c(-4, 4)
)

points(
  df_subset$Rank, log10(df_subset$FOT_ave), pch = 16, col = 'red'
)

df_subset_rank_order = order(df_subset$Rank)
text(
  df_subset$Rank, 
  log10(df_subset$FOT_ave), 
  labels = as.vector(df_subset$Symbol),
  pos = 4, cex = 0.5, offset = 0.5
)

text(
  df_subset$Rank[df_subset_rank_order[7:10]], 
  log10(df_subset$FOT_ave)[df_subset_rank_order[7:10]], 
  labels = as.vector(df_subset$Symbol[df_subset_rank_order[7:10]]),
  pos = 4, cex = 0.5, offset = 1
)

