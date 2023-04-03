#### set work directory ####
BASE_DIR = normalizePath(
  "D:/MAQC/20200508/Analysis/code_library", mustWork = F
)
setwd(BASE_DIR)

#### library ####
library('ggpubr')
library('ggsci')
library(ggforce)
library(ggplot2)
library(webr)
library(moonBook)
library('GGally')
library(networkD3)
library(viridis)

#### source ####
custom_lib = normalizePath(
  file.path(BASE_DIR, 'maqc_analysis_library.R'), mustWork = F
)
source(custom_lib, echo = TRUE)

pir_chart_lib = normalizePath(
  file.path(BASE_DIR, 'plot_piechart.R'), mustWork = F
)
source(pir_chart_lib, echo = TRUE)

#### library ####

#### input_dir ####
INPUT_DIR = normalizePath(
  file.path(BASE_DIR, '..', 'raw_input'), mustWork = F
)
fot_path = normalizePath(
  # file.path(INPUT_DIR,'maqc_proteome_profile_FOT_us=1.xlsx'),
  file.path(INPUT_DIR,'maqc_proteome_464_FOT_us=1.xlsx'),
  mustWork = F
)
metadata_path =  normalizePath(
  file.path(INPUT_DIR,'maqc_proteome_metadata_uniform.xlsx'), mustWork = F
)


#### import data ####
fot_data = read_excel(fot_path)
# fot_data_noqc = read_excel(fot_noqc_path)
meta_data = read_excel(metadata_path, sheet = 2)
# meta_data = read_excel(metadata_path, sheet = 2)

# parameters
p_colors = c(
  # '#FF8080', '#FFFF80', '#80FF40', '#8080FF'
  '#FF8080', '#FFC65D', '#80FF40', '#8080FF'
)
p_flags = c(paste('P', seq(5,8), sep = ''))

color.D5 = '#4CC3D9' #blue
color.D6 = '#7BC8A4' #green
color.F7 = '#FFC65D' #yellow
color.M8 = '#F16745' # red
p_colors = c(
  color.D5, color.D6, color.F7, color.M8
)







if(F){
  # person 
  # meta_data
  meta_data$UPerson = meta_data$Person
  p5.index = which(meta_data$Person == 'P5')
  p6.index = which(meta_data$Person == 'P6')
  p7.index = which(meta_data$Person == 'P7')
  p8.index = which(meta_data$Person == 'P8')
  meta_data$UPerson[p5.index] = 'D5'
  meta_data$UPerson[p6.index] = 'D6'
  meta_data$UPerson[p7.index] = 'F7'
  meta_data$UPerson[p8.index] = 'M8'
  
  # Labcode
  # meta_data
  # meta_data_labcode_map
  meta_data$Labcode = meta_data$Affiliation
  label_code_count = nrow(meta_data_labcode_map)
  for(i in 1:label_code_count){
    i_affiliation = meta_data_labcode_map$Affiliation[i]
    i_affiliation.index = which(meta_data$Affiliation == i_affiliation)
    meta_data$Labcode[i_affiliation.index] = meta_data_labcode_map$ThreeLetterLabCode[i]
  }
  write.csv(
    meta_data, 'maqc_proteome_metadata_uniform.csv', row.names = F
  )
}


# meta_data.replicates = meta_data[which(meta_data$Purpose=='Replications'),]
# meta_data.replicates.fot = get_fot_df_by_expnames(fot_data, meta_data.replicates$Exp_code)
# 
# write.csv(meta_data.replicates, 'maqc_meta_data.replicates_for_yy_20200813.csv', row.names = F)
# write.csv(meta_data.replicates.fot, 'maqc_fot.replicates_for_yy_20200813.csv', row.names = F)



if(T){
  get_fot_ref_df <- function(df){
    # qe.fot
    df.row_mean = rowMeans(df[,-1])
    df.row_expr_label = get_expr_label(df.row_mean)
    df.df = data.frame(
      GeneSymbol = df$GeneSymbol,
      RowID = apply(df[,-1], 1, function(x){length(which(x>0))}),
      MeanFOT = df.row_mean,
      ExprLab = df.row_expr_label
    )
    return(df.df)
  }
  
  # Single PCA
  # qe
  plot_single_pca <- function(fot_df, fot_df.ref, main_lab = 'G1', p_flags, p_colors, cex=2){
    # fot_df = qe.fot
    # fot_df.ref = qe.fot.df
    
    fot_df.mat = as.matrix(fot_df[,-1])
    fot_df.mat.scale = t(scale(t(fot_df.mat)))
    fot_group = rep(p_flags, each = 3)
    par(mfrow = c(2,2))
    snr_list = NULL
    if(T){
      # Low
      mat = fot_df.mat.scale[fot_df.ref$ExprLab=='Low',]
      snr = round(snr_function(mat, fot_group),2)
      colors = rep(p_colors, each=3)
      pca <- prcomp(((t(mat))), center = F, scale = F)
      importance <- summary(pca)$importance
      pc1 = paste('PC1 (', round(importance[2,1]*100, 2), '%)', sep = '')
      pc2 = paste('PC2 (', round(importance[2,2]*100, 2), '%)', sep = '')
      pca_predict <- predict(pca)
      pca_predict.2d = pca_predict[,c(1,2)]
      main = paste(main_lab, '; Low-intensity', ' Proteins (', nrow(mat), ')\nSNR = ', snr, sep = '') 
      plot(
        pca_predict.2d, 
        t = 'n', 
        main = main,
        xlab = pc1, ylab = pc2, 
        xlim = get_xlim(pca_predict.2d), 
        ylim = get_ylim(pca_predict.2d),
      )
      points(pca_predict.2d, pch = 16, col = colors, cex = cex)
      snr_list = c(snr_list, snr)
    }
    
    
    
    if(T){
      # Medium
      mat = fot_df.mat.scale[fot_df.ref$ExprLab=='Medium',]
      snr = round(snr_function(mat, fot_group),2)
      colors = rep(p_colors, each=3)
      pca <- prcomp(((t(mat))), center = F, scale = F)
      importance <- summary(pca)$importance
      pc1 = paste('PC1 (', round(importance[2,1]*100, 2), '%)', sep = '')
      pc2 = paste('PC2 (', round(importance[2,2]*100, 2), '%)', sep = '')
      pca_predict <- predict(pca)
      pca_predict.2d = pca_predict[,c(1,2)]
      main = paste(main_lab, '; Medium-intensity', ' Proteins (', nrow(mat), ')\nSNR = ', snr, sep = '') 
      plot(
        pca_predict.2d, 
        t = 'n', 
        main = main, 
        xlab = pc1, ylab = pc2, 
        xlim = get_xlim(pca_predict.2d), 
        ylim = get_ylim(pca_predict.2d),
      )
      points(pca_predict.2d, pch = 16, col = colors, cex = cex)
      snr_list = c(snr_list, snr)
    }
    
    
    if(T){
      # High
      mat = fot_df.mat.scale[fot_df.ref$ExprLab=='High',]
      snr = round(snr_function(mat, fot_group),2)
      colors = rep(p_colors, each=3)
      pca <- prcomp(((t(mat))), center = F, scale = F)
      importance <- summary(pca)$importance
      pc1 = paste('PC1 (', round(importance[2,1]*100, 2), '%)', sep = '')
      pc2 = paste('PC2 (', round(importance[2,2]*100, 2), '%)', sep = '')
      pca_predict <- predict(pca)
      pca_predict.2d = pca_predict[,c(1,2)]
      main = paste(main_lab, '; High-intensity', ' Proteins (', nrow(mat), ')\nSNR = ', snr, sep = '') 
      plot(
        pca_predict.2d, 
        t = 'n', 
        main = main, 
        xlab = pc1, ylab = pc2, 
        xlim = get_xlim(pca_predict.2d), 
        ylim = get_ylim(pca_predict.2d),
      )
      points(pca_predict.2d, pch = 16, col = colors, cex = cex)
      snr_list = c(snr_list, snr)
    }
    
    
    if(T){
      # Global
      mat = fot_df.mat.scale #[fot_df.ref$ExprLab=='Low',]
      snr = round(snr_function(mat, fot_group),2)
      colors = rep(p_colors, each=3)
      pca <- prcomp(((t(mat))), center = F, scale = F)
      importance <- summary(pca)$importance
      pc1 = paste('PC1 (', round(importance[2,1]*100, 2), '%)', sep = '')
      pc2 = paste('PC2 (', round(importance[2,2]*100, 2), '%)', sep = '')
      pca_predict <- predict(pca)
      pca_predict.2d = pca_predict[,c(1,2)]
      main = paste(main_lab, '; All', ' Proteins (', nrow(mat), ')\nSNR = ', snr, sep = '') 
      plot(
        pca_predict.2d, 
        t = 'n', 
        main = main, 
        xlab = pc1, ylab = pc2, 
        xlim = get_xlim(pca_predict.2d), 
        ylim = get_ylim(pca_predict.2d),
      )
      points(pca_predict.2d, pch = 16, col = colors, cex = cex)
      snr_list = c(snr_list, snr)
    }
    names(snr_list) = c('Low', 'Medium', 'High', 'Global')
    return(snr_list)
  }
}


