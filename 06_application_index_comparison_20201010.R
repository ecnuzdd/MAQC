# Identification
# correlation
# reproducibility
# CV
# SNR


if(T){
  meta_data.cp = meta_data[which(meta_data$Purpose == 'Replications'),]
  meta_data.cp.uids = paste(
    meta_data.cp$Labcode, 
    meta_data.cp$SInstrument, 
    meta_data.cp$Batch, 
    sep = '_'
  )
  meta_data.cp$UID = meta_data.cp.uids
  
  uids = sort(unique(meta_data.cp.uids))
  uids_code = c(
    paste('D0', seq(1,9), sep = ''),
    paste('D1', seq(0,9), sep = ''),
    paste('D2', seq(0,6), sep = '')
  )
  meta_data.cp$Code = NA
  for (i in 1:length(uids)) {
    index = which(meta_data.cp.uids==uids[i])
    meta_data.cp$Code[index] = uids_code[i]
  }
  meta_data.cp$Code = factor(meta_data.cp$Code, levels = uids_code)
  meta_data.cp = meta_data.cp[order(meta_data.cp$Code),]
  
  a = data.frame(
    unique(meta_data.cp$UID),
    unique(meta_data.cp$Code)
  )
  
}





get_pcorr <- function(df, method = 'pearson'){
  df = i.px.fot_data
  mat = as.matrix(i.px.fot_data[,-1])
  non_zero_min = min(mat[mat>0])/10
  # mat[mat==0] = non_zero_min
  mat.col = ncol(mat)
  corvector = NULL
  for (i in 1:mat.col) {
    x = mat[,i]
    for (j in 1:mat.col) {
      if(i!=j){
        y = mat[,j]
        xy.index = which(x>0|y>0)
        xy.cor = cor(x[xy.index], y[xy.index], method = method)
        corvector = c(corvector, xy.cor)
      }
    }
  }
  corvector = unique(corvector)
  cormedian = round(median(corvector), 2)
  return(cormedian)
}

get_reproducibility <- function(df){
  # df = i.px.fot_data
  mat = as.matrix(i.px.fot_data[,-1])
  non_zero_min = min(mat[mat>0])/10
  # mat[mat==0] = non_zero_min
  mat.col = ncol(mat)
  reprovector = NULL
  for (i in 1:mat.col) {
    x = mat[,i]
    for (j in 1:mat.col) {
      if(i!=j){
        y = mat[,j]
        xy.index = which(x>0 & y>0)
        xy.repro = length(xy.index)/length(x[x>0])
        reprovector = c(reprovector, xy.repro)
      }
    }
  }
  reprovector = unique(reprovector)
  repromedian = round(median(reprovector), 2)
  return(repromedian)
}


ref_names = read_excel('26_SNR_Figures_20201015.xlsx')


#### ID/PCorr/Repro/CV #### 
library(ggrepel)
p_list = list()
px_list = c('D5', 'D6', 'F7', 'M8')
# if(T){
for (pxi in 1:4) {
  
  # D5
  
  px = px_list[pxi]
  
  id_list = list()
  uids_code_count = length(uids_code)
  for (i in 1:uids_code_count) {
    i.px.metadata = meta_data.cp[which(meta_data.cp$Code==uids_code[i]),]
    i.px.metadata = i.px.metadata[order(i.px.metadata$Person),]
    i.px.metadata = i.px.metadata[i.px.metadata$UPerson==px,]
    i.px.fot_data = get_fot_df_by_expnames(fot_data, i.px.metadata$Exp_code)
    i.df = data.frame(
      ID = floor(mean(apply(i.px.fot_data[,-1], 2, function(x){length(which(x>0))}))),
      PCorr = get_pcorr(i.px.fot_data),
      Repro = get_reproducibility(i.px.fot_data),
      CV = round(median(apply(i.px.fot_data[,-1], 1, function(x){sd(x)/mean(x)})), 2)
    )
    id_list[[i]] = i.df
    print(i)
  }
  names(id_list) = uids
  print(unlist(lapply(id_list, function(x){x$CV})))
  
  df = data.frame(
    Labcode = uids_code,
    UID = uids,
    ID = unlist(lapply(id_list, function(x){x$ID})),
    PCorr = unlist(lapply(id_list, function(x){x$PCorr})),
    Repro = unlist(lapply(id_list, function(x){x$Repro})),
    CV = unlist(lapply(id_list, function(x){x$CV}))
  )
  mindex = match(df$UID, ref_names$UID)
  
  df$Labcode = ref_names$UID_new
  df$Labcode = factor(df$Labcode, levels = df$Labcode)
  
  
  
  
  p <- ggplot(data=df, aes(x=PCorr, y=Repro, size = ID, color = CV))
  p <- p + geom_point(alpha=0.7) 
  # p <- p + scale_fill_gradientn(colours=magma(100)[10:100])
  # p + scale_size(range = c(0.21, 0.76), name="Population (M)")
  p <- p + scale_color_viridis(discrete=F) 
  p <- p + geom_text_repel(
    aes(PCorr, Repro, label= Labcode),
    color="black",
    box.padding=unit(0.35, "lines"), point.padding=unit(0.5, "lines"),
    segment.colour = "grey50", size = 4
  )
  p <- p + scale_x_continuous(breaks=seq(0.4, 1, 0.1), limits = c(0.4, 1)) 
  # p <- p + scale_y_continuous(breaks=seq(0.8, 0.95, 0.05), limits = c(0.8, 0.95))
  p <- p + xlab('Pearson correlation coefficient') + ylab('Reproducibility')
  p <- p + ggtitle(paste(
    '4D evaluation from three replicates of single sample (', px, ')', sep = ''
  ))
  p <- p + theme_bw() 
  p_list[[pxi]] = p
}

ggarrange(#11*9
  p_list[[1]], p_list[[2]], p_list[[3]], p_list[[4]],
  nrow = 2, ncol = 2,
  labels = c('A', 'B', 'C', 'D')
)


temp <- venn.diagram(
  x = list(
    R1 = i.px.fot_data$GeneSymbol[i.px.fot_data$Exp036599>0],
    R2 = i.px.fot_data$GeneSymbol[i.px.fot_data$Exp036600>0],
    R3 = i.px.fot_data$GeneSymbol[i.px.fot_data$Exp036601>0]
  ),
  filename = NULL,
  main = 'JNU_Lumos_2',
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
  # col = jco_col,
  # fill = jco_col,
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
grid.draw(temp)
