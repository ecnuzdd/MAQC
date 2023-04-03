source('D:/MAQC/20200508/Analysis/code_library/00_initialization.R')

# total
total.id.no_qc = nrow(fot_data_noqc)
total.id = nrow(fot_data)


# P5
p5.expnames = meta_data$Exp_code[meta_data$Person == 'P5']
p5.df.noqc = get_fot_df_by_expnames_with_na(fot_data_noqc, p5.expnames)
p5.df = get_fot_df_by_expnames(fot_data, p5.expnames)
p5.id.noqc = nrow(p5.df.noqc)
p5.id = nrow(p5.df)

# P6
p6.expnames = meta_data$Exp_code[meta_data$Person == 'P6']
p6.df.noqc = get_fot_df_by_expnames_with_na(fot_data_noqc, p6.expnames)
p6.df = get_fot_df_by_expnames(fot_data, p6.expnames)
p6.id.noqc = nrow(p6.df.noqc)
p6.id = nrow(p6.df)

# P7
p7.expnames = meta_data$Exp_code[meta_data$Person == 'P7']
p7.df.noqc = get_fot_df_by_expnames_with_na(fot_data_noqc, p7.expnames)
p7.df = get_fot_df_by_expnames(fot_data, p7.expnames)
p7.id.noqc = nrow(p7.df.noqc)
p7.id = nrow(p7.df)

# P8
p8.expnames = meta_data$Exp_code[meta_data$Person == 'P8']
p8.df.noqc = get_fot_df_by_expnames_with_na(fot_data_noqc, p8.expnames)
p8.df = get_fot_df_by_expnames(fot_data, p8.expnames)
p8.id.noqc = nrow(p8.df.noqc)
p8.id = nrow(p8.df)


bp_df = cbind(
  c(p5.id, p5.id.noqc),
  c(p6.id, p6.id.noqc),
  c(p7.id, p7.id.noqc),
  c(p8.id, p8.id.noqc),
  c(total.id, total.id.no_qc)
)
colnames(bp_df) = c('P5', 'P6', 'P7', 'P8', 'Total')
rownames(bp_df) = c('After', 'Before')
bp_mat = as.matrix(bp_df)
bp_mat[2,] = bp_mat[2,] - bp_mat[1,]
# plot
if(T){
  # jco_col = show_col(pal_jco("default", alpha = 0.3)(10))
  library(RColorBrewer)
  display.brewer.pal(n = 8, name = 'Paired')
  pair_col = brewer.pal(n = 8, name = 'Paired')
  bp = barplot(
    bp_mat,
    main = 'Protein coverage in Quartet',
    ylim = c(0, 20000),
    # xlim = c(0, 5),
    yaxt = 'n',
    ylab = 'Protein indetification',
    # width = 0.2,
    col = pair_col[2:1],# c('firebrick', 'grey'),
    space = 0.6,
    # width = 1,
    names.arg = c('D5 (110)', 'D6 (110)', 'F7 (110)', 'M8 (110)', 'Total (440)')
  )
  axis(2, at = seq(0, 20000, 5000), labels = c(0, '5K', '10K', '15K', '20K'))
  abline(h = 19979, lty = 3, lwd = 2)
  x = bp
  y1 = bp_df[1,]
  text(x, y1, labels = y1, pos = 3, col = pair_col[2])
  y2 = bp_df[2,]
  text(x, y2, labels = y2, pos = 3, col = pair_col[1])
  y3 = 19979
  text(
    x[3], y3,
    labels = 'NCBI human RefSeq protein database \n(released on 04-07-2013, 32,015 entries mapping to 19,979 Genes) )', 
    pos = 1, col = 'black', cex = 0.7
  )
}




#### plot venn ####
library(VennDiagram)
library(tidyverse)
# library(hrbrthemes)
# library(tm)
# library(proustr)
alpha_per = 0.3
venn.diagram(
  x = list(
    P5 = as.vector(p5.df$GeneSymbol),
    P6 = as.vector(p6.df$GeneSymbol),
    P7 = as.vector(p7.df$GeneSymbol),
    P8 = as.vector(p8.df$GeneSymbol)
  ),
  category.names = c("D5 (10,197)" , "D6 (10,266)" , "F7 (10,273)", "M8 (10,346)"),
  filename = 'Venn/venn.png',
  output = TRUE ,
  imagetype="png" ,
  height = 480 , 
  width = 480 , 
  resolution = 300,
  compression = "lzw",
  lwd = 1,
  col = p_colors,
  fill = c(
    alpha(p_colors[1],alpha_per),
    alpha(p_colors[2],alpha_per),
    alpha(p_colors[3],alpha_per),
    alpha(p_colors[4],alpha_per)
  ),
  cex = 0.5,
  fontfamily = "sans",
  cat.cex = 0.3,
  cat.default.pos = "outer",
  # cat.pos = c(-2, 27, -135, 135),
  # cat.dist = c(0.055, 0.055, 0.085),
  cat.fontfamily = "sans"
  # cat.col = p_colors
  # rotation = 1
)


intersect_genes = intersect(
  as.vector(p5.df$GeneSymbol), intersect(
    as.vector(p6.df$GeneSymbol), intersect(
      as.vector(p7.df$GeneSymbol), as.vector(p8.df$GeneSymbol)
    )
  )
)
union_genes = unique(
  c(
    as.vector(p5.df$GeneSymbol), as.vector(p6.df$GeneSymbol),
    as.vector(p7.df$GeneSymbol), as.vector(p8.df$GeneSymbol)
  )
)
diff_genes = setdiff(
  union_genes, intersect_genes
)

library(xlsx)
write.xlsx(
  data.frame(union_genes = union_genes),
  'Quartet_genes.xlsx', col.names = F, sheetName = 'Union'
)
write.xlsx(
  data.frame(intersect_genes = intersect_genes),
  'Quartet_genes.xlsx', col.names = F, sheetName = 'Intersect'
)
write.xlsx(
  data.frame(diff_genes = diff_genes),
  'Quartet_genes.xlsx', col.names = F, sheetName = 'Diff'
)


  
  