source('D:/MAQC/20200508/Analysis/code_library/00_initialization.R')

# P5
p5.expnames = meta_data$Exp_code[meta_data$Person == 'P5']
p5.df = get_fot_df_by_expnames(fot_data, p5.expnames)
p5.deep.expnames = meta_data$Exp_code[which(
  meta_data$Person == 'P5' & meta_data$Purpose == 'DeepCoverage'
)]
p5.id_vector = get_cumulative_curve_vector(p5.df)
p5.deep.index = match(p5.deep.expnames, names(p5.id_vector))

# P6
p6.expnames = meta_data$Exp_code[meta_data$Person == 'P6']
p6.df = get_fot_df_by_expnames(fot_data, p6.expnames)
p6.deep.expnames = meta_data$Exp_code[which(
  meta_data$Person == 'P6' & meta_data$Purpose == 'DeepCoverage'
)]
p6.id_vector = get_cumulative_curve_vector(p6.df)
p6.deep.index = match(p6.deep.expnames, names(p6.id_vector))

# P7
p7.expnames = meta_data$Exp_code[meta_data$Person == 'P7']
p7.df = get_fot_df_by_expnames(fot_data, p7.expnames)
p7.deep.expnames = meta_data$Exp_code[which(
  meta_data$Person == 'P7' & meta_data$Purpose == 'DeepCoverage'
)]
p7.id_vector = get_cumulative_curve_vector(p7.df)
p7.deep.index = match(p7.deep.expnames, names(p7.id_vector))


# P8
p8.expnames = meta_data$Exp_code[meta_data$Person == 'P8']
p8.df = get_fot_df_by_expnames(fot_data, p8.expnames)
p8.deep.expnames = meta_data$Exp_code[which(
  meta_data$Person == 'P8' & meta_data$Purpose == 'DeepCoverage'
)]
p8.id_vector = get_cumulative_curve_vector(p8.df)
p8.deep.index = match(p8.deep.expnames, names(p8.id_vector))


if(T){
  x_seq = seq(1,ncol(p5.df)-1)
  plot(
    NA,
    xlab = 'MS experiments', ylab = 'GPs',
    xlim = c(0, 120), ylim = c(0, 13000), 
    yaxt = 'n', xaxt = 'n'
  )
  axis(1, at = seq(0,120,20), tcl = -0.5, cex.axis = 1)
  axis(1, at = seq(10,120,20), labels = rep('', 6), tcl = -0.3, cex.axis = 1)
  axis(
    2, at = seq(0,13000,2000), 
    labels = c('0', '2K', '4K', '6K', '8K', '10K', '12K'), 
    tcl = -0.5, cex.axis = 1
  )
  # axis(
  #   2, at = seq(1000,13000,2000), 
  #   labels = rep('', 7), 
  #   tcl = -0.3, cex.axis = 1
  # )
  
  # P5
  alpha_per = 1
  points(x_seq[-p5.deep.index], p5.id_vector[-p5.deep.index], pch = 16, col = alpha(p_colors[1],alpha_per))
  points(x_seq[p5.deep.index], p5.id_vector[p5.deep.index], pch = 17, col = alpha(p_colors[1],alpha_per))
  # P6
  points(x_seq[-p6.deep.index], p6.id_vector[-p6.deep.index], pch = 16, col = alpha(p_colors[2],alpha_per))
  points(x_seq[p6.deep.index], p6.id_vector[p6.deep.index], pch = 17, col = alpha(p_colors[2],alpha_per))
  
  # P7
  points(x_seq[-p7.deep.index], p7.id_vector[-p7.deep.index], pch = 16, col = alpha(p_colors[3],alpha_per))
  points(x_seq[p7.deep.index], p7.id_vector[p7.deep.index], pch = 17, col = alpha(p_colors[3],alpha_per))
  
  # P8
  points(x_seq[-p8.deep.index], p8.id_vector[-p8.deep.index], pch = 16, col = alpha(p_colors[4],alpha_per))
  points(x_seq[p8.deep.index], p8.id_vector[p8.deep.index], pch = 17, col = alpha(p_colors[4],alpha_per))
  
  
}
legend(
  5, 12500,
  bty = 'n',
  title = 'Single shot',
  legend = c('D5', 'D6', 'F7', 'M8'),
  # border = 'white',
  ncol = 1,
  # pt.cex = 2,
  cex = 0.8,
  title.adj = 0.5,
  x.intersp = 1.5,
  pt.bg = 'black',
  pt.cex = 1.2,
  # pt.lwd  = 2,
  col = p_colors,
  box.lwd = 0,
  pch = rep(16,4)
)
legend(
  30, 12500,
  bty = 'n',
  title = 'Multiple fractions',
  legend = c('D5', 'D6', 'F7', 'M8'),
  # border = 'white',
  ncol = 1,
  # pt.cex = 2,
  cex = 0.8,
  title.adj = 0.5,
  x.intersp = 1.5,
  pt.bg = 'black',
  pt.cex = 1.2,
  # pt.lwd  = 2,
  col = p_colors,
  box.lwd = 0,
  pch = 17
)







