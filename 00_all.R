source('D:/MAQC/20200508/Analysis/code_library/00_initialization.R')

# Pall
df.pall = get_fot_df_by_expnames(fot_data, meta_data$Exp_code)
df.pall.row_id = apply(df.pall[,-1], 1, function(x){
  length(which(x>0))
})
pall.hc_index = which(df.pall.row_id == 440)
pall.gps = as.vector(df.pall$GeneSymbol[pall.hc_index])
write.xlsx2(
  pall.gps, 'MAQC_HC_GPs=679.xlsx', row.names = F
)


# P5
df.p5 = get_fot_df_by_expnames(fot_data, meta_data$Exp_code[which(meta_data$Person=='P5')])
df.p5.row_id = apply(df.p5[,-1], 1, function(x){
  length(which(x>0))
})
p5.hc_index = which(df.p5.row_id == 110)
p5.gps = as.vector(df.p5$GeneSymbol[p5.hc_index])

# P6
df.p6 = get_fot_df_by_expnames(fot_data, meta_data$Exp_code[which(meta_data$Person=='P6')])
df.p6.row_id = apply(df.p6[,-1], 1, function(x){
  length(which(x>0))
})
p6.hc_index = which(df.p6.row_id == 110)
p6.gps = as.vector(df.p6$GeneSymbol[p6.hc_index])

# P7
df.p7 = get_fot_df_by_expnames(fot_data, meta_data$Exp_code[which(meta_data$Person=='P7')])
df.p7.row_id = apply(df.p7[,-1], 1, function(x){
  length(which(x>0))
})
p7.hc_index = which(df.p7.row_id == 110)
p7.gps = as.vector(df.p7$GeneSymbol[p7.hc_index])

# P8
df.p8 = get_fot_df_by_expnames(fot_data, meta_data$Exp_code[which(meta_data$Person=='P8')])
df.p8.row_id = apply(df.p8[,-1], 1, function(x){
  length(which(x>0))
})
p8.hc_index = which(df.p8.row_id == 110)
p8.gps = as.vector(df.p8$GeneSymbol[p8.hc_index])


hc_gps = unique(
  c(p5.gps, p6.gps, p7.gps, p8.gps)
)

write.xlsx(
  data.frame(HC_GPs = hc_gps),
  'MAQC_HC_GPs=1181.xlsx',
  row.names = F
)




