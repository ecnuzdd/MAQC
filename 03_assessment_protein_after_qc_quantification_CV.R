
source('00_initialization.R')

protein_fot_path = normalizePath(
  file.path(PROT_OUTPUT_DIR, 'merge_ifot_df.csv'), mustWork = F
)
protein_fot_df = read.csv(protein_fot_path)

summary_table_path = CP_METADATA_PATH

cp_site = CP_SITES[1]
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
UID_set_count = length(UID_set)
UID_person_cv_list_list = list()
person_flag = 'P5'
col_v = NULL
x_max_v = NULL
y_min_v = NULL
y_max_v = NULL
loading_v = NULL
for(i in 1:UID_set_count){
  UID_i = UID_set[i]
  cp_summary_table_data_subset = cp_summary_table_data[cp_summary_table_data$UID==UID_i,]
  loading = unique(cp_summary_table_data_subset$loading)
  loading_v = c(loading_v, loading)
  UID_df = get_ifot_df_by_expnames(protein_fot_df, cp_summary_table_data_subset$name)
  UID_df_person_expnames = cp_summary_table_data_subset$name[cp_summary_table_data_subset$person==person_flag]
  UID_cv_list = get_cv_list(UID_df, UID_df_person_expnames)
  UID_person_cv_list_list[[i]] = UID_cv_list
  if(length(which(cp_summary_table_data_subset$loading=='200ng' & cp_summary_table_data_subset$person==person_flag))>0){
    col_i = '#FA9F99'
  }
  if(length(which(cp_summary_table_data_subset$loading=='500ng' & cp_summary_table_data_subset$person==person_flag))>0){
    col_i = '#A3C64C'
  }
  if(length(which(cp_summary_table_data_subset$loading=='660ng' & cp_summary_table_data_subset$person==person_flag))>0){
    col_i = '#4CD2D6'
  }
  if(length(which(cp_summary_table_data_subset$loading=='1ug' & cp_summary_table_data_subset$person==person_flag))>0){
    col_i = '#D8A3FF'
  }
  col_v = c(col_v, col_i)
  y_min_v = c(y_min_v, min(UID_cv_list$all))
  y_max_v = c(y_max_v, max(UID_cv_list$all))
}

order_index = order(loading_v)
order_index_new = c(order_index[11:17], order_index[1:10])


ylim = c(floor(min(y_min_v)), ceiling(max(y_max_v)))

par(mfrow = c(2,1))

boxplot(NA, xlim = c(0, 18), ylim = c(0, 2), xaxt = 'n', yaxt = 'n', main = paste('CV distribution of ', person_flag, sep = ''),  ylab = 'CV')
axis(2, at = seq(0, 2, 0.2))
axis(1, at = seq(1,17), labels = seq(1,17))
abline(v = line_at, lty = 3)
at = 0
for(i in order_index_new){
  at = at + 1
  boxplot(
    UID_person_cv_list_list[[i]][2:4], at = at_df[at,], add = T, 
    col = c('#7A023C', '#FF5983', '#EBEDF4'),
    boxwex = 0.2, xaxt = 'n', yaxt = 'n'
  )
  
}

legend_df = data.frame(
  UID = UID_set[order_index_new],
  Code = seq(1,17)
)
abline(v = 0.5, lty = 3)


boxplot(NA, xlim = c(0, 18), ylim = c(0, 2), xaxt = 'n', yaxt = 'n', ylab = 'CV')
axis(2, at = seq(0, 2, 0.2))
at_df = NULL
at = c(0.25, 0.5, 0.75, 1)
line_at = NULL
for(i in 1:17){
  at = c(i-0.25, i, i+0.25)
  at_df = rbind(at_df, at)
  line_at_i = at[3] + 0.25
  line_at = c(line_at, line_at_i)
}
axis(1, at = seq(1,17), labels = seq(1,17))
# abline(v = line_at, lty = 3)



at = 0
for(i in order_index_new){
  at = at + 1
  boxplot(UID_person_cv_list_list[[i]][1], at = at, add = T, col = col_v[i], boxwex = 1.5, xaxt = 'n', yaxt = 'n')
  
}





