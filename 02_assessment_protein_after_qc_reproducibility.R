
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
UID_person_pairwise_repro_df = NULL
person_flag = 'P8'
for(i in 1:UID_set_count){
  UID_i = UID_set[i]
  cp_summary_table_data_subset = cp_summary_table_data[cp_summary_table_data$UID==UID_i,]
  UID_df = get_ifot_df_by_expnames(protein_fot_df, cp_summary_table_data_subset$name)
  UID_df_person_expnames = cp_summary_table_data_subset$name[cp_summary_table_data_subset$person==person_flag]
  UID_df_person_pairwise_repro_v = get_pariwise_reproducibility(UID_df, UID_df_person_expnames)
  UID_person_pairwise_repro_df = rbind(UID_person_pairwise_repro_df, UID_df_person_pairwise_repro_v)
}
rownames(UID_person_pairwise_repro_df) = UID_set
UID_mean = apply(UID_person_pairwise_repro_df, 1, mean)
UID_sd = apply(UID_person_pairwise_repro_df, 1, sd)
UID_group = UID_set
instrument = apply(data.frame(UID_group), 1, function(x){
  strsplit(x, split = '_')[[1]][2]
})
loading = apply(data.frame(UID_group), 1, function(x){
  strsplit(x, split = '_')[[1]][3]
})

df = data.frame(
  id = seq(1, length(UID_set)),
  UID = UID_group,
  mean = UID_mean,
  sd = UID_sd,
  loading = loading,
  instrument = instrument
)
df1 = df[df$loading=='200ng',]
df1 = df1[order(df1$instrument),]

df2 = df[df$loading=='500ng',]
df2 = df2[order(df2$instrument),]

df3 = df[df$loading=='660ng',]
df3 = df3[order(df3$instrument),]

df4 = df[df$loading=='1ug',]
df4 = df4[order(df4$instrument),]

df5 = rbind.data.frame(df1, df2, df3, df4)
df5$group = c(rep('200ng', nrow(df1)), rep('500ng', nrow(df2)), rep('660ng', nrow(df3)), rep('1ug', nrow(df4)))
df5$group = factor(df5$group, levels = c('200ng', '500ng', '660ng', '1ug'))

df5$id = factor(df5$id, levels = df5$id)

library(ggplot2)
# Default bar plot
p <- ggplot(df5, aes(x=id, y=mean)) + 
  geom_bar(aes(fill = group), stat="identity", position=position_dodge(), alpha=0.7, width = 0.8) +
  geom_errorbar(aes(x=id, ymin=mean-sd, ymax=mean+sd), width=.4, alpha=0.9, size=0.8, colour="orange") 
p <- p + labs(title=paste("Pairwise reproducibility of ", person_flag, ' (US>=1 & S>=2)', sep = ''), x="Center_ID", y="Jaccard index") + theme_classic()
p <- p + theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
p <- p + annotate(geom = 'text', x = (df5$id), y = 0, label = (df5$UID), srt = 90, hjust = -0.05, size=3)
p1 <- p + facet_grid(.~ df5$group, scales = 'free', space = 'free') + scale_y_continuous(breaks = seq(0,1,0.2)) + geom_hline(aes(yintercept=0.8), colour="#990000", linetype="dashed")

print(p1)












