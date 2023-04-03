# e480.px.col_gps, e480.df
e480.repro = get_px_reproducibility(e480.px.col_gps, ref_df = NULL, ref_lab = NULL)
e480.repro.h = get_px_reproducibility(e480.px.col_gps, ref_df = e480.df, ref_lab = 'High')
e480.repro.m = get_px_reproducibility(e480.px.col_gps, ref_df = e480.df, ref_lab = 'Medium')
e480.repro.l = get_px_reproducibility(e480.px.col_gps, ref_df = e480.df, ref_lab = 'Low')
df.e480 = data.frame(
  Intensity = rep(c('Low', 'Medium', 'High', 'Global'), each = 6),
  Reproducibity = c(e480.repro.l, e480.repro.m, e480.repro.h, e480.repro)
)
df.e480$Intensity = factor(df.e480$Intensity, levels = c('Low', 'Medium', 'High', 'Global'))