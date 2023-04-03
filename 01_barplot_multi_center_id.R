source('D:/MAQC/20200508/Analysis/code_library/00_initialization.R')

meta_data.replicates = meta_data[meta_data$Purpose=='Replications',]

# P5
meta_data.replicates.p5 = meta_data.replicates[meta_data.replicates$Person=='P5',]
meta_data.replicates.p5.uids = paste(
  meta_data.replicates.p5$Affiliation, 
  meta_data.replicates.p5$SInstrument,
  meta_data.replicates.p5$Batch, sep = '_'
)
p5.mean.id = tapply(meta_data.replicates.p5$GPs, meta_data.replicates.p5.uids, mean)


# P6
meta_data.replicates.p6 = meta_data.replicates[meta_data.replicates$Person=='P6',]
meta_data.replicates.p6.uids = paste(
  meta_data.replicates.p6$Affiliation, 
  meta_data.replicates.p6$SInstrument,
  meta_data.replicates.p6$Batch, sep = '_'
)
p6.mean.id = tapply(meta_data.replicates.p6$GPs, meta_data.replicates.p6.uids, mean)


# P7
meta_data.replicates.p7 = meta_data.replicates[meta_data.replicates$Person=='P7',]
meta_data.replicates.p7.uids = paste(
  meta_data.replicates.p7$Affiliation, 
  meta_data.replicates.p7$SInstrument,
  meta_data.replicates.p7$Batch, sep = '_'
)
p7.mean.id = tapply(meta_data.replicates.p7$GPs, meta_data.replicates.p7.uids, mean)


# P8
meta_data.replicates.p8 = meta_data.replicates[meta_data.replicates$Person=='P8',]
meta_data.replicates.p8.uids = paste(
  meta_data.replicates.p8$Affiliation, 
  meta_data.replicates.p8$SInstrument,
  meta_data.replicates.p8$Batch, sep = '_'
)
p8.mean.id = tapply(meta_data.replicates.p8$GPs, meta_data.replicates.p8.uids, mean)



# plot
main = "D5"
xlab = "The average GPs (X100)"
ylab = "Test runs (Total=24)"
par(mfrow = c(2,2))
p5.h = hist(
  p5.mean.id, 
  xlab = xlab, ylab = ylab, ylim = c(0,10), xlim = c(1500, 5000),
  main = "D5 (Single-shot)", col = pair_col[1], #'#EDEDED',
  xaxt = 'n'
)
axis(1, at = seq(1000, 5000, 500), labels = seq(1000, 5000, 500)/100)


hist(
  p6.mean.id, 
  xlab = xlab, ylab = ylab, ylim = c(0,10), xlim = c(1500, 5000),
  main = "D6 (Single-shot)", col = pair_col[1],#'#EDEDED',
  xaxt = 'n'
)
axis(1, at = seq(1000, 5000, 500), labels = seq(1000, 5000, 500)/100)

hist(
  p7.mean.id, 
  xlab = xlab, ylab = ylab, ylim = c(0,10), 
  main = "F7 (Single-shot)", col = pair_col[1], #'#EDEDED', 
  xlim = c(1500, 5000),
  xaxt = 'n'
)
axis(1, at = seq(1000, 5000, 500), labels = seq(1000, 5000, 500)/100)

hist(
  p8.mean.id, 
  xlab = xlab, ylab = ylab, ylim = c(0,10), 
  main = "M8 (Single-shot)", col = pair_col[1], # '#EDEDED', 
  xlim = c(1500, 5000),
  xaxt = 'n'
)
axis(1, at = seq(1000, 5000, 500), labels = seq(1000, 5000, 500)/100)

########################################################

meta_data.deepcoverage = meta_data[meta_data$Purpose=='DeepCoverage',]
# P5
meta_data.deepcoverage.p5 = meta_data.deepcoverage[meta_data.deepcoverage$Person=='P5',]
meta_data.deepcoverage.p6 = meta_data.deepcoverage[meta_data.deepcoverage$Person=='P6',]
meta_data.deepcoverage.p7 = meta_data.deepcoverage[meta_data.deepcoverage$Person=='P7',]
meta_data.deepcoverage.p8 = meta_data.deepcoverage[meta_data.deepcoverage$Person=='P8',]

par(mfrow = c(2,2))
xlab = "GPs (X100)"
ylab = "Test runs (Total=8)"
hist(
  meta_data.deepcoverage.p5$GPs, 
  xlab = xlab, ylab = ylab, ylim = c(0,5), xlim = c(4000, 9000),
  main = "D5 (Multiple fractions)", col = pair_col[2], 
  xaxt = 'n'
)
axis(1, at = seq(4000, 9000, 1000), labels = seq(4000, 9000, 1000)/100)

hist(
  meta_data.deepcoverage.p6$GPs, 
  xlab = xlab, ylab = ylab, ylim = c(0,5), xlim = c(4000, 9000),
  main = "D6 (Multiple fractions)", col = pair_col[2], 
  xaxt = 'n'
)
axis(1, at = seq(4000, 9000, 1000), labels = seq(4000, 9000, 1000)/100)

hist(
  meta_data.deepcoverage.p7$GPs, 
  xlab = xlab, ylab = ylab, ylim = c(0,5), xlim = c(4000, 9000),
  main = "F7 (Multiple fractions)", col = pair_col[2], 
  xaxt = 'n'
)
axis(1, at = seq(4000, 9000, 1000), labels = seq(4000, 9000, 1000)/100)

hist(
  meta_data.deepcoverage.p8$GPs, 
  xlab = xlab, ylab = ylab, ylim = c(0,5), xlim = c(4000, 9000),
  main = "M8 (Multiple fractions)",col = pair_col[2], 
  xaxt = 'n'
)
axis(1, at = seq(4000, 9000, 1000), labels = seq(4000, 9000, 1000)/100)
















