#### plot venn ####
library(VennDiagram)
library(tidyverse)
# library(hrbrthemes)
# library(tm)
# library(proustr)
alpha_per = 0.3
venn.diagram(
  x = list(
    Fusion = as.vector(fusion.px.fot$GeneSymbol),
    Lumos = as.vector(lumos.px.fot$GeneSymbol)
  ),
  # category.names = c("P5 (10,197)" , "P6 (10,266)" , "P7 (10,273)", "P8 (10,346)"),
  filename = 'Venn/03_venn.png',
  output = TRUE ,
  imagetype="png" ,
  height = 480 , 
  width = 480 , 
  resolution = 300,
  compression = "lzw",
  lwd = 1,
  # col = p_colors,
  # fill = c(
  #   alpha(p_colors[1],alpha_per),
  #   alpha(p_colors[2],alpha_per),
  #   alpha(p_colors[3],alpha_per),
  #   alpha(p_colors[4],alpha_per)
  # ),
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
