
sankey_plot_input_path = normalizePath(
  file.path(BASE_DIR, 'sankey_plot_input.xlsx'), mustWork = F
)
df = read_excel(sankey_plot_input_path, sheet = 'Sankey3')
colnames(df) = c('source', 'target', 'value')

links_df = as.data.frame(df)
links = links_df[order(links_df$source),]
# From these flows we need to create a node data frame: it lists every entities involved in the flow
nodes <- data.frame(
  name=c(as.character(links$source), 
         as.character(links$target)) %>% unique()
)
# nodes$name = factor(nodes$name, levels = c('Fusion-UP', 'Fusion-DOWN', 'Fusion-NA', 'Lumos-UP', 'Lumos-DOWN', 'Lumos-NA'))

# With networkD3, connection must be provided using id, not using real name like in the links dataframe.. So we need to reformat it.
links$IDsource <- match(links$source, nodes$name)-1 
links$IDtarget <- match(links$target, nodes$name)-1

# Make the Network
# '#33C860', # 绿色
# '#81B0FF', # 蓝色
# '#F9918A'  # 红色
# my_color <- 'd3.scaleOrdinal() .domain(["Fusion-UP", "Fusion-DOWN", "Fusion-NA", "Lumos-UP", "Lumos-DOWN", "Lumos-NA"]) .range(["#F9918A",  "#33C860", "grey", "#F9918A",  "#33C860", "grey"])'

p <- sankeyNetwork(
  Links = links, 
  Nodes = nodes,
  Source = "IDsource", Target = "IDtarget",
  Value = "value", 
  NodeID = "name", 
  fontSize = 15,
  # colourScale = my_color,
  # LinkGroup="group",
  sinksRight = FALSE
)
print(p)


# meta_data.cp = meta_data[which(meta_data$Purpose == 'Replications'),]

rf = c('Peptide', 'Protein')[2]
meta_data.cp = meta_data[which(meta_data$Preservation == rf & meta_data$Person=='P5'),]
a = data.frame(
  table(
    meta_data.cp$Labcode
  )
  
)
a = a[a$Freq>0,]
write.csv(a, 'a.csv', row.names = F)




meta_data.cp = meta_data[which(meta_data$Preservation == rf & meta_data$Person=='P5'),]
a = data.frame(
  table(
    meta_data.cp$Labcode,
    meta_data.cp$SInstrument
    
  )
  
)
a = a[a$Freq>0,]
write.csv(a, 'a1.csv', row.names = F)




























