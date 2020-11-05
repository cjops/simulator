library(ggplot2)
library(ggmuller)

M_plot <- function(){
  tree_file <- Sys.glob("*_tree.csv")[1]
  short_name <- strsplit(tree_file, "_tree.csv")[1]
  pop_file <- paste(short_name, "_pop.csv", sep="")
  output_file <- paste(short_name, ".png", sep="")
  name <- read.table(paste(short_name, "_name.txt", sep=""))[1,1]
  
  #png(output_file, width=1200, height=600, pointsize = 14)
  plot_title = name
  tree <- read.csv(tree_file, colClasses = c("numeric", "character", "character"))
  pop <- read.csv(pop_file, colClasses = c("numeric", "character", "numeric"))
  Muller_df <- get_Muller_df(tree, pop)
  genotype_order <- unique(pop[pop$Population > 0, 'Identity'])
  Muller_df$Genotype <- factor(Muller_df$Identity, levels = genotype_order, ordered = TRUE)
  m <- Muller_plot(Muller_df, colour_by = 'Genotype', xlab = 'Timestep', add_legend = TRUE) + labs(title=plot_title)
  #dev.off()
  
  file.rename(tree_file, paste(short_name, "_tree_.csv", sep=""))
  
  return(m)
}