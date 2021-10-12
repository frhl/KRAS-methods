# libraries (all are strictly required to run the following code)
devtools::load_all()
library(igraph)
library(qgraph)
library(RColorBrewer)
library(genoppi) # can be downloaded from github. Ask frederik.

# data path
path <- 'inst/extdata/210216_SAINTexpress_KRAS.xlsx'
sheets <- excel_sheets(path)
in_sheets <- sheets[grepl('(WT)|(G12)|(G13)|(Q61H)',sheets)]

# get data
data <- read_workbook(path, in_sheets, rename_saint_sheets(in_sheets))
data <- set_saintscore_workbook(data)
data <- filter_workbook(data, SaintScore >= 0.8 & LogFC > 4)

#  get protein complexes data
data("all_corum_complexes")
corum_human <- all_corum_complexes[all_corum_complexes$organism == 'Human',]
corum_names <- unique(corum_human$complex)
corum <-lapply(corum_names, function(x) corum_human$genes[corum_human$complex == x])

names(corum) <- corum_names

# get inweb data
inweb <- get_inweb_list('KRAS')
inweb <- inweb$gene[inweb$significant]

# setup contrasts
cell <- unlist(lapply(strsplit(names(data), split = '_'), function(x) x[1]))
condition <- unlist(lapply(strsplit(names(data), split = '_'), function(x) x[2]))
contrasts <- data.frame(cell = cell, condition = condition)

lst_tables <- list()
pdf('derived/210330_network_contrast_condition_corum4_ss08_logfc4.pdf', width = 14, height = 14)
for (i in seq(1, nrow(contrasts), by = 2)){
  
  # use contrats to get right files
  cell <- contrasts$cell[i]
  conds <- contrasts$condition[contrasts$cell %in% cell]
  cond1 <- paste0(cell,'_',conds[1])
  cond2 <- paste0(cell,'_',conds[2])
  d1 <- data[[ names(data)[grepl(cond1, names(data))] ]]
  d2 <- data[[ names(data)[grepl(cond2, names(data))] ]]
  
  # title
  title = paste(cond1, 'vs', cond2, ' (CORUM pathways with  >60% partners present)\n SaintScore >= 0.8 & LogFC > 4')

  # find overlaps
  d_intersect <- data.frame(Bait = "KRAS", Prey = intersect(d1$Prey, d2$Prey), set = 'overlap', color = 'cyan')
  d1_only <- data.frame(Bait = "KRAS", Prey = d1$Prey[! d1$Prey %in% d2$Prey], set = cond1, color = 'blue')
  d2_only <- data.frame(Bait = "KRAS", Prey = d2$Prey[! d2$Prey %in% d1$Prey], set = cond2, color = 'green')
  
  # combine into single data.frame.
  d <- rbind(d_intersect, d1_only, d2_only)
  stopifnot(!any(duplicated(d$Prey)))
  
  # only include corum complexes if at least 501% proteins are present
  corum_ok <- unlist(lapply(corum, function(x) sum(x %in% d$Prey)/length(x) > 0.6))
  corum_final <- corum
  corum_numbers <- unlist(lapply(corum_final, function(x) paste0('[',sum(x %in% d$Prey),'/',length(x),']') ))
  names(corum_final) <- paste(names(corum_final), corum_numbers)
  
  corum_final <- corum_final[corum_ok]
  corum_final <- corum_final[!duplicated(corum_final)]
  names(corum_final) <- gsub("\\s*\\([^\\)]+\\)","",names(corum_final))
  corum_final <- corum_final[!duplicated(names(corum_final))]
  
  # make graph
  setup_group_graph(d, corum_final, known_interactors = inweb, main = title,
                    size.group = 3, size.bait = 5)
  
  # setup legend
  colors <- c('red',unique(d$color))
  names(colors) <- c('Bait (KRAS)',unique(d$set))
  legend('topright', legend = names(colors), col = colors,
         pch = 20, bty = "n",  pt.cex = 1.5, cex = 0.8, 
         text.col = "black", horiz = FALSE)
  

  
  # setup document 
  corum_stack <- stack(corum)
  colnames(corum_stack) <- c('Prey', 'Pathway')
  corum_merge <- merge(corum_stack, d, all.x = T, all.y = T)
  corum_merge$color <- NULL
  corum_merge <- corum_merge[order(corum_merge$Pathway, corum_merge$set),]
  
  outtable <- do.call(rbind, lapply(unique(as.character(na.omit(corum_merge$Pathway))), function(pathway){
    
    # get summary stats
    proteins <- corum_merge[corum_merge$Pathway %in% pathway,]
    total <- nrow(proteins)
    overlap <- proteins$Prey[proteins$set %in% 'overlap']
    set1 <- proteins$Prey[proteins$set %in% cond1]
    set2 <- proteins$Prey[proteins$set %in% cond2]
    unknown <- proteins$Prey[!proteins$Prey %in% c(overlap, set1, set2)]
    stopifnot(length(set1) + length(set2) + length(overlap) + length(unknown) == total)
    
    # make data.frame
    data.frame(
      set1_name = cond1, 
      set2_name = cond2,
      pathway = pathway,
      in_sample = length(c(overlap, set1, set2)),
      in_database = total,
      pct_in_sample = round(length(c(overlap, set1, set2)) / total, 4) * 100,
      inset1_n = length(set1),
      inset2_n = length(set2),
      inoverlap_n = length(overlap),
      set1_genes = paste(set1, collapse = ';'),
      set2_genes = paste(set2, collapse = ';'),
      overlap_genes = paste(overlap, collapse = ';'),
      missing_genes = paste(unknown, collapse = ';')
    )
    
  }))
  
  # save table
  outtable <- outtable[outtable$in_sample != 0,]
  outtable <-outtable[rev(order(outtable$in_sample)), ]
  lst_tables[[paste0(cond1,'_',cond2)]] <- outtable
  
}
graphics.off()

#write_xlsx(result_conditional,'derived/tables/210319_goa_bp_table_conditional_enrichment_logfc4.xlsx')
