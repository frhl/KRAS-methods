devtools::load_all()
library(UpSetR)
library(cowplot)

# data path
path <- 'inst/extdata/KRAS/210216_SAINTexpress_KRAS.xlsx'
sheets <- excel_sheets(path)
in_sheets <- sheets[grepl('(WT)|(G12)|(G13)|(Q61H)',sheets)]
#in_sheets <- in_sheets[grepl('starv', tolower(in_sheets))]
data <- read_workbook(path, in_sheets, rename_saint_sheets(in_sheets))
names(data)

pdf('derived/KRAS/plots/210827_upset_plots.pdf', width = 10, height = 8)
upset_data_subset <- lapply(data, function(x) unique(x$Prey[x$SaintScore >= 1]))
upset(fromList(upset_data_subset), 8, order.by = 'freq')
write.csv(list_to_matrix(upset_data_subset),'derived/KRAS/tables/210827_upset_plots_tables.csv')

#upset_data_all <- lapply(data, function(x) unique(x$Prey[x$SaintScore >= 0]))
#upset(fromList(upset_data_all), 8, order.by = 'freq')
graphics.off()


# 
#wd = '~/Desktop/to_andreas'
#dir.create(wd)
#write.csv(goa_bp_table, file.path(wd, 'goa_bp_table.csv'), quote = F)
#write.csv(goa_cc_table, file.path(wd, 'goa_cc_table.csv'), quote = F)
#write.csv(goa_mf_table, file.path(wd, 'goa_mf_table.csv'), quote = F)


