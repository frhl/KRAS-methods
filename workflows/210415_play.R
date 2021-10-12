# play around with the data

# load the scripts in the package
devtools::load_all()
library(genoppi)

# check the paths and the data
path <- 'inst/extdata/210216_SAINTexpress_KRAS.xlsx'
sheets <- excel_sheets(path)
print(sheets) # sheet names

# subset sheet names by WT and starcation condition
in_sheets <- sheets[grepl('(WT)|(G12)|(G13)|(Q61H)',sheets)]
in_sheets <- in_sheets[grepl('(starv)', tolower(in_sheets))]
print(in_sheets)

# load the data into memory and rename sheets
newnames <- rename_saint_sheets(in_sheets)
data <- read_workbook(path, in_sheets, newnames)
print(data)

# filter the data
data <- set_saintscore_workbook(data) # only extract relevant columns
data <- filter_workbook(data, SaintScore > 0.8 & LogFC > 2) # filter rows
print(data)

# investigate one dataset
names(data)

current_data <- data[["G12_Starvation"]]

# let's find some complexes
complexes <- find_corum_simple(current_data$Prey)
print(complexes)

# let's subset to the ones we like
complexes <- complexes[names(complexes) %in% c("PI4K2A-WASH complex", "Large Drosha complex")]
print(complexes)

# make network
make_group_graph(current_data, geneset = complexes) # draw full graph
 
# make simplified network
make_group_graph(current_data, geneset = complexes, simple_plot = T) # remove baits

# type in '?make_group_graph' to see the options for this function!!
?make_group_graph



