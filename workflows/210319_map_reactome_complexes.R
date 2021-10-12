# reactome complexes
complex_reactome <- fread('inst/extdata/ReactomeHumanComplexes.txt')

# hmm, missing column?

# also, manual mapping of uniprot to HGNC?


# make reactome based on uniprot
uniprot <- gsub('uniprot\\:','',complex_reactome$participants)
uniprot <- lapply(strsplit(uniprot, split = '\\|'), function(x) x[!grepl('chebi',x)])
names(uniprot) <- complex_reactome$identifier
reactome <- stack(uniprot)
colnames(reactome) <- c('uniprot','identifier')
reactome <- merge(reactome, complex_reactome[,c(2,3,5)])

# map uniprot to hgnc
#write.table(data.frame(uniprot = unique(reactome$uniprot)), 
#            '~/Desktop/210604_uniprot_to_hgnc_mapping.tsv',
#            row.names = F, quote = F) 

# uploaded to uniprot
# retrvied here:
mapping <- read.delim('https://www.uniprot.org/mapping/M20210407A94466D2655679D1FD8953E075198DA81CBDA7G.tab')
colnames(mapping) <- c('uniprot', 'uniprot_to_gene')
reactome_mapped <- merge(mapping,reactome, all.x = T, all.y = T)
#reactome_mapped[order(reactome_mapped$name),]


## map ensemble ids

# assign column  to ensembl gene names
reactome_mapped$ensgid <- NA
newnames <- unlist(lapply(strsplit(reactome_mapped$uniprot[grepl('ensembl', reactome_mapped$uniprot)], split = ':'), function(x) x[2]))
reactome_mapped$ensgid[grepl('ensembl', reactome_mapped$uniprot)] <- newnames

# use biomart to do mapping
library(biomaRt)
ensid <- unlist(lapply(strsplit(reactome$uniprot[grepl('ensembl',reactome$uniprot)], split = ':'), function(x) x[2]))
ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")
att <- listAttributes(ensembl)
ensgid_to_hgnc <- getBM(attributes=c("ensembl_gene_id" ,"hgnc_symbol"), mart = ensembl)
ensgid_to_hgnc <- ensgid_to_hgnc[ensgid_to_hgnc$ensembl_gene_id %in% ensid,]
colnames(ensgid_to_hgnc) <- c('ensgid','ensgid_to_gene')
reactome_mapped_ensgid <- merge(ensgid_to_hgnc, reactome_mapped, all.x = T, all.y = T)

# finally map to a new gene column
d <- reactome_mapped_ensgid
d$gene <- d$uniprot_to_gene
d$gene[is.na(d$gene)] <- d$ensgid_to_gene[is.na(d$gene)]

# still 4% missing. fine for now. Some are mouse genes and ensembl transcripts.
colnames(d)
d <- d[,c(8,3,5,6)]
d <- d[!is.na(d$gene),]
d <- d[order(d$name),]
all_reactome_complexes <- d
colnames(all_reactome_complexes) <- c('genes', 'id', 'complex_id','complex')
save(file = 'data/all_reactome_complexes.rda', all_reactome_complexes, compress = 'xz')


