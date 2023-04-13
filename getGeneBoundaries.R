library(optparse)

option_list = list(make_option(c("-c", "--chromosome"), type='character', default=NULL, help='chromosome of interest', metavar='character'), 
		   make_option(c("-i", "--genepred_cds_path"), type='character', default=NULL, help='GTF file path containing the CDS coordinates in chromosome of interest', metavar='character'),
		   make_option(c("-o", "--output_folder"), type='character', default=NULL, help='output folder path (end path with slash)', metavar='character')
)


print('Parsing arguments')
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

chr = opt$chromosome
genepred_cds_path = opt$genepred_cds_path
output_folder = opt$output_folder

#genepred_exons_path = 'output/CDS-coordinates/mm10.ncbiRefSeq.coding.chrY.gtf'
#chr = 'chrY'
#output_folder = 'CDS-information/'

genepred_cds = read.delim(genepred_cds_path, header=F)

### get gene id info
gene_id_info = genepred_cds[,4]
gene_id_unique = unique(gene_id_info)

gene_cds_information = list()
gene_boundaries_information = NULL
for (i in 1:length(gene_id_unique)){
	print(paste("gene", i, "/", length(gene_id_unique)))
	ind_gene_id = which(gene_id_info == gene_id_unique[i])
	genepred_gene_i = genepred_cds[ind_gene_id,]
	genepred_gene_i[,2] = format(genepred_gene_i[,2], scientific=FALSE)
	genepred_gene_i[,3] = format(genepred_gene_i[,3], scientific=FALSE)
	
	genepred_merged_i = genepred_gene_i
	
	##############
	#############
	
	gene_cds_information[[i]] = genepred_merged_i
	strand = unique(genepred_merged_i[,6])
	if (length(strand)>1){
		strand="NA"
	}
	
	gene_boundaries_i = data.frame("chr"=chr, "start"=min(as.numeric(genepred_merged_i[,2])), "end"=max(as.numeric(genepred_merged_i[,3])), "strand"=strand)
	gene_boundaries_information = rbind(gene_boundaries_information, gene_boundaries_i)
}

names(gene_cds_information) = gene_id_unique
gene_boundaries_information$gene_id = gene_id_unique

saveRDS(gene_cds_information, paste(output_folder, "gene_cds_information_", chr, ".RDS", sep=""))
saveRDS(gene_boundaries_information, paste(output_folder, "gene_boundaries_information_", chr, ".RDS", sep=""))

