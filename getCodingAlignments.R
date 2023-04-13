library(optparse)
require("rphast")

option_list = list(make_option(c('-c', '--chromosome'), type='character', default=NULL, help='chromosome of interest', metavar='character'),
		   make_option(c('-r', '--refseq'), type='character', default=NULL, help='reference sequence', metavar='character'),
		   make_option(c('-i', '--cds_info_folder'), type='character', default=NULL, help='CDS information folder', metavar='character'),
		   make_option(c('-o', '--output_folder'), type='character', default=NULL, help='output folder for coding region alignments', metavar='character'),
		   make_option(c('-a', '--mafFolderPath'), type='character', default=NULL, help='path to folder containing chromosome MAF fragments (end path with slash)', metavar='character'),
		   make_option(c('-p', '--prefix'), type='character', default=NULL, help='prefix of MAF fragments', metavar='character')
)

print('Parsing arguments')
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

chr = opt$chromosome
refseq = opt$refseq
output_folder = opt$output_folder
cds_info_folder = opt$cds_info_folder
mafFolderPath = opt$mafFolderPath
prefix = opt$prefix

#chr= 'chrY'
#refseq = 'mm10'
#output_folder = 'output/coding-region-alignment/chrY/'
#cds_info_folder = 'output/CDS-information/'
#mafFolderPath = 'data/alignment/chrY/'
#prefix = 'chrYsplit'

gene_boundaries_information = readRDS(paste0(cds_info_folder, "gene_boundaries_information_", chr, ".RDS"))
gene_cds_information = readRDS(paste0(cds_info_folder, "gene_cds_information_", chr, ".RDS"))

gene_boundaries_sorted = gene_boundaries_information[order(as.numeric(gene_boundaries_information$start), decreasing=F),]




`%notin%` = Negate(`%in%`)

# sort the mafs
maf_files = list.files(mafFolderPath)

#splitnum = sub(chr, "", sub("split", "", sub(".maf", "", maf_files)))
splitnum = sub(prefix, "", sub(".maf", "", maf_files))
splitnum_sort = sort(as.integer(splitnum), decreasing=F)

gene_count = 0

for (i in 1:length(splitnum_sort)){
	mafPath = paste0(mafFolderPath, prefix, splitnum_sort[i], ".maf")
	maf = read.msa(mafPath)
	coord_range = coord.range.msa(maf)
	maf_start = coord_range[1]
	maf_end = coord_range[2]
	
	ind_gene_start_in_maf = intersect(which(as.numeric(gene_boundaries_sorted[,2]) >= maf_start), which(as.numeric(gene_boundaries_sorted[,2]) <= maf_end))
	if (length(ind_gene_start_in_maf) > 0){
		ind_gene_end_in_maf = intersect(which(as.numeric(gene_boundaries_sorted[,3]) >= maf_start), which(as.numeric(gene_boundaries_sorted[,3]) <= maf_end))
		ind_gene_complete_in_maf = intersect(ind_gene_start_in_maf, ind_gene_end_in_maf)
		ind_gene_partial_in_maf = ind_gene_start_in_maf[which(ind_gene_start_in_maf %notin% ind_gene_end_in_maf)]
		
		if (length(ind_gene_complete_in_maf) > 0){
			complete_gene_in_maf = gene_boundaries_sorted[ind_gene_complete_in_maf,]
			
			for (j in 1:nrow(complete_gene_in_maf)){
				gene_j = complete_gene_in_maf$gene_id[j]
				strand = complete_gene_in_maf$strand[j]
				exons_in_maf = gene_cds_information[gene_j][[1]]
				
				for (k in 1:nrow(exons_in_maf)){
					if (k == 1){
						gene_alignment = sub.msa(maf, start.col=as.numeric(exons_in_maf[k,2]), end.col=as.numeric(exons_in_maf[k,3]), refseq=refseq)
					} else {
						exon_alignment = sub.msa(maf, start.col=as.numeric(exons_in_maf[k,2]), end.col=as.numeric(exons_in_maf[k,3]), refseq=refseq)
						gene_alignment = concat.msa(list(gene_alignment, exon_alignment))
					}
				}

				if (strand == "-"){
					gene_alignment = reverse.complement.msa(gene_alignment)
				}


				msa_name = paste0(output_folder, toupper(gene_j), ".fa")
				write.msa(gene_alignment, file=msa_name, format="FASTA")

				gene_count = gene_count + 1
				print(paste("Completed", gene_count, "/", nrow(gene_boundaries_sorted), "genes"))
			}
		}
		
		
		if (length(ind_gene_partial_in_maf) > 0){
			partial_gene_in_maf = gene_boundaries_sorted[ind_gene_partial_in_maf,]
			
			for (j in 1:nrow(partial_gene_in_maf)){
				gene_j = partial_gene_in_maf$gene_id[j]
				strand = partial_gene_in_maf$strand[j]
				gene_exons_j = gene_cds_information[gene_j][[1]]

				### tidy up exons that lie across two mafs
				ind_split_exon = intersect(intersect(which(as.numeric(gene_exons_j[,2]) >= maf_start), which(as.numeric(gene_exons_j[,2]) <= maf_end)), which(as.numeric(gene_exons_j[,3]) > maf_end))
				if (length(ind_split_exon) > 0){
					split_exon = gene_exons_j[ind_split_exon,]
					if (ind_split_exon == nrow(gene_exons_j)){
						corrected_exon_starts = c(gene_exons_j[1:ind_split_exon,2] , maf_end+1)
					} else {
						corrected_exon_starts = c(gene_exons_j[1:ind_split_exon,2], maf_end+1, gene_exons_j[(ind_split_exon+1):nrow(gene_exons_j),2])
					}

					if (ind_split_exon == 1){
						corrected_exon_ends = c(maf_end, gene_exons_j[ind_split_exon:nrow(gene_exons_j),3])
					} else {
						corrected_exon_ends = c(gene_exons_j[1:(ind_split_exon-1),3], maf_end, gene_exons_j[ind_split_exon:nrow(gene_exons_j),3])
					}
					gene_exons_j = data.frame("chr"=rep(chr, length(corrected_exon_starts)), "start"=corrected_exon_starts, "end"=corrected_exon_ends)
				}

				### find which exons lie in this current maf
				ind_exons_in_maf = which(as.numeric(gene_exons_j[,3]) <= maf_end)
				exons_in_maf = gene_exons_j[ind_exons_in_maf,]

				for (k in 1:nrow(exons_in_maf)){
					if (k == 1){
						gene_alignment = sub.msa(maf, start.col=as.numeric(exons_in_maf[k,2]), end.col=as.numeric(exons_in_maf[k,3]), refseq=refseq)
					} else {
						exon_alignment = sub.msa(maf, start.col=as.numeric(exons_in_maf[k,2]), end.col=as.numeric(exons_in_maf[k,3]), refseq=refseq)
						gene_alignment = concat.msa(list(gene_alignment, exon_alignment))
					}
				}

				######## find exons lying in other mafs
				ind_exons_notin_maf = which(as.numeric(gene_exons_j[,3]) > maf_end)
				exons_notin_maf = gene_exons_j[ind_exons_notin_maf,]

				remaining_exons_count = nrow(exons_notin_maf)
				additional_maf_count = 1
				while (remaining_exons_count > 0){
					mafPlusPath = paste0(mafFolderPath, prefix, splitnum_sort[i+additional_maf_count], ".maf")
					mafplus = read.msa(mafPlusPath)

					coord_range_plus = coord.range.msa(mafplus)
					maf_start_plus = coord_range_plus[1]
					maf_end_plus = coord_range_plus[2]

					### tidy up exons that lie across two mafs
					ind_split_exon = intersect(intersect(which(as.numeric(exons_notin_maf[,2]) >= maf_start_plus), which(as.numeric(exons_notin_maf[,2]) <= maf_end_plus)),
								   which(as.numeric(exons_notin_maf[,3]) > maf_end_plus))
					
					if (length(ind_split_exon) > 0){
						split_exon = exons_notin_maf[ind_split_exon,]
						if (ind_split_exon == nrow(exons_notin_maf)){
							corrected_exon_starts = c(exons_notin_maf[1:ind_split_exon,2], maf_end_plus+1)
						} else {
							corrected_exon_starts = c(exons_notin_maf[1:ind_split_exon,2], maf_end_plus+1, exons_notin_maf[(ind_split_exon+1):nrow(exons_notin_maf),2])
						}

						if (ind_split_exon == 1){
							corrected_exon_ends = c(maf_end_plus, exons_notin_maf[ind_split_exon:nrow(exons_notin_maf),3])
						} else {
							corrected_exon_ends = c(exons_notin_maf[1:(ind_split_exon-1),3], maf_end_plus, exons_notin_maf[ind_split_exon:nrow(exons_notin_maf),3])
						}
						exons_notin_maf = data.frame("chr"=rep(chr, length(corrected_exon_starts)), "start"=corrected_exon_starts, "end"=corrected_exon_ends)
					}

					### find which exons lie in this current maf
					ind_exons_in_mafplus = intersect(which(as.numeric(exons_notin_maf[,2]) >= maf_start_plus),which(as.numeric(exons_notin_maf[,3]) <= maf_end_plus))

					if (length(ind_exons_in_mafplus) > 0){
						exons_in_mafplus = exons_notin_maf[ind_exons_in_mafplus,]

						for (k in 1:nrow(exons_in_mafplus)){
							exon_alignment = sub.msa(mafplus, start.col=as.numeric(exons_in_mafplus[k,2]), end.col=as.numeric(exons_in_mafplus[k,3]), refseq=refseq)
							gene_alignment = concat.msa(list(gene_alignment, exon_alignment))
							remaining_exons_count = remaining_exons_count - 1
							if (remaining_exons_count > 0 && k == nrow(exons_in_mafplus)){
								additional_maf_count = additional_maf_count + 1
							}
						}
					} else {
						additional_maf_count = additional_maf_count + 1
					}
				}

				if (strand == "-"){
					gene_alignment = reverse.complement.msa(gene_alignment)
				}
				
				msa_name = paste0(output_folder,toupper(gene_j), ".fa", sep="")
				write.msa(gene_alignment, file=msa_name, format="FASTA")

				gene_count = gene_count + 1
				print(paste("Completed", gene_count, "/", nrow(gene_boundaries_sorted), "genes"))
			}
		}
	}

}
