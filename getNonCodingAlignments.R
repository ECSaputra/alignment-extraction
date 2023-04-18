require("rphast")
require(optparse)

option_list = list(make_option(c('-c', '--chromosome'), type='character', default=NULL, help='chromosome of interest', metavar='character'),
		   make_option(c('-r', '--refseq'), type='character', default=NULL, help='reference species', metavar='character'),
		   make_option(c('-i', '--elementBedPath'), type='character', default=NULL, help='path to BED file containing the coordinates of regions of interest', metavar='character'),
		   make_option(c('-o', '--output_folder'), type='character', default=NULL, help='output folder for non-coding region alignments', metavar='character'),
		   make_option(c('-a', '--mafFolderPath'), type='character', default=NULL, help='path to folder containing chromosome MAF fragments (end path with slash)', metavar='character'),
		   make_option(c('-p', '--prefix'), type='character', default=NULL, help='prefix of MAF fragments', metavar='character'))

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

chromosome = opt$chromosome
elementBedPath = opt$elementBed
mafFolderPath = opt$mafFolderPath
refseq = opt$refseq
output_folder = opt$output_folder
split_alignment_prefix = opt$prefix

maf_files = list.files(mafFolderPath)
regions = read.delim(elementBedPath, header=F)
if (ncol(regions)==3){
	        regions[,4] = paste0('region_', seq(1,nrow(regions),1))
}

regions.df = as.data.frame(regions)
regions.df[,2] = as.integer(regions.df[,2])
regions.df[,3] = as.integer(regions.df[,3])

regions_sorted = regions.df[order(regions.df[,2], decreasing=F),]

`%notin%` = Negate(`%in%`)

# sort the mafs
splitnum = sub(split_alignment_prefix, '', sub('.maf', '', maf_files))
splitnum_sort = sort(as.integer(splitnum), decreasing=F)


CNE_count = 0

for (i in 1:length(splitnum_sort)){
	print(paste('Scanning fragment',i, '/', length(splitnum_sort)))
	mafPath = paste0(mafFolderPath, split_alignment_prefix, splitnum_sort[i], '.maf')
	
	maf = read.msa(mafPath, pointer.only=T)
	coord_range = coord.range.msa(maf)
	maf_start = coord_range[1]
	maf_end = coord_range[2]
	
	ind_CNE_start_in_maf = intersect(which(regions_sorted[,2] >= maf_start), which(regions_sorted[,2] <= maf_end))
	
	if (length(ind_CNE_start_in_maf)==0){
		print('No CNE in this fragment')
	}
	if (length(ind_CNE_start_in_maf) > 0){
		ind_CNE_end_in_maf = intersect(which(regions_sorted[,3] >= maf_start), which(regions_sorted[,3] <= maf_end))
		ind_CNE_complete_in_maf = intersect(ind_CNE_start_in_maf, ind_CNE_end_in_maf)
		ind_CNE_partial_in_maf = ind_CNE_start_in_maf[which(ind_CNE_start_in_maf %notin% ind_CNE_end_in_maf)]
		
		if (length(ind_CNE_complete_in_maf) > 0){
			complete_CNE_in_maf = regions_sorted[ind_CNE_complete_in_maf,]
			
			for (j in 1:nrow(complete_CNE_in_maf)){
				CNE_alignment = sub.msa(maf, start.col=complete_CNE_in_maf[j,2], end.col=complete_CNE_in_maf[j,3], refseq=refseq)
				msa_name = paste0(output_folder, complete_CNE_in_maf[j,4], '.fa')
				write.msa(CNE_alignment, file=msa_name, format="FASTA")
				
				CNE_count = CNE_count + 1
				print(paste("--- Completed", CNE_count, "/", nrow(regions_sorted), "CNEs"))
			}
		}

		if (length(ind_CNE_partial_in_maf) > 0){
			partial_CNE_in_maf = regions_sorted[ind_CNE_partial_in_maf,]
			
			for (j in 1:nrow(partial_CNE_in_maf)){
				### tidy up CNE that lie across two mafs
				corrected_CNE_starts = c(partial_CNE_in_maf[j,2], maf_end+1)
				corrected_CNE_ends = c(maf_end, partial_CNE_in_maf[j,3])
				corrected_CNE = data.frame("chr"=c(chr, chr), "start"=corrected_CNE_starts, "end"=corrected_CNE_ends)
				CNE_alignment = sub.msa(maf, start.col=corrected_CNE[1,2], end.col=corrected_CNE[1,3], refseq=refseq)
				mafPlusPath = paste(mafFolderPath, split_alignment_prefix, splitnum_sort[i+1], ".maf", sep="")
				mafplus = read.msa(mafPlusPath)

				fragment_alignment = sub.msa(mafplus, start.col=corrected_CNE[2,2], end.col=corrected_CNE[2,3], refseq=refseq)
				CNE_alignment = concat.msa(list(CNE_alignment, fragment_alignment))
				msa_name = paste0(output_folder, partial_CNE_in_maf[j,4], '.fa')
				write.msa(CNE_alignment, file=msa_name, format="FASTA")

				CNE_count = CNE_count + 1
				print(paste("--- Completed", CNE_count, "/", nrow(regions_sorted), "CNEs"))
			}
		}
	}
	if (CNE_count == nrow(regions_sorted)){
		break
	}
}

