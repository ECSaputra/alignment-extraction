library(optparse)

option_list = list(make_option(c("-i", "--inputpath"), type='character', default=NULL, help='path to input GTF file', metavar='character'),
		   make_option(c("-o", "--outputpath"), type='character', default=NULL, help='path to output GTF file', metavar='character')
)

print('Parsing arguments')
opt_parser = OptionParser(option_list=option_list)
opt=parse_args(opt_parser)

gtf_input_path = opt$inputpath
outPath = opt$outputpath
#gtf_input_path = 'data/genepred/mm10.ncbiRefSeq.chrY.gtf'
#outPath = 'output/CDS-coordinates/mm10.ncbiRefSeq.coding.chrY.gtf'


all_gtf = read.delim(gtf_input_path, header=F)

### Isolate CDS and stop codons
retain_cats = c("CDS", "stop_codon")
ind_retain = which(all_gtf[,3] %in% retain_cats)
coding_gtf = all_gtf[ind_retain,]

### get gene names
region_info = strsplit(coding_gtf[,9], ";")
gene_name = NULL
gene_id = NULL
transcript_id = NULL
for (i in 1:length(region_info)){
	ind_gene_name = which(grepl("gene_name", region_info[[i]]))
	ind_gene_id = which(grepl("gene_id", region_info[[i]]))
	ind_transcript_id = which(grepl("transcript_id", region_info[[i]]))
	
	gene_name = c(gene_name, sub(" gene_name ", "", region_info[[i]][ind_gene_name]))
	gene_id = c(gene_id, sub("gene_id ", "", region_info[[i]][ind_gene_id]))
	transcript_id = c(transcript_id, sub(" transcript_id ", "", region_info[[i]][ind_transcript_id]))
}

### construct the CDS table

unq_gene_ids = unique(gene_id)

exon_bed = NULL
for (i in 1:length(unq_gene_ids)){
	print(paste('gene', i, "/", length(unq_gene_ids)))
	ind_gene = which(gene_id == unq_gene_ids[i])
	gene_gtf = coding_gtf[ind_gene,]
	transcripts_i = transcript_id[ind_gene]
	
	unq_transcripts_id = unique(transcripts_i)
	ind_NM = which(grepl("NM", unq_transcripts_id))
	
	if (length(ind_NM) > 1){
		NM_transcripts_id = unq_transcripts_id[ind_NM]
		ind_transcript_id = which(transcripts_i %in% NM_transcripts_id)
		transcripts_gtf = gene_gtf[ind_transcript_id,]
		
		transcript_lengths = rep(0, length(NM_transcripts_id))
		for (j in 1:length(NM_transcripts_id)){
			ind_transcript_j = which(transcripts_i == NM_transcripts_id[j])
			transcript_j = gene_gtf[ind_transcript_j,]
			CDS_lengths = as.numeric(transcript_j[,5]) - as.numeric(transcript_j[,4]) + 1
			transcript_lengths[j] = sum(CDS_lengths) - 3
		}
		use_transcript_id = unq_transcripts_id[which(transcript_lengths == max(transcript_lengths))]
		ind_transcript_id = which(transcripts_i == use_transcript_id[1])
		gene_gtf = gene_gtf[ind_transcript_id,]
	} else if (length(ind_NM) == 1){
		ind_transcript_id = which(transcripts_i == unq_transcripts_id[ind_NM])
		gene_gtf = gene_gtf[ind_transcript_id,]
	} else if (length(ind_NM == 0)){
		transcript_lengths = rep(0, length(unq_transcripts_id))
		for (j in 1:length(unq_transcripts_id)){
			ind_transcript_j = which(transcripts_i == unq_transcripts_id[j])
			transcript_j = gene_gtf[ind_transcript_j,]
			CDS_lengths = as.numeric(transcript_j[,5]) - as.numeric(transcript_j[,4]) + 1
			transcript_lengths[j] = sum(CDS_lengths) - 3
		}
		use_transcript_id = unq_transcripts_id[which(transcript_lengths == max(transcript_lengths))]
		ind_transcript_id = which(transcripts_i == use_transcript_id[3])
		gene_gtf = gene_gtf[ind_transcript_id,]
	}
	
	gene_df = data.frame("chr"=gene_gtf[,1], "start"=as.numeric(gene_gtf[,4]), "end"=as.numeric(gene_gtf[,5]), "gene"=rep(toupper(unq_gene_ids[i]), nrow(gene_gtf)), "score"=gene_gtf[,8], "strand"=gene_gtf[,7])
	gene_df = gene_df[order(gene_df$start, decreasing=F),]

	exon_bed = rbind(exon_bed, gene_df)
}

write.table(exon_bed, outPath, quote=F, row.names=F, col.names=F, sep="\t")

