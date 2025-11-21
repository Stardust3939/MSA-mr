library(rtracklayer)
library(GenomicRanges)
library(vroom)

path_gwas <- "I:\\data\\GCST90406925\\"
gwas_hg19 <- vroom(file = paste0(path_gwas, "GCST90406925.tsv"))
head(gwas_hg19)
gwas_hg19['chr_str'] <- paste0("chr", gwas_hg19$chromosome)
head(gwas_hg19)

gwas_hg19_sel <- gwas_hg19

# IMPORTANT NOTE: We saw that sometimes if you did not add "strand" column info on your data, some bp_hg38 positions will add with +1 bp in the created hg38 file !!! So, we added this column to the data as well.
gwas_hg19_sel["strand"] <- "+"

# Check data types of gwas_hg19_sel
str(gwas_hg19_sel)

# Use this section ONLY if you have "NA" data in "bp_hg19" column
# Check for missing values in "bp_hg19" column. If so, it then removes them
#sum(is.na(gwas_hg19_sel$bp_hg19))
#gwas_hg19_sel <- gwas_hg19_sel[!is.na(gwas_hg19_sel$bp_hg19),]
#sum(is.na(gwas_hg19_sel$bp_hg19))

# It was downloaded from https://hgdownload.cse.ucsc.edu/goldenpath/hg38/liftOver/
chain <- import.chain("C:\\Users\\Jerry\\Downloads\\hg38ToHg19.over.chain")

gr <- makeGRangesFromDataFrame(gwas_hg19_sel, ignore.strand = TRUE, seqnames.field = "chr_str", start.field = "base_pair_location", end.field = "base_pair_location")

hg38 <- liftOver(gr, chain)
hg38_df <- as.data.frame(hg38)

gwas_hg19_sel$rownames_col <- as.numeric(rownames(gwas_hg19_sel))
head(gwas_hg19_sel)
match_pos <- match(gwas_hg19_sel$rownames_col, hg38_df$group)

# Add the corresponding information from hg38_df to gwas_hg19_sel
gwas_hg19_sel$group[!is.na(match_pos)] <- hg38_df$group[match_pos[!is.na(match_pos)]]
gwas_hg19_sel$seqnames[!is.na(match_pos)] <- hg38_df$seqnames[match_pos[!is.na(match_pos)]]
gwas_hg19_sel$bp_hg38[!is.na(match_pos)] <- hg38_df$start[match_pos[!is.na(match_pos)]]
head(gwas_hg19_sel)

# Check for missing values in the "bp_hg38" column of "gwas_hg19_sel" GWAS file, then remove them
sum(is.na(gwas_hg19_sel$bp_hg38))
gwas_hg19_sel <- gwas_hg19_sel[!is.na(gwas_hg19_sel$bp_hg38),]
sum(is.na(gwas_hg19_sel$bp_hg38))
# rename column bp_hg38 to bp_hg19
colnames(gwas_hg19_sel)[colnames(gwas_hg19_sel) == "bp_hg38"] <- "bp_hg19"
# Write the output to a file
vroom_write(gwas_hg19_sel, file = "I:\\data\\GCST90406925\\GCST90406925_hg19.txt", delim = "\t", quote = "none", col_names = TRUE)