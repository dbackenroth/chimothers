library(data.table)
library(stringr)
library(purrr)
library(dplyr)

ProcessPileups <- function() {
  pileups <- Sys.glob("*.pileup")
  read_data <- map_dfr(pileups, function(x){ProcessPileup(x) %>% mutate(name = x)})
  write.csv(read_data, "all_read_data.csv", row.names = F, quote = F)
  return(read_data)
}

ProcessPileup <- function(pileup.file, min_mq){
  p.info <- fread(pileup.file)
  #p.info <- p.info[1:50, ]
  setnames(p.info, c("CHROM", "POS", "REF", "NUMSPANNINGREADS", "PILEUP", "MQ"))
  ins.pattern <- "\\+[0-9]+[ACGTNacgtn]+"
  del.pattern <- "-[0-9]+[ACGTNacgtn]+"
  
  p.info[, `:=`(Alleles=gsub("\\^.", "", PILEUP))]     # start of read marker
  p.info[, `:=`(Alleles=gsub("\\$", "", Alleles))] # end of read marker
  p.info[, `:=`(Alleles=gsub(ins.pattern, "", Alleles))] # insertions
  p.info[, `:=`(Alleles=gsub(del.pattern, "", Alleles))] # deletions
  p.info[, `:=`(Alleles=gsub("*", "", Alleles, fixed=T))] # continuation of deletion
  p.info[, `:=`(Alleles=gsub("a|A", "A", Alleles))] # alternate alleles
  p.info[, `:=`(Alleles=gsub("c|C", "C", Alleles))] # alternate alleles
  p.info[, `:=`(Alleles=gsub("g|G", "G", Alleles))] # alternate alleles
  p.info[, `:=`(Alleles=gsub("t|T", "T", Alleles))] # alternate alleles
  p.info[, `:=`(Alleles=gsub(",|\\.", "R", Alleles))] # reference alleles
  n_snps <- nrow(p.info)
  mqs <- map(asc(p.info$MQ), function(x){x - 33})
  filtered_alleles <- vector(mode = "character", n_snps)
  for (i in 1:n_snps) {
    pass_filter <- which(mqs[[i]] >= min_mq)
    alleles_split <- strsplit(p.info$Alleles[i], "")[[1]]
    alleles_kept <- alleles_split[pass_filter]
    filtered_alleles[[i]] <- paste0(alleles_kept, collapse = "")
  }
  p.info$filtered_alleles <- filtered_alleles
  p.info[, `:=`(R=str_count(filtered_alleles, pattern="R"), 
                A=str_count(filtered_alleles, pattern="A"), 
                C=str_count(filtered_alleles, pattern="C"), 
                G=str_count(filtered_alleles, pattern="G"),
                T=str_count(filtered_alleles, pattern="T"))]
  p.info <- p.info[, .(CHROM, POS, REF, R, A, C, G, T)]
  return(p.info)
}