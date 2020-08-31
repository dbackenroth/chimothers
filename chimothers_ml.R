source("chimothers.R")

TRIO_RELS <- c("Father", "Mother", "Child")

if (F) {
  vcf_file <- "~/Desktop/chimamas/merged_filtered.vcf"
  reads_file <- "~/Desktop/chimamas//all_read_counts.csv"
  processed_vcf <- ProcessVCF(vcf_file)
  processed_reads <- ProcessReadCounts(reads_file, processed_vcf)
}

GetTrios <- function(vcf) {

  vcf_trios <- vcf %>%
    filter(rel %in% TRIO_RELS)
  
  family_stats <- vcf_trios %>%
    group_by(family) %>%
    summarise(num_family_members = length(unique(sample)), 
              has_m_f_c = all(TRIO_RELS %in% rel))
  
  has_three <- family_stats %>%
    filter(has_m_f_c, 
           num_family_members == 3)
  
  return(has_three$family)
}

FilterToHomHomHet <- function(read_data, trios) {
  filtered1 <- FilterToHighConfidence(read_data) %>%
    group_by(CHROM, POS, family) %>%
    filter(n() == 3,
           family %in% trios,
           rel %in% TRIO_RELS)
  filtered2 <- filtered1 %>%
    filter(all(GT[rel %in% c("Mother", "Father")] %in% c("0/0", "1/1")))
  filtered3 <- filtered2 %>%
    filter(!GT[rel == "Mother"] == GT[rel == "Father"], 
           GT[rel == "Child"] == "0/1") %>%
    ungroup() %>%
    arrange(POS, name)
}

OffsetProbs <- function(ref, alt, gt, offset, dat) {
  stopifnot(offset >= 0)
  gg <- filter(dat, 
               REF == ref, ALT == alt, GT == gt)
  ref_prob_col <- paste0("prob_", ref)
  alt_prob_col <- paste0("prob_", alt)
  if (gt == "0/0"){
    gg[, ref_prob_col] <- gg[, ref_prob_col] - offset
    gg[, alt_prob_col] <- gg[, alt_prob_col] + offset
  } else {
    gg[, ref_prob_col] <- gg[, ref_prob_col] + offset
    gg[, alt_prob_col] <- gg[, alt_prob_col] - offset
  }
  return(gg)
}

CalculateMultinomialLogProbs <- function(offset, data, error) {
  with_error <- inner_join(data, error, by = c("CHROM", "POS", "GT", "REF", "ALT"))
  
  prob <- as.data.frame(with_error)[, c("sum_A", "sum_C", "sum_G", "sum_T")] %>%
    as.matrix() %>%
    `+`(1) %>%
    apply(1, function(x){x / sum(x)}) %>%
    t() %>%
    as.data.frame()
  colnames(prob) <- paste0("prob_", c("A", "C", "G", "T"))
  with_error <- cbind(with_error, prob) %>%
    arrange(CHROM, POS, family, GT)
  splits <- with_error %>%
    dplyr::select(ref = REF, alt = ALT, gt = GT) %>%
    distinct()
  offset_data <- pmap_dfr(splits, OffsetProbs, offset = offset, dat = with_error) %>%
    arrange(CHROM, POS, family, GT)
  x <- as.data.frame(offset_data)[, c("AA", "CC", "GG", "TT")] %>%
    as.matrix()
  offset_probs <- as.data.frame(offset_data)[, c("prob_A", "prob_C", "prob_G", "prob_T")] %>%
    as.matrix()
  lp <- dmultinomial(x = x, prob = offset_probs, log = T)
  offset_data$lp <- lp
  offset_data$offset <- offset
  return(offset_data)
}

Analysis <- function(read_data, vcf) {
  error_distribution <- GetErrorDistribution(read_data, c("Father", "Mother"))
  
  trios <- GetTrios(vcf)
  homhom <- FilterToHomHomHet(read_data, trios) %>%
    filter(!rel == "Child")
  homhom_sel <- homhom %>%
    transmute(CHROM, POS, family, rel, REF, ALT, GT, AA, CC, GG, TT, GT) 
  
  sample_error_rates <- GetSampleErrorRates(read_data %>%
                                              anti_join(homhom_sel %>%
                                                          dplyr::select(CHROM, POS, family, rel)))
  offset_lps <- map_dfr(seq(0, 0.005, by = 0.0001), 
                CalculateMultinomialLogProbs, data = homhom_sel, error = error_distribution)
  by_fam_summ <- offset_lps %>%
    group_by(offset, rel, family) %>%
    summarise(sum_lp = sum(lp))
  by_fam_p <- ggplot(by_fam_summ, aes(x = offset, y = sum_lp, col=rel, group = paste(family, rel))) + 
    geom_text(aes(label = family)) + facet_wrap(~rel) + xlim(0, 0.002)
  ggsave("Results/aug_31_by_family_ml_result.pdf")
  overall_summ <- offset_lps %>%
    group_by(offset, rel) %>%
    summarise(sum_lp = sum(lp))
  overall_p <- ggplot(overall_summ, aes(x = offset, y = sum_lp, col=rel)) + 
    geom_line()
  ggsave("Results/aug_31_overall_ml_result.pdf")
  parent_error_rates <- sample_error_rates %>%
    ungroup() %>%
    filter(family %in% by_fam_summ$family) %>%
    transmute(family, rel, GT, n, perc_right, perc_wrong, perc_oth) %>%
    as.data.frame()
  write.csv(parent_error_rates %>% 
              filter(rel == "Father") %>%
              arrange(desc(perc_wrong)), "Results/aug_31_father_error_rates.csv", 
            row.names = F)
  write.csv(parent_error_rates %>% 
              filter(rel == "Mother") %>%
              arrange(desc(perc_wrong)), "Results/aug_31_mother_error_rates.csv", 
            row.names = F)
}

GetSampleErrorRates <- function(read_data) {
  hc <- FilterToHighConfidence(read_data)
  errors <- hc %>%
    filter(GT %in% c("0/0", "1/1")) %>%
    mutate(num_right = if_else(GT == "0/0", num_ref, num_alt), 
           num_wrong = if_else(GT == "1/1", num_ref, num_alt)) %>%
    group_by(name, GT, family, rel) %>%
    summarise(n = n(), 
              tot_right = sum(num_right), 
              tot_wrong = sum(num_wrong), 
              tot_other = sum(AA + CC + GG + TT - num_right - num_wrong)) %>%
    mutate(total = tot_right + tot_wrong + tot_other,
           perc_right = tot_right / total, 
           perc_wrong = tot_wrong / total, 
           perc_oth = tot_other / total)
}

GetErrorDistribution <- function(read_data, rels_exclude = NULL) {
  hc <- FilterToHighConfidence(read_data)
  excluded <- hc %>%
    filter(GT %in% c("0/0", "1/1"), 
           !rel %in% rels_exclude)
  errors <- excluded %>%
    group_by(CHROM, POS, GT, REF, ALT) %>%
    summarise(sum_A = sum(AA), 
              sum_C = sum(CC), 
              sum_G = sum(GG), 
              sum_T = sum(TT))
}
