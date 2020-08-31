library(data.table)
library(tidyverse)
library(broom)

FilterToHighConfidence <- function(reads_dat) {
  filtered <- reads_dat %>%
    filter(num_ref_alt >= 200) %>%
    filter(case_when(GT == "0/1" ~ num_ref / num_ref_alt > 0.3 & num_alt / num_ref_alt > 0.3, 
                     GT == "0/0" ~ num_alt / num_ref_alt <= 0.025 & num_alt <= 20,
                     GT == "1/1" ~ num_ref / num_ref_alt <= 0.025 & num_ref <= 20))
  return(filtered)
}

ProcessReadCounts <- function(reads_file, processed_vcf) {
  alt_ref <- processed_vcf %>%
    select(CHROM, POS, VCF_REF = REF, VCF_ALT = ALT) %>%
    as.data.frame() %>%
    distinct()
  
  reads_dat <- fread(reads_file) %>%
    dplyr::select(-R) %>%
    mutate(name = gsub(".sam.sorted.bam.base.counts.txt", "", name, fixed = T), 
           REF = toupper(REF)) %>%
    rename(AA=A, CC=C, GG=G, TT=T) %>%
    as.data.frame()
  
  with_alt_ref <- reads_dat %>%
    inner_join(alt_ref, by = c("CHROM", "POS")) %>%
    mutate(num_ref = case_when(REF == "A" ~ AA, 
                               REF == "C" ~ CC,
                               REF == "G" ~ GG,
                               REF == "T" ~ TT),
           num_alt = case_when(VCF_ALT == "A" ~ AA, 
                               VCF_ALT == "C" ~ CC, 
                               VCF_ALT == "G" ~ GG, 
                               VCF_ALT == "T" ~ TT))
  joined_dat <- with_alt_ref %>%
    select(CHROM, POS, name, num_ref, num_alt, AA, CC, GG, TT) %>%
    inner_join(processed_vcf, by = c("CHROM", "POS", name = "sample")) %>%
    select(-FORMAT, -info) %>%
    mutate(num_ref_alt = num_ref + num_alt)
  return(joined_dat)
}

ProcessVCF <- function(vcf_file) {
  BASES <- c("A", "C", "G", "T")
  
  # ignore deletions, insertions etc...
  vcf <- fread(vcf_file, skip = "#CHROM") %>%
    filter(REF %in% BASES, 
           ALT %in% BASES) %>%
    rename(CHROM = `#CHROM`) %>%
    mutate(CHROM = paste0("chr", CHROM))
  
  vcf_long <- vcf %>%
    select(-ID, -QUAL, -FILTER, -INFO) %>%
    gather(sample, info, -CHROM, -POS, -REF, -ALT, -FORMAT) %>%
    as.data.table()
  
  REL_DICT <- c(ChA = "Child", CHA = "Child", ChC = "Child", 
                FA = "Father", FC = "Father", jMC = "Mother", 
                GFFC = "GP", GFFH = "GP", GMFC = "GP", 
                kGFMC = "GP", kGFMH = "GP", lGMMC = "GP", lGMMH = "GP")
  
  vcf_long2 <- copy(vcf_long)
  
  vcf_long2[FORMAT == "GT:AD:DP:GQ:PL", c("GT", "AD", "DP", "GQ", "PL") := tstrsplit(info, ":", fixed = T)]
  vcf_long2[FORMAT == "GT:DP:AD:GQ:PL", c("GT", "DP", "AD", "GQ", "PL") := tstrsplit(info, ":", fixed = T)]
  vcf_long2[DP == '.', `:=`(DP = '0')]
  vcf_long2[, `:=`(DP = as.numeric(DP))]
  vcf_long2[, c("family", "rel2") := tstrsplit(sample, "-", keep = c(1,2))]
  vcf_long2[, rel := REL_DICT[rel2]]
  return(vcf_long2)
}




# when see both alleles in child, should be more likely to see both alleles in mother than father

# GetCounts <- function() {
#   
#   BASES <- c("A", "C", "G", "T")
#   
#   vcf <- fread("merged_filtered.vcf", skip = "#CHROM") %>%
#     filter(REF %in% BASES, 
#            ALT %in% BASES, 
#            `#CHROM` %in% 7) %>%
#     rename(CHROM = `#CHROM`) %>%
#     mutate(CHROM = paste0("chr", CHROM))
#   
#   alt_ref <- vcf %>%
#     select(CHROM, POS, VCF_REF = REF, VCF_ALT = ALT)
#   
#   dat <- fread("all_read_data.csv") %>%
#     rename(a=A, c=C, g=G, t=T) %>%
#     mutate(name = gsub(".sam.sorted.bam.pileup", "", name, fixed = T), 
#            REF = toupper(REF)) %>%
#     inner_join(alt_ref, by = c("CHROM", "POS")) %>%
#     mutate(num_ref = R, 
#            num_alt = case_when(VCF_ALT == "A" ~ a, 
#                                VCF_ALT == "C" ~ c, 
#                                VCF_ALT == "G" ~ g, 
#                                VCF_ALT == "T" ~ t), 
#            num_not_alt = case_when(VCF_ALT == "A" ~ c+g+t, 
#                                    VCF_ALT == "C" ~ a+g+t,
#                                    VCF_ALT == "G" ~ a+c+t, 
#                                    VCF_ALT == "T" ~ a+c+g),
#            num_tot = num_ref + num_alt,
#            tot_ref_other = paste0(num_tot, "_", num_ref, "_", num_not_alt))
#   dat <- as.data.table(dat)[, c("family", "rel") := tstrsplit(name, "-", keep = c(1,2))] %>%
#     as.data.frame() %>%
#     mutate(rels = case_when(rel %in% c("ChA", "CHA", "ChC") ~ "Child", 
#                             rel %in% c("FA", "FC") ~ "Father", 
#                             rel %in% c("jMC") ~ "Mother", 
#                             TRUE ~ "Grandparent"))
#   table(dat$family, dat$rels)   
#   #11, #15, #16, #18, #19, #25, #26, #3, #4, #5, #7, #8, #9
#   #c(3,4,5,7,8,9,11,15,16,18,19,25,26)
#   spread <- dat %>% 
#     select(CHROM, POS, name, tot_ref_other) %>%
#     spread(name, tot_ref_other)
#   
#   vcf_processed <- vcf %>%
#     select(-ID, -REF, -ALT, -QUAL, -FILTER, -INFO, -FORMAT) %>%
#     mutate_at(vars(-c("CHROM", "POS")), ~substring(., 1, 3))
#   
#   merged <- inner_join(spread, vcf_processed, by = c("CHROM", "POS"))
# }
# 
# DoFamily <- function(family_number, counts) {
#   this_family <- counts[, c("CHROM", "POS", 
#                             colnames(counts)[grepl(paste0("^", family_number, "-"), colnames(counts))])]
#   cols <- colnames(this_family)
#   child_cols <- grep("Ch|CH", cols)
#   mother_cols <- grep("MA|MC", cols)
#   father_cols <- grep("FA|FC", cols)
#   
#   stopifnot(all(c(3,6) %in% child_cols) &
#             all(c(4,7) %in% father_cols) &
#             all(c(5,8) %in% mother_cols))
#   colnames(this_family) <- c("CHROM", "POS", "cc", "fc", "mc", "cg", "fg", "mg")
#   
#   this_family1 <- this_family %>%
#     mutate(m_err = (cg == "1/1" & (mg == "0/0" | fg == "0/0")) |
#              (cg == "0/0" & (mg == "1/1" | fg == "1/1")) |
#              (cg == "0/1" & ((mg == "0/0" & fg == "0/0") | (mg == "1/1") & fg == "1/1")))
#   keep <- this_family1 %>%
#     filter(!m_err, 
#            !cg == "./.", 
#            !mg == "./.", 
#            !fg == "./.") %>%
#     separate(cc, into = c("cs", "cr", "co"), sep = "_", convert = T) %>%
#     separate(fc, into = c("fs", "fr", "fo"), sep = "_", convert = T) %>%
#     separate(mc, into = c("ms", "mr", "mo"), sep = "_", convert = T) %>%
#     filter(cs >= 1000, fs >= 1000, ms >= 1000)
# 
#   child_both <- keep %>%
#     filter(cr < cs & cr > 0) %>%
#     mutate(mother_both = mr < ms & mr > 0, 
#            father_both = fr < fs & fr > 0)
#   browser()
#   summ <- child_both %>%
#     summarise(mother_both = mean(mother_both), 
#               father_both = mean(father_both))
#   return(summ)
#   mother_data <- keep %>%
#     dplyr::select(CHROM, POS, sum = ms, ref = mr, cg, pg = mg) %>%
#     mutate(parent = "mother")
#   father_data <- keep %>%
#     dplyr::select(CHROM, POS, sum = fs, ref = fr, cg, pg = fg) %>%
#     mutate(parent = "father")
#   
#   both <- bind_rows(mother_data, father_data) %>%
#     filter(pg %in% c("0/0", "1/1")) %>%
#     mutate(child_other = pg == "0/0" & !cg == "0/0" |
#                          pg == "1/1" & !cg == "1/1",
#            other_count = if_else(pg == "0/0", sum - ref, ref), 
#            same_count = if_else(pg == "0/0", ref, sum - ref),
#            reference_bias = pg == "1/1")
#   
#   father_mod <- glm(cbind(other_count, same_count) ~ reference_bias + child_other, 
#              both %>% filter(parent == "father", !POS == 118387469), 
#              family = binomial) %>%
#     tidy(exponentiate = T) %>%
#     mutate(parent = "father")
#   browser()
#   mother_mod <- glm(other_count>0 ~ reference_bias + child_other, 
#                     both %>% filter(parent == "mother"), 
#                     family = binomial) %>%
#     tidy(exponentiate = T) %>%
#     mutate(parent = "mother")
#   return(bind_rows(father_mod, mother_mod) %>% mutate(family = family_number) %>%
#            filter(term == "child_otherTRUE"))
#   
#   mdiff <- keep %>%
#     filter(mg == "0/0" & !cg == "0/0" |
#              mg == "1/1" & !cg == "1/1") %>%
#     mutate(type = "mother")
#   
#   fdiff <- keep %>%
#     filter(fg == "0/0" & !cg == "0/0" |
#              fg == "1/1" & !cg == "1/1") %>%
#     mutate(type = "father")
#   
#   both <- bind_rows(mdiff, fdiff) %>%
#     mutate(parent_genotype = if_else(type == "mother", mg, fg), 
#            total_count = if_else(type == "mother", ms, fs), 
#            ref_count = if_else(type == "mother", mr, fr), 
#            error_count = if_else(type == "mother", mo, fo),
#            other_count = if_else(parent_genotype == "0/0", total_count - ref_count, ref_count), 
#            reference_bias = parent_genotype == "1/1")
#   
#   in_both <- intersect(mdiff$POS, fdiff$POS)
#   if (length(in_both) == 0) return(NULL)
#   summ <- both %>% 
#     group_by(type) %>%
#     summarise(n = n(), 
#               mean_child_allele_prop = mean(other_count / (total_count + error_count)), 
#               mean_error_allele_prop = mean(error_count / (total_count + error_count))) %>%
#     mutate(family = family_number) %>%
#     ungroup()
#   
#   print(family_number)
#   mod <- glm(cbind(other_count, total_count - other_count) ~ type + reference_bias, 
#              both, 
#              family = binomial) %>%
#     tidy()
#   mod <- filter(mod, term == "typemother")
#   summ$OR <- exp(mod$estimate[1])
#   summ$p.value <- mod$p.value[1]
#   return(summ)
# }
# # logistic regression 
# # control for reference bias etc...
# 
# #p <- ggplot(both, aes(x = type, y = q)) + geom_violin()
# #print(p)
# 
# #counts <- GetCounts()
# # 7 has none of the right SNPs
# family_numbers <- c(3,4,5,7,8,9,11,15,16,18,19,25,26)
# 
# res <- map_dfr(family_numbers, DoFamily, counts = counts)
# if (F) {
# ggplot(res, aes(x=type, y=mean_child_allele_prop)) + 
#   geom_violin() + 
#   ylab("Mean child allele proportion") + 
#   theme_bw()
# 
# res %>%
#   dplyr::select(n, family, OR, p.value) %>%
#   group_by(family) %>%
#   mutate(total_n_SNPs = sum(n)) %>%
#   dplyr::select(-n) %>%
#   distinct() %>%
#   kable(format = "rst", digits = 3) %>%
#   print()
# 
# ggplot(res, aes(x=mean_error_allele_prop, mean_child_allele_prop, col = type)) + 
#   geom_point() + 
#   geom_abline(slope=0.5, intercept =0) + 
#   #geom_smooth(method = "rlm") + 
#   theme_bw()
# }