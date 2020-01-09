
# combine ms and 1k genome data into ms_1k 
ms <- fread("data/GCST005531-independently-associated-loci.tsv") %>% as_tibble()
ms <- separate(ms, 'Lead Variant', into = c("chr","pos","A1","A2"), sep = "_")

af <- fread("data/1kg.polymorphic.SNPs.hg19.tsv") %>%
  as_tibble() %>%
  select(ID,CHROM,POS,REF,ALT,EUR_AF) %>%
  rename(chr=CHROM,
         pos=POS) %>%
  filter(chr != "X") %>% # remove X chromosome SNPs
  filter(chr != "Y") # remove Y chromosome SNPs

af$chr <- as.integer(af$chr)
ms$chr <- as.integer(ms$chr)
ms$pos <- as.integer(ms$pos)
ms_1k <- left_join(ms, af, by = c("chr", "pos"))
ms_1k <- ms_1k %>% rename(snp = rsID, a1 = A1, a2 = A2, eaf = EUR_AF)
MDD_sig <- MDD_sig %>% rename(snp = SNP, a1 = effect_allele, a2 = other_allele)
 
#ms_1k$rsID[!ms_1k$rsID%in%MDD_sig$SNP]
#ms_1k[order(ms_1k$rsID),]
