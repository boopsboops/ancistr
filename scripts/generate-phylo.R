#!/usr/bin/env Rscript

# load libs
library("ape")
library("tidyverse")
library("magrittr")
library("ips")
library("phangorn")
library("rentrez")
library("bold")
library("traits")
library("ggtree")
source("https://raw.githubusercontent.com/legalLab/protocols-scripts/master/scripts/tab2fas.R")

# search genbank
ancistrus.results <- entrez_search(db="nuccore", term="(Ancistrus[ORGN] AND COI) OR (Lasiancistrus[ORGN] AND COI)", retmax=99999999, use_history=FALSE)

# get metadata
ancistrus.tab.ncbi <- purrr::map(ancistrus.results$ids,ncbi_byid) %>% bind_rows() %>% as_tibble()

# clean genbank data
ancistrus.tab.ncbi %<>% filter(gi_no!="NCBI_GENOMES") %>% 
    distinct(gi_no, .keep_all=TRUE) %>% 
    mutate(acc_no=str_replace_all(acc_no,"\\.[0-9]",""), source="GENBANK") %>%
    # fix the lat_lon into decimal
    mutate(lat=paste(str_split_fixed(lat_lon, " ", 4)[,1], str_split_fixed(lat_lon, " ", 4)[,2]), lon=paste(str_split_fixed(lat_lon, " ", 4)[,3], str_split_fixed(lat_lon, " ", 4)[,4])) %>%
    mutate(lat=if_else(grepl(" N",lat), true=str_replace_all(lat," N",""), false=if_else(grepl(" S",lat), true=paste0("-",str_replace_all(lat," S","")), false=lat))) %>%
    mutate(lon=if_else(grepl(" E",lon), true=str_replace_all(lon," E",""), false=if_else(grepl(" W",lon), true=paste0("-",str_replace_all(lon," W","")), false=lon))) %>% 
    mutate(lat=str_replace_all(lat,"^ ", NA_character_), lon=str_replace_all(lon,"^ ", NA_character_)) %>%
    mutate(lat=as.numeric(lat), lon=as.numeric(lon)) %>% 
    # tidy up
    select(-taxonomy,-organelle,-keyword,-lat_lon) %>% 
    rename(sciNameOrig=taxon,notesGenBank=gene_desc,dbid=gi_no,gbAccession=acc_no,catalogNumber=specimen_voucher,
    publishedAs=paper_title,publishedIn=journal,publishedBy=first_author,date=uploaded_date,decimalLatitude=lat,decimalLongitude=lon,nucleotides=sequence)

# get bold data
ancistrus.tab.bold <- as_tibble(bold_seqspec("Ancistrus",format="tsv",sepfasta=FALSE,response=FALSE))

# clean up bold data
ancistrus.tab.bold %<>%
    mutate_all(na_if,"") %>%
    mutate(nucleotides=str_replace_all(nucleotides,"-",""), nucleotides=str_replace_all(nucleotides,"N",""), num_bases=nchar(nucleotides)) %>% 
    filter(num_bases > 0) %>%
    filter(institution_storing!="Mined from GenBank, NCBI") %>% 
    mutate(processidUniq=paste(processid,markercode,sep=".")) %>% 
    distinct(processidUniq, .keep_all=TRUE)

# remove duplicates
dups <- intersect(ancistrus.tab.ncbi$gbAccession,ancistrus.tab.bold$genbank_accession)

# remove dups and clean
ancistrus.tab.bold %<>% filter(!genbank_accession %in% dups) %>% 
    mutate(source="BOLD",nucleotides=str_to_lower(nucleotides), length=as.character(str_length(nucleotides))) %>% 
    select(source,processidUniq,genbank_accession,species_name,lat,lon,country,institution_storing,catalognum,nucleotides,length) %>%
    rename(dbid=processidUniq,gbAccession=genbank_accession,sciNameOrig=species_name,decimalLatitude=lat,decimalLongitude=lon,institutionCode=institution_storing,catalogNumber=catalognum)

# merge with genbank
ancistrus.tab.merged <- bind_rows(ancistrus.tab.ncbi,ancistrus.tab.bold,tibble(dbid=c("RC","YPW"),sciNameOrig=c("Common bristlenose","Trinidad bristlenose"),country=c("New Zealand (trade)","Trinidad"),nucleotides=c("ctttacctagtgtttggtgcctgagccggaatggttggtacagccctcagtctcttaattcgagctgagttaagccaacccggttctctattaggtgatgaccagatttataatgtcatcgttaccgcacatgctttcgtaataattttctttatagtcatgccaatcataattgggggctttggaaattgactagttccactaatgattggggcacccgatatagccttcccacgaataaataacatgagcttctgactactgcccccctcattccttcttctactggcctcttcaggggttgaagcgggagctgggacaggttgaactgtatacccacccctcgccggaaacctggcccacgcaggagcttccgttgacctgactattttttcactacacctggctggtgtttcttcaattctgggggcaattaacttcattaccacaatcattaacataaagcccccggctatttcacaataccaaacccccctatttgtgtgagccgtacttgttacagcggtcctactcctgctttccttgcccgttctggccgccggcattacaatactgctcacagatcgaaatctaaacaccacattctttgaccctgcgggcggtggagatcctatcctttatcaacactta","cccggttctctattaggtgatgaccaaatttataatgtcatcgttaccgcacatgctttcgtaataattttctttatagttataccaattatgattgggggctttggaaattgactagtaccactaatgatcggagcacctgacatagccttcccacgaataaataacataagcttctgactactccccccttcattccttctcttgctagcctcatcaggagtcgaggcaggggcagggacaggttgaactgtatacccaccccttgccggaaacttagcccacgcaggggcctcagttgacctaactatcttttcactacacttggctggtgtatcctcaattctaggtgctattaacttcattaccacaattattaacataaaacccccagctatttcacaataccaaacccccttatttgtgtgggccgtacttgtaacagcggttctgctcctgctttccttacccgttctagctgccggtattacaatactactaacagatcgaaacctaaataccacattctttgaccctgcaggcggaggggatcctatcctttatcagcacttattt")))

# make some labels
ancistrus.tab.merged %<>% mutate(country=replace_na(country,"No locality")) %>%
    mutate(labels=paste(dbid,sciNameOrig,str_trunc(country,width=30,side="right",ellipsis="...")))

# convert to fasta
ancistrus.fas <- tab2fas(ancistrus.tab.merged,"nucleotides","labels")

# align
ancistrus.aligned <- ips::mafft(ancistrus.fas,exec="mafft")

# get starting tree
ancistrus.aligned.pars <- acctran(multi2di(pratchet(as.phyDat(ancistrus.aligned))),as.phyDat(ancistrus.aligned))

# optimise ML tree
ancistrus.tree <- optim.pml(pml(tree=ancistrus.aligned.pars, data=as.phyDat(ancistrus.aligned),model="TrN",k=4,shape=0.5,inv=0),model="TrN",optInv=FALSE,optNni=TRUE,optGamma=TRUE)
ancistrus.tr <- midpoint(ancistrus.tree$tree)

# plot tree
p <- ggtree(ancistrus.tr, ladderize=TRUE, color="darkorange2", size=1) + xlim(0,0.4)
p <- p + geom_tiplab(aes(size=3))
ggsave(plot=p, filename="results/ancistrus.png", width=10, height=20, bg="transparent", limitsize=FALSE)
rm(p)
