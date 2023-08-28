#!/usr/bin/env Rscript

# load libs
library("here")
library("ape")
library("tidyverse")
library("ips")
library("phangorn")
library("rentrez")
library("bold")
library("traits")
library("ggtree")
library("rmarkdown")
library("knitr")
library("yaml")
library("BiocManager")
library("BiocVersion")

# load funs
source("https://raw.githubusercontent.com/boopsboops/UTILITIES/main/RScripts/tab2fas.R")
source("https://raw.githubusercontent.com/boopsboops/UTILITIES/main/RScripts/hap_collapse_df.R")

# CLEAN NCBI
clean_ncbi <- function(df) {
    df.clean <- df %>% filter(gi_no!="NCBI_GENOMES") %>% 
    distinct(gi_no, .keep_all=TRUE) %>% 
    mutate(acc_no=str_replace_all(acc_no,"\\.[0-9]",""), source="GENBANK") %>%
    # fix the lat_lon into decimal
    mutate(lat=paste(str_split_fixed(lat_lon, " ", 4)[,1], str_split_fixed(lat_lon, " ", 4)[,2]), lon=paste(str_split_fixed(lat_lon, " ", 4)[,3], str_split_fixed(lat_lon, " ", 4)[,4])) %>%
    mutate(lat=if_else(grepl(" N",lat), true=str_replace_all(lat," N",""), false=if_else(grepl(" S",lat), true=paste0("-",str_replace_all(lat," S","")), false=lat))) %>%
    mutate(lon=if_else(grepl(" E",lon), true=str_replace_all(lon," E",""), false=if_else(grepl(" W",lon), true=paste0("-",str_replace_all(lon," W","")), false=lon))) %>% 
    mutate(lat=str_replace_all(lat,"^ ", NA_character_), lon=str_replace_all(lon,"^ ", NA_character_)) %>%
    mutate(lat=suppressWarnings(as.numeric(lat)), lon=suppressWarnings(as.numeric(lon))) %>% 
    # tidy up
    select(-taxonomy,-organelle,-keyword,-lat_lon) %>% 
    rename(sciNameOrig=taxon,notesGenBank=gene_desc,dbid=gi_no,gbAccession=acc_no,catalogNumber=specimen_voucher,
    publishedAs=paper_title,publishedIn=journal,publishedBy=first_author,date=uploaded_date,decimalLatitude=lat,decimalLongitude=lon,nucleotides=sequence)
    return(df.clean)
}


# CLEAN BOLD
clean_bold <- function(df) {
    df.clean <- df %>% mutate(across(where(is.character), ~na_if(.x,""))) %>%
    mutate(nucleotides=str_replace_all(nucleotides,"-",""), nucleotides=str_replace_all(nucleotides,"N",""), num_bases=nchar(nucleotides)) %>% 
    filter(num_bases > 0) %>%
    filter(institution_storing!="Mined from GenBank, NCBI") %>% 
    mutate(processidUniq=paste(processid,markercode,sep=".")) %>% 
    distinct(processidUniq, .keep_all=TRUE) %>%
    mutate(source="BOLD",nucleotides=str_to_lower(nucleotides), length=as.character(str_length(nucleotides))) %>% 
    mutate(lat=suppressWarnings(as.numeric(lat)), lon=suppressWarnings(as.numeric(lon))) %>%
    select(source,processidUniq,genbank_accession,species_name,lat,lon,country,institution_storing,catalognum,nucleotides,length) %>%
    rename(dbid=processidUniq,gbAccession=genbank_accession,sciNameOrig=species_name,decimalLatitude=lat,decimalLongitude=lon,institutionCode=institution_storing,catalogNumber=catalognum)
    return(df.clean)
}


# COUNT GAP CHARACTERS
prop_gaps <- function(x) {
    lx <- length(x)
    gx <- length(which(x == "-"))
    rx <- gx/lx
    return(rx)
}
