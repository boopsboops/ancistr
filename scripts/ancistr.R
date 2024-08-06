#!/usr/bin/env Rscript

### HOUSEKEEPING ###

# load up libs and funs
source(here::here("scripts/load-libs-funs.R"))


### RETRIEVE DATA ###

# search ncbi genbank
ancistrus.results <- rentrez::entrez_search(db="nuccore", term="(Ancistrus[ORGN] AND COI) OR (Ancistrus[ORGN] AND COX1) OR (Lasiancistrus[ORGN] AND COI) AND (490:1700[SLEN])", retmax=9999, use_history=FALSE)
#ancistrus.results <- entrez_search(db="nuccore", term="(Ancistrus[ORGN] AND cytb) OR (Lasiancistrus[ORGN] AND cytb)", retmax=99999, use_history=FALSE)

# get ncbi metadata - may need to chunk these in future if these are more than ~200 hits
ancistrus.tab.ncbi <- ncbi_byid(ancistrus.results$ids) %>% as_tibble()
#ancistrus.tab.ncbi <- purrr::map(ancistrus.results$ids,ncbi_byid) %>% bind_rows() %>% as_tibble()

# get bold data
ancistrus.tab.bold <- bold::bold_seqspec("Ancistrus",format="tsv",sepfasta=FALSE,response=FALSE) %>% as_tibble()

# check data
if(is.data.frame(ancistrus.tab.ncbi)!=TRUE) {stop(writeLines("\nError retrieving data from NCBI. Please try again.\n"))}
if(length(ancistrus.results$ids) != nrow(ancistrus.tab.ncbi)) {stop(writeLines("\nError retrieving data from NCBI. Please try again.\n"))}
if(is.data.frame(ancistrus.tab.bold)!=TRUE) {stop(writeLines("\nError retrieving data from BOLD. Please try again.\n"))}


### CLEAN UP AND FORMAT DATA ###

# clean genbank data
ancistrus.tab.ncbi.clean <- clean_ncbi(df=ancistrus.tab.ncbi)

# clean up bold data
ancistrus.tab.bold.clean <- clean_bold(df=ancistrus.tab.bold)

# remove duplicates
dups <- intersect(ancistrus.tab.ncbi.clean$gbAccession,ancistrus.tab.bold.clean$gbAccession)
ancistrus.tab.bold.clean <- ancistrus.tab.bold.clean %>% filter(!gbAccession %in% dups) 

# merge genbank and bold
ancistrus.tab.merged <- dplyr::bind_rows(ancistrus.tab.ncbi.clean,ancistrus.tab.bold.clean)

# tidy names of undescribed spp.
ancistrus.tab.merged.tidy <- ancistrus.tab.merged %>%
    mutate(species=str_replace_all(sciNameOrig,"sp\\. .+","sp.")) %>% 
    mutate(species=str_replace_all(species,"\\. ",".")) %>% 
    mutate(species=paste(str_split_fixed(species," ",3)[,1],str_split_fixed(species," ",3)[,2],sep=" ")) %>%
    mutate(delimitation=paste(species,publishedAs))

# dereplicate
ancistrus.tab.merged.derep <- ancistrus.tab.merged.tidy %>% 
    group_by(delimitation) %>% 
    dplyr::group_map(~hap_collapse_df(df=.x,lengthcol="length",nuccol="nucleotides",cores=1),.keep=TRUE) %>%
    dplyr::bind_rows()

# load in-house data
house.ancistrus <- readr::read_csv(file=here("assets/ancistrus-coi.csv"),show_col_types=FALSE)

# merge with genbank
ancistrus.tab.merged.derep.house <- dplyr::bind_rows(ancistrus.tab.merged.derep,house.ancistrus)

# make some labels
ancistrus.tab.merged.derep.house.labels <- ancistrus.tab.merged.derep.house %>% 
    mutate(country=str_split_fixed(country,":",2)[,1]) %>%
    mutate(country=replace_na(country,"No locality")) %>%
    mutate(species=if_else(is.na(species),sciNameOrig,species)) %>%
    mutate(acc=if_else(is.na(gbAccession),dbid,gbAccession)) %>%
    mutate(tiplab=paste(acc,species,country,sep=" | ")) %>%
    mutate(tipcol=if_else(dbid=="RC","domestic","other")) %>%
    relocate(dbid,.before=1)

# make stats
today.long <- format(lubridate::today(),"%d %B %Y")
n.ancistrus <- ancistrus.tab.merged.derep.house.labels %>% 
    filter(!grepl("sp\\.",species) & !grepl("Lasiancistrus",species)) %>% 
    distinct(species) %>% 
    pull() %>% 
    length()
gb.version <- read.table("https://ftp.ncbi.nih.gov/genbank/GB_Release_Number")$V1

### PHYLOGENETIC ANALYSIS ###

# convert to fasta
ancistrus.fas <- tab2fas(ancistrus.tab.merged.derep.house.labels,"nucleotides","dbid")

# align
ancistrus.aligned <- ips::mafft(ancistrus.fas,exec="mafft")

# trim alignment at 80% gaps
ancistrus.aligned.trimmed <- ancistrus.aligned[,which(apply(as.character(ancistrus.aligned),2,prop_gaps) < 0.8)]

# get starting tree
ancistrus.aligned.trimmed.pars <- phangorn::acctran(multi2di(pratchet(as.phyDat(ancistrus.aligned.trimmed))),as.phyDat(ancistrus.aligned.trimmed))

# optimise ML tree
ancistrus.tree <- phangorn::optim.pml(pml(tree=ancistrus.aligned.trimmed.pars, data=as.phyDat(ancistrus.aligned.trimmed),model="TrN",k=4,shape=0.5,inv=0),model="TrN",optInv=FALSE,optNni=TRUE,optGamma=TRUE,rearrangement="stochastic")#NNI
#ancistrus.tr <- phangorn::midpoint(ancistrus.tree$tree)

# root tree and drop Lasiancistrus
lasi.tips <- ancistrus.tab.merged.derep.house.labels %>% filter(grepl("Lasiancistrus",tiplab)) %>% pull(dbid)
ancistrus.tr <- ape::root(ancistrus.tree$tree,node=ape::getMRCA(phy=ancistrus.tree$tree,tip=lasi.tips),resolve.root=TRUE)
ancistrus.tr.rooted <- ape::drop.tip(phy=ancistrus.tr,tip=lasi.tips)


### PLOT TREE ###

# ggtree plot
p <- ggtree(ancistrus.tr.rooted,ladderize=TRUE,right=TRUE,size=1,color="gray30") %<+% ancistrus.tab.merged.derep.house.labels
p <- p + geom_tiplab(aes(label=tiplab,color=tipcol),offset=0.001) + 
    scale_color_brewer(palette="Set1") +
    geom_tippoint(size=2,color="gray30") +
    xlim(0,0.3) +
    theme(legend.position="none")
#print(p)


### MAKE DISTANCE MATRIX ###

ancistrus.aligned.trimmed.dist <- ape::dist.dna(ancistrus.aligned.trimmed,model="raw",pairwise.deletion=TRUE,as.matrix=TRUE)

# convert to data frame in two cols
ancistrus.aligned.trimmed.dist.rc <- as_tibble(ancistrus.aligned.trimmed.dist,rownames="comp1") %>% 
    pivot_longer(cols=-comp1,names_to="comp2",values_to="dist") 

# filter and clean up
ancistrus.aligned.trimmed.dist.rc.filt <- ancistrus.aligned.trimmed.dist.rc %>% 
    filter(comp1=="RC" & comp2!="RC") %>%
    arrange(dist) %>%
    slice_head(n=10) %>%
    select(-comp1) %>%
    mutate(dist=(1-dist)*100) %>%
    rename(dbid=comp2)

# add info
ancistrus.aligned.trimmed.dist.rc.filt.tidy <- ancistrus.aligned.trimmed.dist.rc.filt %>% 
    left_join(select(ancistrus.tab.merged.derep.house.labels,dbid,acc,tiplab)) %>%
    select(tiplab,dist) %>%
    separate_wider_delim(cols=tiplab,names=c("Accession","Species","Locality"),delim=" | ") %>%
    rename(`Sequence similarity (%)`=dist)

# convert to table
#ancistrus.aligned.trimmed.dist.rc.filt.tidy.table <- ancistrus.aligned.trimmed.dist.rc.filt.tidy %>% kable(digits=1)


### WRITE OUT FILES ###

# make date dir for results
today.dir <- here("results",paste0("Results_",Sys.Date()))
if(!dir.exists(today.dir)) {dir.create(today.dir,recursive=TRUE)}

# write files
readr::write_csv(ancistrus.tab.merged,file=here(today.dir,"ncbi-bold-ancistrus.csv"))
readr::write_csv(ancistrus.aligned.trimmed.dist.rc,file=here(today.dir,"ancistrus-distances.csv"))
ape::write.FASTA(ancistrus.aligned.trimmed,file=here(today.dir,"ncbi-bold-ancistrus.fasta"))
ape::write.tree(ancistrus.tr.rooted,file=here(today.dir,"ncbi-bold-ancistrus.nwk"))
ggsave(plot=p, filename=here(today.dir,"ncbi-bold-ancistrus.pdf"), width=11, height=12, bg="transparent", limitsize=FALSE)
