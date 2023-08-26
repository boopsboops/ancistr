#!/usr/bin/env Rscript

### HOUSEKEEPING ###

# load up libs and funs
source(here::here("scripts/load-libs-funs.R"))


### RETRIEVE DATA ###

# search ncbi genbank
ancistrus.results <- rentrez::entrez_search(db="nuccore", term="(Ancistrus[ORGN] AND COI) OR (Lasiancistrus[ORGN] AND COI)", retmax=99999999, use_history=FALSE)
#ancistrus.results <- entrez_search(db="nuccore", term="(Ancistrus[ORGN] AND cytb) OR (Lasiancistrus[ORGN] AND cytb)", retmax=99999999, use_history=FALSE)

# get ncbi metadata
ancistrus.tab.ncbi <- ncbi_byid(ancistrus.results$ids) %>% as_tibble()
#ancistrus.tab.ncbi <- purrr::map(ancistrus.results$ids,ncbi_byid) %>% bind_rows() %>% as_tibble()

# get bold data
ancistrus.tab.bold <- bold::bold_seqspec("Ancistrus",format="tsv",sepfasta=FALSE,response=FALSE) %>% as_tibble()

# check data
if(is.data.frame(ancistrus.tab.ncbi)!=TRUE) {stop(writeLines("\nError retrieving data from NCBI. Please try again.\n"))}
if(is.data.frame(ancistrus.tab.bold)!=TRUE) {stop(writeLines("\nError retrieving data from BOLD. Please try again.\n"))}


### CLEAN UP AND FORMAT DATA ###

# clean genbank data
ancistrus.tab.ncbi.clean <- clean_ncbi(ancistrus.tab.ncbi)

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

# load in house data
house.ancistrus <- readr::read_csv(file=here("assets/ancistrus-coi.csv"),show_col_types=FALSE)

# merge with genbank
ancistrus.tab.merged.derep.house <- dplyr::bind_rows(ancistrus.tab.merged.derep,house.ancistrus)
#ancistrus.tab.merged <- ancistrus.tab.ncbi 

# make some labels
ancistrus.tab.merged.derep.house.labels <- ancistrus.tab.merged.derep.house %>% 
    mutate(country=str_split_fixed(country,":",2)[,1]) %>%
    mutate(country=replace_na(country,"No locality")) %>%
    mutate(species=if_else(is.na(species),sciNameOrig,species)) %>%
    mutate(acc=if_else(is.na(gbAccession),dbid,gbAccession)) %>%
    mutate(tiplab=paste(acc,species,country,sep=" | ")) %>%
    mutate(tipcol=if_else(dbid=="RC","domestic","other")) %>%
    relocate(dbid,.before=1)

# convert to fasta
ancistrus.fas <- tab2fas(ancistrus.tab.merged.derep.house.labels,"nucleotides","dbid")


### PHYLOGENETIC ANALYSIS ###

# align
ancistrus.aligned <- ips::mafft(ancistrus.fas,exec="mafft")

# get starting tree
ancistrus.aligned.pars <- phangorn::acctran(multi2di(pratchet(as.phyDat(ancistrus.aligned))),as.phyDat(ancistrus.aligned))

# optimise ML tree
ancistrus.tree <- phangorn::optim.pml(pml(tree=ancistrus.aligned.pars, data=as.phyDat(ancistrus.aligned),model="TrN",k=4,shape=0.5,inv=0),model="TrN",optInv=FALSE,optNni=TRUE,optGamma=TRUE)
ancistrus.tr <- phangorn::midpoint(ancistrus.tree$tree)


### PLOT TREE ###

# ggtree plot
p <- ggtree(ancistrus.tr,ladderize=TRUE,right=TRUE,size=1,color="gray30") %<+% ancistrus.tab.merged.derep.house.labels
p <- p + geom_tiplab(aes(label=tiplab,color=tipcol),offset=0.001) + 
    scale_color_brewer(palette="Set1") +# values = c("red","grey20")
    geom_tippoint(size=2,color="gray30") + #aes(color=species),
    #geom_hilight(node=ape::getMRCA(phy=ancistrus.tr,tip=c("RC","1108491107")),fill="yellow",alpha=0.5) +
    xlim(0,0.3) +
    theme(legend.position="none")
print(p)


### MAKE DISTANCE MATRIX ###

ancistrus.aligned.dist <- ape::dist.dna(ancistrus.aligned,model="raw",pairwise.deletion=TRUE,as.matrix=TRUE)

# convert to data frame in two cols
ancistrus.aligned.dist.rc <- as_tibble(ancistrus.aligned.dist,rownames="comp1") %>% 
    pivot_longer(cols=-comp1,names_to="comp2",values_to="dist") 

# filter and clean up
ancistrus.aligned.dist.rc.filt <- ancistrus.aligned.dist.rc %>% 
    filter(comp1=="RC" & comp2!="RC") %>%
    arrange(dist) %>%
    slice_head(n=10) %>%
    select(-comp1) %>%
    mutate(dist=(1-dist)*100) %>%
    rename(dbid=comp2)

# add info
ancistrus.aligned.dist.rc.filt.tidy <- ancistrus.aligned.dist.rc.filt %>% 
    left_join(select(ancistrus.tab.merged.derep.house.labels,dbid,acc,tiplab)) %>%
    select(tiplab,dist) %>%
    separate_wider_delim(cols=tiplab,names=c("Accession","Species","Locality"),delim=" | ") %>%
    rename(`Genetic similarity (%)`=dist)

# convert to table
#ancistrus.aligned.dist.rc.filt.tidy.table <- ancistrus.aligned.dist.rc.filt.tidy %>% kable(digits=1)


### WRITE OUT FILES ###

# make date dir for results
today.dir <- here("results",paste0("Results_",Sys.Date()))
if(!dir.exists(today.dir)) {dir.create(today.dir,recursive=TRUE)}

# write files
readr::write_csv(ancistrus.tab.merged,file=here(today.dir,"ncbi-bold-ancistrus.csv"))
readr::write_csv(ancistrus.aligned.dist.rc,file=here(today.dir,"ancistrus-distances.csv"))
ape::write.FASTA(ancistrus.aligned,file=here(today.dir,"ncbi-bold-ancistrus.fasta"))
ape::write.tree(ancistrus.tr,file=here(today.dir,"ncbi-bold-ancistrus.nwk"))
ggsave(plot=p, filename=here(today.dir,"ncbi-bold-ancistrus.pdf"), width=11, height=12, bg="transparent", limitsize=FALSE)
