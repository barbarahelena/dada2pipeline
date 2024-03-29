library(dada2)
require(tidyverse)
packageVersion('dada2')

# raw reads (fastq.gz files) for should be in separate folders, by sequencing run: 'run1' ,'run2', etc 

##############################
forward_pattern <- '_R1.fastq.gz'
reverse_pattern <- '_R2.fastq.gz'
##############################

# get runs
runs <- list.dirs('raw')[grep('run', list.dirs('raw'))]
runs
length(runs)

# remove primers with seqtk

# >>> install seqtk with conda/mamba: "conda install -c bioconda seqtk"

# check if fwd_primer is present at the start of fwd reads!
# check if rev_primers (or reverse complement of rev_primers) is present at the start of reverse reads!
# use grep and wildcard

fwd_primer <- "CCTACGGGAGGCAGCAG"
rev_primer <- "TACNVGGGTATCTAAKCC"
# fwd_primer <- ""
# rev_primer <- ""
seqtk <- '/Users/barbaraverhaar/miniconda3/bin/seqtk' # change as needed!

dir.create('trimmed')
# removing 16S rRNA primers

for (i in runs) {
    run_name <-  str_remove(str_remove(i, "[[:punct:]]"), 'raw')
    dir.create(file.path('trimmed', run_name))
    cat('\n\n***** Trimming 16S rRNA primers from', run_name,'...\n')
    l <- list.files(i)
    ll <- file.path(i, l)
    fwd_reads <- ll[grep(forward_pattern, ll)]
    rev_reads <- ll[grep(reverse_pattern, ll)]
    cat("\n***", run_name,'has', length(fwd_reads), 'samples.\n')
    samples_fwd <- stringr::str_remove(fwd_reads, forward_pattern)
    samples_fwd <- stringr::str_remove(samples_fwd, paste0(i,'/'))
    samples_rev <- stringr::str_remove(rev_reads, reverse_pattern)
    samples_rev <- stringr::str_remove(samples_rev, paste0(i,'/'))
    df <- data.frame(fwd_reads, rev_reads, samples_fwd, samples_rev)
    for (j in 1:nrow(df)) {
        if (df$samples_fwd[j] == df$samples_rev[j]) {
            cat('\nTrimming 16S rRNA primers from sample', df$samples_fwd[j],'...\n')
            fwd_trimmed <- str_replace(df$fwd_reads[j], 'raw', 'trimmed')
            rev_trimmed <- str_replace(df$rev_reads[j], 'raw', 'trimmed')
            trim_fwd <- paste(seqtk, "trimfq -b", nchar(fwd_primer), df$fwd_reads[j], '| gzip >', fwd_trimmed, sep = ' ')
            trim_rev <- paste(seqtk, "trimfq -b", nchar(rev_primer), df$rev_reads[j], '| gzip >', rev_trimmed, sep = ' ')
            system(trim_fwd)
            system(trim_rev)
        } else {
            cat("Forwards and reverse fastq files for sample", df$samples_fwd[j], "do not match!\n\n")
        }
    }
}


# QC 
qc_path <- 'qc_trimmed'
dir.create(qc_path)
for (i in runs) {
    run_name <-  str_remove(str_remove(i, "[[:punct:]]"), 'raw')
    i <- str_replace(i, 'raw', 'trimmed')
    l <- list.files(i)
    ll <- file.path(i, l)
    fwd_reads <- ll[grep(forward_pattern, ll)]
    rev_reads <- ll[grep(reverse_pattern, ll)]
    samples_fwd <- stringr::str_remove(fwd_reads, forward_pattern)
    samples_fwd <- stringr::str_remove(samples_fwd, paste0(i,'/'))
    samples_rev <- stringr::str_remove(rev_reads, reverse_pattern)
    samples_rev <- stringr::str_remove(samples_rev, paste0(i,'/'))
    if (!all(samples_fwd == samples_rev)) {
        cat('Some samples are missing forwards or reverse reads!\n\n Quitting analysis...')
    } else {
        samples <- samples_fwd
    }
    cat('\n\n***** Getting quality profiles for trimmed', run_name, "...\n")
    cat('\n***', run_name, 'has', length(samples), 'samples: \n')
    
    # make an aggregate plot for all fwd / reverse files from the run
    cat('\nPlotting forward reads..\n')
    pl1 <- plotQualityProfile(fwd_reads, aggregate = T)
    ggsave(pl1, path = qc_path, filename = paste0('QC_plot_', run_name, '_forward_reads.pdf'), width = 6, height = 5, device = 'pdf')
    cat('\nPlotting reverse reads..\n')
    pl2 <- plotQualityProfile(rev_reads, aggregate = T)
    ggsave(pl2, path = qc_path, filename = paste0('QC_plot_', run_name, '_reverse_reads.pdf'), width = 6, height = 5, device = 'pdf')
}

# examine trimmed reads plots and decide on trimming thresholds!
    
###### parameters trimming and filtering ###

trunq <- 4
maxEEfwd <- 4
maxEErev <- 4
trimming_fwd <- 270
trimming_rev <- 230
maxN <- 0

nthreads <- 6 

######### accepted ASV length range #####

min_length <- 380
max_length <- 480

### pooling ###

pooling <- F # T is a bit more sensitive for low abundance ASVs, but a lot slower and a bit more false positives
# pooling <- T
# pooling <- 'pseudo' # pseudo-pooling gives intermediary speed and sensitivity

#########################################

# assign taxonomy with the dada2 assignTaxonomy function OR with the IdTaxa algorithm from the DECIPHER package
# taxonomy_assignment <- "IdTaxa"
taxonomy_assignment <- "dada2"

#########################################

dir.create("filtered")
dir.create('results')

# for each run: filter, learn error rates, infer ASVs, and get ASV table
getN <- function(x) sum(getUniques(x))

# for each run, trim and filter reads, learn error rates, infer ASVs, and make ASV table
for (i in runs) {
    run_name <-  str_remove(str_remove(i, "[[:punct:]]"), 'raw')
    fastqFs <- sort(list.files('trimmed', pattern = forward_pattern, recursive = T))
    fastqRs <- sort(list.files('trimmed', pattern = reverse_pattern, recursive = T))
    fastqFs <- fastqFs[grep(run_name, fastqFs)]
    fastqRs <- fastqRs[grep(run_name, fastqRs)]
    if (length(fastqFs) != length(fastqRs)) stop("Forward and reverse files do not match.")
    
    # trim and filter reads
    out <- filterAndTrim(fwd = file.path('trimmed', fastqFs), filt = file.path('filtered', fastqFs),
                         rev = file.path('trimmed', fastqRs), filt.rev = file.path('filtered', fastqRs),
                         truncLen = c(trimming_fwd, trimming_rev), 
                         # maxEE = c(maxEEfwd, maxEErev), 
                         truncQ = trunq, 
                         maxN = maxN, 
                         rm.phix = T,
                         compress = T, # set to 'True' to save space at the cost a bit of speed 
                         verbose = T, 
                         multithread = nthreads)
    
    filts <- file.path('filtered', run_name, list.files(file.path('filtered', run_name)))
    filtFs <- filts[grep(forward_pattern, filts)]
 
    # Learn error rates
    cat('\nLearning error rates for', run_name,'...\n\n')
    
    cat('\nLearning forward error rates...\n\n')
    errF <- learnErrors(filtFs, multithread = nthreads, verbose = T)
    pl1 <- plotErrors(errF)
    ggsave(pl1, path = qc_path, filename = paste0('error_plot_', run_name, '_forward_reads.pdf'), width = 6, height = 5, device = 'pdf')
    
    filtRs <- filts[grep(reverse_pattern, filts)]
    cat('\nLearning reverse error rates...\n\n')
    errR <- learnErrors(filtRs, multithread = nthreads, verbose = T)
    pl2 <- plotErrors(errF)
    ggsave(pl2, path = qc_path, filename = paste0('error_plot_', run_name, '_reverse_reads.pdf'), width = 6, height = 5, device = 'pdf')
    
    # infer ASVs
    cat('\nInferring forward ASVs...\n\n')
    dadaFs <- dada(filtFs, err = errF, multithread = nthreads)
    print(dadaFs[[1]])
    
    cat('\nInferring reverse ASVs...\n\n')
    dadaRs <- dada(filtRs, err = errR, multithread = nthreads, pool = pooling)
    print(dadaRs[[1]])
    
    # merge forward and reverse ASVs
    mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose = T)
    
    # make ASV table
    seqtab <- makeSequenceTable(mergers)
    table(nchar(getSequences(seqtab)))
    
    # save run tab
    saveRDS(seqtab, file.path('results', paste0(run_name,'.RDS')))
    
    # track reads
    track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN))
    colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged")
    sample.names <- str_remove(sapply(str_split(fastqFs,'/'),'[',2),forward_pattern)
    rownames(track) <- sample.names
    track <- as.data.frame(track)
    track$sample <- rownames(track)

    # save read tracking
    write.csv(track, file.path('results', paste0(run_name,'_track.csv')), row.names = F)
}

# integrate read tracking so far from multiple runs
ff <- list.files('results')
tracks <- file.path('results', ff[grep('track.csv', ff)])
tt <- matrix(nrow = 0, ncol = 6)
colnames(tt) <- c('input', 'filtered', 'denoisedF', 'denoisedR', 'merged', 'sample')

# import and append read tracking from all runs
for (i in 1:length(runs)) {
    tr <- read.csv(tracks[i])
    tt <- rbind(tt, tr)
}
tracks <- as.data.frame(tt)
tracks

# merge ASV tables from multiple runs
tables <- file.path('results', ff[grep('.RDS', ff)])
tables <- tables[grep("run..RDS", tables)]
st.all <- mergeSequenceTables(tables = tables)
rownames(st.all) <- str_remove(rownames(st.all), forward_pattern)

# remove ASVs outside accepted ASV sequence length range
seqtab <- st.all[, nchar(colnames(st.all)) %in% min_length:max_length]
no_removed <- ncol(st.all) - ncol(seqtab)
no_removed 
table(nchar(getSequences(seqtab)))
track_in_len <- rowSums(st.all)
tracks$in_len_range <- track_in_len[match(tracks$sample, names(track_in_len))]
cat('\nProportion of counts remaining after ASV length range filter:\n')
cat(sum(seqtab)/sum(st.all), "\n\n")

# remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method = "consensus", multithread = nthreads, verbose = T)
rownames(seqtab.nochim ) <- str_remove(rownames(seqtab.nochim), forward_pattern)
table(nchar(getSequences(seqtab.nochim)))
track_non_chimera <- rowSums(seqtab.nochim)
tracks$non_chimera <- track_non_chimera[match(tracks$sample, names(track_non_chimera))]
cat('Proportion of counts remaining after chimera filter:\n')
cat(sum(seqtab.nochim)/sum(st.all), "\n\n")
s <- seqtab.nochim

# assign taxonomy (SILVA v138)

if (taxonomy_assignment == 'dada2') {
    # check if tax db exists, else download
    dada2_tax_db <- 'tax_db/silva_nr99_v138_train_set.fa.gz'
    dada2_tax_db_link <- 'https://zenodo.org/record/3986799/files/silva_nr99_v138_train_set.fa.gz\\?download\\=1'
    dir.create('tax_db')
    if (!file.exists(dada2_tax_db)) {
        cat('Downloading SILVA v138 taxonomy db\n\n')
        get_tax <- paste0("wget -c ", dada2_tax_db_link, " -O ", dada2_tax_db)
        system(get_tax)
    } else {
        cat('SILVA v138 taxonomy db is already downloaded!\n\n')
    }
    
    # check if species-level tax db exists, else download
    dada2_tax_db_species <- 'tax_db/silva_species_assignment_v138.fa.gz'
    dada2_tax_db_species_link <- 'https://zenodo.org/record/3986799/files/silva_species_assignment_v138.fa.gz\\?download\\=1'
    if (!file.exists(dada2_tax_db_species)) {
        cat('Downloading SILVA v138 taxonomy species db\n\n')
        get_tax_species <- paste0("wget -c ", dada2_tax_db_species_link, " -O ", dada2_tax_db_species)
        system(get_tax_species)
    } else {
        cat('SILVA v138 taxonomy species db is already downloaded!\n\n')
    }
    
    # Assign taxonomy using the DADA2 function
    cat("Assigning taxonomy...\n\n")
    s.tax <- assignTaxonomy(s, dada2_tax_db, multithread = nthreads, minBoot = 80)
    # Assign species taxonomy
    cat("Assigning species-level taxonomy...\n\n")
    taxa <- addSpecies(s.tax, dada2_tax_db_species, verbose = T, allowMultiple = 3, tryRC = T)
    library(Biostrings)
    dna <- DNAStringSet(getSequences(s)) # Create a DNAStringSet from the ASVs
    names(dna) <- paste0('ASV',1:length(dna))
    writeXStringSet(dna, 'results/ASVs.fasta', width = 600) # write ASVs to fasta file
} else if (taxonomy_assignment == 'IdTaxa') {
    # OR, assign taxonomy using the IDTAXA algorithm
    library(DECIPHER); packageVersion("DECIPHER")
    
    # check if IDTAXA tax db exists, else download
    idtaxa_tax_db <- 'tax_db/SILVA_SSU_r138_2019.RData'
    id_taxa_link <- 'http://www2.decipher.codes/Classification/TrainingSets/SILVA_SSU_r138_2019.RData'
    if (!file.exists(idtaxa_tax_db)) {
        cat('Downloading IDTAXA SILVA v138 taxonomy db\n\n')
        get_idtaxa <- paste0("wget -c ", id_taxa_link," -O ", idtaxa_tax_db)
        system(get_idtaxa)
    } else {
        cat('IDTAXA SILVA v138 taxonomy db is already downloaded!\n\n')
    }
    
    library(Biostrings)
    dna <- DNAStringSet(getSequences(s)) # Create a DNAStringSet from the ASVs
    names(dna) <- paste0('ASV',1:length(dna))
    writeXStringSet(dna, 'results/ASVs.fasta', width = 600) # write ASVs to fasta file
    load(idtaxa_tax_db) # load IdTaxa databases
    # assign taxonomy using IdTaxa
    ids <- IdTaxa(dna, trainingSet, strand = "both", processors = NULL, 
                  verbose = T, threshold = 60) # use all processors
    ranks <- c("domain", "phylum", "class", "order", "family", "genus", "species") # ranks of interest
    # Convert the output object of class "Taxa" to a matrix analogous to the output from assignTaxonomy
    taxid <- t(sapply(ids, function(x) {
        m <- match(ranks, x$rank)
        taxa <- x$taxon[m]
        taxa[startsWith(taxa, "unclassified_")] <- NA
        taxa
    }))
    colnames(taxid) <- ranks; rownames(taxid) <- getSequences(s)
    taxa <- taxid
}

# make multiple sequence alignment with MAFFT binary available in the 'bin' folder
# READ: https://mafft.cbrc.jp/alignment/software/

mafft <- '/usr/local/bin/mafft'

make_multiple_alignment <- paste0(mafft, " --auto --thread ", nthreads, " results/ASVs.fasta > results/ASVs.msa")
make_multiple_alignment
system(make_multiple_alignment)

# make phylogenetic tree with FastTreeDbl (Double Precision) --> needs FastTreeDbl available in the 'bin' folder
# READ: http://www.microbesonline.org/fasttree/#BranchLen (read for MAC OS)

# check if FastreeDbl is available
fastree <- '/Users/barbaraverhaar/miniconda3/bin/fasttree'

make_tree <- paste0(fastree, " -nt -gtr < results/ASVs.msa > results/ASVs.tree")
system(make_tree)

# make phyloseq object
library(phyloseq); packageVersion("phyloseq")

# replace ASV sequences with ASV names (matching those in fasta)
asv_seqs <- as.character(dna)
asvs_names <- names(asv_seqs)

# set ASV names in ASV table
colnames(s) <- asvs_names[match(colnames(s), asv_seqs)]
# edit sample names in ASV table
rownames(s) <- str_remove(rownames(s), forward_pattern)

# edit sample names in taxonomy table
rownames(taxa) <- asvs_names[match(rownames(taxa), asv_seqs)]

# import and root tree

tree <- ape::read.tree('results/ASVs.tree')
tree_rooted <- phytools::midpoint.root(tree)

ps <- phyloseq(otu_table(s, taxa_are_rows = FALSE), tax_table(taxa), seq, tree_rooted)
ps@refseq <- dna

ps
sample_sums(ps)
saveRDS(ps, 'results/phyloseq_object.RDS')

tracks
write.csv(tracks, 'results/read_tracking.csv', row.names = F)
