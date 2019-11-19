library(rtracklayer)
library(ShortRead)
library(Biostrings)
library(BSgenome)

## you know what I give up on Bioc and try to go tidy for bio too...
library(biobroom)
library(tidyverse)


## to create the example data
EfalGenome <- readDNAStringSet("/home/ele/ToxoDB-46_EfalciformisBayerHaberkorn1970_Genome.fasta")
EfalAnnotation <- import.gff("/home/ele/ToxoDB-46_EfalciformisBayerHaberkorn1970.gff")

## EfalAnnotationTxDb <- makeTxDbFromGFF("/home/ele/ToxoDB-46_EfalciformisBayerHaberkorn1970.gff")

save(EfalGenome, file="data/EfalGenome.rda", compress='xz')
save(EfalAnnotation, file="data/EfalAnnotation.rda", compress='xz')


## From the original transcriptome work on A. crassus
## if(!exists("VCF")){
##     source("scripts/3_SNP_analysis.R")
## }



EfalGenomeNamed <- EfalGenome
names(EfalGenomeNamed) <- gsub("(.*?) .*", "\\1", names(EfalGenomeNamed))

EfalCDS <- subset(EfalAnnotation, type%in%"CDS")

EfalCDSseq <- getSeq(EfalGenomeNamed, EfalCDS)
EfalCDSParented <- as.factor(unlist(EfalCDS@elementMetadata$Parent))
EfalCDSStrand <- as.factor(EfalCDS@strand)

EfalCDSTranscriptStrands <- by(EfalCDSStrand, EfalCDSParented, function(x) unique(as.character(x)))

EfalCDSTranscripts <- by(EfalCDSseq, EfalCDSParented, function(x) do.call(c, x))

transONT <- lapply(seq_along(EfalCDSTranscriptStrands), function(i) {
    if(EfalCDSTranscriptStrands[[i]]=="-"){
       transcript <- reverseComplement(EfalCDSTranscripts[[i]])
    } else {
       transcript <- EfalCDSTranscripts[[i]]
    }
    tran <- as.character(transcript)
    codons <- substring(tran,
                        seq(1, nchar(tran), by=3),
                        seq(3, nchar(tran), by=3))
    ont <- base.ontology.encode(tolower(codons))
    otl <- paste0(ont, collapse="")
    if(EfalCDSTranscriptStrands[[i]]=="-"){
        reverse(otl)
    } else {
        otl
    }
})


## That's a complete mess...
## we'd now need to split again by cds then replace the cds in the genome...
## by(EfalCDSseq, EfalCDSParented, function (x) cumsum(width(x)))

cat("foo")

cat("bar")

## ## BASE ONTOLOGY
## #
## # Code from Mark Blaxter, modified by John Davey,
## # Translated from Perl to R by Emanuel Heitlinger:

## # A phase 1 any change is nonsynonymous
## # B phase 2 any change is nonsynonymous
## # C phase 3 any change is nonsynonymous
## # D phase 1 change to CT is nonsynonymous
## # E phase 2 change to CT is nonsynonymous
## # F phase 3 change to CT is nonsynonymous
## # G phase 1 change to AG is nonsynonymous
## # H phase 2 change to AG is nonsynonymous
## # I phase 3 change to AG is nonsense
## # K phase 1 change to GT is nonsynonymous
## # L phase 2 change to A is nonsense, to anything else is nonsynonymous
## # J phase 3 change to G is nonsynonymous
## # M phase 3 change to G is nonsense, to A is nonsynonymous
## # N phase 3 any change synonymous
## # O phase 1 change to T nonsense, others nonsynonymous
## # P phase 3 change to AG is nonsynonymous
## # Q phase 1 change to T nonsense, to G nonsynonymous
## # R phase 2 change to AG nonsense, others nonsynonymous
## # S phase 3 change to A nonsense, others nonsynonymous
## # T phase 3 change to A nonsense, G nonsynonymous

## # W all changes are unknown # EH added 08/23/2011

## #        a           g           c           t
## #
## # a     aaa K OBF   aga R QBF   aca T ABN   ata I ABJ
## #       aag K OBF   agg R KBF   acg T ABN   atg M ABC
## #       aac N ABP   agc S ABP   acc T ABN   atc I ABJ
## #       aat N ABP   agt S ABP   act T ABN   att I ABJ
## #
## # g     gaa E OBF   gga G OBN   gca A ABN   gta V ABN
## #       gag E OBF   ggg G ABN   gcg A ABN   gtg V ABN
## #       gac D ABP   ggc G ABN   gcc A ABN   gtc V ABN
## #       gat D ABP   ggt G ABN   gct A ABN   gtt V ABN
## #
## # c     caa Q OBF   cga R QBN   cca P ABN   cta L GBN
## #       cag Q OBF   cgg R KBN   ccg P ABN   ctg L GBN
## #       cac H ABP   cgc R ABN   ccc P ABN   ctc L ABN
## #       cat H ABP   cgt R ABN   cct P ABN   ctt L ABN
## #
## # t     taa * AEF   tga * AEC   tca S ARN   tta L GRF
## #       tag * ABF   tgg W ALS   tcg S ALN   ttg L GLF
## #       tac Y ABI   tgc C ABT   tcc S ABN   ttc F ABP
## #       tat Y ABI   tgt C ABT   tct S ABN   ttt F ABP

base.ontology.encode <- function(x){  
    ## Set up Base Ontology vector
    base.ontology.encode.string = c(
        "aaa" = "OBF",
        "aag" = "OBF",
        "aac" = "ABP",
        "aat" = "ABP",
        "aga" = "QBF",
        "agg" = "KBF",
        "agc" = "ABP",
        "agt" = "ABP",
        "aca" = "ABN",
        "acg" = "ABN",
        "acc" = "ABN",
        "act" = "ABN",
        "ata" = "ABJ",
        "atg" = "ABC",
        "atc" = "ABJ",
        "att" = "ABJ",
        ##
        "gaa" = "OBF",
        "gag" = "OBF",
        "gac" = "ABP",
        "gat" = "ABP",
        "gga" = "OBN",
        "ggg" = "ABN",
        "ggc" = "ABN",
        "ggt" = "ABN",
        "gca" = "ABN",
        "gcg" = "ABN",
        "gcc" = "ABN",
        "gct" = "ABN",
        "gta" = "ABN",
        "gtg" = "ABN",
        "gtc" = "ABN",
        "gtt" = "ABN",
        ##
        "caa" = "OBF",
        "cag" = "OBF",
        "cac" = "ABP",
        "cat" = "ABP",
        "cga" = "QBN",
        "cgg" = "KBN",
        "cgc" = "ABN",
        "cgt" = "ABN",
        "cca" = "ABN",
        "ccg" = "ABN",
        "ccc" = "ABN",
        "cct" = "ABN",
        "cta" = "GBN",
        "ctg" = "GBN",
        "ctc" = "ABN",
        "ctt" = "ABN",
        ##
        "taa" = "AEF",
        "tag" = "ABF",
        "tac" = "ABI",
        "tat" = "ABI",
        "tga" = "AEC",
        "tgg" = "ALS",
        "tgc" = "ABT",
        "tgt" = "ABT",
        "tca" = "ARN",
        "tcg" = "ALN",
        "tcc" = "ABN",
        "tct" = "ABN",
        "tta" = "GRF",
        "ttg" = "GLF",
        "ttc" = "ABP",
        "ttt" = "ABP",
        ##
        ## bad strings
        "x"="W",
        "xx"="WW",
        "xxx"="WWW"
    );
    
    ## a replace for iupac and other "bad" bases
    tr <- function (y) paste0(rep("x", times=nchar(y)), collapse="")
    offenders <- which(nchar(gsub("[^agct]", "", x))!=3)
    x[offenders] <- sapply(x[offenders], tr)
    ## do the magic
    x <- unlist(x)
    base.ontology.encode.string[x]
}
  
base.ontology.decode = list(
  "A" = c("A" = "Nonsynonymous",  "C" = "Nonsynonymous",
    "G" = "Nonsynonymous",  "T" = "Nonsynonymous" ),
  "B" = c("A" = "Nonsynonymous",  "C" = "Nonsynonymous",
    "G" = "Nonsynonymous",  "T" = "Nonsynonymous" ),
  "C" = c( "A" = "Nonsynonymous",  "C" = "Nonsynonymous",
    "G" = "Nonsynonymous",  "T" = "Nonsynonymous" ),
  "D" = c( "A" = "Synonymous", "C" = "Nonsynonymous",
    "G" = "Synonymous", "T" = "Nonsynonymous" ),
  "E" = c( "A" = "Synonymous", "C" = "Nonsynonymous",
    "G" = "Synonymous", "T" = "Nonsynonymous" ),
  "F" = c( "A" = "Synonymous", "C" = "Nonsynonymous",
    "G" = "Synonymous", "T" = "Nonsynonymous" ),
  "G" = c( "A" = "Nonsynonymous",  "C" = "Synonymous",
    "G" = "Nonsynonymous",  "T" = "Synonymous" ),
  "H" = c( "A" = "Nonsynonymous",  "C" = "Synonymous",
    "G" = "Nonsynonymous",  "T" = "Synonymous" ),
  "I" = c( "A" = "Nonsense",  "C" = "Synonymous",
    "G" = "Nonsense",  "T" = "Synonymous" ),
  "J" = c( "A" = "Synonymous", "C" = "Synonymous",
    "G" = "Nonsynonymous",  "T" = "Synonymous" ),
  "K" = c( "A" = "Synonymous", "C" = "Synonymous",
    "G" = "Nonsynonymous",  "T" = "Nonsynonymous" ),
  "L" = c( "A" = "Nonsense",  "C" = "Nonsynonymous",
    "G" = "Nonsynonymous",  "T" = "Nonsynonymous" ),
  "M" = c( "A" = "Nonsynonymous",  "C" = "Synonymous",
    "G" = "Nonsense",  "T" = "Synonymous" ),
  "N" = c( "A" = "Synonymous", "C" = "Synonymous",
    "G" = "Synonymous", "T" = "Synonymous" ),
  "O" = c( "A" = "Nonsynonymous",  "C" = "Nonsynonymous",
    "G" = "Nonsynonymous",  "T" = "Nonsense" ),
  "P" = c( "A" = "Nonsynonymous",  "C" = "Synonymous",
    "G" = "Nonsynonymous",  "T" = "Synonymous" ),
  "Q" = c( "A" = "Synonymous", "C" = "Synonymous",
    "G" = "Nonsynonymous",  "T" = "Nonsense" ),
  "R" = c( "A" = "Nonsense",  "C" = "Nonsynonymous",
    "G" = "Nonsense",  "T" = "Nonsynonymous" ),
  "S" = c( "A" = "Nonsense",  "C" = "Nonsynonymous",
    "G" = "Nonsynonymous",  "T" = "Nonsynonymous" ),
  "T" = c( "A" = "Nonsense",  "C" = "Synonymous",
    "G" = "Nonsynonymous",  "T" = "Synonymous" ),
  "X" = c( "A" = "Nonsense",  "C" = "Nonsense",
    "G" = "Nonsense",  "T" = "Nonsense" ),
  "W" = c("A" = NA,  "C" = NA,
    "G" = NA, "T" = NA , "R" = NA,
    "Y" = NA, "S" = NA, "W" = NA,
    "K" = NA, "M" = NA, "B" = NA,
    "D" = NA, "H" = NA, "V" = NA,
    "N" = NA, "X" = NA),
  "Y" = c("A" = "low coverage",  "C" = "low coverage",
    "G" = "low coverage", "T" = "low coverage" , "R" = "low coverage",
    "Y" = "low coverage", "S" = "low coverage", "W" = "low coverage",
    "K" = "low coverage", "M" = "low coverage", "B" = "low coverage",
    "D" = "low coverage", "H" = "low coverage", "V" = "low coverage",
      "N" = "low coverage", "X" = "low coverage"),
  "Z" = c("A" = "outside ORF",  "C" = "outside ORF",
    "G" = "outside ORF", "T" = "outside ORF" , "R" = "outside ORF",
    "Y" = "outside ORF", "S" = "outside ORF", "W" = "outside ORF",
    "K" = "outside ORF", "M" = "outside ORF", "B" = "outside ORF",
    "D" = "outside ORF", "H" = "outside ORF", "V" = "outside ORF",
    "N" = "outside ORF", "X" = "outside ORF")
  );


get.coding.seq <- function(transcript, start, end, strand){
    coding <- tolower(substr(transcript, as.numeric(start),
                             as.numeric(end)))
    if(nchar(coding)%%3!=0){
        warning("coding region has not multiple of 3 lenght:\n",
                transcript, "\tstart:", start,
                "\tend:", end, "\tstrand:", strand)
    }
    if(strand%in%"-"){
        coding <- revcom(coding) # revcom  on the minus strand
    }
    return(coding)
}

get.ontology <- function(transcript, start, end, strand) {
    utr1 <- substr(transcript, 1,
                    as.numeric(start) - 1)
    utr2 <- substr(transcript,
                    as.numeric(end)+1,
                   nchar(transcript))
    coding <- get.coding.seq(transcript, start, end, strand)
    ## split the cds by 3 and get the base ontology for each codon
    codons <- substring(coding,
                        seq(1, nchar(coding), by=3),
                        seq(3,nchar(coding), by=3))
    ont <- lapply(codons, base.ontology.encode)    
    ontology <- paste(ont, collapse="")
    utr1 <- gsub("\\w", "Z", utr1)
    utr2 <- gsub("\\w", "Z", utr2)
    if(strand%in%"+"){
        return(paste(utr1, ontology, utr2, sep=""))
    }
    ### Need to strReverse the wrong way round cds
    if(strand%in%"-"){
        return(paste(utr1, strReverse(ontology), utr2, sep=""))
    }
}

cds.gff.df$base.ontology <-
    as.character(apply(cds.gff.df, 1, function (x) {
        get.ontology(x["transcript"], x["start"],
                     x["end"], x["strand"])}))

## all ontologies have the right length
table(nchar(as.character(cds.gff.df$transcript))==
      nchar(as.character(cds.gff.df$base.ontology)))

VCF$SNP <- rownames(VCF)
VAR <- merge(cds.gff.df, VCF, by.x = "seqnames", by.y = "V1")

VAR$transcript <- as.character(VAR$transcript)

transversion.transition <- function (from, to){
  get.trans <- function (x) {
    trans <- summary(as.factor(from):as.factor(to))
    trans <- trans[trans!=0]
  }
  transf <- get.trans(VARobj)
  get.vers <- function (x){
    transitions <- c("A:G", "G:A", "C:T", "T:C")
    sition <- sum(x[names(x)%in%transitions])
    version <- sum(x[!names(x)%in%transitions])
    res <- cbind(transitions=sition, transversions=version, ratioTS.TV=sition/version)
    return(as.data.frame(res))
  }
  get.vers(transf)
}

transversion.transition(VAR$V4, VAR$V5)


get.effect <- function (ontology, base, to, strand){
    if(strand%in%"-"){
        to <- revcom(to)
    }
    code <- unlist(strsplit(ontology, ""))[[base]]
    base.ontology.decode[[code]][to]
}

VAR$effect <- apply(VAR, 1, function(x){
                        get.effect(as.character(x["base.ontology"]),
                                   as.numeric(x["V2"]),
                                   as.character(x["V5"]),
                                   as.character(x["strand"]))
})

tapply(VAR$effect, as.character(VAR$strand), table)


get.sites <- function (ontology, transcript){
  s.sites <- sapply(1:length(ontology), function (i) {
     split.ont <- unlist(strsplit(ontology[[i]], ""))
     split.cod <- unlist(strsplit(transcript[[i]], ""))
     decoded <- lapply(split.ont, function (a){
       base.ontology.decode[[a]]})
     reduced.decoded <- lapply(1:length(decoded), function (x) {
       subset(decoded[[x]], names(decoded[[x]])!=toupper(split.cod[[x]]))})
     s <- lapply(reduced.decoded, function (w) {
       length(w[grepl("Syn", w)])/length(w[grepl("Syn|Non", w)])})
     n <- lapply(reduced.decoded, function (w) {
       length(w[grepl("Non", w)])/length(w[grepl("Syn|Non", w)])})
     nsyn.sites <- sum(unlist(n), na.rm=TRUE)
     syn.sites <- sum(unlist(s), na.rm=TRUE)
     cbind(nsyn.sites, syn.sites)
   })
  data.frame(t(s.sites))
}

sites <- get.sites(VAR[!duplicated(VAR$seqnames), "base.ontology"],
                   VAR[!duplicated(VAR$seqnames), "transcript"])

rownames(sites) <-  VAR$seqnames[!duplicated(VAR$seqnames)]
names(sites) <- c("nsyn.sites", "syn.sites")

VAR <- merge(VAR, sites, by.x = "seqnames", by.y = 0)

get.pn.ps <- function(VARobj){
    npSNP <- nrow(VARobj)
    uni <- !duplicated(as.character(VARobj[,"seqnames" ]))
    pN <- nrow(VARobj[grepl("Non*", VARobj$effect),])/
        sum(as.numeric(VARobj[uni, "nsyn.sites"]), na.rm=T)
    pS <- nrow(VARobj[VARobj$effect=="Synonymous",])/
        sum(as.numeric(VARobj[uni, "syn.sites"]), na.rm=T)
    pNpS <- pN/pS
    return(list(npSNP=npSNP, pN=pN, pS=pS, pNpS=pNpS))
}

pn.ps.overall <- get.pn.ps(VAR)

## use only SNPs which are ther in at least 3 genotypes
VAR.div <- VAR[VAR$SNP%in%rownames(GT[rowSums(GT,na.rm = TRUE)>2, ]),]
get.pn.ps(VAR.div)

pn.ps.list <- by(VAR, VAR$seqnames, function (x){
    list(get.pn.ps(x))
})

pn.ps.list.div <- by(VAR.div, VAR.div$seqnames, function (x){
    list(get.pn.ps(x))
})

pn.ps.frame <- do.call(rbind, pn.ps.list)
pn.ps.frame <- as.data.frame(apply(pn.ps.frame, 2, unlist))

## pn.ps.frame.div <- do.call(rbind, pn.ps.list.div)
## pn.ps.frame.div <- as.data.frame(apply(pn.ps.frame.div, 2, unlist))

pn.ps.frame$pNpS[is.infinite(pn.ps.frame$pNpS)] <- 7
pn.ps.frame$pNpS[is.nan(pn.ps.frame$pNpS)] <- 0
pn.ps.frame <- pn.ps.frame[!is.na(pn.ps.frame$npSNP),]


pn.ps.transset <- rownames(subset(pn.ps.frame, pNpS>1 & npSNP>3))
pn.ps.geneset <- gsub("_seq\\d+$", "", pn.ps.transset)

## pn.ps.transset.div <- rownames(subset(pn.ps.frame.div, pNpS>1&npSNP>3))
## pn.ps.geneset.div <- gsub("_seq\\d+$", "", pn.ps.transset.div)



