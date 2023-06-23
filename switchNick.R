library("Biostrings")
library("rmelting")
library("crisprVerse")
library("crisprBwa")
library("dplyr")

#dna <- readDNAStringSet("/home/files/168H2and746H13.txt")
fnm <- file.choose()
dna <- readDNAStringSet(fnm)
lowlen <- 10
highlen <- 14
lowmt <- 45
highmt <- 62
idealloc <- 100000
rangeopt <- 30000
highloc <- idealloc+rangeopt
lowloc <- idealloc-rangeopt


findNickPairs <- function(dna, lowlen, highlen, lowmt, highmt, lowloc, highloc) {
  methods <- c('-+'="pamout",'+-'="pamin","--"="-inline","++"="+inline")
  seq <- dna[[1]]
  data(SpCas9, package="crisprBase")
  allsites <- findSpacers(dna, crisprNuclease=SpCas9)
  print(paste0("All spacers in sequence: ",length(allsites)))
  sallsites <- subset(allsites, between(as.numeric(cutSites(allsites)), lowloc, highloc))
  print(paste0("Spacers within range: ",length(sallsites)))
  end <- length(sallsites)-1
  arange <- 1:(end-1)
  pairs <- data.frame(matrix(ncol = 16, nrow = 0))
  colnames(pairs) <- c("site1","site1_spacer","site1_cut","site1_strand",
                       "site2","site2_spacer","site2_cut","site2_strand",
                       "start","end","insequence","length","meltingtemp",
                       "pam_arrangement","percent_GC","greater_sequence")
  usedt <- c()
  for (n in arange) {
    site <- sallsites[n]
    cut <- as.numeric(cutSites(site))
    strand <- as.character(strand(site))
    i <- n +1
    while(TRUE){
        i <- i + 1
        if (i >= end) {break}
        fsite <- sallsites[i]
        fcut <- as.numeric(cutSites(fsite))
        len <- fcut-cut
        if (len < lowlen) {next}
        if (len > highlen) {break}
        fstrand <- as.character(strand(fsite))
        pstart <- cut + (((-1*as.numeric(paste0(strand,1)))+1)/2)
        pend <- fcut + (((-1*as.numeric(paste0(fstrand,1)))-1)/2)
        inseq <- seq[pstart:pend]
        meltt <- melting(sequence = as.character(inseq), nucleic.acid.conc = 0.00001,hybridisation.type = "dnadna", Na.conc = 0.1)
        if (meltt[[3]][[5]] < lowmt){next}
        if (meltt[[3]][[5]] > highmt){break}
        mseq <- as.character(seq[(pstart-20):(pend+20)])
        pgc <- letterFrequency(inseq,letters="GC",as.prob=TRUE)
        pairs <- rbind(pairs, c(names(site),as.character(site$protospacer),cut, strand,
                              names(fsite),as.character(fsite$protospacer),fcut, fstrand,
                              pstart,pend,as.character(inseq),length(inseq),meltt[[3]][[5]],
                              methods[(paste0(strand,fstrand))],pgc,mseq))
        usedt <- c(usedt, c(n,i))
    }
  }
  subsites <- sallsites[unique(usedt)]
  subsites <- addSpacerAlignments(subsites, aligner="biostrings", custom_seq=seq, n_mismatches=4, all_alignments=TRUE)
  ita <- names(offTargets(subsites))
  subsites <- sallsites[unique(usedt)]
  subsubsites <- subset(subsites, !(names(subsites) %in% ita))
  subsubsites <- addSequenceFeatures(subsubsites)
  subpairs <- subset(pairs, (!(pairs[,1] %in% ita)) & (!(pairs[,5] %in% ita)))
  colnames(subpairs) <- c("site1","site1_spacer","site1_cut","site1_strand",
                          "site2","site2_spacer","site2_cut","site2_strand",
                          "start","end","insequence","length","melting_temp",
                          "pam_arrangement","percent_GC","greater_sequence")
  print(paste0("Pairs without off-targets: ",nrow(subpairs)))
  subsubsites <- addOnTargetScores(subsubsites,methods=c("deephf", "deepspcas9","ruleset3"))
  subsubsites <- rankSpacers(subsubsites)
  subsubsites$name <- names(subsubsites)
  subsubsites$score_mean <- (subsubsites$score_deephf+ subsubsites$score_deepspcas9 + ((subsubsites$score_ruleset3+3)/6))/3
  subpairs <- transform(subpairs,site1_score_deephf=(subsubsites[site1]$score_deephf))
  subpairs <- transform(subpairs,site1_score_deepspcas9=(subsubsites[site1]$score_deepspcas9))
  subpairs <- transform(subpairs,site1_score_ruleset3=(subsubsites[site1]$score_ruleset3))
  subpairs <- transform(subpairs,site1_score_composite=(subsubsites[site1]$score_composite))
  subpairs <- transform(subpairs,site1_score_mean=(subsubsites[site1]$score_mean))
  subpairs <- transform(subpairs,site1_percentGC=(subsubsites[site1]$percentGC))
  subpairs <- transform(subpairs,site2_score_deephf=(subsubsites[site2]$score_deephf))
  subpairs <- transform(subpairs,site2_score_deepspcas9=(subsubsites[site2]$score_deepspcas9))
  subpairs <- transform(subpairs,site2_score_ruleset3=(subsubsites[site2]$score_ruleset3))
  subpairs <- transform(subpairs,site2_score_composite=(subsubsites[site2]$score_composite))
  subpairs <- transform(subpairs,site2_score_mean=(subsubsites[site2]$score_mean))
  subpairs <- transform(subpairs,site2_percentGC=(subsubsites[site2]$percentGC))
  subpairs <- transform(subpairs,cum_score=(subsubsites[site1]$score_composite + subsubsites[site2]$score_composite))
  subpairs <- transform(subpairs,cum_means=(subsubsites[site1]$score_mean + subsubsites[site2]$score_mean))
  subpairs <- subpairs[order(-subpairs$cum_score),]
  name <- paste0(lowloc,"-",highloc)
  write.csv(subsubsites, paste0("output/",name,"_sites.csv"), row.names=FALSE, quote=FALSE)
  write.csv(subpairs, paste0("output/",name,"_pairs.csv"), row.names=FALSE, quote=FALSE)
  return(list(subsubsites, subpairs))
}

findNickPairs(dna,lowlen, highlen, lowmt,highmt,lowloc,highloc)
