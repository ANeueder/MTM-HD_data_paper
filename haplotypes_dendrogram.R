
  # set bases
    A <- 'A'
    T <- 'T'
    G <- 'G'
    C <- 'C'

    # rs149109767
      R <- 'AGAG'
      D <- 'A'

  # rs IDs
    rs <- c('rs2857845','rs2471347','rs2798296',
      'rs3856973','rs2285086','rs10015979','rs2071655',
      'rs363082','rs363066','rs6855981','rs11731237',
      'rs363096','rs2298969','rs363092','rs916171',
      'rs82333','rs110501','rs149109767','rs362272',
      'rs2269499','rs3095073')

  # haplotype matrix
    haps.old <- data.frame(
      Hap01 = c(A,A,G,G,A,G,T,T,T,G,T,T,A,C,C,A,T,D,G,C,G),
      Hap10 = c(A,A,G,G,A,G,T,C,T,G,T,T,A,C,C,A,T,D,G,C,G),
      Hap05 = c(T,A,G,G,A,G,T,T,T,G,T,T,A,C,C,A,T,D,G,C,G),
      Hap11 = c(A,A,G,G,A,G,T,T,T,G,T,T,A,C,C,A,T,R,G,C,A),
      Hap12 = c(A,A,G,G,A,G,T,T,T,G,T,T,A,C,C,A,T,R,G,C,G),
      Hap03 = c(A,G,A,G,A,G,T,T,T,G,T,T,A,C,C,A,T,R,G,C,G),
      Hap15 = c(A,A,G,G,A,G,T,T,T,G,C,T,A,C,C,A,T,R,G,C,G),
      Hap06 = c(A,A,G,G,A,A,T,T,T,G,C,C,A,C,C,A,T,R,G,C,G),
      Hap02 = c(T,A,G,G,A,A,T,T,T,G,C,C,A,C,C,A,T,R,G,C,G),
      Hap07 = c(T,A,G,G,A,A,T,T,T,G,C,C,A,C,C,A,T,R,G,C,A),
      Hap09 = c(T,A,G,G,A,A,T,T,T,G,C,C,A,C,C,A,T,R,A,T,A),
      Hap13 = c(A,A,A,G,A,A,T,T,T,G,C,T,G,C,C,A,T,R,G,C,G),
      Hap04 = c(T,A,G,A,G,A,G,C,T,G,C,C,G,A,G,G,C,R,G,C,G),
      Hap16 = c(A,A,G,A,G,A,G,C,T,G,C,C,G,A,G,G,C,R,G,C,G),
      Hap08 = c(A,G,A,A,G,A,G,T,G,A,C,C,G,A,G,G,C,R,A,T,A),
      Hap14 = c(A,A,A,A,G,A,G,T,T,A,C,C,G,A,G,G,C,R,G,C,G),
      stringsAsFactors = FALSE, row.names = rs)

    haps <- read.table('library/hap.definition.170115.for.KC.csv',
      sep = ',', header = TRUE, row.names = 1, stringsAsFactors = FALSE)
    all(colnames(haps) == rs)
    haps$rs149109767 <- sub('I', 'AGAG', haps$rs149109767)
    haps <- t(haps)

  # now determine haplotype in test sample
    suppressPackageStartupMessages(library(circlize))
    suppressPackageStartupMessages(library(dendextend))

  # read in predicted mHTT phase (corrected / reviewed by Jian)
    htt_hap <- read.table('library/predicted_mHTT_phase.tsv', header = FALSE, stringsAsFactors = FALSE)
    htt_hap$HTT <- ifelse(htt_hap[,3] == 'HP1', 'HP2',
      ifelse(htt_hap[,3] == 'HP2', 'HP1',
        ifelse(htt_hap[,3] == 'inconclusive', NA, NA)))
    htt_hap <- htt_hap[,-2]
    colnames(htt_hap) <- c('ID', 'mHTT', 'HTT')
    htt_hap$mHTT[htt_hap$mHTT == 'inconclusive'] <- NA

  # function to determine nearest haplotype match if non-perfect match
    nearestmatch <- function(HP, haps) {
      tab <- apply((HP == haps), 2, function(x) table(x)[c('TRUE', 'FALSE')])
      tab[is.na(tab)] <- 0
      tab <- tab[,order(tab[1,], decreasing = TRUE)]
      max <- tab[1,1]
      #ret <- paste(colnames(tab[,which(tab[1,] == max)]), collapse = ';')
      ret <- colnames(tab)[1]
      return(ret)
    }

  # function to determine individual genotype matches (or not)
    matchgenotypes <- function(HP, haps, match) {
      seq <- data.frame(HP, haps[,match])
      seq[is.na(seq)] <- 'na'
      seq <- apply(seq, 1, function(x) {
        ifelse(x[1] == x[2],
          paste(x, collapse = '=='),
          paste(x, collapse = ' != '))})
      return(seq)
    }

  # loop through HTT genotypes and infer haplotype
    genotypes <- list.files('haplotyping/', pattern = '*ChaoGenotype*', recursive = TRUE, full.name = TRUE)
    ids <- unlist(lapply(strsplit(genotypes, '/'), function(x) x[3]))

    write.table(paste(c('Sample', 'mHTT', 'HTT',
      'HP1 exact match', 'HP2 exact match',
      'HP1 nearest match', 'HP2 nearest match'), collapse = ','),
      'PhaseHaplotypesv1.csv',
      row.names = FALSE, quote = FALSE, col.names = FALSE)

    write.table(paste(c('Sample', 'Phase', rs), collapse = ','),
      'GenotypeMatches.csv',
      row.names = FALSE, quote = FALSE, col.names = FALSE)

    res <- htt_hap
    res$HP1 <- rep(NA, nrow(res))
    res$HP2 <- rep(NA, nrow(res))
    res$HP1_NearestMatch <- rep(NA, nrow(res))
    res$HP2_NearestMatch <- rep(NA, nrow(res))
    for (i in 1:length(genotypes)) {
      message('--genotypes file is:\t', genotypes[i])
      message('--sample ID is:\t', ids[i])

      sam <- read.table(genotypes[i], header = TRUE, sep = '\t', stringsAsFactors = FALSE)

      bIncomplete <- FALSE

      if (nrow(sam<21)) bIncomplete <- TRUE

      sam <- sam[match(rs, sam$rs),]

      hp1_match <- names(which(apply(sam$HP1 == haps, 2, all) == TRUE))
      hp2_match <- names(which(apply(sam$HP2 == haps, 2, all) == TRUE))
      hp1_nearestmatch <- NA
      hp2_nearestmatch <- NA
      hp1_matchgenotypes <- matchgenotypes(sam$HP1, haps, hp1_match)
      hp2_matchgenotypes <- matchgenotypes(sam$HP2, haps, hp2_match)

      if (length(hp1_match) == 0 & bIncomplete) {
        hp1_match <- 'NA'
        hp1_nearestmatch <- nearestmatch(sam$HP1, haps)
        hp1_matchgenotypes <- matchgenotypes(sam$HP1, haps,
          unlist(lapply(strsplit(hp1_nearestmatch, ';'), function(x) x[1])))
      } else if (length(hp1_match) == 0 & !bIncomplete) {
        hp1_match <- 'NA'
        hp1_nearestmatch <- nearestmatch(sam$HP1, haps)
        hp1_matchgenotypes <- matchgenotypes(sam$HP1, haps,
          unlist(lapply(strsplit(hp1_nearestmatch, ';'), function(x) x[1])))
      }

      if (length(hp2_match) == 0 & bIncomplete) {
        hp2_match <- 'NA'
        hp2_nearestmatch <- nearestmatch(sam$HP2, haps)
        hp2_matchgenotypes <- matchgenotypes(sam$HP2, haps,
          unlist(lapply(strsplit(hp2_nearestmatch, ';'), function(x) x[1])))
      } else if (length(hp2_match) == 0 & !bIncomplete) {
        hp2_match <- 'NA'
        hp2_nearestmatch <- nearestmatch(sam$HP2, haps)
        hp2_matchgenotypes <- matchgenotypes(sam$HP2, haps,
          unlist(lapply(strsplit(hp2_nearestmatch, ';'), function(x) x[1])))
      }

      message('--inferred haplotypes:\t', hp1_match, '|', hp2_match)
      res$HP1[i] <- hp1_match
      res$HP2[i] <- hp2_match

      message('--inferred nearest haplotypes:\t', hp1_nearestmatch, '|', hp2_nearestmatch)
      res$HP1_NearestMatch[i] <- hp1_nearestmatch
      res$HP2_NearestMatch[i] <- hp2_nearestmatch

      message('--genotypes aligned for output:\t',
        all(c(names(hp1_matchgenotypes) == rs, names(hp2_matchgenotypes) == rs)),
        '\n')

      write.table(paste(c(ids[i], 'HP1', hp1_matchgenotypes), collapse = ','),
        'GenotypeMatches.csv',
        row.names = FALSE, quote = FALSE, col.names = FALSE, append = TRUE)
      write.table(paste(c(ids[i], 'HP2', hp2_matchgenotypes), collapse = ','),
        'GenotypeMatches.csv',
        row.names = FALSE, quote = FALSE, col.names = FALSE, append = TRUE)
    }

  # write out results  
    write.table(res, 'PhaseHaplotypesv1.csv',
      row.names = FALSE, quote = FALSE, col.names = FALSE, sep = ',', append = TRUE)

    a <- subset(res, mHTT == 'HP1')[,c('ID','mHTT','HP1','HP1_NearestMatch')]
    b <- subset(res, HTT == 'HP1')[,c('ID','HTT','HP1','HP1_NearestMatch')]
    a$Allele <- 'mHTT'
    b$Allele <- 'HTT'
    colnames(a) <- c('Sample', 'Phase', 'Exact Haplotype', 'Nearest Match Haplotype', 'Allele')
    colnames(b) <- c('Sample', 'Phase', 'Exact Haplotype', 'Nearest Match Haplotype', 'Allele')

    c <- subset(res, mHTT == 'HP2')[,c('ID','mHTT','HP2','HP2_NearestMatch')]
    d <- subset(res, HTT == 'HP2')[,c('ID','HTT','HP2','HP2_NearestMatch')]
    c$Allele <- 'mHTT'
    d$Allele <- 'HTT'
    colnames(c) <- c('Sample', 'Phase', 'Exact Haplotype', 'Nearest Match Haplotype', 'Allele')
    colnames(d) <- c('Sample', 'Phase', 'Exact Haplotype', 'Nearest Match Haplotype', 'Allele')

    res2 <- rbind(a,b,c,d)[,c('Sample', 'Phase', 'Allele', 'Exact Haplotype', 'Nearest Match Haplotype')]
    res2 <- res2[order(res2$Sample, res2$Phase),]

    write.table(res2, 'PhaseHaplotypesv2.csv',
      row.names = FALSE, quote = FALSE, col.names = TRUE, sep = ',')

  # generate cluster dendrograms
    # a binary matrix of 0 (GRCh38 ref) and 1 (variant) will be produced
    sam <- read.table(genotypes[2], header = TRUE, sep = '\t', stringsAsFactors = FALSE)
    binarymat <- ifelse(sam$GRCh38 == haps, 0, 1)
    for (i in 1:length(genotypes)) {
      message('--genotypes file is:\t', genotypes[i])
      message('--sample ID is:\t', ids[i])

      sam <- read.table(genotypes[i], header = TRUE, sep = '\t', stringsAsFactors = FALSE)

      sam <- sam[match(rs, sam$rs),]

      binarymat <- cbind(binarymat, ifelse(sam$GRCh38 == sam$HP1, 0, 1))
      binarymat <- cbind(binarymat, ifelse(sam$GRCh38 == sam$HP2, 0, 1))
      colnames(binarymat)[(ncol(binarymat)-1):ncol(binarymat)] <- c(paste0(ids[i], '.HP1'), paste0(ids[i], '.HP2'))
    }

    # cluster via binary distance
      binarymat[is.na(binarymat)] <- 0
      d <- dist(t(data.frame(binarymat)), method = 'binary')
      hc <- hclust(d, method = 'ward.D2')
      dend <- as.dendrogram(hc)

      pdf('Dendrogram.pdf', width = 23, height = 10)
        par(mfrow=c(1,1))

        # Get the heights for each branch
          heights <- round(get_branches_heights(dend, sort = FALSE), 1)

        # Get max height
          maxHeight <- max(heights)

        # Set label and dendrogram height for circular dendrogram
          labelHeight <- 0.1
          dendHeight <- 0.8

        # Draw the circular dendrogram
          dend.col <- dend %>%
            set('branches_k_color', k = 3,
              value = c('red2', 'forestgreen', 'royalblue')) %>% 
            set('branches_lwd', 2) %>%
            set('labels_colors', k = 3, value = c('red2', 'forestgreen', 'royalblue')) %>% 
            set('labels_cex', 0.5)
          circlize_dendrogram(dend.col,
            facing = 'outside',
            labels = TRUE,
            labels_track_height = labelHeight,
            dend_track_height = dendHeight)

         # Create tick co-ordinates and values for the new axis
         # We have to ensure that we don't overlap the label plot region
         #   (height specified by labelHeight), nor the central region of the
         #   plot (1-(dendHeight+labelHeight))
           ticks <- seq(from = (1-(dendHeight+labelHeight)),
             to = (1-labelHeight), length.out=5)
           values <- round(rev(seq(from=0, to=maxHeight, length.out=5)), 1)

         # Add the new axis
           suppressPackageStartupMessages(library(plotrix))
           ablineclip(h = 0, v = ticks, col = 'black',
             x1 = 1-(dendHeight+labelHeight),
             x2 = 1-labelHeight,
             y1 = 0,
             y2 = 0.04,
             lwd = 2)
           text(ticks, 0+0.08, values, cex = 0.8)
           text(
             (1-labelHeight)-(((1-labelHeight)-(1-(dendHeight+labelHeight)))/2),
             0+0.2,
             'Distance', cex = 1.2)

        # regular dendrogram
          labels_cex(dend) <- 0.4
          plot(color_labels(
            dend,
            col = ifelse(grepl('HP1|HP2', labels(dend)), 2, 1)),
            main = '',
            sub = 'Distance & linkage methods: binary distance with Ward\'s linkage\nTree cut for 3 clusters (first page)',
            cex.sub = 1.2)
      dev.off()

  # generate cluster dendrograms: wild-type HTT
    # original IDs
      # a binary matrix of 0 (GRCh38 ref) and 1 (variant) will be produced
      sam <- read.table(genotypes[2], header = TRUE, sep = '\t', stringsAsFactors = FALSE)
      binarymat <- ifelse(sam$GRCh38 == haps, 0, 1)
      for (i in 1:length(genotypes)) {
        message('--genotypes file is:\t', genotypes[i])
        message('--sample ID is:\t', ids[i])

        sam <- read.table(genotypes[i], header = TRUE, sep = '\t', stringsAsFactors = FALSE)

        sam <- sam[match(rs, sam$rs),]

        binarymat <- cbind(binarymat, ifelse(sam$GRCh38 == sam$HP1, 0, 1))
        binarymat <- cbind(binarymat, ifelse(sam$GRCh38 == sam$HP2, 0, 1))
        colnames(binarymat)[(ncol(binarymat)-1):ncol(binarymat)] <- c(paste0(ids[i], '.HP1'), paste0(ids[i], '.HP2'))
      }

      # filter out just wild-type HTT
        idx <- c(
          grep('hap', colnames(binarymat)),
          which(colnames(binarymat) %in% apply(subset(res2, Allele == 'HTT')[,1:2], 1, function(x) paste(x, collapse = '.'))))
        binarymat <- binarymat[,idx]

      # cluster via binary distance
        binarymat[is.na(binarymat)] <- 0
        d <- dist(t(data.frame(binarymat)), method = 'binary')
        hc <- hclust(d, method = 'ward.D2')
        dend <- as.dendrogram(hc)

        pdf('Dendrogram_HTT.pdf', width = 21, height = 9.5)
          par(mfrow=c(1,1))

          # Get the heights for each branch
            heights <- round(get_branches_heights(dend, sort = FALSE), 1)

          # Get max height
            maxHeight <- max(heights)

          # Set label and dendrogram height for circular dendrogram
            labelHeight <- 0.1
            dendHeight <- 0.8

          # Draw the circular dendrogram
            dend.col <- dend %>%
              set('branches_k_color', k = 3,
                value = c('red2', 'forestgreen', 'royalblue')) %>% 
              set('branches_lwd', 2) %>%
              set('labels_colors', k = 3, value = c('red2', 'forestgreen', 'royalblue')) %>% 
              set('labels_cex', 0.5)
            circlize_dendrogram(dend.col,
              facing = 'outside',
              labels = TRUE,
              labels_track_height = labelHeight,
              dend_track_height = dendHeight)
  
           # Create tick co-ordinates and values for the new axis
           #   We have to ensure that we don't overlap the label plot region
           #     (height specified by labelHeight), nor the central region of the
           #   plot (1-(dendHeight+labelHeight))
             ticks <- seq(from = (1-(dendHeight+labelHeight)),
               to = (1-labelHeight), length.out=5)
             values <- round(rev(seq(from=0, to=maxHeight, length.out=5)), 1)

           # Add the new axis
             suppressPackageStartupMessages(library(plotrix))
             ablineclip(h = 0, v = ticks, col = 'black',
               x1 = 1-(dendHeight+labelHeight),
               x2 = 1-labelHeight,
               y1 = 0,
               y2 = 0.04,
               lwd = 2)
               text(ticks, 0+0.08, values, cex = 0.8)
             text(
               (1-labelHeight)-(((1-labelHeight)-(1-(dendHeight+labelHeight)))/2),
               0+0.2,
               'Distance', cex = 1.2)

          # regular dendrogram
            labels_cex(dend) <- 0.4
            plot(color_labels(
              dend,
              col = ifelse(grepl('HP1|HP2', labels(dend)), 2, 1)),
              main = 'Wild-type HTT',
              sub = 'Distance & linkage methods: binary distance with Ward\'s linkage\nTree cut for 3 clusters (first page)',
              cex.sub = 1.2)
        dev.off()

    # new IDs
      lookup <- data.frame(readxl::read_xlsx(
        'library/Cryopreserved Lines at Roslin 01-02-2020.xlsx')[,c('Foundation ID', 'Repository ID')])

      # a binary matrix of 0 (GRCh38 ref) and 1 (variant) will be produced
      sam <- read.table(genotypes[2], header = TRUE, sep = '\t', stringsAsFactors = FALSE)
      binarymat <- ifelse(sam$GRCh38 == haps, 0, 1)
      for (i in 1:length(genotypes)) {
        message('--genotypes file is:\t', genotypes[i])
        message('--sample ID is:\t', ids[i])

        sam <- read.table(genotypes[i], header = TRUE, sep = '\t', stringsAsFactors = FALSE)

        sam <- sam[match(rs, sam$rs),]

        binarymat <- cbind(binarymat, ifelse(sam$GRCh38 == sam$HP1, 0, 1))
        binarymat <- cbind(binarymat, ifelse(sam$GRCh38 == sam$HP2, 0, 1))
        colnames(binarymat)[(ncol(binarymat)-1):ncol(binarymat)] <- c(paste0(ids[i], '.HP1'), paste0(ids[i], '.HP2'))
      }

      # filter out just wild-type HTT
        idx <- c(
          grep('hap', colnames(binarymat)),
          which(colnames(binarymat) %in% apply(subset(res2, Allele == 'HTT')[,1:2], 1, function(x) paste(x, collapse = '.'))))
        binarymat <- binarymat[,idx]

      # lookup up new IDs
        colnames(binarymat) <- sub('_pilot\\.HP[12]$|\\.HP[12]', '', colnames(binarymat))
        colnames(binarymat) <- unlist(lapply(colnames(binarymat), function(x) {
          ifelse(any(grepl(x, lookup$Repository.ID)),
            lookup[which(lookup$Repository.ID == x),'Foundation.ID'],
            x)}))

      # cluster via binary distance
        binarymat[is.na(binarymat)] <- 0
        d <- dist(t(data.frame(binarymat)), method = 'binary')
        hc <- hclust(d, method = 'ward.D2')
        dend <- as.dendrogram(hc)
        labels(dend) <- sub('CHDI\\.', 'CHDI\\-', labels(dend))

        tiff('Dendrogram_HTT.dpi300.tiff', units = 'in',
          res = 300, width = 21, height = 9.5)
          par(mfrow=c(1,1))

          # Get the heights for each branch
            heights <- round(get_branches_heights(dend, sort = FALSE), 1)

          # Get max height
            maxHeight <- max(heights)

          # Set label and dendrogram height for circular dendrogram
            labelHeight <- 0.1
            dendHeight <- 0.8

          # Draw the circular dendrogram
            dend.col <- dend %>%
              set('branches_k_color', k = 3,
                value = c('red2', 'forestgreen', 'royalblue')) %>% 
              set('branches_lwd', 3) %>%
              set('labels_colors', k = 3, value = c('red2', 'forestgreen', 'royalblue')) %>% 
              set('labels_cex', 0.6)
            circlize_dendrogram(dend.col,
              facing = 'outside',
              labels = TRUE,
              labels_track_height = labelHeight,
              dend_track_height = dendHeight)
  
           # Create tick co-ordinates and values for the new axis
             ticks <- seq(from = 0.09,
               to = 0.86, length.out=5, lwd = 4)
             values <- round(rev(seq(from=0, to=maxHeight, length.out=5)), 1)

           # Add the new axis
             suppressPackageStartupMessages(library(plotrix))
             ablineclip(h = 0, v = ticks, col = 'black',
               x1 = 0.09,
               x2 = 0.86,
               y1 = 0,
               y2 = 0.04,
               lwd = 4)
               text(ticks, 0+0.08, values, cex = 1.4)
             text(
               (1-labelHeight)-(((1-labelHeight)-(1-(dendHeight+labelHeight)))/2),
               0+0.175,
               'Distance', cex = 1.4)
        dev.off()

        tiff('Dendrogram_HTT.dpi100.tiff', units = 'in',
          res = 100, width = 21, height = 9.5)
          par(mfrow=c(1,1))

          # Get the heights for each branch
            heights <- round(get_branches_heights(dend, sort = FALSE), 1)

          # Get max height
            maxHeight <- max(heights)

          # Set label and dendrogram height for circular dendrogram
            labelHeight <- 0.1
            dendHeight <- 0.8

          # Draw the circular dendrogram
            dend.col <- dend %>%
              set('branches_k_color', k = 3,
                value = c('red2', 'forestgreen', 'royalblue')) %>% 
              set('branches_lwd', 3) %>%
              set('labels_colors', k = 3, value = c('red2', 'forestgreen', 'royalblue')) %>% 
              set('labels_cex', 0.6)
            circlize_dendrogram(dend.col,
              facing = 'outside',
              labels = TRUE,
              labels_track_height = labelHeight,
              dend_track_height = dendHeight)
  
           # Create tick co-ordinates and values for the new axis
             ticks <- seq(from = 0.09,
               to = 0.86, length.out=5, lwd = 4)
             values <- round(rev(seq(from=0, to=maxHeight, length.out=5)), 1)

           # Add the new axis
             suppressPackageStartupMessages(library(plotrix))
             ablineclip(h = 0, v = ticks, col = 'black',
               x1 = 0.09,
               x2 = 0.86,
               y1 = 0,
               y2 = 0.04,
               lwd = 4)
               text(ticks, 0+0.08, values, cex = 1.4)
             text(
               (1-labelHeight)-(((1-labelHeight)-(1-(dendHeight+labelHeight)))/2),
               0+0.175,
               'Distance', cex = 1.4)
        dev.off()



  # generate cluster dendrograms: mutant HTT
    # original IDs
      # a binary matrix of 0 (GRCh38 ref) and 1 (variant) will be produced
      sam <- read.table(genotypes[2], header = TRUE, sep = '\t', stringsAsFactors = FALSE)
      binarymat <- ifelse(sam$GRCh38 == haps, 0, 1)
      for (i in 1:length(genotypes)) {
        message('--genotypes file is:\t', genotypes[i])
        message('--sample ID is:\t', ids[i])

        sam <- read.table(genotypes[i], header = TRUE, sep = '\t', stringsAsFactors = FALSE)

        sam <- sam[match(rs, sam$rs),]

        binarymat <- cbind(binarymat, ifelse(sam$GRCh38 == sam$HP1, 0, 1))
        binarymat <- cbind(binarymat, ifelse(sam$GRCh38 == sam$HP2, 0, 1))
        colnames(binarymat)[(ncol(binarymat)-1):ncol(binarymat)] <- c(paste0(ids[i], '.HP1'), paste0(ids[i], '.HP2'))
      }

      # filter out just mutant HTT
        idx <- c(
          grep('hap', colnames(binarymat)),
          which(colnames(binarymat) %in% apply(subset(res2, Allele == 'mHTT')[,1:2], 1, function(x) paste(x, collapse = '.'))))
        binarymat <- binarymat[,idx]

      # cluster via binary distance
        binarymat[is.na(binarymat)] <- 0
        d <- dist(t(data.frame(binarymat)), method = 'binary')
        hc <- hclust(d, method = 'ward.D2')
        dend <- as.dendrogram(hc)

        pdf('Dendrogram_mHTT.pdf', width = 21, height = 9.5)
          par(mfrow=c(1,1))

          # Get the heights for each branch
            heights <- round(get_branches_heights(dend, sort = FALSE), 1)

          # Get max height
            maxHeight <- max(heights)

          # Set label and dendrogram height for circular dendrogram
            labelHeight <- 0.1
            dendHeight <- 0.8

          # Draw the circular dendrogram
            dend.col <- dend %>%
              set('branches_k_color', k = 3,
                value = c('red2', 'forestgreen', 'royalblue')) %>% 
              set('branches_lwd', 2) %>%
              set('labels_colors', k = 3, value = c('red2', 'forestgreen', 'royalblue')) %>% 
              set('labels_cex', 0.5)
            circlize_dendrogram(dend.col,
              facing = 'outside',
              labels = TRUE,
              labels_track_height = labelHeight,
              dend_track_height = dendHeight)

            # Create tick co-ordinates and values for the new axis
            # We have to ensure that we don't overlap the label plot region
            #   (height specified by labelHeight), nor the central region of the
            #   plot (1-(dendHeight+labelHeight))
              ticks <- seq(from = (1-(dendHeight+labelHeight)),
                to = (1-labelHeight), length.out=5)
              values <- round(rev(seq(from=0, to=maxHeight, length.out=5)), 1)

            # Add the new axis
              suppressPackageStartupMessages(library(plotrix))
              ablineclip(h = 0, v = ticks, col = 'black',
                x1 = 1-(dendHeight+labelHeight),
                x2 = 1-labelHeight,
                y1 = 0,
                y2 = 0.04,
                lwd = 2)
              text(ticks, 0+0.08, values, cex = 0.8)
              text(
                (1-labelHeight)-(((1-labelHeight)-(1-(dendHeight+labelHeight)))/2),
                0+0.2,
                'Distance', cex = 1.2)

            # regular dendrogram
              labels_cex(dend) <- 0.4
              plot(color_labels(
              dend,
              col = ifelse(grepl('HP1|HP2', labels(dend)), 2, 1)),
              main = 'Mutant HTT',
              sub = 'Distance & linkage methods: binary distance with Ward\'s linkage\nTree cut for 3 clusters (first page)',
              cex.sub = 1.2)
        dev.off()

    # new IDs
      lookup <- data.frame(readxl::read_xlsx(
        'library/Cryopreserved Lines at Roslin 01-02-2020.xlsx')[,c('Foundation ID', 'Repository ID')])

      # a binary matrix of 0 (GRCh38 ref) and 1 (variant) will be produced
      sam <- read.table(genotypes[2], header = TRUE, sep = '\t', stringsAsFactors = FALSE)
      binarymat <- ifelse(sam$GRCh38 == haps, 0, 1)
      for (i in 1:length(genotypes)) {
        message('--genotypes file is:\t', genotypes[i])
        message('--sample ID is:\t', ids[i])

        sam <- read.table(genotypes[i], header = TRUE, sep = '\t', stringsAsFactors = FALSE)

        sam <- sam[match(rs, sam$rs),]

        binarymat <- cbind(binarymat, ifelse(sam$GRCh38 == sam$HP1, 0, 1))
        binarymat <- cbind(binarymat, ifelse(sam$GRCh38 == sam$HP2, 0, 1))
        colnames(binarymat)[(ncol(binarymat)-1):ncol(binarymat)] <- c(paste0(ids[i], '.HP1'), paste0(ids[i], '.HP2'))
      }

      # filter out just mutant HTT
        idx <- c(
          grep('hap', colnames(binarymat)),
          which(colnames(binarymat) %in% apply(subset(res2, Allele == 'mHTT')[,1:2], 1, function(x) paste(x, collapse = '.'))))
        binarymat <- binarymat[,idx]

      # lookup up new IDs
        colnames(binarymat) <- sub('_pilot\\.HP[12]$|\\.HP[12]', '', colnames(binarymat))
        colnames(binarymat) <- unlist(lapply(colnames(binarymat), function(x) {
          ifelse(any(grepl(x, lookup$Repository.ID)),
            lookup[which(lookup$Repository.ID == x),'Foundation.ID'],
            x)}))

      # cluster via binary distance
        binarymat[is.na(binarymat)] <- 0
        d <- dist(t(data.frame(binarymat)), method = 'binary')
        hc <- hclust(d, method = 'ward.D2')
        dend <- as.dendrogram(hc)
        labels(dend) <- sub('CHDI\\.', 'CHDI\\-', labels(dend))

        tiff('Dendrogram_mHTT.dpi300.tiff', units = 'in',
          res = 300, width = 21, height = 9.5)
          par(mfrow=c(1,1))

          # Get the heights for each branch
            heights <- round(get_branches_heights(dend, sort = FALSE), 1)

          # Get max height
            maxHeight <- max(heights)

          # Set label and dendrogram height for circular dendrogram
            labelHeight <- 0.1
            dendHeight <- 0.8

          # Draw the circular dendrogram
            dend.col <- dend %>%
              set('branches_k_color', k = 3,
                value = c('red2', 'forestgreen', 'royalblue')) %>% 
              set('branches_lwd', 3) %>%
              set('labels_colors', k = 3, value = c('red2', 'forestgreen', 'royalblue')) %>% 
              set('labels_cex', 0.6)
            circlize_dendrogram(dend.col,
              facing = 'outside',
              labels = TRUE,
              labels_track_height = labelHeight,
              dend_track_height = dendHeight)
  
           # Create tick co-ordinates and values for the new axis
             ticks <- seq(from = 0.09,
               to = 0.86, length.out=5, lwd = 4)
             values <- round(rev(seq(from=0, to=maxHeight, length.out=5)), 1)

           # Add the new axis
             suppressPackageStartupMessages(library(plotrix))
             ablineclip(h = 0, v = ticks, col = 'black',
               x1 = 0.09,
               x2 = 0.86,
               y1 = 0,
               y2 = 0.04,
               lwd = 4)
               text(ticks, 0+0.08, values, cex = 1.4)
             text(
               (1-labelHeight)-(((1-labelHeight)-(1-(dendHeight+labelHeight)))/2),
               0+0.175,
               'Distance', cex = 1.4)
        dev.off()

        tiff('Dendrogram_mHTT.dpi100.tiff', units = 'in',
          res = 100, width = 21, height = 9.5)
          par(mfrow=c(1,1))

          # Get the heights for each branch
            heights <- round(get_branches_heights(dend, sort = FALSE), 1)

          # Get max height
            maxHeight <- max(heights)

          # Set label and dendrogram height for circular dendrogram
            labelHeight <- 0.1
            dendHeight <- 0.8

          # Draw the circular dendrogram
            dend.col <- dend %>%
              set('branches_k_color', k = 3,
                value = c('red2', 'forestgreen', 'royalblue')) %>% 
              set('branches_lwd', 3) %>%
              set('labels_colors', k = 3, value = c('red2', 'forestgreen', 'royalblue')) %>% 
              set('labels_cex', 0.6)
            circlize_dendrogram(dend.col,
              facing = 'outside',
              labels = TRUE,
              labels_track_height = labelHeight,
              dend_track_height = dendHeight)
  
           # Create tick co-ordinates and values for the new axis
             ticks <- seq(from = 0.09,
               to = 0.86, length.out=5, lwd = 4)
             values <- round(rev(seq(from=0, to=maxHeight, length.out=5)), 1)

           # Add the new axis
             suppressPackageStartupMessages(library(plotrix))
             ablineclip(h = 0, v = ticks, col = 'black',
               x1 = 0.09,
               x2 = 0.86,
               y1 = 0,
               y2 = 0.04,
               lwd = 4)
               text(ticks, 0+0.08, values, cex = 1.4)
             text(
               (1-labelHeight)-(((1-labelHeight)-(1-(dendHeight+labelHeight)))/2),
               0+0.175,
               'Distance', cex = 1.4)
        dev.off()

