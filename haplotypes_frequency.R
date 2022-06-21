
  # generate barplots to show the number of wtHTT and mHTT alleles per binned haplotype
    options(scipen = 9)

  # wtHTT: exact match only
    df <- read.table('main/haplotyping/PhaseHaplotypesv2.csv',
      header = TRUE, sep = ',', stringsAsFactors = FALSE)

    # remove CH01855_pilot so as to not have duplicate samples
      df <- df[-which(df$Sample == 'CH01855_pilot'),]

    # wild-type only
      df <- subset(df, Allele == 'HTT')

    # summarise haps into bins
      haps <- as.numeric(sub('hap\\.', '', df$Exact.Haplotype))
      haps[is.na(haps)] <- 9999
      ggdata <- data.frame(table(
        ifelse(haps <= 6, '1-6',
          ifelse(haps <= 12, '7-12',
            ifelse(haps <= 18, '13-18',
              ifelse(haps <= 24, '19-24',
                ifelse(haps <= 30, '25-30',
                  ifelse(haps <= 36, '31-36',
                    ifelse(haps <= 42, '37-42',
                      ifelse(haps <= 48, '43-48',
                        ifelse(haps <= 54, '49-54',
                          ifelse(haps <= 60, '55-60',
                            ifelse(haps <= 66, '61-66',
                              ifelse(haps <= 72, '67-72',
                                ifelse(haps <= 78, '73-78',
                                  ifelse(haps <= 84, '79-84',
                                    ifelse(haps <= 90, '85-90',
                                      ifelse(haps <= 96, '91-96',
                                        ifelse(haps <= 102, '97-102',
                                          ifelse(haps <= 108, '103-108',
                                            NA))))))))))))))))))))
      levels <- c('1-6','7-12','13-18','19-24','25-30','31-36','37-42',
        '43-48','49-54','55-60','61-66','67-72','73-78','79-84',
        '85-90','91-96','97-102','103-108')
      ggdata <- ggdata[match(levels, ggdata$Var1),]
      ggdata$Var1 <- levels
      ggdata$Freq[is.na(ggdata$Freq)] <- 0
      ggdata$Var1 <- factor(ggdata$Var1,
        levels = c('1-6','7-12','13-18','19-24','25-30','31-36','37-42',
          '43-48','49-54','55-60','61-66','67-72','73-78','79-84',
          '85-90','91-96','97-102','103-108'))

    # plot
      require(ggplot2)
      require(ggrepel)
      mytheme <- theme_bw(base_size = 24) + theme(
        legend.position = 'none',
        legend.background = element_rect(),
        plot.title = element_text(angle = 0, size = 20, face = 'bold', vjust = 1),
        plot.subtitle = element_text(angle = 0, size = 18, face = 'plain', vjust = 1),
        plot.caption = element_text(angle = 0, size = 16, face = 'plain', hjust = 0, vjust = 1),
        axis.text.x = element_text(angle = 0, size = 16, face = 'bold', hjust = 0.5),
        axis.text.y = element_text(angle = 0, size = 16, face = 'bold', vjust = 0.5),
        axis.title = element_text(size = 18, face = 'bold'),
        legend.key = element_blank(),
        legend.key.size = unit(1, 'cm'),
        legend.text = element_text(size = 16),
        title = element_text(size = 16))

      p1 <- ggplot(ggdata, aes(x = Var1, y = Freq, label = Freq)) +
        geom_bar(stat = 'identity', fill = 'royalblue') +
        xlab('Haplotype ID') +
        ylab('Allele counts') +
        ylim(0, 20) +
        labs(
          title = 'HTT Allele Haplotype Distribution',
          subtitle = NULL,
          caption = NULL) +
        #scale_y_continuous(breaks = round(seq(0, 20, by = 5), 0)) +
        mytheme +
        guides(colour = guide_legend(override.aes = list(size = 2.5))) +
         #coord_flip() +
        geom_label_repel(
          data = subset(ggdata, Freq > 0),
          min.segment.length = 0,
          box.padding = 0.5,
          size = 6,
          nudge_x = 0.3,
          nudge_y = 0.9)



  # mHTT: exact match only
    df <- read.table('main/haplotyping/PhaseHaplotypesv2.csv',
      header = TRUE, sep = ',', stringsAsFactors = FALSE)

    # remove CH01855_pilot so as to not have duplicate samples
      df <- df[-which(df$Sample == 'CH01855_pilot'),]

    # mutant only
      df <- subset(df, Allele == 'mHTT')

    # summarise haps into bins
      haps <- as.numeric(sub('hap\\.', '', df$Exact.Haplotype))
      haps[is.na(haps)] <- 9999
      ggdata <- data.frame(table(
        ifelse(haps <= 6, '1-6',
          ifelse(haps <= 12, '7-12',
            ifelse(haps <= 18, '13-18',
              ifelse(haps <= 24, '19-24',
                ifelse(haps <= 30, '25-30',
                  ifelse(haps <= 36, '31-36',
                    ifelse(haps <= 42, '37-42',
                      ifelse(haps <= 48, '43-48',
                        ifelse(haps <= 54, '49-54',
                          ifelse(haps <= 60, '55-60',
                            ifelse(haps <= 66, '61-66',
                              ifelse(haps <= 72, '67-72',
                                ifelse(haps <= 78, '73-78',
                                  ifelse(haps <= 84, '79-84',
                                    ifelse(haps <= 90, '85-90',
                                      ifelse(haps <= 96, '91-96',
                                        ifelse(haps <= 102, '97-102',
                                          ifelse(haps <= 108, '103-108',
                                            NA))))))))))))))))))))
      levels <- c('1-6','7-12','13-18','19-24','25-30','31-36','37-42',
        '43-48','49-54','55-60','61-66','67-72','73-78','79-84',
        '85-90','91-96','97-102','103-108')
      ggdata <- ggdata[match(levels, ggdata$Var1),]
      ggdata$Var1 <- levels
      ggdata$Freq[is.na(ggdata$Freq)] <- 0
      ggdata$Var1 <- factor(ggdata$Var1,
        levels = c('1-6','7-12','13-18','19-24','25-30','31-36','37-42',
          '43-48','49-54','55-60','61-66','67-72','73-78','79-84',
          '85-90','91-96','97-102','103-108'))

    # plot
      p2 <- ggplot(ggdata, aes(x = Var1, y = Freq, label = Freq)) +
        geom_bar(stat = 'identity', fill = 'orange') +
        xlab('Haplotype ID') +
        ylab('Allele counts') +
        ylim(0, 20) +
        labs(
          title = 'mHTT Allele Haplotype Distribution',
          subtitle = NULL,
          caption = 'Only exact haplotype matches considered') +
        #scale_y_continuous(breaks = round(seq(0, 20, by = 5), 0)) +
        mytheme +
        guides(colour = guide_legend(override.aes = list(size = 2.5))) +
        #coord_flip() +
        geom_label_repel(
          data = subset(ggdata, Freq > 0),
          min.segment.length = 0,
          box.padding = 0.5,
          size = 6,
          nudge_x = 0.3,
          nudge_y = 0.9)



  # export
    tiff('main/haplotyping/PhaseHaplotypes.Freqs.Exact.dpi300.tiff', units = 'in',
      res = 300, width = 16, height = 9)
      cowplot::plot_grid(p1,p2,nrow=2)
    dev.off()
    tiff('main/haplotyping/PhaseHaplotypes.Freqs.Exact.dpi100.tiff', units = 'in',
      res = 100, width = 16, height = 9)
      cowplot::plot_grid(p1,p2,nrow=2)
    dev.off()



  # wtHTT: exact + nearest matches
    df <- read.table('main/haplotyping/PhaseHaplotypesv2.csv',
      header = TRUE, sep = ',', stringsAsFactors = FALSE)

    # remove CH01855_pilot so as to not have duplicate samples
      df <- df[-which(df$Sample == 'CH01855_pilot'),]

    # wild-type only
      df <- subset(df, Allele == 'HTT')

    # if no exact match, use nearest
      df$Exact.Haplotype <- ifelse(is.na(df$Exact.Haplotype),
        df$Nearest.Match.Haplotype, df$Exact.Haplotype)

    # summarise haps into bins
      haps <- as.numeric(sub('hap\\.', '', df$Exact.Haplotype))
      haps[is.na(haps)] <- 9999
      ggdata <- data.frame(table(
        ifelse(haps <= 6, '1-6',
          ifelse(haps <= 12, '7-12',
            ifelse(haps <= 18, '13-18',
              ifelse(haps <= 24, '19-24',
                ifelse(haps <= 30, '25-30',
                  ifelse(haps <= 36, '31-36',
                    ifelse(haps <= 42, '37-42',
                      ifelse(haps <= 48, '43-48',
                        ifelse(haps <= 54, '49-54',
                          ifelse(haps <= 60, '55-60',
                            ifelse(haps <= 66, '61-66',
                              ifelse(haps <= 72, '67-72',
                                ifelse(haps <= 78, '73-78',
                                  ifelse(haps <= 84, '79-84',
                                    ifelse(haps <= 90, '85-90',
                                      ifelse(haps <= 96, '91-96',
                                        ifelse(haps <= 102, '97-102',
                                          ifelse(haps <= 108, '103-108',
                                            NA))))))))))))))))))))
      levels <- c('1-6','7-12','13-18','19-24','25-30','31-36','37-42',
        '43-48','49-54','55-60','61-66','67-72','73-78','79-84',
        '85-90','91-96','97-102','103-108')
      ggdata <- ggdata[match(levels, ggdata$Var1),]
      ggdata$Var1 <- levels
      ggdata$Freq[is.na(ggdata$Freq)] <- 0
      ggdata$Var1 <- factor(ggdata$Var1,
        levels = c('1-6','7-12','13-18','19-24','25-30','31-36','37-42',
          '43-48','49-54','55-60','61-66','67-72','73-78','79-84',
          '85-90','91-96','97-102','103-108'))

    # plot
      require(ggplot2)
      require(ggrepel)
      mytheme <- theme_bw(base_size = 24) + theme(
        legend.position = 'none',
        legend.background = element_rect(),
        plot.title = element_text(angle = 0, size = 20, face = 'bold', vjust = 1),
        plot.subtitle = element_text(angle = 0, size = 18, face = 'plain', vjust = 1),
        plot.caption = element_text(angle = 0, size = 16, face = 'plain', hjust = 0, vjust = 1),
        axis.text.x = element_text(angle = 0, size = 16, face = 'bold', hjust = 0.5),
        axis.text.y = element_text(angle = 0, size = 16, face = 'bold', vjust = 0.5),
        axis.title = element_text(size = 18, face = 'bold'),
        legend.key = element_blank(),
        legend.key.size = unit(1, 'cm'),
        legend.text = element_text(size = 16),
        title = element_text(size = 16))

      p1 <- ggplot(ggdata, aes(x = Var1, y = Freq, label = Freq)) +
        geom_bar(stat = 'identity', fill = 'royalblue') +
        xlab('Haplotype ID') +
        ylab('Allele counts') +
        ylim(0, 35) +
        labs(
          title = 'HTT Allele Haplotype Distribution',
          subtitle = NULL,
          caption = NULL) +
        #scale_y_continuous(breaks = round(seq(0, 20, by = 5), 0)) +
        mytheme +
        guides(colour = guide_legend(override.aes = list(size = 2.5))) +
         #coord_flip() +
        geom_label_repel(
          data = subset(ggdata, Freq > 0),
          min.segment.length = 0,
          box.padding = 0.5,
          size = 6,
          nudge_x = 0.3,
          nudge_y = 0.9)



  # mHTT: exact + nearest matches
    df <- read.table('main/haplotyping/PhaseHaplotypesv2.csv',
      header = TRUE, sep = ',', stringsAsFactors = FALSE)

    # remove CH01855_pilot so as to not have duplicate samples
      df <- df[-which(df$Sample == 'CH01855_pilot'),]

    # mutant only
      df <- subset(df, Allele == 'mHTT')

    # if no exact match, use nearest
      df$Exact.Haplotype <- ifelse(is.na(df$Exact.Haplotype),
        df$Nearest.Match.Haplotype, df$Exact.Haplotype)

    # summarise haps into bins
      haps <- as.numeric(sub('hap\\.', '', df$Exact.Haplotype))
      haps[is.na(haps)] <- 9999
      ggdata <- data.frame(table(
        ifelse(haps <= 6, '1-6',
          ifelse(haps <= 12, '7-12',
            ifelse(haps <= 18, '13-18',
              ifelse(haps <= 24, '19-24',
                ifelse(haps <= 30, '25-30',
                  ifelse(haps <= 36, '31-36',
                    ifelse(haps <= 42, '37-42',
                      ifelse(haps <= 48, '43-48',
                        ifelse(haps <= 54, '49-54',
                          ifelse(haps <= 60, '55-60',
                            ifelse(haps <= 66, '61-66',
                              ifelse(haps <= 72, '67-72',
                                ifelse(haps <= 78, '73-78',
                                  ifelse(haps <= 84, '79-84',
                                    ifelse(haps <= 90, '85-90',
                                      ifelse(haps <= 96, '91-96',
                                        ifelse(haps <= 102, '97-102',
                                          ifelse(haps <= 108, '103-108',
                                            NA))))))))))))))))))))
      levels <- c('1-6','7-12','13-18','19-24','25-30','31-36','37-42',
        '43-48','49-54','55-60','61-66','67-72','73-78','79-84',
        '85-90','91-96','97-102','103-108')
      ggdata <- ggdata[match(levels, ggdata$Var1),]
      ggdata$Var1 <- levels
      ggdata$Freq[is.na(ggdata$Freq)] <- 0
      ggdata$Var1 <- factor(ggdata$Var1,
        levels = c('1-6','7-12','13-18','19-24','25-30','31-36','37-42',
          '43-48','49-54','55-60','61-66','67-72','73-78','79-84',
          '85-90','91-96','97-102','103-108'))

    # plot
      p2 <- ggplot(ggdata, aes(x = Var1, y = Freq, label = Freq)) +
        geom_bar(stat = 'identity', fill = 'orange') +
        xlab('Haplotype ID') +
        ylab('Allele counts') +
        ylim(0, 35) +
        labs(
          title = 'mHTT Allele Haplotype Distribution',
          subtitle = NULL,
          caption = 'Exact & nearest haplotype matches considered') +
        #scale_y_continuous(breaks = round(seq(0, 20, by = 5), 0)) +
        mytheme +
        guides(colour = guide_legend(override.aes = list(size = 2.5))) +
        #coord_flip() +
        geom_label_repel(
          data = subset(ggdata, Freq > 0),
          min.segment.length = 0,
          box.padding = 0.5,
          size = 6,
          nudge_x = 0.3,
          nudge_y = 0.9)



  # export
    tiff('main/haplotyping/PhaseHaplotypes.Freqs.AnyMatch.dpi300.tiff', units = 'in',
      res = 300, width = 16, height = 9)
      cowplot::plot_grid(p1,p2,nrow=2)
    dev.off()
    tiff('main/haplotyping/PhaseHaplotypes.Freqs.AnyMatch.dpi100.tiff', units = 'in',
      res = 100, width = 16, height = 9)
      cowplot::plot_grid(p1,p2,nrow=2)
    dev.off()

