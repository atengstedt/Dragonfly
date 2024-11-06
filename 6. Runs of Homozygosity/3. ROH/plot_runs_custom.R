#custom version of the plot_Runs function from the detectRUNS R package: https://github.com/bioinformatics-ptp/detectRUNS

plot_Runs_custom <- function(runs, suppressInds=FALSE, savePlots=FALSE, separatePlots=FALSE, outputName=NULL) {

  # avoid notes
  chrom <- NULL ; from <- NULL ; to <- NULL ; group <- NULL

  chr_order <- c((0:10000),"X","Y","XY","MT","Z","W")
  list_chr=unique(runs$chrom)
  new_list_chr=as.vector(sort(factor(list_chr,levels=chr_order, ordered=TRUE)))

  plot_list <- list()
  for (chromosome in new_list_chr) {

    #subset by chromosome
    krom <- subset(runs,chrom==chromosome)

    #rearrange subset
    teilsatz <- krom[,c(5,6,2,1)]
    teilsatz <- teilsatz[order(teilsatz$group),]

    #new progressive numerical ID
    newID <- seq(1,length(unique(teilsatz$id)))
    id <- unique(teilsatz$id)
    teilsatz$NEWID=newID[match(teilsatz$id,id)]

    optionen <- ggplot2::scale_y_discrete("IDs",limits=unique(teilsatz$id))
    alfa <- 1
    grosse <- 1

    if (length(id) > 50) {

      optionen <- ggplot2::theme(axis.text.y=element_blank(), axis.title.y=element_blank(),axis.ticks.y=element_blank())
      alfa <- 0.75
      grosse <- 1
    }

    if (suppressInds) optionen <- ggplot2::theme(axis.text.y=element_blank(),
                                                 axis.title.y=element_blank(),axis.ticks.y=element_blank())

    #size in mb
    teilsatz$from <- (teilsatz$from/(10^6))
    teilsatz$to <- (teilsatz$to/(10^6))

    row.names(teilsatz) <- NULL

    teilsatz$id <- as.factor(teilsatz$id)
    teilsatz$id <- factor(teilsatz$id, levels = unique(teilsatz$id[order(teilsatz$NEWID)]))

    p <- ggplot2::ggplot(teilsatz)
    p <- p + ggplot2::geom_segment(data=teilsatz,aes(x = from, y = id, xend = to,
                                                     yend = id,colour=as.factor(group)),alpha=alfa, size=grosse)
    p <- p + ggplot2::scale_color_manual(values=c("NWJ"="#9b5fe0","EJ"="#16a4d8","WJ"="#60dbe8","SJ"="#8bd346","NEZ"="#efdf48","SZ"="#f9a52c","B"="#d64e12"))
    p <- p + ggplot2::xlim(0, max(teilsatz$to)) + ggplot2::ggtitle(paste('Scaffold ',chromosome))
    p <- p + ggplot2::guides(colour=guide_legend(title="Population")) + ggplot2::xlab("Mbps")
    p <- p + theme(plot.title = element_text(hjust = 0.5))
    p <- p + optionen


    if(savePlots & separatePlots) {
      if (! is.null(outputName)) {
        fileNameOutput <- paste(outputName, "Chr", chromosome, '.pdf', sep="_")
      } else { fileNameOutput <- paste("Chr", chromosome, '.pdf',sep="_") }
      ggsave(filename = fileNameOutput , plot = p, device = "pdf")
    } else if (savePlots) { plot_list[[chromosome]] <- p
    } else { print(p)}

  }

  if(savePlots & !separatePlots) {

    if (! is.null(outputName)) {
      fileNameOutput <- paste(outputName, "AllChromosomes.pdf", sep="_")
    } else {
      fileNameOutput <- "Runs_AllChromosome.pdf"
    }
      plot_list_final <- gridExtra::marrangeGrob(plot_list, nrow=1, ncol=1, top=NULL)
      ggsave(filename = fileNameOutput , plot = plot_list_final, device = "pdf", units="cm", height=21, width=15.92, dpi=300)
  }
}
