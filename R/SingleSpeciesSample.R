SingleSpeciesSample <- setClass(Class = "SingleSpeciesSample",
                                slots = c(species1 = "character",
                                          cells = "vector",
                                          genes = "vector",
                                          dge = "data.frame"))

setMethod("initialize",
          "SingleSpeciesSample",
          function (.Object, species1 = "", cells = c(), genes = c(), dge = data.frame()) {
            .Object@species1 = species1
            .Object@cells = names(dge)
            .Object@genes = rownames(dge)
            .Object@dge = dge
            .Object
          })

#' Compute genes per cell
#'
#' @param object A Single species sample.
#' @return A \code{data.frame} with cells, gene counts and species.
setGeneric("computeGenesPerCell",
           function(object, min.reads=2, ...) {standardGeneric("computeGenesPerCell")})
setMethod("computeGenesPerCell",
          "SingleSpeciesSample",
          function(object, min.reads) {
            genes <- data.frame("cells" = names(colSums(object@dge >= min.reads)),
                                "counts" = as.numeric(colSums(object@dge >= min.reads)),
                                "species" = object@species1)
            return (genes)
          })

#' Compute transcripts per cell
#'
#' @param object A Single species sample.
#' @return A \code{data.frame} with cells, transcript counts and species.
setGeneric("computeTranscriptsPerCell",
           function(object, ...) {standardGeneric("computeTranscriptsPerCell")})
setMethod("computeTranscriptsPerCell",
          "SingleSpeciesSample",
          function(object) {
            transcripts <- data.frame("cells" = object@cells,
                                      "counts" = as.numeric(colSums(object@dge)),
                                      "species" = object@species1)
            return (transcripts)
          })

setMethod("removeLowQualityCells",
          "SingleSpeciesSample",
          function(object, min.genes) {
            return (new("SingleSpeciesSample", species1=object@species1, dge=removeLowQualityCells(object@dge, min.genes)))
          })

setMethod("keepBestCells",
          "SingleSpeciesSample",
          function(object, num.cells, min.num.trans) {
            return (new("SingleSpeciesSample", species1=object@species1,
                        dge=keepBestCells(object@dge, num.cells, min.num.trans)))
          })

setMethod("removeLowQualityGenes",
          "SingleSpeciesSample",
          function(object, min.cells) {
            return (new("SingleSpeciesSample", species1=object@species1, dge=removeLowQualityGenes(object@dge, min.cells)))
          })

setMethod("geneExpressionMean",
          "SingleSpeciesSample",
          function(object) {
            return (geneExpressionMean(object@dge))
          })

setMethod("geneExpressionDispersion",
          "SingleSpeciesSample",
          function(object) {
            return (geneExpressionDispersion(object@dge))
          })

setMethod("geneExpressionVariability",
          "SingleSpeciesSample",
          function(object, bins, low, do.plot) {
            return (geneExpressionVariability(object@dge, bins, low, do.plot))
          })

#' List cells that are candidates for collapsing.
#'
#' Identify and list cells which share 11 bases in their barcodes and only the last
#' one is different. The cells are marked as candidates if and only if they're classified
#' as belonging to the same species.
#' @param A \code{SingleSpeciesSample} object.
#' @return A list of pairs od candidate cells marked for collapsing.
setGeneric("listCellsToCollapse",
           function(object, ...) {
             standardGeneric("listCellsToCollapse")
           })
setMethod("listCellsToCollapse",
          "SingleSpeciesSample",
          function (object) {
            theListOfCellPairs <- list()
            cells <- sort(object@cells)
            for (cell in 1:(length(cells)-1)) {
              if (substr(cells[cell], 1, 11) == substr(cells[cell+1], 1, 11)) {
                if (substr(cells[cell], nchar(cells[cell]), nchar(cells[cell])) == "N" |
                    substr(cells[cell+1], nchar(cells[cell+1]), nchar(cells[cell+1])) == "N") {
                  theListOfCellPairs <- c(theListOfCellPairs, list(c(cells[cell], cells[cell+1])))
                }
              }
            }
            return(theListOfCellPairs)
          })

#' Collapse cells by barcodes similarity
#'
#' Collapse cells which have barcodes differing by one mutation on the last base.
setGeneric("collapseCellsByBarcode",
           function(object, ...) {
             standardGeneric("collapseCellsByBarcode")
           })
setMethod("collapseCellsByBarcode",
          "SingleSpeciesSample",
          function(object) {
            listOfCells <- listCellsToCollapse(object)
            if (length(listOfCells) == 0) {
              return(object)
            }

            for (index in 1:length(listOfCells)) {
              object@dge <- cbind(object@dge, rowSums(object@dge[, listOfCells[[index]]]))
            }
            object@dge <- object@dge[, !names(object@dge) %in% unlist(listOfCells)]

            names(object@dge)[(length(names(object@dge)) -
                                     length(listOfCells) + 1):length(names(object@dge))] <- unlist(listOfCells)[seq(1, length(unlist(listOfCells)), 2)]
            object@cells <- names(object@dge)
            return (object)
          })

#' Compare gene expression levels between two single species samples
#'
#' @param object1 A \code{SingleSpeciesSample} object.
#' @param object2 A \code{SingleSpeciesSample} object.
#' @param name1 The name of the first sample.
#' @param name2 The name of the second sample.
#' @return A plot comparing the gene expression levels.
setGeneric("compareGeneExpressionLevels",
           function(object1, object2, name1="sample1", name2="sample2",
                    col="steelblue", method="pearson", ...) {
             standardGeneric("compareGeneExpressionLevels")
           })
setMethod("compareGeneExpressionLevels",
          "SingleSpeciesSample",
          function(object1, object2, name1, name2, col, method) {
            common.genes <- intersect(object1@genes, object2@genes)

            object1 <- object1@dge[common.genes, ]
            object2 <- object2@dge[common.genes, ]

            big.df <- data.frame("genes" = common.genes,
                                 "sample1" = log(as.numeric(rowSums(object1))+1, 2),
                                 "sample2" = log(as.numeric(rowSums(object2))+1, 2))

            cor <- signif(cor(big.df$sample1, big.df$sample2, method=method), 2)

            comp.plot <- (ggplot(data = big.df, aes(x = sample1, y = sample2))
                          + ggtitle(paste0("R= ", cor))
                          + xlab(name1) + ylab(name2)
                          + geom_point(col=col, alpha=0.5, size=2)
                          + scale_y_continuous(expand = c(0.02, 0.02)) + scale_x_continuous(expand = c(0.02, 0.02))
                          + theme_minimal() + plotCommonGrid + plotCommonTheme
                          + theme(axis.ticks.y = element_blank()))

            return (comp.plot)
          })

setMethod("computeCorrelationSingleCellsVersusBulk",
          "SingleSpeciesSample",
          function(single.cells, bulk.data, measure, method) {
            corr = computeCorrelationSingleCellsVersusBulk(single.cells@dge, bulk.data, measure, method)[[1]]
            df = computeCorrelationSingleCellsVersusBulk(single.cells@dge, bulk.data, measure, method)[[2]]
            return (list(corr, df))
          })

setMethod("computeCorrelationCellToCellVersusBulk",
          "SingleSpeciesSample",
          function(single.cells, bulk.data, measure, method) {
            return (computeCorrelationCellToCellVersusBulk(single.cells@dge, bulk.data, measure, method))
          })

setMethod("compareSingleCellsAgainstBulk",
          "SingleSpeciesSample",
          function(single.cells, bulk.data, measure, method, ylab, col) {
            return (compareSingleCellsAgainstBulk(single.cells@dge, bulk.data,
                                                  measure, method, ylab, col))
            })

setMethod("computeCellGeneFilteringFromBulk",
          "SingleSpeciesSample",
          function(single.cells, bulk.data, log.space, method,
                   min.cells, max.cells, iteration.steps) {
            computeCellGeneFilteringFromBulk(single.cells@dge, bulk.data, log.space, method,
                                             min.cells, max.cells, iteration.steps)
            })

setGeneric("assignCellCyclePhases",
           function(object, gene.file="~/Desktop/things/git/dropseq/data/cell_cycle_genes.xlsx", do.plot=T, ...) {
             standardGeneric("assignCellCyclePhases")
           })
setMethod("assignCellCyclePhases",
          "SingleSpeciesSample",
          function(object, gene.file, do.plot, cycleOrder=c("G1","G1.S","S","G2","G2.M","M")) {
            require(xlsx)
            require(plyr)
            file.ext = gsub(".*\\.(.*)$","\\1", gene.file)
            
			if (file.ext == "xlsx") {
				if (object@species1 == "human") {
					cc_genes <- read.xlsx(gene.file, sheetIndex = 2, stringsAsFactors = F)
				}

				if (object@species1 == "mouse") {
					cc_genes <- read.xlsx(gene.file, sheetIndex = 3, stringsAsFactors = F)
				}

			} else if (file.ext == "csv") {
			    cc_genes = read.table(gene.file, header=T, stringsAsFactors=F, sep=",")
				cc_genes$Phase = gsub("/", ".", cc_genes$Phase)
            }  

			cc_genes$Symbol = gsub(" ", "", cc_genes$Symbol)
			
			ScorePhasesPerCell = function(dge, genes) {
					return ( mean(log(dge[genes]+1, 2), na.rm=T) )
			}
			
			ScorePhase = function(genes, dge=object@dge) {
					return( apply(dge, 2, ScorePhasesPerCell, genes) )
			}
			
			
			# Use with/by to subset dynamically the gene symbols by phase, pass it to the 
			# wrapper function ScorePhase (takes a DGE matrix and iterates cells) 
			# and ScorePhasesPerCell (scores each cell and returns a vector to the wrapper
			phase_score = do.call(cbind,with(cc_genes,by(Symbol,Phase,ScorePhase)))
			phase_score = phase_score[,cycleOrder]
			
			# Normalize and reannotate rows
            phase_score_norm <- data.frame(
						apply(phase_score, 2, scale),row.names=rownames(phase_score)
					)
            
            # Normalize 2 and reannotate columns
            phase_score_norm2 <- setNames(
					data.frame(t(apply(phase_score_norm, 1, scale))),
					colnames(phase_score_norm)
				)
            
            
            
            # Take the maximum and annotate the phase
            phase_score_norm2$Phase = colnames(phase_score)[apply(phase_score_norm2, 1, which.max)]

			# Order columns using a reference column i
			# Then, rearrange the rest to get -x,-x around column i
			# including columns
			OrderColumns = function(i,v) {
					# v with columns
					# i == reference column
					c.i = which(v == i)
					c.o = c(c.i)
					# Number of iterations (floor)
 					n.it = floor((length(v) - 1)/2)
 					for (it in seq(1,n.it,1))
 					{	c.o = c(c.o,
								ifelse(c.i - it < 1, length(v)+(c.i-it), c.i-it),
								ifelse(c.i + it > length(v), length(v) - c.i+it, c.i+it)
							)
 					}
 					
 					# If the length is uneven, add last value to the vector
 					if ( (length(v)-1) %% 2 == 1 ) {
						it = n.it + 1
						c.o = c(c.o,
							ifelse(c.i - it < 1, length(v)+(c.i-it), c.i-it)
							)
					}
						return(c.o)
			}
				
			
			RoundAndOrder = function(i,df) {
 					# Subset by Phase
 					df=phase_score_norm2
 					df <- df[df$Phase == i,grep("Phase", colnames(df), invert=T)]
 					
 					# Order and round values for the column i
 					df <- round(df[order(df[[ i ]], decreasing = T), ], 2)
 					
 					# Create rearrenged column vector		
 					column.order = colnames(df)[OrderColumns(i,colnames(df))]
					
 					# Order df
 					df = df[do.call("order", c(df[,column.order], decreasing=T)), ]
 					df$Phase = i
					return( df )
			}
			
			
			phase_score_norm2 = do.call(rbind, lapply(colnames(phase_score_norm2)[seq(1,length(phase_score_norm2)-1,1)],
				   RoundAndOrder, df=phase_score_norm2)
			)


            heatmap_palette <- colorRampPalette(c("#3794bf", "#FFFFFF", "#cc4140"))
            if (do.plot) {
              heatmap.2(as.matrix(t(phase_score_norm2[,seq(1, ncol(phase_score_norm2)-1,1)])), 	
						trace='none', Rowv=F, Colv=F, dendrogram='none', 
						labCol=F, col=heatmap_palette(10))
            }

#             return()
			  return(list(Scaled=phase_score_norm2,
						  Raw=phase_score))
          })
