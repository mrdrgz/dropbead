MixedSpeciesSample <- setClass(Class = "MixedSpeciesLabelledSample",
                         slots = c(T2C = "data.frame",
								   T2CExpanded = "data.frame"),
                         contains = "MixedSpeciesSample"
                         )

setMethod("initialize",
          "MixedSpeciesLabelledSample",
          function (.Object, species1="", species2="", cells=c(), genes=c(), dge=data.frame()) {
            .Object@species1 = species1
            .Object@species2 = species2
            .Object@cells = names(dge)
            .Object@genes = rownames(dge)
            .Object@dge = dge
            .Object@T2C = NA
            .Object@T2CExpanded = NA
            .Object
          })

#' Add T2C information
#'
#' @param object A \code{MixedSpeciesSample} object.
#' @return A list of \code{data.frames} corresponding to the genes of the two species.
setGeneric(name = "addT2CStats",
           def = function(object, data.file = NA) {standardGeneric("addT2CStats")})
setMethod(f = "addT2CStats",
          signature = "MixedSpeciesSample",
          function(object, data.file) {
            library(data.table)
            T2C.table = data.frame(fread(paste('zcat <', data.file, "| grep INTERGENIC", sep=" "), header=F, stringsAsFactors=F,colClasses=c("character","numeric","character","character", "character","character","character","character","character","character"), col.names=c("Chr","Position","Transition","Alignment","Read","Nucleotide","Cell_barcode","Molecular_barcode","Gene",
"Region", "NH")))
	    object@T2CExpanded = T2C.table
	    return(object)            

          })
