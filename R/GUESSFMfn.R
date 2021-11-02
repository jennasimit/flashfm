## These are extracted from the GUESSFM package by Chris Wallace
## generics from S3
setGeneric("length")
setGeneric("summary")


#' Create a union of groups, snppicker or tags objects
#'
#' First, objects are converted to class groups.  Then any groups
#' which overlap are identified, and replaced by their union.  Groups
#' which do not overlap are unchanged.  The combined set of groups is
#' returned.
#' @rdname union
#' @author Chris Wallace
#' @param x object of class \code{groups}, \code{snppicker} or \code{tags}
#' @param y object of same class as x
#' @return object of class groups
setGeneric("union", function(x,y) standardGeneric("union"))

#' overlap
#' 
#' overlap can be used to examine the overlap between groups before the union is created
#' @author Chris Wallace
#' @param x object of class \code{groups}, \code{snppicker} or \code{tags}
#' @param y object of same class as x
#' @return object of class groups
setGeneric("overlap", function(x,y) standardGeneric("overlap"))

#' Check whether a snp is in a snppicker, groups or tags object
#'
#' @rdname snpin
#' @author Chris Wallace
#' @param x character vector of SNPs to be checked
#' @param y object of class snppicker, groups or tags
#' @return logical matrix of nrow equal to length(x) and ncol equal to number of groups in y
setGeneric("snpin", function(x,y) standardGeneric("snpin"))
#setGeneric("snpdrop", function(x,y) standardGeneric("drop"))

##' Convert from old definitions of groups, tags classes to new
##' 
##' DON'T USE THIS FUNCTION UNLESS YOU HAVE OBJECTS STORED FROM A
##' PREVIOUS PRE-DEVELOPMENT VERSION!
##' 
#' @rdname conversion
##' @param object GUESSFM object from pre-development version
##' @return new S4 structure
setGeneric("convert",function(object) standardGeneric("convert"))

##' Accessor functions
##'
##' \code{snps} shows list of snps in an object of class groups and
##' returns a list of character vectors
##'
##' @title Accessors for groups objects
##' @param object object from which items should be extracted
##' @return a list of character vectors
##' @author Chris Wallace
##' @rdname accessors
setGeneric("snps",function(object) standardGeneric("snps"))

##' @details \code{tags} shows tags from an object of class groups and
##' returns a character vector of tag SNPs
##' @rdname accessors
setGeneric("tags",function(object) standardGeneric("tags"))

##' \code{tagsof} shows tags for a named character vector of SNPs
##' @rdname groups-subset
##' @param object tags or groups object
##' @return data.frame of tags and their tagged SNPs
setGeneric("tagsof", function(object,i) standardGeneric("tagsof"))

##' \code{taggedby} shows SNPs tagged by a named character vector of tag SNPs
##' @rdname groups-subset
setGeneric("taggedby", function(object,i) standardGeneric("taggedby"))


#' Class to hold data relating to multiple models fitted to SNP data
#'
#' @slot snps data.frame containing marginal probabilities of inclusion for the SNPs
#' @slot models data.frame containing summaries for each model
#' @slot model.snps list containing the SNPs for each model. May be removed.
#'
#' new("snpmod")
setClass("snpmod",
         slots=c(snps="data.frame",
                        models="data.frame",
                        model.snps="list"),
         validity=function(object) {
           if(nrow(object@models)!=length(object@model.snps))
             stop("Model summary should contain same number of models as model.snps decodes")
         })
#' Group focused class for holding information about sets of SNPs
#' defined by their mutual LD
#'
#' \code{groups} and \code{tags} are two structures for holding the
#' same information, depending on whether your focus is on the sets of
#' SNPs or their index members.  It is easy to convert one to another
#' and perhaps, in future, one class may be deprecated.
#'
#' @slot tags character vector giving tag SNPs.  Each tag indexes one group of SNPs
#' @slot .Data list of character vectors giving the SNP membership of each group
#' 
#' new("groups")
setClass("groups",
         slots=c(tags="character"),
         contains="list",
         validity=function(object) {
           if(length(object@tags)!=length(object@.Data)) {
             stop("groups must be named by their tag")
           }
         })
#' Tags focused class for holding information about sets of SNPs
#' defined by their mutual LD
#' @slot tags character vector giving tag SNPs, one per SNP in
#' in \code{.Data}, repeated as necessary
#' @slot .Data character vector giving SNPs included in this tags object
#' new("tags")
setClass("tags",
         slots=c(tags="character"),
         contains="character",
         validity=function(object) {
           if(length(object@tags)!=length(object@.Data))
             stop("tags must be in tags vector, tagging themselves")
         })
         
#' Class to hold results of snp.picker algorithm
#'
#' \code{snp.picker} groups SNPs according to LD and model inclusion
#' and outputs objects of class snppicker.
#'
#' @slot groups list of data.frames describing SNPs in each group
#' @slot plotsdata list of additional data relating to the snp.picker
#' process that allows a summary of that process to be plotted via
#' \code{plot}.
#' new("snppicker")
setClass("snppicker",
         slots=c(plotsdata="list",groups="list"),
         validity=function(object) {
           if(length(object@groups)!=length(object@plotsdata))
             stop("groups and plotsdata should be lists of equal length")
         })

#' Class to hold results of pp.nsnp
#'
#' \code{pp.nsnp} summarises prior and posterior support for the
#' number of SNPs required to model a trait or several traits
#'
#' @slot .Data contains a named list of numeric vectors giving the support for each number
#' @slot plot contains a ggplot object
#' @slot traits names of traits, corresponding to items in .Data
#'
#'
#' new("ppnsnp")
setClass("ppnsnp",
         slots=c(.Data="list",plot="ANY",traits="character"),
         validity=function(object) {
           if(length(object@.Data)!=length(object@traits))
             stop("traits should be same length as .Data")
           for(i in seq_along(object@.Data)) {
             if(!is(object[[i]],"array"))
               stop("ppnsnp should contain a list of arrays")
           }})

######################################################################

##                          display                                 ##

######################################################################
##' Show
##'
##' Methods to briefly summarise objects to screen. \code{show} is
##' called implicitly when you type an object's name at the command
##' line.  Use \code{summary}, where available, to get more details.
##' @param object the thing to show
##' @export
##' @rdname show-methods 
setMethod("show", signature="snpmod",
          function(object) {
            nmod <- nrow(object@models)
            nsnp <- nrow(object@snps)
            maxpp <- max(object@models$PP)
            minpp <- min(object@models$PP)
            spp <- sum(object@models$PP)
            message("snpmod object, containing information on ",nmod," models / ",nsnp," SNPs.")
            message(sprintf("PP ranges from %4.3f-%4.3f (sum: %4.3f).",minpp,maxpp,spp))
          })

##' @rdname show-methods
setMethod("show", signature="snppicker",
          function(object) {
            ngroup <- length(object@groups)
            nsnps <- if(ngroup==0) { 0 } else { sapply(object@groups,nrow) }
            message("snppicker object, containing ",sum(nsnps)," SNPs grouped into ",ngroup," groups.")
          })
##' @rdname show-methods
setMethod("show", signature="ppnsnp",
          function(object) {
            L <- object@.Data
            names(L) <- object@traits
            show(L)
          })

##' @rdname show-methods
setMethod("show", signature="tags",
          function(object) {
            ntags <- length(unique(object@tags))
            nsnps <- length(object@.Data)
            message("tags object, containing ",nsnps," SNPs in ",ntags," groups.")
          })

##' @rdname show-methods
setMethod("show", signature="groups",
          function(object) {
            ntags <- length(object@tags)
            nsnps <- length(unlist(object@.Data))
            message("groups object, containing ",nsnps," SNPs in ",ntags," groups.")
          })

######################################################################

##          detailed summary                                        ##

######################################################################


##' Summaries
##' 
##' Print summary of an object
##'
##' @param object the thing to summarise
##' @rdname summary
setMethod("summary",signature="snppicker",
          function(object){
            ngroups <- length(object@groups)
            cmpi <- unlist(lapply(object@groups, function(x) max(x$cmpi)))
            nsnp <- unlist(lapply(object@groups, nrow))
            index <- unlist(lapply(object@groups, function(x) x[1,"var"]))
            maxr2 <- unlist(lapply(object@groups, function(x) x[nrow(x),"r2"]))
            data.frame(SNP.index=index,
                       SNP.count=nsnp,
                       min.R2=1-maxr2,
                       gMMPI=cmpi)
          })
##' @rdname summary
setMethod("summary",signature="groups",
          function(object) {
            data.frame(tag=object@tags,
                       nsnp=sapply(object@.Data,length))
          })


################################################################################

## convert between snppicker, groups and tags

################################################################################

setAs("groups", "tags",      
      def=function(from) {
        if(!length(from))
          return(new("tags")) # return empty object
        new("tags",unlist(from@.Data),tags=rep(from@tags,times=sapply(from@.Data,length)))
      })
setAs("tags", "groups",
      def=function(from) {
        if(!length(from@tags))
          return(new("groups"))
        new("groups",split(from@.Data,from@tags),tags=unique(from@tags))
      })
setAs("snppicker","groups",
      def=function(from) {
        if(!length(from@groups))
          return(new("groups")) # return empty object
        new("groups",
            lapply(from@groups, "[[", "var"),
            tags=unlist(lapply(from@groups,function(x) x[1,"var"])))
      })
setAs("snppicker","tags",
      def=function(from) {
        if(!length(from@groups))
          return(new("tags")) # return empty object
        g <- new("groups",
                 lapply(from@groups, "[[", "var"),
                 tags=unlist(lapply(from@groups,function(x) x[1,"var"])))
        as(g,"tags")
      })

################################################################################

## subsets

################################################################################

##' @rdname groups-subset
setMethod("taggedby",signature=c(object="tags",i="character"),
          function(object,i) {
            wh <- which(object@tags %in% i)
            if(!length(wh))
              stop("tags not found")
            ret <- data.frame(tag=object@tags[wh],snps=object@.Data[wh])
            return(ret[order(ret$tag),])  
          })
##' @rdname groups-subset
setMethod("tagsof",signature=c(object="tags",i="character"),
          function(object,i) {
            wh <- which(object@.Data %in% i)
             if(!length(wh))
              stop("SNPs not found")
           ret <- data.frame(tag=object@tags[wh],snps=object@.Data[wh])
            return(ret[order(ret$tag),])  
          })
##' @rdname groups-subset
setMethod("taggedby",signature=c(object="groups",i="character"),
          function(object,i) {
            taggedby(as(object,"tags"),i)
          })
##' @rdname groups-subset
setMethod("tagsof",signature=c(object="groups",i="character"),
          function(object,i) {
            tagsof(as(object,"tags"),i)
          })

##' Subset groups or tags objects
##'
##' '[' will extract another object of the same class.  '[[' will extract a single element.
##' @param x groups or tags object
##' @param i numeric, logical or character vector to index SNPs or tags
##' @return subsetted groups or tags object
##' @rdname groups-subset
setMethod("[",signature=c(x="groups",i="character",j="missing",drop="missing"),
          function(x,i) {
            wh <- which(x@tags %in% i)
            new("groups",x@.Data[wh],tags=x@tags[wh])
          })
##' @rdname groups-subset
setMethod("[",signature=c(x="groups",i="numeric",j="missing",drop="missing"),
          function(x,i) {
            new("groups",x@.Data[i],tags=x@tags[i])
          })
##' @rdname groups-subset
setMethod("[",signature=c(x="groups",i="logical",j="missing",drop="missing"),
          function(x,i) {
            new("groups",x@.Data[i],tags=x@tags[i])
          })
##' @rdname groups-subset
setMethod("[",signature=c(x="tags",i="character",j="missing",drop="missing"),
          function(x,i) {
            wh <- sapply(i,function(ii) which(snps(x)==ii))
            tags(tags)[wh]
          })
##' @rdname groups-subset
setMethod("[[",signature=c(x="groups",i="numeric"),
          function(x,i) {
            x@.Data[[i]]
          })
##' @rdname groups-subset
setMethod("[[",signature=c(x="groups",i="logical"),
          function(x,i) {
            x@.Data[[i]]
          })
##' @rdname groups-subset
setMethod("[[",signature=c(x="groups",i="character"),
          function(x,i) {
            wh <- which(x@tags %in% i)
            x@.Data[[wh]]
          })

######################################################################

##                         conversion                               ##

######################################################################

#' @rdname conversion
setMethod("convert",signature=c(object="groups"), function(object) {
  object@.Data=object@groups
  return(object)
})
#' @rdname conversion
setMethod("convert",signature=c(object="tags"), function(object) {
  object@.Data=object@snps
  return(object)
})



#' @rdname snpin
setMethod("snpin",signature(x="character",y="snppicker"),definition=function(x,y) {
  snpin(x,as(y,"groups"))
})

#' @rdname snpin
setMethod("snpin",signature(x="character",y="tags"),definition=function(x,y) {
  snpin(x,as(y,"groups"))
})

#' @rdname snpin
setMethod("snpin",signature(x="character",y="groups"),definition=function(x,y) {
  if(!length(y@.Data) || !length(x))
    return(NULL)
  names(y@.Data) <- paste0("group",seq_along(y@.Data))
  ret <- sapply(y@.Data,function(yg) x %in% yg)
  if(is.null(dim(ret)))
    ret <- matrix(ret,nrow=1)
  rownames(ret) <- x
  return(ret)
})

######################################################################

##                        concatenate                               ##

######################################################################
## two groups in x (8, 11) match to one group in y (1) with dups
## x.7 also matches to y.1

#' @rdname union
setMethod("union",signature(x="snppicker",y="snppicker"),definition=function(x,y) {

  ##   xgr <- as(x,"groups")
  ## ygr <- as(y,"groups")
  ## int <- .group.intersection(xgr,ygr)
    M <- new("snppicker",
             plotsdata=c(x@plotsdata,y@plotsdata),
             groups=c(x@groups,y@groups))
    summppi <- sapply(M@groups, function(k) sum(k$Marg_Prob_Incl))    
    M <- M[order(summppi,decreasing=TRUE)]
    wh <- .group.intersection(as(M,"groups"),as(M,"groups"),as.which=TRUE)
    if(!nrow(wh))
        return(M)
    print(wh)
    Mkeep <- rep(TRUE,length(M@groups))
    for(i in 1:nrow(wh)) { # overlap exists
        ix <- wh[i,1]
        iy <- wh[i,2]
        ##    iy.M <- iy + length(x@groups)
        gx <- M@groups[[ ix ]]
        gy <- M@groups[[ iy ]]
        gint <- intersect(rownames(gx),rownames(gy))
        ## complete subset
        if(all(rownames(gy) %in% rownames(gx))) {
            Mkeep[[ iy ]] <- FALSE
            next
        }
        if(all(rownames(gx) %in% rownames(gy))) {
            Mkeep[[ ix ]] <- FALSE
            next
        }
        pintx <- sum(gx[gint,"Marg_Prob_Incl"])
        pinty <- sum(gy[gint,"Marg_Prob_Incl"])
        pnintx <- sum(gx[setdiff(rownames(gx),gint),"Marg_Prob_Incl"])
        pninty <- sum(gy[setdiff(rownames(gy),gint),"Marg_Prob_Incl"])
        if(pnintx > pintx || pninty > pinty) { # unmerge
            if(pintx/(pintx + pnintx) > pinty/(pinty + pninty)) { # keep int in x
                M@groups[[iy]] <- M@groups[[iy]][ !(rownames(M@groups[[iy]]) %in% gint), ]
            } else { # keep int in y
                M@groups[[ix]] <- M@groups[[ix]][ !(rownames(M@groups[[ix]]) %in% gint), ]
            }
        } else { # merge
            tmp <- rbind(M@groups[[ix]],
                         M@groups[[iy]][ !(rownames(M@groups[[iy]]) %in% gint), ])
            M@groups[[ix]] <- tmp[ order(tmp$Marg_Prob_Incl, decreasing=TRUE), ]
            Mkeep[[ iy ]] <- FALSE
        }
    }

    M <- new("snppicker",
        plotsdata=M@plotsdata[which(Mkeep)],
        groups=M@groups[which(Mkeep)])
})

#' @rdname union
setMethod("union",signature(x="tags",y="tags"),definition=function(x,y) {
  ugr <- union(as(x,"groups"),as(y,"groups"))
  as(ugr,"tags") })

.group.intersection <- function(x,y,drop=numeric(0),as.which=FALSE) {
  int <- matrix(FALSE,length(x),length(y))
  for(i in seq_along(x)) {
    for(j in seq_along(y)) {      
      int[i,j] <- sum(x[[i]] %in% y[[j]])
    }
  }
  if(as.which) {
      wh <- which(int>0,arr.ind=TRUE)
      wh <- cbind(wh, int[wh])
      wh <- wh[ wh[,1] < wh[,2], ,drop=FALSE ] # remove diagonals and only count pairs once
      wh <- wh[order(wh[,3],decreasing = TRUE),,drop=FALSE]
      if(length(drop))
          wh <- wh[ !(wh[,1] %in% drop | wh[,2] %in% drop), , drop=FALSE ]
      return(wh)
  }
  return(int)
}

#' @rdname union
setMethod("union",signature(x="groups",y="groups"),definition=function(x,y) {
  if(!length(x))
    return(y)
  if(!length(y))
    return(x)
  ## find intersecting groups, and make unions
  wh <- .group.intersection(x,y,as.which=TRUE)
  if(nrow(wh)) {
    ## merge y into x
    for(i in 1:nrow(wh)) {
      xi <- wh[i,1]
      yi <- wh[i,2]
      x[[xi]] <- unique(c(x[[xi]],y[[yi]]))
    }
    ## check - duplicate y's mean x groups must be merged
    dups <- wh[,2][ duplicated(wh[,2]) ]
    if(length(dups)) {
      ## do the merge
      for(ydup in dups) {
        xi <- wh[ wh[,2]==ydup, 1]
        x[[ xi[1] ]] <- unique(unlist(x[xi]))
      }
      ## drop NULLs now we don't need wh any more
      to.drop <- vector("list",length(dups))
      for(i in seq_along(dups)) {
        yi <- dups[[i]]
        to.drop[[i]] <- (wh[ wh[,2]==ydup, 1])[-1]
      }
      to.drop <- unique(unlist(to.drop))
      x@.Data <- x@.Data[ -to.drop ]
      x@tags <- x@tags[ -to.drop ]
    }
    
    add <- setdiff(seq_along(y),unique(wh[,2]))
    if(!length(add))
      return(x)
    return(new("groups",c(x,y[add]),tags=c(x@tags,y@tags[add])))
  }
  ## no overlap: simply concatenate
  new("groups",c(x,y=y),tags=c(x@tags,y@tags))
})

######################################################################

## split

######################################################################

## # @rdname union-methods
## # @aliases union,snppicker,snppicker-method
## setMethod("union",signature(x="snppicker",y="snppicker"),definition=function(x,y) {
##   union(as(x,"groups"),as(y,"groups")) })
## # @rdname union-methods
## # @aliases union,tags,tags-method
## setMethod("union",signature(x="tags",y="tags"),definition=function(x,y) {
##   as(union(as(x,"groups"),as(y,"groups")),"tags") })
## # @rdname split-methods
## # @aliases split,groups,groups-method
## setMethod("split",signature(x="groups",y="groups"),definition=function(x,y) {


######################################################################

##                      accessors                                   ##

######################################################################

##' @rdname accessors
setMethod("snps",signature(object="groups"), function(object) { object@.Data })
##' @rdname accessors
setMethod("tags",signature(object="groups"), function(object) { object@tags })
##' @rdname accessors
setMethod("snps",signature(object="tags"), function(object) { object@.Data })
##' @rdname accessors
setMethod("tags",signature(object="tags"), function(object) { object@tags })

## setMethod("snpdrop",signature(x="snpmod",y="character"),
##           function(x,y) {
##             wh <- which(sapply(x@model.snps,function(x) any(x %in% y)))
##             message(length(wh), " / ", length(x@model.snps), " models (",
##                     format.pval(100*length(wh)/length(x@model.snps)), "%) will be dropped.")
##             if(!length(wh))
##               return(d)
##             x@models <- x@models[-wh,]
##             x@model.snps <- x@model.snps[-wh]  
##             return(marg.snps(x))  
##           })

## setMethod("+",signature(e1="snpmod",e2="snpmod"), function(e1,e2) snpmod.add(e1,e2))

