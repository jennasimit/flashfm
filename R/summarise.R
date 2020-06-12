# internal function for PPmodGroups
 gmodK.fn <- function(modsnps,snpgroups) {
 tmp <- as.character(modsnps)
 msnps <- unlist(strsplit(tmp,"%"))

 gsnps <- groupIDs.fn(snpgroups,msnps)
 out <- paste(gsnps[order(gsnps)],collapse="%")

 return(out)

 }



PPmodGroups <- function(modPP,snpgroups,minPP=0.05) {

 snpmods <- names(modPP)
 pp <- modPP
 Gmods <- apply(as.matrix(snpmods),1,gmodK.fn,snpgroups)
 modpp <- data.frame(model=Gmods,PP=pp)
 mods <- unique(Gmods)
 GPP <- apply(as.matrix(mods),1,function(x) sum(modpp$PP[modpp$model==x]))
 out <- data.frame(model=mods,PP=GPP)
 out <- out[order(out$PP,decreasing=TRUE),]
 Gout <- out[out$PP>minPP,]
 rownames(Gout) <- Gout$model
 Sout <- data.frame(SNPmod=snpmods,Gmod=Gmods,PP=pp)
 return(list(snp=Sout,group=Gout))
 }



MPPmodGroups <- function(PP1) {

    mnames <- rownames(PP1)
    msep <- apply(matrix(1:length(mnames), ncol = 1), 1, MFM:::sep.fn, 
        mnames)
    gnames <- unique(unlist(msep))
    mpp1 <- NULL   
        tmp1 <- apply(matrix(1:length(mnames), ncol = 1), 1, 
            MFM:::check.fn, msep, PP1$PP, gnames)
        mpp1 <- apply(as.matrix(tmp1), 1, sum)   
    mpp1 <- data.frame(mpp1, row.names = gnames)  
    return(mpp1)
}


#' @title Summarise PP and MPP results from single-trait fine-mapping and flashfm, by SNP and by SNP group
#' @param fm.multi output from flashfm function
#' @param snpGroups  list of two sets of snp groups; output from makeSNPgroups2
#' @param minPP a single value such that output consists of snp/group models where PP > minPP for at least one set of fine-mapping results; default is 0.05 
#' @return a list of 4 objects: MPP lists trait-specific PP of SNP inclusion in a model; MPPg lists trait-specific PP of SNP group inclusion in a model;  
#'			PP lists trait-specific model PP; PPg lists trait-specific model PP in terms of  SNP group
#' @author Jenn Asimit
#' @export
PPsummarise <- function (fm.multi, snpGroups, minPP=0.05) 
# for one TOdds setting in fm.multi, snpGroups is a list of 2 snpgroups
{

	MPP <- fm.multi$MPP
	pp <- fm.multi$PP	
	M <- length(pp)
	snpgroups <- snpGroups[[1]]
	snpgroups2 <- snpGroups[[2]]
	tnames <- names(fm.multi$PP)
	
	PPout <- PPout2 <- outPPg <- outPP <- vector("list",M)
	MPPout <- MPPout2 <- outMPPg <- outMPP <- vector("list",M)
	
	for(i in 1:M) {
	 PPout[[i]] <- PPmodGroups(fm.multi$PP[[i]][,1],snpgroups, minPP)
	 PPout2[[i]] <- PPmodGroups(fm.multi$PP[[i]][,2],snpgroups2, minPP)
	 pplist <- list(t(data.frame(PPout[[i]]$group[,2],row.names=rownames(PPout[[i]]$group))),t(data.frame(PPout2[[i]]$group[,2],row.names=rownames(PPout2[[i]]$group))))
	 outPPg[[i]] <- t(do.call("smartbind",c(pplist,fill=0)))
	 colnames(outPPg[[i]]) <- fm.multi$sharing
	 pplist <- list(t(data.frame(PPout[[i]]$snp[,3],row.names=rownames(PPout[[i]]$snp))),t(data.frame(PPout2[[i]]$snp[,3],row.names=rownames(PPout2[[i]]$snp))))
	 outPP[[i]] <- t(do.call("smartbind",c(pplist,fill=0)))
	 colnames(outPP[[i]]) <- fm.multi$sharing
	 
	 MPPout[[i]] <- MPPmodGroups(PPout[[i]]$snp)
	 MPPout2[[i]] <- MPPmodGroups(PPout2[[i]]$snp)
	 pplist <- list(t(data.frame(MPPout[[i]][,1],row.names=rownames(MPPout[[i]]))),t(data.frame(MPPout2[[i]][,1],row.names=rownames(MPPout2[[i]]))))
	 outMPP[[i]] <- t(do.call("smartbind",c(pplist,fill=0)))
	 colnames(outMPP[[i]]) <- fm.multi$sharing
	 
	  MPPout[[i]] <- MPPmodGroups(PPout[[i]]$group)
	  MPPout2[[i]] <- MPPmodGroups(PPout2[[i]]$group)
	  pplist <- list(t(data.frame(MPPout[[i]][,1],row.names=rownames(MPPout[[i]]))),t(data.frame(MPPout2[[i]][,1],row.names=rownames(MPPout2[[i]]))))
	  outMPPg[[i]] <- t(do.call("smartbind",c(pplist,fill=0)))
	  colnames(outMPPg[[i]]) <- fm.multi$sharing
	  }
	  
	  
	 names(outMPPg) <- names(outMPP) <- names(outPP)  <- names(outPPg) <- tnames
	 
	return(list(MPP=outMPP, MPPg=outMPPg, PP=outPP, PPg=outPPg)) 
	}


 



