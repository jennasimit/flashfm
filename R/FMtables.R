#' @title Formats flashfm output to table format
#' @param PP MFM:::MPP.PP.groups.fn output from 'flashfm' in Rdata format.
#' @param SW stepwise list of results from Stepwise model. By def. = NULL
#' @param stepwise TRUE if list of results from sw is provided. By def. = FALSE
#' @param regname Region name for table caption
#' @param path.input path where the input files are located.
#' @param path.output path to save the output latex table.
#' @param trait.id id. number if the traits.
#' @param trait.names trait names. If not provided a vector is constructed. 'Trait_1', 'Trait_2', ...
#' @return Table results as txt file.
#' @author Nico Hernandez
#' @export
FMtables <- function(PP,SW, stepwise=F, regname, path.input, path.output, trait.id, trait.names=NULL){
         
load(paste0(path.input,PP,'.Rdata'))
load(paste0(path.input,SW,'.Rdata'))
mpp.pp<-get(PP)
sw<-get(SW)
  
# TRAITS
qt <- trait.id
if (is.null(trait.names)) {
    ts.aux <- rep('trait_',length(qt));
    ts<-c()
    for (i in 1:length(qt)){ts[i]<-paste0(ts.aux[i],i)}
} else {ts <- trait.names}

M = length(ts)
aux<-matrix(c(1:length(ts),qt),nrow=length(ts),ncol=2)

RES_FM=RES.sw=ST<-list()
for (i in 1:M){
    
#### STOCHASTIC FineMap 
a<-as.data.frame(mpp.pp$gPP[[i]])
a<-a[which(a[,1]>=0.05|a[,2]>=0.05),]
model<-rownames(a)
result<-data.frame()
for (j in 1:dim(a)[1]){
  result[j,1]<-as.character(model[j])
  result[j,2]<-round(a[j,1],3)
  result[j,3]<-as.character(model[j])
  result[j,4]<-round(a[j,2],3)
}
    
res.aux1<-result[,1:2];res.aux2<-result[,3:4]
result2<-cbind(res.aux1[order(-res.aux1$V2),],res.aux2[order(-res.aux2$V4),])
result2$V1[which(result2$V2<0.05)]<-0;result2$V3[which(result2$V4<0.05)]<-0
result2[result2 < 0.05] <- 0
    
aux<-c()
for (l in 1:nrow(result2)){
 aux[l]<-!all(result2[l,]==0)
}
result2<-result2[aux,]
result2[result2 ==0] <- '--'
RES_FM[[i]]<-result2

#### STEPWISE FineMap

if (stepwise==TRUE){
    
    chr.name<-strsplit(as.character(sw[[i]][,1]), split="\\,")
    snps.aux<-c()
    if (length(chr.name)==1){ snps.aux[1]<-chr.name[[1]]
    } else {
      snps.aux[1]<-chr.name[[1]]
      for (m in 2:length(chr.name)){
        snps.aux[m]<-chr.name[[m]][which(!(chr.name[[m]]%in%chr.name[[m-1]]))]
      }
    }
    
    snps.aux <- gsub('_',':',snps.aux)
    sw.model<-c()
    for (h in 1:dim(sw[[i]])[1]){
      if ( grepl("\\d", sw[[i]][h,2]) | is.na(sw[[i]][h,2]) ){  sw.model[h]<-snps.aux[h] }  
      else {sw.model[h]<-paste0(snps.aux[h],'/',sw[[i]][h,2])}
    }
    
    aux.res.sw<-as.data.frame(cbind(sw.model, sw[[i]][,3]))
    aux.res.sw[,2]<-signif(as.numeric(levels(aux.res.sw$V2))[aux.res.sw$V2],digits=3)
    RES.sw[[i]]<-aux.res.sw
    
    ##### MERGING RESUTLS
    ST[[i]]<-RES_FM[[i]]
    ST[[i]][is.na(ST[[i]])] <- '--'
        
    if(dim(ST[[i]])[1]==1){
      ST[[i]]<-cbind(RES.sw[[i]],ST[[i]])
        } else if (dim(ST[[i]])[1]>dim(RES.sw[[i]])[1]) {
                comp<-as.data.frame(matrix(rep('--',(dim(ST[[i]])[1]-dim(RES.sw[[i]])[1])*dim(RES.sw[[i]])[2]),nrow = dim(ST[[i]])[1]-dim(RES.sw[[i]])[1]))
                colnames(comp)<-colnames(RES.sw[[i]])
                RES.sw[[i]]<-rbind(RES.sw[[i]],comp)
                ST[[i]]<-cbind(RES.sw[[i]],ST[[i]])
        } else if (dim(ST[[i]])[1]<dim(RES.sw[[i]])[1]) {
          comp<-as.data.frame(matrix(rep('--',(dim(RES.sw[[i]])[1]-dim(ST[[i]])[1])*dim(ST[[i]])[2]),nrow = dim(RES.sw[[i]])[1]-dim(ST[[i]])[1]))
          colnames(comp)<-colnames(ST[[i]])
          ST[[i]]<-rbind(ST[[i]],comp)
          ST[[i]]<-cbind(RES.sw[[i]],ST[[i]])
        } else { ST[[i]]<-cbind(RES.sw[[i]],ST[[i]]) }
    
    } else {
      
      ST[[i]]<-RES_FM[[i]]
      ST[[i]][is.na(ST[[i]])] <- '--'
  
    }
        
} # END LOOP
       
# PREPARING XTABLE
if (stepwise==TRUE){
names(ST)<-c(ts)
ST<-data.table::rbindlist(ST,idcol = T)
colnames(ST)<-c('Traits','SNP/Model','P-value','Model','PP','Model', 'PP')
ST$Traits[which(duplicated(ST$Traits))]<-''
cols <- colnames(ST)
addtorow <- list()
addtorow$pos <- list(0,0)
addtorow$command <- c(paste0("&\\multicolumn{2}{c}{Stepwise} & \\multicolumn{2}{c}{FineMap} & \\multicolumn{2}{c}{FM-FlashFM}\\\\\n"), 
                        paste(paste(cols, collapse=" & "), "\\\\\n") )
TABLE<-print(xtable::xtable(ST, caption = regname,
                             align = c("l","l","c","c","c","c","c","c")), add.to.row=addtorow, include.colnames=F, include.rownames = F,
                             NA.string="-", booktabs = F)
} else {
  
  names(ST)<-c(ts)
  ST<-data.table::rbindlist(ST,idcol = T)
  colnames(ST)<-c('Traits','Model','PP','Model', 'PP')
  ST$Traits[which(duplicated(ST$Traits))]<-''
  cols <- colnames(ST)
  addtorow <- list()
  addtorow$pos <- list(0,0)
  addtorow$command <- c(paste0("&\\multicolumn{2}{c}{FineMap} & \\multicolumn{2}{c}{FM-FlashFM}\\\\\n"), 
                        paste(paste(cols, collapse=" & "), "\\\\\n") )
  TABLE<-print(xtable::xtable(ST, caption = regname,
                      align = c("l","l","c","c","c","c")), add.to.row=addtorow, include.colnames=F, include.rownames = F,
               NA.string="-", booktabs = F)
}

write.table(TABLE, paste0(path.output,'TABLE_',regname,'.txt'), col.names = F, row.names = F)
return(TABLE)

} 
