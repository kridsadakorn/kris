#' Read binary PLINK files (BED, BIM and FAM)
#'
#' @description Reads a SNP dataset stored in the binary PLINK format. All
#' three files are needed: BED (genotypes), BIM (SNP information) and FAM
#' (individual information).
#'
#' @param bed Path to the BED file.
#' @param bim Path to the BIM file.
#' @param fam Path to the FAM file.
#' @param only.snp If TRUE, return only the SNP matrix. Default FALSE.
#'
#' @return A list with:
#' \itemize{
#' \item \code{$snp}: a matrix of SNPs in additive coding (0, 1, 2), with
#' individuals in rows and SNPs in columns.
#' \item \code{$snp.info}: a data.frame of SNP information from the BIM file.
#' \item \code{$ind.info}: a data.frame of individual information from the FAM
#' file.
#' }
#' If \code{only.snp = TRUE}, the list holds only \code{$snp}. If a file is
#' missing or the BED file is not valid, the function prints a message and
#' returns NA.
#'
#' @details The file format is described at
#' \url{https://www.cog-genomics.org/plink/1.9/formats#bed}.
#'
#' @export
#' @md
#'
#' @importFrom utils read.table
#'
#' @seealso \code{\link{write.bed}}
#'
#' @examples
#'
#' #Use the example files embedded in the package.
#' bed <- system.file("extdata", "example_SNP.bed", package="KRIS")
#' bim <- system.file("extdata", "example_SNP.bim", package="KRIS")
#' fam <- system.file("extdata", "example_SNP.fam", package="KRIS")
#' snp <- read.bed(bed, bim, fam )
#'
#' #Check the objects inside 'snp'
#' ls(snp)
#'
#' #Preview $snp
#' print(snp$snp[1:10, 1:10])
#'
#' #Preview $snp.info
#' head(snp$snp.info)
#'
#' #Preview $ind.info
#' head(snp$ind.info)

read.bed <- function(bed, bim, fam, only.snp = FALSE){

  ret = NA

  if (!file.exists(bed)){
    cat(paste0("Error: BED file doesn't exist: ",bed,"\n"))
    return(ret)
  }
  if (!file.exists(bim)){
    cat(paste0("Error: BIM file doesn't exist: ",bim,"\n"))
    return(ret)
  }
  if (!file.exists(fam)){
    cat(paste0("Error: FAM file doesn't exist: ",fam,"\n"))
    return(ret)
  }
  #Read BIM file
  snp.info <- read.table(bim,header=FALSE, colClasses = c("factor","factor","factor","factor","factor","factor"))
  colnames(snp.info) = c("chr","ID","GD","position","allele1","allele2")

  #Read FAM file
  ind.info = read.table(fam,header=FALSE, colClasses = c("factor","factor","factor","factor","factor","factor"))
  colnames(ind.info) = c("FamID","IndID","PatID","MatID","sex","phenotype")

  no.ind = dim(ind.info)[1]
  no.snp = dim(snp.info)[1]

  if (only.snp == TRUE){
    snp.info = NA
    ind.info = NA
  }

  #Read BED file
  fh = file(bed, 'rb')
  #Read the first three bytes to check file format
  buff = readBin(fh, what="raw",n=3)
  if (sum(buff[1:2] == c('6c','1b')) != 2){
    cat(paste0("Error: BED file is not in a correct format: ",bed,"\n"))
    return(ret)
  }
  snp.major = (buff[3] == '01')
  if (snp.major){
    # Each SNP is a block of ceiling(no.ind/4) bytes
    no.byte.to.read = ceiling(no.ind/4.0)
    no.loop = no.snp
    no.keep = no.ind
  }else{
    # Individual-major: each individual is a block of ceiling(no.snp/4) bytes
    no.byte.to.read = ceiling(no.snp/4.0)
    no.loop = no.ind
    no.keep = no.snp
  }

  snp = readBin(fh, what="raw",n=(no.byte.to.read*no.loop))
  close(fh)
  if (length(snp) != no.byte.to.read*no.loop){
    cat(paste0("Error: BED file is shorter than expected from the BIM and FAM files: ",bed,"\n"))
    return(ret)
  }

  # Decode whole blocks at a time, so padding bits at the end of each block
  # are dropped in the right place. Two bits per genotype, low bit first:
  # 00 = 0, 10 = 1 (heterozygous), 11 = 2, 01 = missing.
  decode.blocks <- function(raw.block, n.block){
    bits = matrix(as.integer(rawToBits(raw.block)),ncol=2,byrow=TRUE)
    geno = bits[,1]+bits[,2]
    geno[which((bits[,2]-bits[,1])<0)] = NA
    matrix(geno, ncol=n.block)[1:no.keep, , drop=FALSE]
  }

  block.per.chunk = max(1, ceiling(no.loop / 10))
  res = vector("list", ceiling(no.loop / block.per.chunk))
  for (i in seq_along(res)){
    b1 = (i-1)*block.per.chunk + 1
    b2 = min(i*block.per.chunk, no.loop)
    idx = ((b1-1)*no.byte.to.read + 1):(b2*no.byte.to.read)
    res[[i]] = decode.blocks(snp[idx], b2-b1+1)
  }
  snp = do.call(cbind, res)

  if (!snp.major){
    snp = t(snp)
  }

  if (only.snp == FALSE){
    return(list("snp"=snp,"snp.info"=snp.info,"ind.info"=ind.info))
  }else{
    return(list("snp"=snp))
  }
}


#' Write a SNP object to binary PLINK files (BED, BIM and FAM)
#'
#' @description Saves a SNP object, such as one returned by
#' \code{\link{read.bed}}, as BED, BIM and FAM files.
#'
#' @param object A list with \code{$snp}, \code{$snp.info} and
#' \code{$ind.info}; see Details.
#' @param file Output path without extension. The function adds .bed, .bim and
#' .fam. Existing files are overwritten.
#'
#' @return \code{NULL}, invisibly.
#'
#' @details \code{object} must contain:
#' \itemize{
#' \item \code{$snp}: a matrix of SNPs in additive coding (0, 1, 2), with
#' individuals in rows and SNPs in columns.
#' \item \code{$snp.info}: a data.frame with the 6 columns of a BIM file.
#' \item \code{$ind.info}: a data.frame with the 6 columns of a FAM file.
#' }
#' The file format is described at
#' \url{https://www.cog-genomics.org/plink/1.9/formats#bed}.
#'
#' @export
#'
#' @importFrom utils write.table
#'
#' @seealso \code{\link{read.bed}}
#'
#' @examples
#'
#'
#' #Load example data
#' data(example_SNP)
#'
#' #Save 'simsnp' to the file as defined in 'save.file'
#' save.file <- file.path(tempdir(),"new_SNP")
#' write.bed(simsnp , save.file)

write.bed <- function(object, file){

  bed = paste0(file,'.bed')
  bim = paste0(file,'.bim')
  fam = paste0(file,'.fam')
  if (file.exists(bed)){
    cat(paste0("Overwriting the existing bed file\n"))
  }
  if (file.exists(bim)){
    cat(paste0("Overwriting the existing bim file\n"))
  }
  if (file.exists(fam)){
    cat(paste0("Overwriting the existing fam file\n"))
  }

  #Write BIM file
  if (is.null(object$snp.info)){
    cat(paste0("Couldn't find object$snp.info: a bim file was not created\n"))
  }else{
    write.table(object$snp.info,file=bim,quote=F,row.names=F,col.names=F)
  }

  #Write FAM file
  if (is.null(object$ind.info)){
    cat(paste0("Couldn't find object$ind.info: a fam file was not created\n"))
  }else{
    write.table(object$ind.info,file=fam,quote=F,row.names=F,col.names=F)
  }

  #Header of BED file
  header = c(0,0,1,1,0,1,1,0,1,1,0,1,1,0,0,0,1,0,0,0,0,0,0,0)

  no.ind = dim(object$ind.info)[1]
  no.snp = dim(object$snp.info)[1]

  fill.up = 4-(no.ind %% 4)
  tmp = NA
  if (fill.up != 4){
    tmp = matrix(rep(0,no.snp*fill.up),ncol=no.snp)
  }

  vec = object$snp
  if (!anyNA(tmp)){
    vec = rbind(vec,tmp)
  }

  vec = as.vector(vec)
  buff = matrix(rep(0,length(vec)*2),ncol=2)
  buff[which(vec == 1),1] = 0
  buff[which(vec == 1),2] = 1
  buff[which(vec == 2),1] = 1
  buff[which(vec == 2),2] = 1
  buff[which(is.na(vec)),1] = 1
  buff[which(is.na(vec)),2] = 0
  vec = NULL
  buff = c(header,as.vector(t(buff)))
  buff=packBits(as.raw(buff))

  fh = file(bed, 'wb')
  writeBin(buff,fh)
  close(fh)

  invisible(NULL)
}




