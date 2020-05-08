#' Get Allele frequencies for SNPs
#'
#' @param variants List of variants (rsids)
#' @param bfile If this is provided then will use the API. Default = NULL
#' @param plink_bin If null and bfile is not null then will detect packaged plink binary for specific OS. Otherwise specify path to plink binary. Default = NULL
#'
#' @export
#' @return matrix with allele frequencies
af_matrix <- function(variants, bfile=NULL, plink_bin = getOption("tools_plink")){
    if(!is.null(bfile)){
        return(af_matrix_local(variants, bfile=bfile, plink_bin=plink_bin))
    }
}




#' Get AF matrix using local plink binary and reference dataset
#'
#' @param variants List of variants (rsids)
#' @param bfile Path to bed/bim/fam ld reference panel
#' @param plink_bin Specify path to plink binary. Default = NULL. See https://github.com/explodecomputer/plinkbinr for convenient access to plink binaries
#'
#' @export
#' @return data frame
af_matrix_local <- function(variants, bfile, plink_bin, with_alleles=TRUE){
    ## Make textfile
    shell <- ifelse(Sys.info()['sysname'] == "Windows", "cmd", "sh")
    fn <- tempfile()
    write.table(data.frame(variants), file=fn, row.names=F, col.names=F, quote=F)
    
    
    fun1 <- paste0(
        shQuote(plink_bin, type=shell),
        " --bfile ", shQuote(bfile, type=shell),
        " --extract ", shQuote(fn, type=shell), 
        " --freq gz ", 
        " --out ", shQuote(fn, type=shell)
    )
    system(fun1)
    res <- read.table(paste0(fn, ".frq.gz"), header=T) %>%
        dplyr::select(SNP, MAF) %>%
        column_to_rownames("SNP") %>%
        as.matrix()
    return(res)
}
