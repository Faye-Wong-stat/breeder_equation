project.name <- "phenomic_pred"
project.dir <- paste0("~/breeder_equation_project/data/interest_of_phenomic_prediction/", 
                      project.name)
data.dir <- paste0(project.dir, "/data")
results.dir <- paste0(project.dir, "/results")
src.dir <- paste0(project.dir,"/src")



suppressPackageStartupMessages(library(readxl))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(apercu))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(signal))
suppressPackageStartupMessages(library(prospectr))
suppressPackageStartupMessages(library(ade4))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(viridis))
suppressPackageStartupMessages(library(VCA))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(rnirs)) # on GitHub
suppressPackageStartupMessages(library(doParallel)) # for parallelization
suppressPackageStartupMessages(library(foreach)) # for parallelization
suppressPackageStartupMessages(library(snow)) # for creating cluster
source(paste0(src.dir,"/PS_functions.R")) # for reading asd files
source(paste0(src.dir,"/ortho.nirs.R"))



t0 <- proc.time()



nirs <- asd_read_dir(directory=paste0(data.dir, "/spectra/p279_wood_2020"))
dim(nirs) # 1896 x 2151
ap(nirs)
length(unique(rownames(nirs))) #1896



p2f <- paste0(data.dir, "/spectra/tableau manip NIRS-bois_final.xlsx")
tools::md5sum(p2f)
dat_geno <- read_xlsx(p2f,sheet = "CC279", na = "x")
dat_geno <- as.data.frame(dat_geno)
dim(dat_geno) # 510 x 12
head(dat_geno)
str(dat_geno)

dat_geno2 <- reshape2::melt(dat_geno,id.vars=c("Pop", "Bloc", "Rang", "Plac",
                                               "Id", "Pop-id", "Bl-Rg-Pl","Remarque"),
                            variable.name="NIRS.origin", value.name="NIRS.nb")
head(dat_geno2)
dat_geno2$Id <- as.factor(dat_geno2$Id)
dat_geno2$Bloc <- as.factor(dat_geno2$Bloc)
dat_geno2$NIRS.nb <- as.numeric(dat_geno2$NIRS.nb)
# remove missing values
dat_geno2 <- dat_geno2[!is.na(dat_geno2$NIRS.nb),]



sample.info <- data.frame("NIRS.id" = rownames(nirs),stringsAsFactors = FALSE)
sample.info$NIRS.nb <- as.numeric(substring(sample.info$NIRS.id, 9))
head(sample.info)

summary(sample.info$NIRS.nb)  # min 1 max 1896
summary(dat_geno2$NIRS.nb) # min 1 max 1896

stopifnot(length(setdiff(sample.info$NIRS.nb, dat_geno2$NIRS.nb)) == 0,
          length(setdiff(dat_geno2$NIRS.nb,sample.info$NIRS.nb)) == 0)

length(unique(sample.info$NIRS.nb)) # 1896
length(unique(dat_geno2$NIRS.nb)) # 1896
stopifnot(!any(duplicated(sample.info$NIRS.nb)),
          !any(duplicated(dat_geno2$NIRS.nb)),
          nrow(dat_geno2) == length(unique(dat_geno2$NIRS.nb)))

sample.info <- merge(sample.info, dat_geno2, by="NIRS.nb", all=FALSE) 
dim(sample.info) # 1896 x 11
identical(rownames(nirs), sample.info$NIRS.id)
sample.info$Id_bloc <- paste0(sample.info$Id, "_", sample.info$Bloc)



subpop <- as.data.frame(read_xlsx(paste0(data.dir, "/genotypic/CC279_GBS/origin_CC279.xlsx")))
sample.info <- merge(sample.info, subpop, by.x="Id", by.y="code.intro", all=FALSE)

infos_geno <- as.data.frame(read_xlsx(paste0(data.dir, "/genotypic/CC279_GBS/ids_279.xlsx"), na = ""))
sample.info <- merge(sample.info, infos_geno, by.x="cultivar.number", by.y="CodeVar_TL")

sample.infoCtl <- sample.info[sample.info$Id %in% c("TEMOIN"),]
sample.infop279 <- sample.info[!sample.info$Id %in% c("TEMOIN"),]

dim(sample.infoCtl)
dim(sample.infop279)

head(sample.infop279)



nirs <- adj_asd(Xi = nirs, iadj  = c(651, 1451))



nirs.meanp279 <- apply(nirs.cut[sample.infop279$NIRS.nb,], 2, function(x) {
  tapply(X = x, INDEX = sample.infop279$Id_bloc, FUN = mean)
})
# For controls, several controls per bloc
nirs.meanCtl <- apply(nirs.cut[sample.infoCtl$NIRS.nb,], 2, function(x) {
  tapply(X = x, INDEX = sample.infoCtl$Id_bloc, FUN = mean)
})
nirs.mean <- rbind(nirs.meanp279, nirs.meanCtl)



dim(nirs.mean)

# create a new data frame with one line per Id_bloc
idx <- !duplicated(sample.info$Id_bloc)
sample.info.mean <- data.frame(Id_bloc=sample.info$Id_bloc[idx],
                               Id=sample.info$Id[idx], NIRS=I(nirs.mean),
                               stringsAsFactors=FALSE)
dim(sample.info.mean)

out <- ortho.nirs(data=sample.info.mean, SpName="NIRS", corr.effect=NULL,
                  max.effect="Id", ncomp=1, plot=TRUE, verbose=1)
head(out$loadings)
dim(out$loadings)
dim(out$corrected)
nirs.mean.corr <- out$corrected

nirs.mean.corr <- nirs.mean # don't apply correction for the moment



nirs.smooth <- t(apply(nirs.mean.corr, 1, function(x) {
  sgolayfilt(x, p = 1, n = 21, m = 0, ts = 1)}))
colnames(nirs.smooth) <- colnames(nirs.mean.corr)



nirs.dt <- prospectr::detrend(X = nirs.smooth,
                              wav = as.integer(colnames(nirs.smooth)))



nirs.der2 <- t(apply(nirs.snv, 1, function(x) {
  sgolayfilt(x, p = 3, n = 81, m = 2, ts = 1)}))

colnames(nirs.der2) <- colnames(nirs.snv)



nirs.list <- list("raw" = nirs.mean, "smooth" = nirs.smooth,
                  "snv" = nirs.snv, "dt" = nirs.dt,
                  "der1" = nirs.der1, "der2" = nirs.der2)
nirs.listp279 <- lapply(nirs.list, function(x) {
  merge(sample.infop279, data.frame("id" = rownames(x), "NIRS" = I(x)),
        by.x="Id_bloc", by.y = "id", all=FALSE)
})
nirs.listp279 <- lapply(nirs.listp279, function(x){
  x[!duplicated(x$Id_bloc),]})
lapply(nirs.listp279, dim)



nirs.listCtl <- lapply(nirs.list, function(x) {
  merge(sample.infoCtl, data.frame("id" = rownames(x), "NIRS" = I(x)),
        by.x="Id_bloc", by.y = "id", all=FALSE)
})
nirs.listCtl <- lapply(nirs.listCtl, function(x){
  x[!duplicated(x$Id_bloc),]})
lapply(nirs.listCtl, dim)



nirs.list <- lapply(nirs.list, function(x) {
  merge(sample.info, data.frame("id" = rownames(x), "NIRS" = I(x)),
        by.x="Id_bloc", by.y = "id", all=FALSE)
})
nirs.list <- lapply(nirs.list, function(x){
  x[!duplicated(x$Id_bloc),]})
lapply(nirs.list, dim)

# save spectra
p2f <- paste0(results.dir,"/spectrum_correction/wood_p279_2020_nirs-list-2-blocks.Rdata")
save(nirs.list, file=p2f)


































