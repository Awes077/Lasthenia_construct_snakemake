

library(geosphere)
library(conStruct)
library(stringr)
library(dplyr)
library(parallel)
library(doParallel)
library(foreach)
args=commandArgs(trailingOnly = T)
if(length(args)!=5){
  stop('This script requires five arguments: input vcf, csv containing coordinates for each sample, k value, ending for ouput, and file containing all individual names')
}else{
  input <- args[1]
  coords <- args[2]
  k_val <- args[3]
  output_name <- args[4]
  ind_file <- args[5]
}


freq_matrix <- readRDS(input)
dup_inds <- which(duplicated(freq_matrix))
ind_f <- read.table(ind_file)
ind_f <- read.table('Desktop/Diascia/clust_vcfs/For_construct/output/inds_diascia.txt')

inds <- ind_f$V1
missing_val <- unlist(str_split(input, pattern='_',n=3))[2]
missing_val <- unlist(str_split(missing_val,pattern='/'))[1]
pwd <- getwd()
output_path <- paste0(pwd,'/output/miss_',missing_val,
                      '/cross_val',
                      '/miss_',missing_val,'_',output_name)

coords <- read.csv(coords)
if(length(dup_inds)>0){
  freq_matrix <- freq_matrix[-dup_inds,]
  coords <- coords[-dup_inds,]
  dropped_inds <- inds[dup_inds]
  write.table(dropped_inds,paste0(output_path,'1_dropped_individuals.txt'), quote = F)
}
setwd('Desktop/Diascia/clust_vcfs/For_construct/output/')
write.table(as.data.frame(inds), 'help.txt', quote = F)

if(as.numeric(missing_val)<1){
  freq_boo <- is.na(freq_matrix)
  sum_uh<-as.data.frame(freq_boo)%>%mutate(sum = rowSums(across(where(is.logical))))
  sum_uh$prop <- sum_uh$sum/(ncol(freq_matrix)-1)
  missing_inds <- which(sum_uh$prop>0.7)
  
  
  if(length(missing_inds)>0){
    freq_matrix <- freq_matrix[-missing_inds,]
    coords <- coords[-missing_inds,]
    distances <- distm(coords[,c(2,1)])
    dropped_inds <- inds[missing_inds]
    write.table(dropped_inds,paste0(output_path,'2_dropped_individuals.txt'), quote = F)
  }else{
    distances <- distm(coords[,c(2,1)])
  }
}else{
  distances <- distm(coords[,c(2,1)])
}
cl <- makeCluster(20, type='FORK')
registerDoParallel(cl)
model <- x.validation(
  train.prop=0.9,
  n.reps=8,
  K=as.integer(1:k_val),
  freqs=freq_matrix,
  geoDist = distances,
  coords = as.matrix(coords),
  prefix = output_path,
  n.chains=1,
  n.iter=60,
  parallel = T,
  cores=8,
  save.files = T,
  make.figs = T,
  control = setNames(list(0.9,15), c('adapt_delta','max_treedepth'))
)
stopCluster(cl)

