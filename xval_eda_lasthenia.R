layers_over_reps <- function(path, k, n_reps){
  block_pattern <- paste0('*_sp*',k,'*block.Robj')
  results_pattern <- paste0('*_sp*',k,'*results.Robj')
  fit_pattern <- paste0('*_sp*',k,'*fit.Robj')
  k <- as.integer(strsplit(k, '')[[1]][2])
  block_fnames <- list.files(path=path, pattern=glob2rx(block_pattern))
  res_fnames <- list.files(path=path, pattern=glob2rx(results_pattern))
  fit_fnames <- list.files(path=path, pattern=glob2rx(fit_pattern))
  #load('miss_0.8/cross_val/miss_0.8_bi_las_nsp_rep1K2_conStruct.results.Robj')
  load(paste0(path,res_fnames[1]))
  #load('miss_0.8/cross_val/miss_0.8_bi_las_sp_rep1K2_data.block.Robj')
  load(paste0(path,block_fnames[1]))
  layers_mat <- matrix(NA, nrow=k, ncol=n_reps)
  layers_mat[,1] <- c(calculate.layer.contribution(conStruct.results[[1]], data.block))
  tmp2 <- conStruct.results[[1]]$MAP$admix.proportions
  
  for(i in 2:n_reps){
    load(paste0(path,res_fnames[i]))
    load(paste0(path,block_fnames[i]))
    tmp2.order <- match.layers.x.runs(tmp2, conStruct.results[[1]]$MAP$admix.proportions)
    layers_mat[,i] <- c(calculate.layer.contribution(conStruct.results = conStruct.results[[1]],
                                                     data.block = data.block,
                                                     layer.order = tmp2.order))
  }
  
  pal <- createPalette(k,c("#ff0000", "#00ff00", "#0000ff"))
  oyy <- barplot(layers_mat,
                 col=pal,
                 ylab='layer contributions')
  return(oyy)
}
huh <- layers_over_reps('miss_1/cross_val/', 'K2', 8)

library(conStruct)
library(Polychrome)
library(ggplot2)
library(tidyr)
library(dplyr)


######## This will require extensive annotation to make clear what it is doing but I think this solution
#works for now. I'm not sure that is the case with the changes that I made to the way I now run construct,
#where we write multiple files. Might though. You never know. Regardless, it works here for Lasthenia.
admix_over_reps <- function(path, k, n_reps, ind_file){
  #create patterns for finding files
  results_pattern <- paste0('*_sp*',k,'*results.Robj')
  missing_pattern <- glob2rx('*dropped_individuals*txt')
  
  #read in indidviduals and make vector
  ind_file <- read.table(ind_file, header = F)
  #ALL INDIVIDUALS SAMPLED FOR ALL RUNS
  inds <- ind_file$V1
  #NUMBER OF INDIVIDUALS
  ninds <- length(inds)
  
  #get list of files containing missing individuals, ie. individuals dropped for missingness or duplication
  #ALL MISSING FILES
  missing_files <- list.files(pattern = missing_pattern, recursive = T, include.dirs = T)
  #NUMBER OF FILES
  num_files <- length(missing_files)
  # if that list actually contains files, then we remove those individuals from the ind file for labels
  if(num_files>0){
    miss_df <- data.frame(inds=NA)
    for(i in 1:num_files){
      tmp <-  read.table(missing_files[i])
      names(tmp) <- 'inds'
      miss_df <- miss_df %>%
        add_row(tmp)
    }
    #MISS DF CONTAINS ALL MISSING INDIVIDUALS OVER ALL MISSINGNESS RUNS
    miss_df <- na.omit(miss_df)
    #ONLY MISSING IS A VECTOR OF THE SAMPLE IDS RATHER THAN A DF
    only_missing <- unique(miss_df$inds)
    #NUMER OF MISSING INDIVIDUALS
    num_miss <- length(only_missing)
    #THE INDICES OF MISSING INDIVIDUALS IN THE TOTAL IND LIST
    miss <- which(inds %in% only_missing)
    #REMOVE THEM - KEEP INDS INTACT SO WE HAVE A FULL LIST OF INDIVIDUALS
    inds_for_filter <- inds[-miss]
  }
  

    #split the integer K we are interested in from the letter supplied
  k <- as.integer(strsplit(k, '')[[1]][2])
  
  #get list of file names for results
  res_fnames <- list.files(path=path, pattern=glob2rx(results_pattern))
  
  #load the first result obj which will act as a reference for orienting all other results admixture props
  #- note that this always loads these files into the environment as 'conStruct.results' 
  load(paste0(path,res_fnames[1]))
  single_missing_fn <- list.files(path=path, pattern=missing_pattern)
  if(length(single_missing_fn!=0)){
    ad_ref <- conStruct.results[[1]]$MAP$admix.proportions
    single_missing <- read.table(paste0(path,single_missing_fn))
    single_missing <- single_missing$x
    sing_miss <- which(inds %in% single_missing)
    rownames(ad_ref) <- inds[-sing_miss]
    ad_prop_1 <- as.data.frame(ad_ref[inds_for_filter,])
  }else{
    ad_ref <- conStruct.results[[1]]$MAP$admix.proportions
    rownames(ad_ref) <- inds
    ad_prop_1 <- as.data.frame(ad_ref[inds_for_filter,])
  }

  #construct column names
  new_names <- rep('pop',k)
  new_nam <- paste0(new_names, c(1:k))
  #give them column names
  names(ad_prop_1) <- new_nam
  #add new column for inds
  ad_prop_1$inds <- inds_for_filter

  #PIVOT
  ad_prop_1_long <- ad_prop_1 %>%
    pivot_longer(cols=c(new_nam[1:k]), names_to='pop', values_to='prop')
  #plot
  ad_1_plot <- ggplot(data = ad_prop_1_long, aes(y=prop,x=inds, fill=pop))+
     geom_bar(stat='identity')+
     theme(axis.text.x = element_text(angle=90, size=5))+
     scale_x_discrete()
  #create plot list and put current plot in first position.
  pl_list <- list()
  pl_list[[1]] <- ad_1_plot
  # set up our reference admixture
   

   for(i in 2:n_reps){
    
     load(paste0(path,res_fnames[i]))
     ad_new <- conStruct.results[[1]]$MAP$admix.proportions
     if(length(single_missing_fn) != 0){

       single_missing <- read.table(paste0(path,single_missing_fn))
       single_missing <- single_missing$x
       sing_miss <- which(inds %in% single_missing)
       rownames(ad_new) <- inds[-sing_miss]

       match_order <- conStruct::match.layers.x.runs(ad_ref[inds_for_filter,],
                                                     ad_new[inds_for_filter,])
       ad_prop <- as.data.frame(ad_new[inds_for_filter,c(match_order)])
     }else{
       print(ad_new)
       rownames(ad_new) <- inds
       ad_prop <- as.data.frame(ad_new[inds_for_filter,])
       match_order <- conStruct::match.layers.x.runs(ad_ref[inds_for_filter,],
                                                     ad_new[inds_for_filter,])
       ad_prop <- as.data.frame(ad_new[inds_for_filter,c(match_order)])
     }

     names(ad_prop) <- new_nam
     ad_prop$inds <- inds_for_filter
     ad_prop <- ad_prop %>%
       pivot_longer(cols=c(new_nam[1:k]), names_to='pop', values_to='prop')

     ad_plot <- ggplot(data = ad_prop, aes(y=prop,x=inds, fill=pop))+
       geom_bar(stat='identity')+
       theme(axis.text.x = element_text(angle=90, size = 5))+
       scale_x_discrete()
      pl_list[[i]] <- ad_plot
   }
   plot_grid(plotlist=pl_list)
   return(pl_list)
  
}
test_list <- admix_over_reps('miss_0.4/cross_val/', 'K3',8, '../both_lanes_high_cov_individuals.txt')

plot_grid(plotlist = test_list)


####### Unresolved ########
setwd('~/Desktop/Lasthenia/for_construct/output/')


block_pattern <- glob2rx('*_sp*K2*block.Robj')
results_pattern <- glob2rx('*_sp*K2*results.Robj')
fit_pattern <- glob2rx('*_sp*K2*fit.Robj')


######## admix over missingness ########

admix_over_missing <- function(k, ind_file){
  #block_pattern <- paste0('*_sp*',k,'*block.Robj')
  results_pattern <- paste0('*_sp*_rep1',k,'*results.Robj')
  #fit_pattern <- paste0('*_sp*',k,'*fit.Robj')
  missing_pattern <- glob2rx('*dropped_individuals*txt')
  
  #all_blocks <- list.files(, pattern = glob2rx(block_pattern), recursive = T, include.dirs = T)
  all_results <- list.files(, pattern=glob2rx(results_pattern), recursive = T, include.dirs = T)
  missing_vals <- strsplit(all_results, '/')
  missing_treats <- length(all_results)

  #all_fits <- list.files(,pattern = glob2rx(fit_pattern), recursive = T, include.dirs = T)
  all_missing <- list.files(, pattern=missing_pattern, recursive = T, include.dirs = T)
  num_mis_file <- length(all_missing)
  miss_list <- list()
  for(i in 1:num_mis_file){
    tmp_mis <- read.table(all_missing[i])
    
    miss_list[[i]] <- tmp_mis$x
  }

  missing_inds <- unlist(miss_list)
  uni_miss <- unique(missing_inds)
  individs <- read.table(ind_file)
  
  inds <- individs$V1
  all_inds <- length(inds)
  ind_match <- which(inds %in% uni_miss)
  inds <- inds[-ind_match]
  load(all_results[1])
  ad_prop_1 <- as.data.frame(conStruct.results[[1]]$MAP$admix.proportions)
  ad_ref <- conStruct.results[[1]]$MAP$admix.proportions
  k <- as.integer(strsplit(k, '')[[1]][2])
  
  new_names <- rep('pop',k)
  new_nam <- paste0(new_names, c(1:k))
  names(ad_prop_1) <- new_nam
  ad_prop_1$inds <- inds

  ad_prop_1_long <- ad_prop_1 %>%
    pivot_longer(cols=c(new_nam[1:k]), names_to='pop', values_to='prop')
  miss_val <- missing_vals[[1]][1]
  ad_1_plot <- ggplot(data = ad_prop_1_long, aes(y=prop,x=inds, fill=pop))+
    geom_bar(stat='identity')+
    ggtitle(miss_val)+
    theme(axis.text.x = element_text(angle=90))+
    scale_x_discrete()
  ad_1_plot
  
  
  pl_list <- list()
  pl_list[[1]] <- ad_1_plot
  
  for(i in 2:missing_treats){
    load(all_results[i])

    num_inds <- nrow(conStruct.results[[1]]$MAP$admix.proportions)
    num_missing <- all_inds - num_inds

    if(num_missing==0){
      match_order <- conStruct::match.layers.x.runs(ad_ref,
                                                    conStruct.results[[1]]$MAP$admix.proportions[-ind_match,])
      ad_prop <- as.data.frame(conStruct.results[[1]]$MAP$admix.proportions[-ind_match,c(match_order)])
      names(ad_prop) <- new_nam
      ad_prop$inds <- inds
      ad_prop <- ad_prop %>%
        pivot_longer(cols=c(new_nam[1:k]), names_to='pop', values_to='prop')
      #print(ad_prop_1_long)
      miss_val_i <- missing_vals[[i]][1]
      ad_plot <- ggplot(data = ad_prop, aes(y=prop,x=inds, fill=pop))+
        geom_bar(stat='identity')+
        theme(axis.text.x = element_text(angle=90))+scale_x_discrete()+
        ggtitle(miss_val_i)

      pl_list[[i]] <- ad_plot
    }else if (num_missing==1){

      match_order <- conStruct::match.layers.x.runs(ad_ref,
                                                    conStruct.results[[1]]$MAP$admix.proportions[-ind_match[2],])
      ad_prop <- as.data.frame(conStruct.results[[1]]$MAP$admix.proportions[-ind_match[2],c(match_order)])
      names(ad_prop) <- new_nam
      ad_prop$inds <- inds
      ad_prop <- ad_prop %>%
        pivot_longer(cols=c(new_nam[1:k]), names_to='pop', values_to='prop')
      #print(ad_prop_1_long)
      miss_val_i2 <- missing_vals[[i]][1]
      ad_plot <- ggplot(data = ad_prop, aes(y=prop,x=inds, fill=pop))+
        geom_bar(stat='identity')+
        theme(axis.text.x = element_text(angle=90))+scale_x_discrete()+
        ggtitle(miss_val_i2)

      pl_list[[i]] <- ad_plot
    }else if(num_missing==2){
      match_order <- conStruct::match.layers.x.runs(ad_ref,
                                                    conStruct.results[[1]]$MAP$admix.proportions)
      ad_prop <- as.data.frame(conStruct.results[[1]]$MAP$admix.proportions[,c(match_order)])
      names(ad_prop) <- new_nam
      ad_prop$inds <- inds
      ad_prop <- ad_prop %>%
        pivot_longer(cols=c(new_nam[1:k]), names_to='pop', values_to='prop')
      #print(ad_prop_1_long)
      miss_val_i3 <- missing_vals[[i]][1]
      ad_plot <- ggplot(data = ad_prop, aes(y=prop,x=inds, fill=pop))+
        geom_bar(stat='identity')+
        theme(axis.text.x = element_text(angle=90))+scale_x_discrete()+
        ggtitle(miss_val_i3)
      pl_list[[i]] <- ad_plot
    }
    

  }
  
  plot_grid(plotlist=pl_list)
  
  return(pl_list)
}


test <- conStruct.results$chain_1$MAP$admix.proportions
indy <- read.table('../both_lanes_high_cov_individuals.txt')
indy$V1
rownames(test) <- indy$V1
test["BTM1_8",]
huh <- admix_over_missing('K3', '../both_lanes_high_cov_individuals.txt')
huh
plot_grid(plotlist=huh)

load("miss_1/cross_val/miss_1_bi_las_sp_rep1K2_conStruct.results.Robj")

ad <- conStruct.results[[1]]$MAP$admix.proportions

conStruct.results$chain_1$MAP$index.iter
setwd('~/Desktop/Lasthenia/for_construct/')
wel <- list.files(, pattern = block_pattern, recursive = T, include.dirs = T)
grep('miss_0.6', wel)
wel
plot_grid(plotlist=mget(paste0("test_plot", 1:8)))

plot_grid(test_list)
