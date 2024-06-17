
library(conStruct)
setwd('~/Desktop/Lasthenia/for_construct/output/')

#### xval comparison missing 1 #####
sp_x_val_results <- read.table('miss_1/cross_val/miss_1_bi_las_sp_xval_results.txt', header=T)
nsp_x_val_results <- read.table('miss_1/cross_val/miss_1_bi_las_nsp_xval_results.txt', header=T)
sp.CIs <-  apply(sp_x_val_results,1,function(x){mean(x) + c(-1.96,1.96) * sd(x)/length(x)})
nsp.CIs <- apply(nsp_x_val_results,1,function(x){mean(x) + c(-1.96,1.96) * sd(x)/length(x)})
plot(rowMeans(nsp_x_val_results),
     pch=19,col="blue",
     ylab="predictive accuracy",xlab="values of K",
     ylim=range(nsp.CIs,sp.CIs),
     main="cross-validation results")
points(rowMeans(sp_x_val_results),col="green",pch=19)

segments(x0 = 1:nrow(nsp_x_val_results),
         y0 = nsp.CIs[1,],
         x1 = 1:nrow(nsp_x_val_results),
         y1 = nsp.CIs[2,],
         col = "blue",lwd=2)




######## x val comparison missing 0.4 ######

sp_x_val_results <- read.table('miss_0.4/cross_val/miss_0.4_bi_las_sp_xval_results.txt', header=T)
nsp_x_val_results <- read.table('miss_0.4/cross_val/miss_0.4_bi_las_nsp_xval_results.txt', header=T)
sp.CIs <-  apply(sp_x_val_results,1,function(x){mean(x) + c(-1.96,1.96) * sd(x)/length(x)})
nsp.CIs <- apply(nsp_x_val_results,1,function(x){mean(x) + c(-1.96,1.96) * sd(x)/length(x)})
plot(rowMeans(sp_x_val_results),
     pch=19,col="blue",
     ylab="predictive accuracy",xlab="values of K",
     ylim=range(sp.CIs),
     main="cross-validation results")
points(rowMeans(sp_x_val_results),col="green",pch=19)

segments(x0 = 1:nrow(sp_x_val_results),
         y0 = sp.CIs[1,],
         x1 = 1:nrow(sp_x_val_results),
         y1 = sp.CIs[2,],
         col = "blue",lwd=2)


###### admix k=4 missing 0.4 ########


load('miss_0.4/cross_val/miss_0.4_bi_las_sp_rep4K4_conStruct.results.Robj')
load('miss_0.4/cross_val/miss_0.4_bi_las_sp_rep4K4_data.block.Robj')
load('miss_0.4/cross_val/miss_0.4_bi_las_sp_rep4K4_model.fit.Robj')
m_4_k_4_results <- conStruct.results
m_4_k_4_data_block <- data.block
m_4_k_4_mod_fit <- model.fit


ind_counts <- read.table('../../both_lanes_high_cov_individuals.txt', header = F)



admi_44 <- as.data.frame(m_4_k_4_results$chain_1$MAP$admix.proportions)

admi_44$V5 <- inds[1:17]

names(admi_44) <- c('pop1', 'pop2', 'pop3', 'pop4','ind')

admi_long_44 <- admi_44 %>%
  pivot_longer(cols=c(pop1, pop2, pop3, pop4), names_to='pop', values_to='prop')

ad_44 <- ggplot(data = admi_long_44, aes(y=prop,x=ind, fill=pop))+
  ylab('Admixture proportion')+xlab('Individual')+
  geom_bar(stat='identity')+
  theme(axis.text.x = element_text(angle=90))+ggtitle('60% Missing Allowed')

ad_44

load('miss_0.4/cross_val/miss_0.4_bi_las_nsp_rep1K3_conStruct.results.Robj')
load('miss_0.4/cross_val/miss_0.4_bi_las_sp_rep1K3_data.block.Robj')

meh <- c(calculate.layer.contribution(conStruct.results[[1]], data.block))
tmp <- conStruct.results[[1]]$MAP$admix.proportions
m4_k3_8rep_layers <- matrix(NA, nrow=3, ncol=8)
m4_k3_8rep_layers[,1] <- meh
for(i in 2:8){
  
  
  load(sprintf('miss_0.4/cross_val/miss_0.4_bi_las_nsp_rep%sK3_conStruct.results.Robj',i))
  load(sprintf('miss_0.4/cross_val/miss_0.4_bi_las_sp_rep%sK3_data.block.Robj',i))
  tmp.order <- match.layers.x.runs(tmp, conStruct.results[[1]]$MAP$admix.proportions)
  m4_k3_8rep_layers[,i] <- c(calculate.layer.contribution(conStruct.results = conStruct.results[[1]],
                                                          data.block = data.block,
                                                          layer.order = tmp.order))
  
}

m4_k3_8rep_layers

barplot(m4_k3_8rep_layers,
        col=c('blue','red','goldenrod1'),
        ylab='layer contributions')

layers_m4_k2 <- matrix(NA, nrow=2, ncol=8)

load('miss_0.4/cross_val/miss_0.4_bi_las_nsp_rep1K2_conStruct.results.Robj')
load('miss_0.4/cross_val/miss_0.4_bi_las_sp_rep1K2_data.block.Robj')
layers_m4_k2[,1] <- c(calculate.layer.contribution(conStruct.results[[1]], data.block))
tmp2 <- conStruct.results[[1]]$MAP$admix.proportions

for(i in 2:8){
  load(sprintf('miss_0.4/cross_val/miss_0.4_bi_las_nsp_rep%sK2_conStruct.results.Robj',i))
  load(sprintf('miss_0.4/cross_val/miss_0.4_bi_las_sp_rep%sK2_data.block.Robj',i))
  tmp2.order <- match.layers.x.runs(tmp2, conStruct.results[[1]]$MAP$admix.proportions)
  layers_m4_k2[,i] <- c(calculate.layer.contribution(conStruct.results = conStruct.results[[1]],
                                                     data.block = data.block,
                                                     layer.order = tmp2.order))
}

barplot(layers_m4_k2,
        col=c('blue','red'),
        ylab='layer contributions')




##### xval comparison missing 0.6 #####

sp_x_val_results <- read.table('miss_0.6/cross_val/miss_0.6_bi_las_sp_xval_results.txt', header=T)
nsp_x_val_results <- read.table('miss_0.6/cross_val/miss_0.6_bi_las_nsp_xval_results.txt', header=T)
sp.CIs <-  apply(sp_x_val_results,1,function(x){mean(x) + c(-1.96,1.96) * sd(x)/length(x)})
nsp.CIs <- apply(nsp_x_val_results,1,function(x){mean(x) + c(-1.96,1.96) * sd(x)/length(x)})
plot(rowMeans(nsp_x_val_results),
     pch=19,col="blue",
     ylab="predictive accuracy",xlab="values of K",
     ylim=range(nsp.CIs,sp.CIs),
     main="cross-validation results")
points(rowMeans(sp_x_val_results),col="green",pch=19)

segments(x0 = 1:nrow(nsp_x_val_results),
         y0 = nsp.CIs[1,],
         x1 = 1:nrow(nsp_x_val_results),
         y1 = nsp.CIs[2,],
         col = "blue",lwd=2)


load('miss_0.6/cross_val/miss_0.6_bi_las_nsp_rep1K3_conStruct.results.Robj')
load('miss_0.6/cross_val/miss_0.6_bi_las_sp_rep1K3_data.block.Robj')

meh <- c(calculate.layer.contribution(conStruct.results[[1]], data.block))
tmp <- conStruct.results[[1]]$MAP$admix.proportions
m6_k3_8rep_layers <- matrix(NA, nrow=3, ncol=8)
m6_k3_8rep_layers[,1] <- meh
for(i in 2:8){
  
  
  load(sprintf('miss_0.6/cross_val/miss_0.6_bi_las_nsp_rep%sK3_conStruct.results.Robj',i))
  load(sprintf('miss_0.6/cross_val/miss_0.6_bi_las_sp_rep%sK3_data.block.Robj',i))
  tmp.order <- match.layers.x.runs(tmp, conStruct.results[[1]]$MAP$admix.proportions)
  m6_k3_8rep_layers[,i] <- c(calculate.layer.contribution(conStruct.results = conStruct.results[[1]],
                                                          data.block = data.block,
                                                          layer.order = tmp.order))
  
}

m6_k3_8rep_layers

barplot(m6_k3_8rep_layers,
        col=c('blue','red','goldenrod1'),
        ylab='layer contributions')

layers_m6_k2 <- matrix(NA, nrow=2, ncol=8)

load('miss_0.6/cross_val/miss_0.6_bi_las_nsp_rep1K2_conStruct.results.Robj')
load('miss_0.6/cross_val/miss_0.6_bi_las_sp_rep1K2_data.block.Robj')
layers_m6_k2[,1] <- c(calculate.layer.contribution(conStruct.results[[1]], data.block))
tmp2 <- conStruct.results[[1]]$MAP$admix.proportions

for(i in 2:8){
  load(sprintf('miss_0.6/cross_val/miss_0.6_bi_las_nsp_rep%sK2_conStruct.results.Robj',i))
  load(sprintf('miss_0.6/cross_val/miss_0.6_bi_las_sp_rep%sK2_data.block.Robj',i))
  tmp2.order <- match.layers.x.runs(tmp2, conStruct.results[[1]]$MAP$admix.proportions)
  layers_m6_k2[,i] <- c(calculate.layer.contribution(conStruct.results = conStruct.results[[1]],
                                                     data.block = data.block,
                                                     layer.order = tmp2.order))
}

barplot(layers_m6_k2,
        col=c('blue','red'),
        ylab='layer contributions')





###### x val comparions missing 0.8 ######

sp_x_val_results <- read.table('miss_0.8/cross_val/miss_0.8_bi_las_sp_xval_results.txt', header=T)
nsp_x_val_results <- read.table('miss_0.8/cross_val/miss_0.8_bi_las_nsp_xval_results.txt', header=T)
sp.CIs <-  apply(sp_x_val_results,1,function(x){mean(x) + c(-1.96,1.96) * sd(x)/length(x)})
nsp.CIs <- apply(nsp_x_val_results,1,function(x){mean(x) + c(-1.96,1.96) * sd(x)/length(x)})
plot(rowMeans(nsp_x_val_results),
     pch=19,col="blue",
     ylab="predictive accuracy",xlab="values of K",
     ylim=range(nsp.CIs,sp.CIs),
     main="cross-validation results")
points(rowMeans(sp_x_val_results),col="green",pch=19)

segments(x0 = 1:nrow(nsp_x_val_results),
         y0 = nsp.CIs[1,],
         x1 = 1:nrow(nsp_x_val_results),
         y1 = nsp.CIs[2,],
         col = "blue",lwd=2)





###### admix plots m = 0.8 k = 3 #######



load('miss_0.8/cross_val/miss_0.8_bi_las_sp_rep3K3_conStruct.results.Robj')
load('miss_0.8/cross_val/miss_0.8_bi_las_sp_rep3K3_data.block.Robj')
load('miss_0.8/cross_val/miss_0.8_bi_las_sp_rep3K3_model.fit.Robj')
m_8_k_3_results <- conStruct.results
m_8_k_3_data_block <- data.block
m_8_k_3_mod_fit <- model.fit



load('miss_0.6/cross_val/miss_0.6_bi_las_sp_rep1K3_conStruct.results.Robj')
load('miss_0.6/cross_val/miss_0.6_bi_las_sp_rep1K3_data.block.Robj')
load('miss_0.6/cross_val/miss_0.6_bi_las_sp_rep1K3_model.fit.Robj')
m_6_k_3_results <- conStruct.results
m_6_k_3_data_block <- data.block
m_6_k_3_mod_fit <- model.fit


load('miss_0.4/cross_val/miss_0.4_bi_las_sp_rep1K3_conStruct.results.Robj')
load('miss_0.4/cross_val/miss_0.4_bi_las_sp_rep1K3_data.block.Robj')
load('miss_0.4/cross_val/miss_0.4_bi_las_sp_rep1K3_model.fit.Robj')
m_4_k_3_results <- conStruct.results
m_4_k_3_data_block <- data.block
m_4_k_3_mod_fit <- model.fit
m_4_k_3_data_block$N
m_6_k_3_data_block$N
m_1_k_3_data_block$N
m_8_k_3_data_block$N
load('miss_1/cross_val/miss_1_bi_las_sp_rep1K3_conStruct.results.Robj')
load('miss_1/cross_val/miss_1_bi_las_sp_rep1K3_data.block.Robj')
load('miss_1/cross_val/miss_1_bi_las_sp_rep1K3_model.fit.Robj')
m_1_k_3_results <- conStruct.results
m_1_k_3_data_block <- data.block
m_1_k_3_mod_fit <- model.fit

ind_counts <- read.table('../../both_lanes_high_cov_individuals.txt', header = F)
inds <- ind_counts$V1
m43_inds <- conStruct::match.layers.x.runs(m_1_k_3_results$chain_1$MAP$admix.proportions,
                                           m_4_k_3_results$chain_1$MAP$admix.proportions)
m63_inds <- conStruct::match.layers.x.runs(m_1_k_3_results$chain_1$MAP$admix.proportions,
                                           m_6_k_3_results$chain_1$MAP$admix.proportions)
m83_inds <- conStruct::match.layers.x.runs(m_1_k_3_results$chain_1$MAP$admix.proportions,
                                           m_8_k_3_results$chain_1$MAP$admix.proportions)
admi_13 <- as.data.frame(m_1_k_3_results$chain_1$MAP$admix.proportions)
admi_43 <- as.data.frame(m_4_k_3_results$chain_1$MAP$admix.proportions[,c(m43_inds)])
admi_63 <- as.data.frame(m_6_k_3_results$chain_1$MAP$admix.proportions[,c(m63_inds)])
admi_83 <-  as.data.frame(m_8_k_3_results$chain_1$MAP$admix.proportions[,c(m83_inds)])
admi_13$V4 <- inds
admi_43$V4 <- inds
admi_63$V4 <- inds
admi_83$V4 <- inds

names(admi_13) <- c('pop1', 'pop2','pop3', 'ind' )
names(admi_43) <- c('pop1', 'pop2','pop3', 'ind' )
names(admi_63) <- c('pop1', 'pop2', 'pop3','ind' )
names(admi_83) <- c('pop1', 'pop2', 'pop3','ind' )

library(tidyr)
admi_long_13 <- admi_13 %>%
  pivot_longer(cols=c(pop1, pop2, pop3), names_to='pop', values_to='prop')
admi_long_43 <- admi_43 %>%
  pivot_longer(cols=c(pop1, pop2, pop3), names_to='pop', values_to='prop')
admi_long_63 <- admi_63 %>%
  pivot_longer(cols=c(pop1, pop2, pop3), names_to='pop', values_to='prop')
admi_long_83 <- admi_83 %>%
  pivot_longer(cols=c(pop1, pop2, pop3), names_to='pop', values_to='prop')
library(ggplot2)
library(cowplot)
ad_13 <- ggplot(data = admi_long_13, aes(y=prop,x=ind, fill=pop))+
  ylab('Admixture proportion')+xlab('Individual')+
  geom_bar(stat='identity')+
  theme(axis.text.x = element_text(angle=90))+ggtitle('No Missing Data')
ad_43 <- ggplot(data = admi_long_43, aes(y=prop,x=ind, fill=pop))+
  ylab('Admixture proportion')+xlab('Individual')+
  geom_bar(stat='identity')+
  theme(axis.text.x = element_text(angle=90))+ggtitle('60% Missing Allowed')
ad_63 <- ggplot(data = admi_long_63, aes(y=prop,x=ind, fill=pop))+
  ylab('Admixture proportion')+xlab('Individual')+
  geom_bar(stat='identity')+
  theme(axis.text.x = element_text(angle=90))+ggtitle('40% Missing Allowed')
ad_83 <- ggplot(data = admi_long_83, aes(y=prop,x=ind, fill=pop))+
  ylab('Admixture proportion')+xlab('Individual')+
  geom_bar(stat='identity')+
  theme(axis.text.x = element_text(angle=90))+ggtitle('20% Missing Allowed')

plot_grid(ad_83, ad_13)



load('miss_0.8/cross_val/miss_0.8_bi_las_nsp_rep1K3_conStruct.results.Robj')
load('miss_0.8/cross_val/miss_0.8_bi_las_sp_rep1K3_data.block.Robj')

meh <- c(calculate.layer.contribution(conStruct.results[[1]], data.block))
tmp <- conStruct.results[[1]]$MAP$admix.proportions
m8_k3_8rep_layers <- matrix(NA, nrow=3, ncol=8)
m8_k3_8rep_layers[,1] <- meh
for(i in 2:8){
  
  
  load(sprintf('miss_0.8/cross_val/miss_0.8_bi_las_nsp_rep%sK3_conStruct.results.Robj',i))
  load(sprintf('miss_0.8/cross_val/miss_0.8_bi_las_sp_rep%sK3_data.block.Robj',i))
  tmp.order <- match.layers.x.runs(tmp, conStruct.results[[1]]$MAP$admix.proportions)
  m8_k3_8rep_layers[,i] <- c(calculate.layer.contribution(conStruct.results = conStruct.results[[1]],
                                                          data.block = data.block,
                                                          layer.order = tmp.order))
  
}

m8_k3_8rep_layers

barplot(m8_k3_8rep_layers,
        col=c('blue','red','goldenrod1'),
        ylab='layer contributions')

layers_m8_k2 <- matrix(NA, nrow=2, ncol=8)

load('miss_0.8/cross_val/miss_0.8_bi_las_nsp_rep1K2_conStruct.results.Robj')
load('miss_0.8/cross_val/miss_0.8_bi_las_sp_rep1K2_data.block.Robj')
layers_m8_k2[,1] <- c(calculate.layer.contribution(conStruct.results[[1]], data.block))
tmp2 <- conStruct.results[[1]]$MAP$admix.proportions

for(i in 2:8){
  load(sprintf('miss_0.8/cross_val/miss_0.8_bi_las_nsp_rep%sK2_conStruct.results.Robj',i))
  load(sprintf('miss_0.8/cross_val/miss_0.8_bi_las_sp_rep%sK2_data.block.Robj',i))
  tmp2.order <- match.layers.x.runs(tmp2, conStruct.results[[1]]$MAP$admix.proportions)
  layers_m8_k2[,i] <- c(calculate.layer.contribution(conStruct.results = conStruct.results[[1]],
                                                     data.block = data.block,
                                                     layer.order = tmp2.order))
}

barplot(layers_m8_k2,
        col=c('blue','red'),
        ylab='layer contributions')

k3_m82_inds <- conStruct::match.layers.x.runs(m_8_k_3_results$chain_1$MAP$admix.proportions,
                                              m_8_k_3_results$chain_2$MAP$admix.proportions)
k3_m83_inds <- conStruct::match.layers.x.runs(m_8_k_3_results$chain_1$MAP$admix.proportions,
                                              m_8_k_3_results$chain_3$MAP$admix.proportions)
k3_m84_inds <- conStruct::match.layers.x.runs(m_8_k_3_results$chain_1$MAP$admix.proportions,
                                              m_8_k_3_results$chain_4$MAP$admix.proportions)

k3_admi_c81 <- as.data.frame(m_8_k_3_results$chain_1$MAP$admix.proportions)
k3_admi_c82 <- as.data.frame(m_8_k_3_results$chain_2$MAP$admix.proportions[,c(k3_m82_inds)])
k3_admi_c83 <- as.data.frame(m_8_k_3_results$chain_3$MAP$admix.proportions[,c(k3_m83_inds)])
k3_admi_c84 <-  as.data.frame(m_8_k_3_results$chain_4$MAP$admix.proportions[,c(k3_m84_inds)])
k3_admi_c81$V4 <- inds
k3_admi_c82$V4 <- inds
k3_admi_c83$V4 <- inds
k3_admi_c84$V4 <- inds

names(k3_admi_c81) <- c('pop1', 'pop2','pop3', 'ind' )
names(k3_admi_c82) <- c('pop1', 'pop2','pop3', 'ind' )
names(k3_admi_c83) <- c('pop1', 'pop2','pop3', 'ind' )
names(k3_admi_c84) <- c('pop1', 'pop2','pop3', 'ind' )

library(tidyr)
k3_admi_long_c81 <- k3_admi_c81 %>%
  pivot_longer(cols=c(pop1, pop2,pop3), names_to='pop', values_to='prop')
k3_admi_long_c82 <- k3_admi_c82 %>%
  pivot_longer(cols=c(pop1, pop2,pop3), names_to='pop', values_to='prop')
k3_admi_long_c83 <- k3_admi_c83 %>%
  pivot_longer(cols=c(pop1, pop2,pop3), names_to='pop', values_to='prop')
k3_admi_long_c84 <- k3_admi_c84 %>%
  pivot_longer(cols=c(pop1, pop2,pop3), names_to='pop', values_to='prop')

k3_ad_c81 <- ggplot(data = k3_admi_long_c81, aes(y=prop,x=ind, fill=pop))+
  geom_bar(stat='identity')+
  theme(axis.text.x = element_text(angle=90))
k3_ad_c82 <- ggplot(data = k3_admi_long_c82, aes(y=prop,x=ind, fill=pop))+
  geom_bar(stat='identity')+
  theme(axis.text.x = element_text(angle=90))
k3_ad_c83 <- ggplot(data = k3_admi_long_c83, aes(y=prop,x=ind, fill=pop))+
  geom_bar(stat='identity')+
  theme(axis.text.x = element_text(angle=90))
k3_ad_c84 <- ggplot(data = k3_admi_long_c84, aes(y=prop,x=ind, fill=pop))+
  geom_bar(stat='identity')+
  theme(axis.text.x = element_text(angle=90))

plot_grid(k3_ad_c81, k3_ad_c82,k3_ad_c83, k3_ad_c84)
library(Polychrome)
library(ggplot2)
library(cowplot)



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

library(ggplot2)
library(tidyr)
library(dplyr)

admix_over_reps <- function(path, k, n_reps, ind_file){
  block_pattern <- paste0('*_sp*',k,'*block.Robj')
  results_pattern <- paste0('*_sp*',k,'*results.Robj')
  fit_pattern <- paste0('*_sp*',k,'*fit.Robj')
  missing_pattern <- glob2rx('*dropped_individuals*txt')

  ind_file <- read.table(ind_file, header = F)
  missing_files <- list.files(path=path, pattern = missing_pattern)
  missing_ind_file <- read.table(paste0(path,missing_files))
  inds <- ind_file$V1
  miss <- which(inds %in% missing_ind_file$x)
  inds <- inds[-miss]
  k <- as.integer(strsplit(k, '')[[1]][2])
  block_fnames <- list.files(path=path, pattern=glob2rx(block_pattern))
  res_fnames <- list.files(path=path, pattern=glob2rx(results_pattern))
  fit_fnames <- list.files(path=path, pattern=glob2rx(fit_pattern))
  
  load(paste0(path,res_fnames[1]))
  #print(res_fnames[1])
  ad_prop_1 <- as.data.frame(conStruct.results[[1]]$MAP$admix.proportions)
  ad_ref <- conStruct.results[[1]]$MAP$admix.proportions
  new_names <- rep('pop',k)
  new_nam <- paste0(new_names, c(1:k))
  names(ad_prop_1) <- new_nam
  ad_prop_1$inds <- inds
  k1 <- k+1
  #print(new_nam[1:k])
  #print(ad_prop_1)
  ad_prop_1_long <- ad_prop_1 %>%
    pivot_longer(cols=c(new_nam[1:k]), names_to='pop', values_to='prop')
  #print(ad_prop_1_long)
  ad_1_plot <- ggplot(data = ad_prop_1_long, aes(y=prop,x=inds, fill=pop))+
     geom_bar(stat='identity')+
     theme(axis.text.x = element_text(angle=90))+
     scale_x_discrete()
   ad_1_plot
   pl_list <- list()
   pl_list[[1]] <- ad_1_plot
   for(i in 2:n_reps){
    
     load(paste0(path,res_fnames[i]))
     
     match_order <- conStruct::match.layers.x.runs(ad_ref,
                                                   conStruct.results[[1]]$MAP$admix.proportions)
     ad_prop <- as.data.frame(conStruct.results[[1]]$MAP$admix.proportions[,c(match_order)])
     names(ad_prop) <- new_nam
     ad_prop$inds <- inds
     ad_prop <- ad_prop %>%
       pivot_longer(cols=c(new_nam[1:k]), names_to='pop', values_to='prop')
     #print(ad_prop_1_long)
     ad_plot <- ggplot(data = ad_prop, aes(y=prop,x=inds, fill=pop))+
       geom_bar(stat='identity')+
       theme(axis.text.x = element_text(angle=90))+
       scale_x_discrete()
      pl_list[[i]] <- ad_plot
   }
   plot_grid(plotlist=pl_list)
   return(pl_list)
  # 
}


####### Unresolved ########
setwd('~/Desktop/Lasthenia/for_construct/output/')
test_list <- admix_over_reps('miss_0.6/cross_val/', 'K3',8, '../both_lanes_high_cov_individuals.txt')

plot_grid(plotlist = test_list)


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
  k1 <- k+1
  #print(new_nam[1:k])
  #print(ad_prop_1)
  ad_prop_1_long <- ad_prop_1 %>%
    pivot_longer(cols=c(new_nam[1:k]), names_to='pop', values_to='prop')
  #print(ad_prop_1_long)
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

huh <- admix_over_missing('K3', '../both_lanes_high_cov_individuals.txt')
huh
plot_grid(plotlist=huh)

load("miss_0.4/cross_val/miss_0.4_bi_las_sp_rep1K2_conStruct.results.Robj")

ad <- conStruct.results[[1]]$MAP$admix.proportions
ads <- ad[]

conStruct.results$chain_1$MAP$index.iter
setwd('~/Desktop/Lasthenia/for_construct/')
wel <- list.files(, pattern = block_pattern, recursive = T, include.dirs = T)
grep('miss_0.6', wel)
wel
plot_grid(plotlist=mget(paste0("test_plot", 1:8)))

plot_grid(test_list)
