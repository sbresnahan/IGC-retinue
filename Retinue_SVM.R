#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly=TRUE)
if(identical(args, character(0))){args = c("TEMP.RData",600,"NA",getwd())}
print(args)

library(e1071)
library(DESeq2)
library(limma)
library(pracma)
library(sva)
library(tidyverse)
library(doParallel)

setwd(args[[4]])
cores = detectCores()
print(paste0(cores," cores detected"))
registerDoParallel(cores=cores)

if(!file.exists(args[[1]])){

  ##### Load data
  counts = read.csv("datExpr.csv",row.names=1)
  metadata = read.csv("metadata.csv")
  counts = as.matrix(counts[metadata$ID])
  
  ##### Correct counts for batch effect of study
  counts_adj <- ComBat_seq(counts, batch=metadata$block, group=metadata$Role)
  
  
  #### Perform initial classification
  ##### Divide data into training (control) set and test (treatment) set
  treatment = metadata[sample(nrow(metadata), 6), "ID"]
  metadata[metadata$ID%in%treatment,"Role"] =
    paste(metadata[metadata$ID%in%treatment,"Role"],"test",sep="_")
  svm.data = metadata[,c("ID","Role")]
  svm.data.train = subset(svm.data, Role %in% c("responsive","unresponsive"))
  svm.data.test = subset(svm.data, !(Role %in% c("responsive","unresponsive")))
  
  ###### best model: radial kernel, gamma = 0.01, cost = 0.5 | best performance: 0.2638667
  
  
  #### Define SVM function
  svm.train = function(readcounts, traindata, testdata = NA, referencelevel = "queen", kerneltype = "radial", crossfold = 5, vstCheck = T){
  
    svm.counts.test=NA
  
    # normalise data
    svm.counts = readcounts
    # perform DESeq's variance stabilizing tranformation, which is preferable to logging for gene expression data
    if(vstCheck){svm.counts.vst = vst(svm.counts)}else{svm.counts.vst = varianceStabilizingTransformation(svm.counts)}
    # scale counts and remove zero-variance features
    svm.counts.vst.quantiles.scale = t(scale(t(svm.counts.vst)))
    svm.counts.vst.quantiles.scale = na.omit(svm.counts.vst.quantiles.scale)
  
    # Divide transcriptomic data into training set (queens and workers from control) and test set (individuals from treatment)
    svm.counts.train = svm.counts.vst.quantiles.scale[,which(colnames(svm.counts.vst.quantiles.scale) %in% traindata$ID)]
    if(length(testdata)>1){
      svm.counts.test = svm.counts.vst.quantiles.scale[,which((colnames(svm.counts.vst.quantiles.scale) %in% testdata$ID))]
    }
  
    # Perform a grid search to optimise SVM parameters
    svm.counts.tuneResult = tune("svm",
                                 train.x = t(svm.counts.train),
                                 train.y = as.numeric(traindata$Role == referencelevel),
                                 probability = TRUE,
                                 scale = FALSE,
                                 kernel = kerneltype,
                                 tunecontrol = tune.control(sampling = "cross",
                                                            cross = crossfold),
                                 ranges = list(gamma = 10^(-3:-1),
                                               cost = 2^(-2:0))
    )
  
    # Final classifier
    svm.counts.classifier = svm.counts.tuneResult$best.model
  
    svm.counts.prediction = NULL
    if(length(testdata)>1){
      # Make predictions for the test data, if test data were provided.
      svm.counts.prediction = predict(svm.counts.classifier,
                                      t(svm.counts.test),
                                      type = "class",
                                      probability = TRUE)
    }
  
    #output prediction for test data and cross-validation error for training data
    svm.result = list("prediction" = svm.counts.prediction,
                      "validation_error" = signif(svm.counts.tuneResult$best.performance,4),
                      "traincounts" = svm.counts.train,
                      "testcounts" = svm.counts.test)
  
    #return results
    return(svm.result)
  }
  
  # apply svm to entire set of genes
  svm.full = svm.train(counts_adj,
                       svm.data.train,
                       svm.data.test,
                       crossfold = 3,
                       vstCheck = F,
                       referencelevel = "responsive")
  print(paste0("Root mean cross-validation error rate for full model: ",svm.full$validation_error))
  
  
  ### Perform feature selection
  # create copy of training data that we can subject to repeated trimming while preserving original frame
  svm.counts.train.iterate = svm.full$traincounts
  #record original number of features
  nfeatures = nrow(svm.counts.train.iterate)
  #target number of features
  nfeatures_target = 100
  traindata = svm.data.train
  #instantiate data frame to hold data on the error of each model
  iterations = data.frame(feature = character(),
                         error_before_removal = numeric())
}

if(file.exists(args[[1]])){
  load(args[[1]])
  print(paste0("Restarting from feature # ", nfeatures))
}

time.start = strptime(Sys.time(), "%Y-%m-%d %H:%M:%S")
print(paste0("Start time: ",time.start))

#iteratively remove features until target number is reached
while(nfeatures > nfeatures_target){
  #Check if elapsed time is less than time limit in seconds (set at beginning of script as args[[1]])
  #If true, and if currently running on the cluster (providing submit script path in args[[3]]),
  #Restart R script from submit script and continue
  time.diff = as.numeric(difftime(strptime(Sys.time(),"%Y-%m-%d %H:%M:%S"),time.start,units="sec"))
  if(time.diff > as.numeric(args[[2]])){
    save.image(args[[1]])
    if(!args[[3]]=="NA"){
      system(paste0("qsub ",args[[3]]))
      quit()
    }
  }
  print(paste(nfeatures,Sys.time(),sep=" "))
  #run repeatedly to account for stochasticity in cross-validation
  tunelist = c()
  tunelist = foreach(i=1:10) %dopar% {
    #Perform a grid search to optimise SVM parameters
    svm.counts.tuneResult = tune("svm",
                                 train.x = t(svm.counts.train.iterate),
                                 train.y =  as.numeric(traindata$Role == "responsive"),
                                 probability = TRUE,
                                 scale = FALSE,
                                 kernel = "radial",
                                 tunecontrol = tune.control(sampling = "cross",
                                                            cross = 3),
                                 ranges = list(gamma = 10^(-3:-1),
                                               cost = 2^(-2:0)))
    #record error
    error = svm.counts.tuneResult$best.performance
    iterlist = c(svm.counts.tuneResult,error=error)
    return(iterlist)
  }
  #sample classifier
  svm.counts.classifier = tunelist[length(tunelist)][[1]][["best.model"]]
  #return mean error value
  error = signif(mean(sapply(tunelist, get, x="error")),4)
  #extract feature weights
  weights = (t(svm.counts.classifier$coefs) %*% svm.counts.classifier$SV)
  #calculate feature with lowest weight (for ties, choose arbitrarily)
  weakfeature = colnames(weights)[which(abs(weights) == min(abs(weights)))[1]]
  #remove lowest-weight feature from data frame
  svm.counts.train.iterate = subset(svm.counts.train.iterate, !(rownames(svm.counts.train.iterate) %in% c(weakfeature)))
  #in a dataframe, store removed feature name and error value before removing that feature
  iterations = rbind(iterations, tibble(feature = weakfeature,
                                        error_before_removal = error))
  #tick down
  nfeatures = (nfeatures-1)
  #output every 100 runs to track progress
  if((nfeatures/100)%%1==0){print(paste0("Features remaining: ",nfeatures))}
}
iterLength = 1:nrow(iterations)
# take moving average to smooth out variation
moving_avg = movavg(iterations$error_before_removal, 100, "s")

save.image("SVM.RData")
load("SVM.RData")

# plot data to ensure we have the expected 'hockeystick' shape 
# note that this will be truncated for the subsetted dataset provided for demoing!
hockeyData = data.frame(num = iterLength, error = moving_avg)
hockeyData_plot = hockeyData
hockeyData_plot$num = abs(iterLength - (max(iterLength)+1))
hockeyData_plot$error <- sqrt(hockeyData_plot$error)

library(ggprism)
g <- ggplot(data=hockeyData_plot, aes(x=num, y=error, group=1)) +
  geom_line(linewidth=1.5) + theme_prism() +
  scale_y_continuous(breaks=seq.int(.15,.55,.1),limits=c(.15,.55)) +
  scale_x_reverse(breaks=c(12000,10000,8000,6000,4000,2000,0)) +
  xlab("Genes") + ylab("RMSE") +
  geom_vline(xintercept=100, linetype="dashed", color = "red") +
  geom_hline(yintercept=0.2093285, linetype="dashed", color = "red")
g
ggsave("figS6.png",g,width=7,height=5,dpi=300)
# get minimum of this curve to find the point at which the error window is at its minimum
optimal_removal = which(moving_avg == min(moving_avg));
# list the features to be removed from the original set of genes
features_to_remove = iterations$feature[1:optimal_removal]
# new dataframe with less-useful features removed
counts_clean_subsample = subset(counts_adj,!(rownames(counts_adj) %in% features_to_remove))
# re-perform support vector classification using the new, optimally caste-separating set of features
svm.optimal = svm.train(counts_clean_subsample, 
                        referencelevel = "responsive",
                        svm.data.train, 
                        svm.data.test,
                        crossfold = 3,
                        vstCheck = F)

print(paste0("Number of genes included in optimised model: ", nrow(counts_clean_subsample)))
print(paste0("Root mean cross-validation error rate for optimised model: ", sqrt(svm.optimal$validation_error)))
write.csv(data.frame(GeneID=row.names(counts_clean_subsample)),
          "SVM_GeneIDs.csv",row.names=F)


fit.full <- list()
fit.optimal <- list()
metadata = read.csv("metadata.csv")
for(i in 1:10){
  treatment = metadata[sample(nrow(metadata), 6), "ID"]
  metadata0 <- metadata
  metadata0[metadata0$ID%in%treatment,"Role"] =
    paste(metadata0[metadata0$ID%in%treatment,"Role"],"test",sep="_")
  svm.data.boot = metadata0[,c("ID","Role")]
  svm.data.train.boot = subset(svm.data.boot, Role %in% c("responsive","unresponsive"))
  svm.data.test.boot = subset(svm.data.boot, !(Role %in% c("responsive","unresponsive")))
  
  fit.full[[i]] <- svm.train(counts_adj,
                       svm.data.train.boot,
                       svm.data.test.boot,
                       crossfold = 3,
                       vstCheck = F,
                       referencelevel = "responsive")[["prediction"]]
  
  fit.optimal[[i]] <- svm.train(counts_clean_subsample, 
                          referencelevel = "responsive",
                          svm.data.train.boot, 
                          svm.data.test.boot,
                          crossfold = 3,
                          vstCheck = F)[["prediction"]]
}






dat.full.fit <- data.frame(value=as.numeric(unlist(fit.full)),
                           ID=names(unlist(fit.full)))
dat.full.fit <- left_join(dat.full.fit,metadata,"ID")
dat.full.fit <- dat.full.fit[,1:3]
dat.full.fit$fit <- "full"


dat.optimal.fit <- data.frame(value=as.numeric(unlist(fit.optimal)),
                              ID=names(unlist(fit.optimal)))
dat.optimal.fit <- left_join(dat.optimal.fit,metadata,"ID")
dat.optimal.fit <- dat.optimal.fit[,1:3]
dat.optimal.fit$fit <- "optimal"

dat.fit <- rbind(dat.optimal.fit,dat.full.fit)
dat.fit$Role <- factor(dat.fit$Role,
                       levels=c("unresponsive","responsive"))
dat.fit$fit <- factor(dat.fit$fit,
                      levels=c("full","optimal"))
levels(dat.fit$fit) <- c("Full","Optimal")

g <- ggplot(dat.fit, aes(x=fit, y=value, color=Role,group=Role)) + 
  geom_point(position = position_jitterdodge(jitter.width = .05, 
                                             jitter.height = .005,
                                             dodge.width=.75),size=3) +
  scale_color_manual(values=alpha(c("blue", "red"), 0.2),
                     guide = guide_legend(override.aes = list(alpha=1))) +
  ylab("Classifier Estimate") + xlab("Model") +
  geom_point(stat="summary",size=3,
             position=position_dodge(0.75),color="black") +
  geom_errorbar(stat="summary",position=position_dodge(0.75),width=.25,
                color="black") +
  scale_y_continuous(breaks=c(0,.25,.5,.75,1),limits=c(0,1)) +
  theme_prism()
g
ggsave("figS4.png",g,width=6,height=5,dpi=300)

sessionInfo()