
#########################################################
##                                                     ##
##   RescueTag: we score sample multiplexing!          ##
##                                                     ##
#########################################################

ScoringTags<-function(path_to_table, x, y, first_tag_number, last_tag_number, scAD){
  sample_tag_table<-read.csv(path_to_table, row.names = 1)
  sample_tag_table<-sample_tag_table[,x:y]
  nCount_Tag<-rowSums(sample_tag_table)
  scAD[["nCount_Tag"]]<-nCount_Tag
  tags<-colnames(sample_tag_table)
  for(i in 1:length(tags)){
    assign(paste(tags[i]), sample_tag_table[,i])
  }
  for(i in 1:length(tags)){
    obj<-get(paste(tags[i]))
    scAD[[paste(tags[i])]]<-obj
  }
  df<-data.frame(matrix(,nrow=nrow(sample_tag_table),ncol=0))
  for(i in 1:length(tags)){
     df<-cbind(df,sample_tag_table[,i])
  }
  tag_names<-paste("Tag", first_tag_number:last_tag_number, sep="_")
  names(df)<-tag_names
  tr_df<-as.data.frame(t(df))
  first_best<-c()
  second_best<-c()
  delta<-c()
  for (i in 1:length(tr_df[1,])){
    x<-tr_df[,i]
    x<-x[order(x, decreasing = TRUE)]
    first_best<-c(first_best, x[1])
    second_best<-c(second_best,x[2])
    delta<-c(delta, x[1]-x[2])
  }
  final_df<-data.frame('First_best'=first_best, 'Second_best'=second_best, 'Delta'=delta)
  scAD[['First_best_tag']]<-final_df$First_best
  scAD[['Delta']]<-final_df$Delta
  df$putative_tags<-colnames(df)[apply(df,1,which.max)]
  scAD[['Putative_tags']]<-df$putative_tags
  for(i in 1:length(tags)){
    scAD[[paste(tags[i])]]<-c()
  }
  return(scAD)
}


NormToDepth<-function(scAD){
  all_reads<-scAD$nCount_Tag+scAD$nCount_RNA
  scAD[['NTD_first_best']]<-scAD$First_best_tag/all_reads
  return(scAD)
}

PreQCFilter<-function(scAD){
  selected_c <- WhichCells(scAD, expression = nFeature_RNA > 200)
  selected_f <- rownames(scAD)[Matrix::rowSums(scAD) > 3]
  scAD <- subset(scAD, features = selected_f, cells = selected_c)
  return(scAD)
}

FilterLowSeqDepth<-function(scAD, distribution="default"){
  if (distribution == "default"){
    scAD<-subset(scAD, subset=NTD_first_best<0.5)
    return(scAD)
  }
  else if (distribution == "normal"){
    x<-scAD$NTD_first_best
    threshold_for_NTD<-round(mean(x)+2*sd(x),2)
    scAD<-subset(scAD, subset=NTD_first_best<threshold_for_NTD)
    return(scAD)
  }
  else if (distribution == "bimodal"){
    # The function was taken from https://rpubs.com/H_Zhu/246450
    x<-scAD$NTD_first_best
    mem <- kmeans(x,2)$cluster
    mu1 <- mean(x[mem==1])
    mu2 <- mean(x[mem==2])
    sigma1 <- sd(x[mem==1])
    sigma2 <- sd(x[mem==2])
    pi1 <- sum(mem==1)/length(mem)
    pi2 <- sum(mem==2)/length(mem)
    
    sum.finite <- function(x) {
      sum(x[is.finite(x)])
    }
    
    Q <- 0
    # starting value of expected value of the log likelihood
    Q[2] <- sum.finite(log(pi1)+log(dnorm(x, mu1, sigma1))) + sum.finite(log(pi2)+log(dnorm(x, mu2, sigma2)))
    
    k <- 2
    while (abs(Q[k]-Q[k-1])>=1e-6) {
      # E step
      comp1 <- pi1 * dnorm(x, mu1, sigma1)
      comp2 <- pi2 * dnorm(x, mu2, sigma2)
      comp.sum <- comp1 + comp2
      
      p1 <- comp1/comp.sum
      p2 <- comp2/comp.sum
      
      # M step
      pi1 <- sum.finite(p1) / length(x)
      pi2 <- sum.finite(p2) / length(x)
      
      mu1 <- sum.finite(p1 * x) / sum.finite(p1)
      mu2 <- sum.finite(p2 * x) / sum.finite(p2)
      
      sigma1 <- sqrt(sum.finite(p1 * (x-mu1)^2) / sum.finite(p1))
      sigma2 <- sqrt(sum.finite(p2 * (x-mu2)^2) / sum.finite(p2))
      
      p1 <- pi1 
      p2 <- pi2
      
      k <- k + 1
      Q[k] <- sum(log(comp.sum))
    }
    p<-c(p1,p2)
    p<-sort(p, decreasing = TRUE)
    mu<-c(mu1,mu2)
    mu<-sort(mu, decreasing = TRUE)
    sg<-c(sigma1.sigma2)
    sg<-sort(sg, decreasing = TRUE)
    gm<-normalmixEM(x,k=3,lambda=c(round(p[1],2),round(p[2],2)),mu=c(round(mu[1],2),round(mu[2],2)),sigma=c(round(sg[1],2),round(sg[2],2)))
    hist(x, prob=T, breaks=32, xlim=c(range(x)[1], range(x)[2]), main='')
    lines(density(x), col="green", lwd=2)
    x1 <- seq(from=range(x)[1], to=range(x)[2], length.out=1000)
    y <- pi1 * dnorm(x1, mean=mu1, sd=sigma1) + pi2 * dnorm(x1, mean=mu2, sd=sigma2)
    lines(x1, y, col="red", lwd=2)
    legend('topright', col=c("green", 'red'), lwd=2, legend=c("kernal", "fitted"))
    print(gm$mu)
    question<-(readline("Which of the above numbers is bigger (1 or 2)?"))
    question<-as.numeric(question)
    if(question==1){
      threshold_for_NTD<-round(gm$mu[1]-2*gm$sigma[1], 2)
      scAD<-subset(scAD, subset=NTD_first_best<threshold_for_NTD)
      return(scAD)
    }
    else if(question==2){
      threshold_for_NTD<-round(gm$mu[2]-2*gm$sigma[2], 2)
      scAD<-subset(scAD, subset=NTD_first_best<threshold_for_NTD)
      return(scAD)
    }
  }
}

NormTagQCParams<-function(scAD){
  scAD[['Postnorm_first_best']]<-scAD$First_best_tag/(scAD$nCount_Tag)
  scAD[['Postnorm_delta']]<-scAD$Delta/(scAD$nCount_Tag)
  return(scAD)
}

RidgeTagQC<-function(scAD){
  x<-c(scAD$Postnorm_first_best)
  x<-c(x, scAD$Postnorm_delta)
  y<-c(rep('First_best', length(scAD$Postnorm_first_best)))
  y<-c(y, rep('Delta', length(scAD$Postnorm_first_best)))
  df<-data.frame('Normalized_values'=x, 'Categories'=y)
  ggplot(df, aes(x=Normalized_values, y=Categories, fill=Categories))+geom_density_ridges(scale=1)+theme_ridges(font_size = 19)+theme(legend.position = "none",panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank())
}

FirstBestThreshold<-function(scAD, distribution){
  first_best<-scAD$Postnorm_first_best
  if(distribution == "normal"){
    threshold_first_best<-mean(first_best)-2*sd(first_best)
    return(threshold_first_best)
  }
  else if (distribution == "bimodal"){
    x<-scAD$Postnorm_first_best
    mem <- kmeans(x,2)$cluster
    mu1 <- mean(x[mem==1])
    mu2 <- mean(x[mem==2])
    sigma1 <- sd(x[mem==1])
    sigma2 <- sd(x[mem==2])
    pi1 <- sum(mem==1)/length(mem)
    pi2 <- sum(mem==2)/length(mem)
    
    sum.finite <- function(x) {
      sum(x[is.finite(x)])
    }
    
    Q <- 0
    # starting value of expected value of the log likelihood
    Q[2] <- sum.finite(log(pi1)+log(dnorm(x, mu1, sigma1))) + sum.finite(log(pi2)+log(dnorm(x, mu2, sigma2)))
    
    k <- 2
    
    while (abs(Q[k]-Q[k-1])>=1e-6) {
      # E step
      comp1 <- pi1 * dnorm(x, mu1, sigma1)
      comp2 <- pi2 * dnorm(x, mu2, sigma2)
      comp.sum <- comp1 + comp2
      
      p1 <- comp1/comp.sum
      p2 <- comp2/comp.sum
      
      # M step
      pi1 <- sum.finite(p1) / length(x)
      pi2 <- sum.finite(p2) / length(x)
      
      mu1 <- sum.finite(p1 * x) / sum.finite(p1)
      mu2 <- sum.finite(p2 * x) / sum.finite(p2)
      
      sigma1 <- sqrt(sum.finite(p1 * (x-mu1)^2) / sum.finite(p1))
      sigma2 <- sqrt(sum.finite(p2 * (x-mu2)^2) / sum.finite(p2))
      
      p1 <- pi1 
      p2 <- pi2
      
      k <- k + 1
      Q[k] <- sum(log(comp.sum))
    }
    p<-c(p1,p2)
    p<-sort(p, decreasing = TRUE)
    mu<-c(mu1,mu2)
    mu<-sort(mu, decreasing = TRUE)
    sg<-c(sigma1,sigma2)
    sg<-sort(sg, decreasing = TRUE)
    gm<-normalmixEM(x,k=2,lambda=c(round(p[1],2),round(p[2],2)),mu=c(round(mu[1],2),round(mu[2],2)),sigma=c(round(sg[1],2),round(sg[2],2)))
    
    hist(x, prob=T, breaks=32, xlim=c(range(x)[1], range(x)[2]), main='')
    lines(density(x), col="green", lwd=2)
    x1 <- seq(from=range(x)[1], to=range(x)[2], length.out=1000)
    y <- pi1 * dnorm(x1, mean=mu1, sd=sigma1) + pi2 * dnorm(x1, mean=mu2, sd=sigma2)
    lines(x1, y, col="red", lwd=2)
    legend('topright', col=c("green", 'red'), lwd=2, legend=c("kernal", "fitted"))
    print(gm$mu)
    question<-(readline("Which of the above numbers is lower (1 or 2)?"))
    question<-as.numeric(question)
    if(question==1){
      threshold_first_best<-round(gm$mu[1]+2*gm$sigma[1], 2)
      return(threshold_first_best)
    }
    else if(question==2){
      threshold_first_best<-round(gm$mu[2]+2*gm$sigma[2], 2)
      return(threshold_first_best)
    }
  }
  else{
    stop("invalid distribution for this function")
  }
}

FirstBestMultiMode<-function(scAD, number_of_peaks){
  if(number_of_peaks == 3){
    x<-scAD$Postnorm_first_best
    mem <- kmeans(x,3)$cluster
    mu1 <- mean(x[mem==1])
    mu2 <- mean(x[mem==2])
    mu3 <- mean(x[mem==3])
    sigma1 <- sd(x[mem==1])
    sigma2 <- sd(x[mem==2])
    sigma3 <- sd(x[mem==3])
    pi1 <- sum(mem==1)/length(mem)
    pi2 <- sum(mem==2)/length(mem)
    pi3 <- sum(mem==3)/length(mem)
    
    sum.finite <- function(x) {
      sum(x[is.finite(x)])
    }
    
    Q <- 0
    # starting value of expected value of the log likelihood
    Q[2] <- sum.finite(log(pi1)+log(dnorm(x, mu1, sigma1))) + sum.finite(log(pi2)+log(dnorm(x, mu2, sigma2))) + sum.finite(log(pi3)+log(dnorm(x, mu3, sigma3)))
    
    k <- 2
    
    while (abs(Q[k]-Q[k-1])>=1e-6) {
      # E step
      comp1 <- pi1 * dnorm(x, mu1, sigma1)
      comp2 <- pi2 * dnorm(x, mu2, sigma2)
      comp3 <- pi3 * dnorm(x, mu3, sigma3)
      comp.sum <- comp1 + comp2 + comp3
      
      p1 <- comp1/comp.sum
      p2 <- comp2/comp.sum
      p3 <- comp3/comp.sum
      
      # M step
      pi1 <- sum.finite(p1) / length(x)
      pi2 <- sum.finite(p2) / length(x)
      pi3 <- sum.finite(p3) / length(x)
      
      mu1 <- sum.finite(p1 * x) / sum.finite(p1)
      mu2 <- sum.finite(p2 * x) / sum.finite(p2)
      mu3 <- sum.finite(p3 * x) / sum.finite(p3)
      
      sigma1 <- sqrt(sum.finite(p1 * (x-mu1)^2) / sum.finite(p1))
      sigma2 <- sqrt(sum.finite(p2 * (x-mu2)^2) / sum.finite(p2))
      sigma3 <- sqrt(sum.finite(p3 * (x-mu3)^2) / sum.finite(p3))
      
      p1 <- pi1 
      p2 <- pi2
      p3 <- pi3
      
      k <- k + 1
      Q[k] <- sum(log(comp.sum))
    }
    p<-c(p1,p2,p3)
    p<-sort(p, decreasing = TRUE)
    mu<-c(mu1,mu2,mu3)
    mu<-sort(mu, decreasing = TRUE)
    sg<-c(sigma1,sigma2,sigma3)
    sg<-sort(sg, decreasing = TRUE)
    gm<-normalmixEM(x,k=3,lambda=c(round(p[1],2),round(p[2],2),round(p[3],3)),mu=c(round(mu[1],2),round(mu[2],2),round(mu[3],2)),sigma=c(round(sg[1],2),round(sg[2],2),round(sg[3],2)))
    
    hist(x, prob=T, breaks=32, xlim=c(range(x)[1], range(x)[2]), main='')
    lines(density(x), col="green", lwd=2)
    x1 <- seq(from=range(x)[1], to=range(x)[2], length.out=1000)
    y <- pi1 * dnorm(x1, mean=mu1, sd=sigma1) + pi2 * dnorm(x1, mean=mu2, sd=sigma2) + pi3 * dnorm(x1, mean=mu3, sd=sigma3)
    lines(x1, y, col="red", lwd=2)
    legend('topright', col=c("green", 'red'), lwd=2, legend=c("kernal", "fitted"))
    print(gm$mu)
    question<-(readline("Choose a number which lies on x-axis closest to the rightmost peak that will mark untagged cells (1,2 or 3)"))
    question<-as.numeric(question)
    if(question==1){
      threshold_first_best<-round(gm$mu[1]+2*gm$sigma[1], 2)
      return(threshold_first_best)
    }
    else if(question==2){
      threshold_first_best<-round(gm$mu[2]+2*gm$sigma[2], 2)
      return(threshold_first_best)
    }
    else if (question==3){
      threshold_first_best<-round(gm$mu[3]+2*gm$sigma[3], 2)
      return(threshold_first_best)
    }
    else{
      stop("Please choose a number from 1 to 3")
    }
  }
  else if (number_of_peaks == 4){
    x<-scAD$Postnorm_first_best
    mem <- kmeans(x,4)$cluster
    mu1 <- mean(x[mem==1])
    mu2 <- mean(x[mem==2])
    mu3 <- mean(x[mem==3])
    mu4 <- mean(x[mem==4])
    sigma1 <- sd(x[mem==1])
    sigma2 <- sd(x[mem==2])
    sigma3 <- sd(x[mem==3])
    sigma4 <- sd(x[mem==4])
    pi1 <- sum(mem==1)/length(mem)
    pi2 <- sum(mem==2)/length(mem)
    pi3 <- sum(mem==3)/length(mem)
    pi4 <- sum(mem==4)/length(mem)
    
    sum.finite <- function(x) {
      sum(x[is.finite(x)])
    }
    
    Q <- 0
    # starting value of expected value of the log likelihood
    Q[2] <- sum.finite(log(pi1)+log(dnorm(x, mu1, sigma1))) + sum.finite(log(pi2)+log(dnorm(x, mu2, sigma2))) + sum.finite(log(pi3)+log(dnorm(x, mu3, sigma3))) + sum.finite(log(pi4)+log(dnorm(x, mu4, sigma4)))
    
    k <- 2
    
    while (abs(Q[k]-Q[k-1])>=1e-6) {
      # E step
      comp1 <- pi1 * dnorm(x, mu1, sigma1)
      comp2 <- pi2 * dnorm(x, mu2, sigma2)
      comp3 <- pi3 * dnorm(x, mu3, sigma3)
      comp4 <- pi4 * dnorm(x, mu4, sigma4)
      comp.sum <- comp1 + comp2 + comp3 + comp4
      
      p1 <- comp1/comp.sum
      p2 <- comp2/comp.sum
      p3 <- comp3/comp.sum
      p4 <- comp4/comp.sum
      
      # M step
      pi1 <- sum.finite(p1) / length(x)
      pi2 <- sum.finite(p2) / length(x)
      pi3 <- sum.finite(p3) / length(x)
      pi4 <- sum.finite(p4) / length(x)
      
      mu1 <- sum.finite(p1 * x) / sum.finite(p1)
      mu2 <- sum.finite(p2 * x) / sum.finite(p2)
      mu3 <- sum.finite(p3 * x) / sum.finite(p3)
      mu4 <- sum.finite(p4 * x) / sum.finite(p4)
      
      sigma1 <- sqrt(sum.finite(p1 * (x-mu1)^2) / sum.finite(p1))
      sigma2 <- sqrt(sum.finite(p2 * (x-mu2)^2) / sum.finite(p2))
      sigma3 <- sqrt(sum.finite(p3 * (x-mu3)^2) / sum.finite(p3))
      sigma4 <- sqrt(sum.finite(p4 * (x-mu4)^2) / sum.finite(p4))
      
      p1 <- pi1 
      p2 <- pi2
      p3 <- pi3
      p4 <- pi4
      
      k <- k + 1
      Q[k] <- sum(log(comp.sum))
    }
    p<-c(p1,p2,p3, p4)
    p<-sort(p, decreasing = TRUE)
    mu<-c(mu1,mu2,mu3, mu4)
    mu<-sort(mu, decreasing = TRUE)
    sg<-c(sigma1.sigma2.sigma3, sigma4)
    sg<-sort(sg, decreasing = TRUE)
    gm<-normalmixEM(x,k=3,lambda=c(round(p[1],2),round(p[2],2),round(p[3],3), round(p[4],2)),mu=c(round(mu[1],2),round(mu[2],2),round(mu[3],2), round(mu[4],2)),sigma=c(round(sg[1],2),round(sg[2],2),round(sg[3],2), round(sg[4],2)))
    
    hist(x, prob=T, breaks=32, xlim=c(range(x)[1], range(x)[2]), main='')
    lines(density(x), col="green", lwd=2)
    x1 <- seq(from=range(x)[1], to=range(x)[2], length.out=1000)
    y <- pi1 * dnorm(x1, mean=mu1, sd=sigma1) + pi2 * dnorm(x1, mean=mu2, sd=sigma2) + pi3 * dnorm(x1, mean=mu3, sd=sigma3) + pi4 * dnorm(x1, mean=mu4, sd=sigma4)
    lines(x1, y, col="red", lwd=2)
    legend('topright', col=c("green", 'red'), lwd=2, legend=c("kernal", "fitted"))
    print(gm$mu)
    question<-(readline("Choose a number which lies on x-axis closest to the rightmost peak that will mark untagged cells (1,2 or 3)"))
    question<-as.numeric(question)
    if(question==1){
      threshold_first_best<-round(gm$mu[1]+2*gm$sigma[1], 2)
      return(threshold_first_best)
    }
    else if(question==2){
      threshold_first_best<-round(gm$mu[2]+2*gm$sigma[2], 2)
      return(threshold_first_best)
    }
    else if (question==3){
      threshold_first_best<-round(gm$mu[3]+2*gm$sigma[3], 2)
      return(threshold_first_best)
    }
    else if (question==4){
      threshold_first_best<-round(gm$mu[4]+2*gm$sigma[4], 2)
      return(threshold_first_best)
    }
    else{
      stop("Please choose a number from 1 to 4")
    }
  }
  else{
    stop("For bimodal distribution use FirstBestThreshold; otherwise the dataset is too noisy to select a confident threshold")
  }
}


DeltaThreshold<-function(scAD, distribution){
  delta<-scAD$Postnorm_delta
  if(distribution == "normal"){
    threshold_delta<-mean(delta)-2*sd(delta)
    return(threshold_delta)
  }
  else if (distribution == "bimodal"){
    x<-scAD$Postnorm_delta
    mem <- kmeans(x,2)$cluster
    mu1 <- mean(x[mem==1])
    mu2 <- mean(x[mem==2])
    sigma1 <- sd(x[mem==1])
    sigma2 <- sd(x[mem==2])
    pi1 <- sum(mem==1)/length(mem)
    pi2 <- sum(mem==2)/length(mem)
    
    sum.finite <- function(x) {
      sum(x[is.finite(x)])
    }
    
    Q <- 0
    # starting value of expected value of the log likelihood
    Q[2] <- sum.finite(log(pi1)+log(dnorm(x, mu1, sigma1))) + sum.finite(log(pi2)+log(dnorm(x, mu2, sigma2)))
    
    k <- 2
    
    while (abs(Q[k]-Q[k-1])>=1e-6) {
      # E step
      comp1 <- pi1 * dnorm(x, mu1, sigma1)
      comp2 <- pi2 * dnorm(x, mu2, sigma2)
      comp.sum <- comp1 + comp2
      
      p1 <- comp1/comp.sum
      p2 <- comp2/comp.sum
      
      # M step
      pi1 <- sum.finite(p1) / length(x)
      pi2 <- sum.finite(p2) / length(x)
      
      mu1 <- sum.finite(p1 * x) / sum.finite(p1)
      mu2 <- sum.finite(p2 * x) / sum.finite(p2)
      
      sigma1 <- sqrt(sum.finite(p1 * (x-mu1)^2) / sum.finite(p1))
      sigma2 <- sqrt(sum.finite(p2 * (x-mu2)^2) / sum.finite(p2))
      
      p1 <- pi1 
      p2 <- pi2
      
      k <- k + 1
      Q[k] <- sum(log(comp.sum))
    }
    p<-c(p1,p2)
    p<-sort(p, decreasing = TRUE)
    mu<-c(mu1,mu2)
    mu<-sort(mu, decreasing = TRUE)
    sg<-c(sigma1,sigma2)
    sg<-sort(sg, decreasing = TRUE)
    gm<-normalmixEM(x,k=2,lambda=c(round(p[1],2),round(p[2],2)),mu=c(round(mu[1],2),round(mu[2],2)),sigma=c(round(sg[1],2),round(sg[2],2)))
    
    hist(x, prob=T, breaks=32, xlim=c(range(x)[1], range(x)[2]), main='')
    lines(density(x), col="green", lwd=2)
    x1 <- seq(from=range(x)[1], to=range(x)[2], length.out=1000)
    y <- pi1 * dnorm(x1, mean=mu1, sd=sigma1) + pi2 * dnorm(x1, mean=mu2, sd=sigma2)
    lines(x1, y, col="red", lwd=2)
    legend('topright', col=c("green", 'red'), lwd=2, legend=c("kernal", "fitted"))
    print(gm$mu)
    question<-(readline("Which of the above numbers is lower (1 or 2)?"))
    question<-as.numeric(question)
    if(question==1){
      threshold_delta<-round(gm$mu[1]+2*gm$sigma[1], 2)
      return(threshold_delta)
    }
    else if(question==2){
      threshold_delta<-round(gm$mu[2]+2*gm$sigma[2], 2)
      return(threshold_delta)
    }
  }
  else{
    stop("invalid distribution for this function")
  }
}

DeltaMultiMode<-function(scAD, number_of_peaks){
  if(number_of_peaks == 3){
    x<-scAD$Postnorm_delta
    mem <- kmeans(x,3)$cluster
    mu1 <- mean(x[mem==1])
    mu2 <- mean(x[mem==2])
    mu3 <- mean(x[mem==3])
    sigma1 <- sd(x[mem==1])
    sigma2 <- sd(x[mem==2])
    sigma3 <- sd(x[mem==3])
    pi1 <- sum(mem==1)/length(mem)
    pi2 <- sum(mem==2)/length(mem)
    pi3 <- sum(mem==3)/length(mem)
    
    sum.finite <- function(x) {
      sum(x[is.finite(x)])
    }
    
    Q <- 0
    # starting value of expected value of the log likelihood
    Q[2] <- sum.finite(log(pi1)+log(dnorm(x, mu1, sigma1))) + sum.finite(log(pi2)+log(dnorm(x, mu2, sigma2))) + sum.finite(log(pi3)+log(dnorm(x, mu3, sigma3)))
    
    k <- 2
    
    while (abs(Q[k]-Q[k-1])>=1e-6) {
      # E step
      comp1 <- pi1 * dnorm(x, mu1, sigma1)
      comp2 <- pi2 * dnorm(x, mu2, sigma2)
      comp3 <- pi3 * dnorm(x, mu3, sigma3)
      comp.sum <- comp1 + comp2 + comp3
      
      p1 <- comp1/comp.sum
      p2 <- comp2/comp.sum
      p3 <- comp3/comp.sum
      
      # M step
      pi1 <- sum.finite(p1) / length(x)
      pi2 <- sum.finite(p2) / length(x)
      pi3 <- sum.finite(p3) / length(x)
      
      mu1 <- sum.finite(p1 * x) / sum.finite(p1)
      mu2 <- sum.finite(p2 * x) / sum.finite(p2)
      mu3 <- sum.finite(p3 * x) / sum.finite(p3)
      
      sigma1 <- sqrt(sum.finite(p1 * (x-mu1)^2) / sum.finite(p1))
      sigma2 <- sqrt(sum.finite(p2 * (x-mu2)^2) / sum.finite(p2))
      sigma3 <- sqrt(sum.finite(p3 * (x-mu3)^2) / sum.finite(p3))
      
      p1 <- pi1 
      p2 <- pi2
      p3 <- pi3
      
      k <- k + 1
      Q[k] <- sum(log(comp.sum))
    }
    p<-c(p1,p2,p3)
    p<-sort(p, decreasing = TRUE)
    mu<-c(mu1,mu2,mu3)
    mu<-sort(mu, decreasing = TRUE)
    sg<-c(sigma1.sigma2.sigma3)
    sg<-sort(sg, decreasing = TRUE)
    gm<-normalmixEM(x,k=3,lambda=c(round(p[1],2),round(p[2],2),round(p[3],3)),mu=c(round(mu[1],2),round(mu[2],2),round(mu[3],2)),sigma=c(round(sg[1],2),round(sg[2],2),round(sg[3],2)))
    
    hist(x, prob=T, breaks=32, xlim=c(range(x)[1], range(x)[2]), main='')
    lines(density(x), col="green", lwd=2)
    x1 <- seq(from=range(x)[1], to=range(x)[2], length.out=1000)
    y <- pi1 * dnorm(x1, mean=mu1, sd=sigma1) + pi2 * dnorm(x1, mean=mu2, sd=sigma2) + pi3 * dnorm(x1, mean=mu3, sd=sigma3)
    lines(x1, y, col="red", lwd=2)
    legend('topright', col=c("green", 'red'), lwd=2, legend=c("kernal", "fitted"))
    print(gm$mu)
    question<-(readline("Choose a number which lies on x-axis closest to the rightmost peak that will mark untagged cells (1,2 or 3)"))
    question<-as.numeric(question)
    if(question==1){
      threshold_delta<-round(gm$mu[1]+2*gm$sigma[1], 2)
      return(threshold_delta)
    }
    else if(question==2){
      threshold_delta<-round(gm$mu[2]+2*gm$sigma[2], 2)
      return(threshold_delta)
    }
    else if (question==3){
      threshold_delta<-round(gm$mu[3]+2*gm$sigma[3], 2)
      return(threshold_delta)
    }
    else{
      stop("Please choose a number from 1 to 3")
    }
  }
  else if (number_of_peaks == 4){
    x<-scAD$Postnorm_delta
    mem <- kmeans(x,4)$cluster
    mu1 <- mean(x[mem==1])
    mu2 <- mean(x[mem==2])
    mu3 <- mean(x[mem==3])
    mu4 <- mean(x[mem==4])
    sigma1 <- sd(x[mem==1])
    sigma2 <- sd(x[mem==2])
    sigma3 <- sd(x[mem==3])
    sigma4 <- sd(x[mem==4])
    pi1 <- sum(mem==1)/length(mem)
    pi2 <- sum(mem==2)/length(mem)
    pi3 <- sum(mem==3)/length(mem)
    pi4 <- sum(mem==4)/length(mem)
    
    sum.finite <- function(x) {
      sum(x[is.finite(x)])
    }
    
    Q <- 0
    # starting value of expected value of the log likelihood
    Q[2] <- sum.finite(log(pi1)+log(dnorm(x, mu1, sigma1))) + sum.finite(log(pi2)+log(dnorm(x, mu2, sigma2))) + sum.finite(log(pi3)+log(dnorm(x, mu3, sigma3))) + sum.finite(log(pi4)+log(dnorm(x, mu4, sigma4)))
    
    k <- 2
    
    while (abs(Q[k]-Q[k-1])>=1e-6) {
      # E step
      comp1 <- pi1 * dnorm(x, mu1, sigma1)
      comp2 <- pi2 * dnorm(x, mu2, sigma2)
      comp3 <- pi3 * dnorm(x, mu3, sigma3)
      comp4 <- pi4 * dnorm(x, mu4, sigma4)
      comp.sum <- comp1 + comp2 + comp3 + comp4
      
      p1 <- comp1/comp.sum
      p2 <- comp2/comp.sum
      p3 <- comp3/comp.sum
      p4 <- comp4/comp.sum
      
      # M step
      pi1 <- sum.finite(p1) / length(x)
      pi2 <- sum.finite(p2) / length(x)
      pi3 <- sum.finite(p3) / length(x)
      pi4 <- sum.finite(p4) / length(x)
      
      mu1 <- sum.finite(p1 * x) / sum.finite(p1)
      mu2 <- sum.finite(p2 * x) / sum.finite(p2)
      mu3 <- sum.finite(p3 * x) / sum.finite(p3)
      mu4 <- sum.finite(p4 * x) / sum.finite(p4)
      
      sigma1 <- sqrt(sum.finite(p1 * (x-mu1)^2) / sum.finite(p1))
      sigma2 <- sqrt(sum.finite(p2 * (x-mu2)^2) / sum.finite(p2))
      sigma3 <- sqrt(sum.finite(p3 * (x-mu3)^2) / sum.finite(p3))
      sigma4 <- sqrt(sum.finite(p4 * (x-mu4)^2) / sum.finite(p4))
      
      p1 <- pi1 
      p2 <- pi2
      p3 <- pi3
      p4 <- pi4
      
      k <- k + 1
      Q[k] <- sum(log(comp.sum))
    }
    p<-c(p1,p2,p3, p4)
    p<-sort(p, decreasing = TRUE)
    mu<-c(mu1,mu2,mu3, mu4)
    mu<-sort(mu, decreasing = TRUE)
    sg<-c(sigma1.sigma2.sigma3, sigma4)
    sg<-sort(sg, decreasing = TRUE)
    gm<-normalmixEM(x,k=3,lambda=c(round(p[1],2),round(p[2],2),round(p[3],3), round(p[4],2)),mu=c(round(mu[1],2),round(mu[2],2),round(mu[3],2), round(mu[4],2)),sigma=c(round(sg[1],2),round(sg[2],2),round(sg[3],2), round(sg[4],2)))
    
    hist(x, prob=T, breaks=32, xlim=c(range(x)[1], range(x)[2]), main='')
    lines(density(x), col="green", lwd=2)
    x1 <- seq(from=range(x)[1], to=range(x)[2], length.out=1000)
    y <- pi1 * dnorm(x1, mean=mu1, sd=sigma1) + pi2 * dnorm(x1, mean=mu2, sd=sigma2) + pi3 * dnorm(x1, mean=mu3, sd=sigma3) + pi4 * dnorm(x1, mean=mu4, sd=sigma4)
    lines(x1, y, col="red", lwd=2)
    legend('topright', col=c("green", 'red'), lwd=2, legend=c("kernal", "fitted"))
    print(gm$mu)
    question<-(readline("Choose a number which lies on x-axis closest to the rightmost peak that will mark untagged cells (1,2 or 3)"))
    question<-as.numeric(question)
    if(question==1){
      threshold_delta<-round(gm$mu[1]+2*gm$sigma[1], 2)
      return(threshold_delta)
    }
    else if(question==2){
      threshold_delta<-round(gm$mu[2]+2*gm$sigma[2], 2)
      return(threshold_delta)
    }
    else if (question==3){
      threshold_delta<-round(gm$mu[3]+2*gm$sigma[3], 2)
      return(threshold_delta)
    }
    else if (question==4){
      threshold_delta<-round(gm$mu[4]+2*gm$sigma[4], 2)
      return(threshold_delta)
    }
    else{
      stop("Please choose a number from 1 to 4")
    }
  }
  else{
    stop("For bimodal distribution use DeltaThreshold; otherwise the dataset is too noisy to select a confident threshold")
  }
}

FinalTagging<-function(scAD,threshold_ratio=0.5){
  if(threshold_ratio > 0.1 & threshold_ratio <= 1){
    final_tags<-scAD$Putative_tags
    future_undetermined<-WhichCells(scAD, expression=Ratio<threshold_ratio)
    for(i in future_undetermined){
      final_tags[i]<-'Undetermined'
    }
    scAD[['Final_tags']]<-final_tags
    return(scAD)
  }
  else if (threshold_ratio <= 0.1){
    stop("Please choose a number between 0.1 (exclusively) and 1")
  }
  else{
    stop("Please choose a number lower or equal to 1")
  }
}

WhyUntagged<-function(scAD, threshold_first_best, threshold_delta){
  undet_seurat<-subset(scAD, subset=Final_tags=='Undetermined')
  flags_for_undet<-c()
  for(i in 1:length(undet_seurat$Postnorm_first_best)){
    if(undet_seurat$Postnorm_first_best[i] < threshold_first_best & undet_seurat$Postnorm_delta[i] < threshold_delta){
      flags_for_undet<-c(flags_for_undet, 'Both failed')
    }
    else if (undet_seurat$Postnorm_first_best[i] < threshold_first_best & undet_seurat$Postnorm_delta[i] >= threshold_delta){
      flags_for_undet<-c(flags_for_undet, 'First best failed')
    }
    else{
      flags_for_undet<-c(flags_for_undet, 'Delta failed')
    }
  }
  undet_seurat[['Flags']]<-flags_for_undet
  df_undet<-data.frame('First_best'=undet_seurat$Postnorm_first_best, 'Delta'=undet_seurat$Postnorm_delta, 'Final-tags'=undet_seurat$Final_tags, 'Flags'=undet_seurat$Flags)
  question<-readline("How do you want to name the csv file with the untagged cells? (Don't forget to put .csv in the end)")
  question<-as.character(question)
  write.csv(df_undet, file=question)
  nums<-c()
  for(i in 1:length(df_undet$Flags)){
    if(df_undet$Flags[i]=='Both failed'){
      nums<-c(nums, length(which(df_undet$Flags=='Both failed')))
    }
    else if(df_undet$Flags[i]=='First best failed'){
      nums<-c(nums, length(which(df_undet$Flags=='First best failed')))
    }
    else{
      nums<-c(nums, length(which(df_undet$Flags=='Delta failed')))
    }
  }
  df_undet$Nums<-nums
  df_undet$Label = paste(df_undet$Flags," (", df_undet$Nums,")", sep = "")
  return(df_undet)
}


#####################################################
##                                                 ## 
##    RescueCluster: we don't filter MGPCs!        ##
##                                                 ##
#####################################################

RefinedDimPlot<-function(unfiltered_seurat, percent.mt_threshold){
  x<-WhichCells(unfiltered_seurat, expression = seurat_clusters == as.character(0) & percent.mt <= percent.mt_threshold)
  l<-list(x)
  alphabet<-c(rev(LETTERS),rev(letters))
  for (i in 1:(length(levels(unfiltered_seurat$seurat_clusters))-1)){
    x<-WhichCells(unfiltered_seurat, expression = seurat_clusters == as.character(i) & percent.mt <= percent.mt_threshold)
    l[[i+1]]<-x
  }
  l<-setNames(l, alphabet[1:(length(levels(unfiltered_seurat$seurat_clusters)))])
  DimPlot(unfiltered_seurat, reduction = "umap",label=T, cells.highlight = l, cols.highlight = hue_pal()(length(levels(unfiltered_seurat$seurat_clusters))), cols= 'grey')+NoLegend()
  
}

ClusterQC<-function(unfiltered_seurat, cluster_feature){
  if(cluster_feature == "cell.ID"){
    clusters<-levels(unfiltered_seurat$cell.ID)
    mean_percent_mt<-c()
    x<-unfiltered_seurat[[c('percent.mt','cell.ID')]]
    for (i in clusters){
      y<-x[x$cell.ID==i,]
      mean_percent_mt<-c(mean_percent_mt,mean(y$percent.mt))
    }
    
    mean_nfeatures<-c()
    x<-unfiltered_seurat[[c('nFeature_RNA','cell.ID')]]
    for (i in clusters){
      z<-x[x$cell.ID==i,]
      mean_nfeatures<-c(mean_nfeatures,mean(z$nFeature_RNA))
    }
    df<-data.frame(percent.mt=mean_percent_mt, nfeatures=mean_nfeatures)
    x1<-c()
    x2<-c()
    y1<-c()
    y2<-c()
    conditions<-c()
    for (i in 1:length(mean_percent_mt)){
      if (df[i,1]<=30 & df[i,2]<=1000){
        conditions<-c(conditions,'low gene count, low mgpc')
        x1<-c(x1,-Inf)
        x2<-c(x2,30)
        y1<-c(y1,-Inf)
        y2<-c(y2,1000)
      } else if (df[i,1]>30 & df[i,2]<=1000){
        conditions<-c(conditions,'low gene count, high mgpc')
        x1<-c(x1,30)
        x2<-c(x2,Inf)
        y1<-c(y1,-Inf)
        y2<-c(y2,1000)
      } else if (df[i,1]<=30 & df[i,2]>1000){
        conditions<-c(conditions,'high gene count, low mgpc')
        x1<-c(x1,-Inf)
        x2<-c(x2,30)
        y1<-c(y1,1000)
        y2<-c(y2,Inf)
      } else {
        conditions<-c(conditions,'high gene count, high mgpc')
        x1<-c(x1,30)
        x2<-c(x2,Inf)
        y1<-c(y1,1000)
        y2<-c(y2,Inf)
      }
    }
    
    df<-data.frame(percent.mt=mean_percent_mt, nfeatures=mean_nfeatures, x1=x1,x2=x2,y1=y1,y2=y2,conditions=conditions)
    ggplot(df, aes(x=percent.mt,y=nfeatures))+geom_rect(data=df,mapping=aes(xmin=x1,xmax=x2,ymin=y1,ymax=y2,fill=conditions), size=0.5,alpha=0.2)+scale_fill_manual(values=c('high gene count, high mgpc'='cyan', 'high gene count, low mgpc'='green','low gene count, high mgpc'='salmon','low gene count, low mgpc'= 'yellow'))+geom_point(size=2.5)+geom_text(label=clusters)+theme_classic()+theme(text = element_text(size=20))
  }
  else if (cluster_feature == "seurat_clusters"){
    clusters<-levels(unfiltered_seurat$seurat_clusters)
    mean_percent_mt<-c()
    x<-unfiltered_seurat[[c('percent.mt','seurat_clusters')]]
    for (i in clusters){
      y<-x[x$seurat_clusters==i,]
      mean_percent_mt<-c(mean_percent_mt,mean(y$percent.mt))
    }
    
    mean_nfeatures<-c()
    x<-unfiltered_seurat[[c('nFeature_RNA','seurat_clusters')]]
    for (i in clusters){
      z<-x[x$seurat_clusters==i,]
      mean_nfeatures<-c(mean_nfeatures,mean(z$nFeature_RNA))
    }
    df<-data.frame(percent.mt=mean_percent_mt, nfeatures=mean_nfeatures)
    x1<-c()
    x2<-c()
    y1<-c()
    y2<-c()
    conditions<-c()
    for (i in 1:length(mean_percent_mt)){
      if (df[i,1]<=30 & df[i,2]<=1000){
        conditions<-c(conditions,'low gene count, low mgpc')
        x1<-c(x1,-Inf)
        x2<-c(x2,30)
        y1<-c(y1,-Inf)
        y2<-c(y2,1000)
      } else if (df[i,1]>30 & df[i,2]<=1000){
        conditions<-c(conditions,'low gene count, high mgpc')
        x1<-c(x1,30)
        x2<-c(x2,Inf)
        y1<-c(y1,-Inf)
        y2<-c(y2,1000)
      } else if (df[i,1]<=30 & df[i,2]>1000){
        conditions<-c(conditions,'high gene count, low mgpc')
        x1<-c(x1,-Inf)
        x2<-c(x2,30)
        y1<-c(y1,1000)
        y2<-c(y2,Inf)
      } else {
        conditions<-c(conditions,'high gene count, high mgpc')
        x1<-c(x1,30)
        x2<-c(x2,Inf)
        y1<-c(y1,1000)
        y2<-c(y2,Inf)
      }
    }
    
    df<-data.frame(percent.mt=mean_percent_mt, nfeatures=mean_nfeatures, x1=x1,x2=x2,y1=y1,y2=y2,conditions=conditions)
    ggplot(df, aes(x=percent.mt,y=nfeatures))+geom_rect(data=df,mapping=aes(xmin=x1,xmax=x2,ymin=y1,ymax=y2,fill=conditions), size=0.5,alpha=0.2)+scale_fill_manual(values=c('high gene count, high mgpc'='cyan', 'high gene count, low mgpc'='green','low gene count, high mgpc'='salmon','low gene count, low mgpc'= 'yellow'))+geom_point(size=2.5)+geom_text(label=clusters)+theme_classic()+theme(text = element_text(size=20))
  }
  else{
    stop("Please use either seurat_clusters or cell.ID as cluster features.")
  }
}
####################################################################################
##                                                                                ## 
##  AutoClusterType: automated annotation by top 3 best annotations in a cluster! ##
##                                                                                ##
####################################################################################

AutoClusterType<-function(reference_file_path, mca_result, unfiltered_seurat){
  references<-read.csv(reference_file_path, header = TRUE)
  x<-as.data.frame(mca_result$scMCA_probility)
  for (i in 0:(length(levels(unfiltered_seurat$seurat_clusters))-1)){
    assign(paste('Cell_type',as.character(i),sep='_'), WhichCells(unfiltered_seurat, expression = seurat_clusters == as.character(i)))
  }
  for (i in 0:(length(levels(unfiltered_seurat$seurat_clusters))-1)){
    obj<-get(paste('Cell_type',as.character(i),sep='_'))
    ann<-data.frame()
    for(j in 1:length(obj)){
      ann<-rbind(ann, x[which(x$Cell==obj[j]),])
    }
    assign(paste('Annotation',as.character(i),sep='_'), ann)
  }
  
  for(i in 0:(length(levels(unfiltered_seurat$seurat_clusters))-1)){
    obj<-get(paste('Annotation',as.character(i),sep='_'))
    inter_first<-data.frame()
    for(j in seq(1, length(obj$Cell),3)){
      inter_first<-rbind(inter_first, obj[j,])
    }
    inter_second<-data.frame()
    for(k in seq(2, length(obj$Cell),3)){
      inter_second<-rbind(inter_second, obj[k,])
    }
    inter_third<-data.frame()
    for(l in seq(3, length(obj$Cell),3)){
      inter_third<-rbind(inter_third, obj[l,])
    }
    final_inter<-cbind('Cell'=inter_first$Cell,'First_best'=as.character(inter_first$`Cell type`), 'First_best_score'=inter_first$Score, 'Second_best'=as.character(inter_second$`Cell type`), 'Second_best_score'=inter_second$Score, 'Third_best'=as.character(inter_third$`Cell type`), 'Third_best_score'=inter_third$`Score`)
    assign(paste('Final_annotation',as.character(i),sep='_'), final_inter)
  }
  
  
  Final_first_best<-c()
  Final_first_best_score<-c()
  Final_second_best<-c()
  Final_second_best_score<-c()
  Final_third_best<-c()
  Final_third_best_score<-c()
  for (i in 0:(length(levels(unfiltered_seurat$seurat_clusters))-1)){
    x<-get(paste('Final_annotation',as.character(i),sep='_'))
    x<-as.data.frame(x)
    one_cell<-c()
    one_org<-c()
    two_cell<-c()
    two_org<-c()
    three_cell<-c()
    three_org<-c()
    for (j in 1:length(x$Cell)){
      one_cell<-c(one_cell,references[references$Cell.types==x$First_best[j],][2][1,])
      one_org<-c(one_org,references[references$Cell.types==x$First_best[j],][3][1,])
      two_cell<-c(two_cell,references[references$Cell.types==x$Second_best[j],][2][1,])
      two_org<-c(two_org,references[references$Cell.types==x$Second_best[j],][3][1,])
      three_cell<-c(three_cell,references[references$Cell.types==x$Third_best[j],][2][1,])
      three_org<-c(three_org,references[references$Cell.types==x$Third_best[j],][3][1,])
    }
    d<-data.frame('First_best_cell_type'=one_cell, 'First_best_organ'=one_org, 'First_best_score'=x$First_best_score, 'Second_best_cell_type'=two_cell, 'Second_best_organ'=two_org, 'Second_best_score'=x$Second_best_score, 'Third_best_cell_type'=three_cell, 'Third_best_organ'=three_org, 'Third_best_score'=x$Third_best_score)
    d$Full_first_best<-paste(d$First_best_cell_type,d$First_best_organ,sep='_')
    d$Full_second_best<-paste(d$Second_best_cell_type,d$Second_best_organ,sep='_')
    d$Full_third_best<-paste(d$Third_best_cell_type,d$Third_best_organ,sep='_')
    assign(paste('Refined_annotation',as.character(i),sep='_'), d)
  }
  
  for (i in 0:(length(levels(unfiltered_seurat$seurat_clusters))-1)){
    obj<-get(paste('Refined_annotation',as.character(i),sep='_'))
    fb<-unique(obj$Full_first_best)
    fbs<-c()
    fbp<-c()
    sb<-unique(obj$Full_second_best)
    sbs<-c()
    sbp<-c()
    tb<-unique(obj$Full_third_best)
    tbs<-c()
    tbp<-c()
    for(j in 1:length(fb)){
      x<-obj[which(obj$Full_first_best==fb[j]),]
      fbs<-c(fbs,mean(as.numeric(x$First_best_score)))
      fbp<-c(fbp,length(x$First_best_cell_type)/length(obj$First_best_cell_type)*100)
    }
    for(j in 1:length(sb)){
      x<-obj[which(obj$Full_second_best==sb[j]),]
      sbs<-c(sbs,mean(as.numeric(x$Second_best_score)))
      sbp<-c(sbp,length(x$Second_best_cell_type)/length(obj$Second_best_cell_type)*100)
    }
    for(j in 1:length(tb)){
      x<-obj[which(obj$Full_third_best==tb[j]),]
      tbs<-c(tbs,mean(as.numeric(x$Third_best_score)))
      tbp<-c(tbp,length(x$Third_best_cell_type)/length(obj$Third_best_cell_type)*100)
    }
    max.len<-max(length(fb),length(fbs),length(fbp),length(sb), length(sbs), length(sbp), length(tb),length(tbs),length(tbp))
    fb<-c(fb, rep(' ', max.len - length(fb)))
    fbs<-c(fbs, rep(' ', max.len - length(fbs)))
    fbp<-c(fbp, rep(' ', max.len - length(fbp)))
    sb<-c(sb, rep(' ', max.len - length(sb)))
    sbs<-c(sbs, rep(' ', max.len - length(sbs)))
    sbp<-c(sbp, rep(' ', max.len - length(sbp)))
    tb<-c(tb, rep(' ', max.len - length(tb)))
    tbs<-c(tbs, rep(' ', max.len - length(tbs)))
    tbp<-c(tbp, rep(' ', max.len - length(tbp)))
    r_1<-data.frame('First_best'=fb,'Average_first_best_score'=fbs,'First_best_percent_in_cluster'=fbp)
    r_1$'First_best_percent_in_cluster'<-as.numeric(as.character(r_1$'First_best_percent_in_cluster'))
    r_1<-r_1[order(r_1$'First_best_percent_in_cluster', decreasing=TRUE),]
    r_2<-data.frame('Second_best'=sb,'Average_second_best_score'=sbs,'Second_best_percent_in_cluster'=sbp)
    r_2$'Second_best_percent_in_cluster'<-as.numeric(as.character(r_2$'Second_best_percent_in_cluster'))
    r_2<-r_2[order(r_2$'Second_best_percent_in_cluster', decreasing = TRUE),]
    r_3<-data.frame('Third_best'=tb,'Average_third_best_score'=tbs,'Third_best_percent_in_cluster'=tbp)
    r_3$'Third_best_percent_in_cluster'<-as.numeric(as.character(r_3$'Third_best_percent_in_cluster'))
    r_3<-r_3[order(r_3$'Third_best_percent_in_cluster', decreasing = TRUE),]
    r<-cbind(r_1,r_2,r_3)
    r[is.na(r)]<-''
    assign(paste('Summary_cluster_',as.character(i),sep=''), r)
    write.csv(r,file=paste('Summary_cluster_',as.character(i),'.csv',sep=''), row.names=FALSE)
  }
  annotation_df<-data.frame()
  for(i in 0:(length(levels(unfiltered_seurat$seurat_clusters))-1)){
    obj<-get(paste('Summary_cluster_',as.character(i),sep=''))
    annotation_df<-rbind(annotation_df, obj[1,])
  }
  rownames(annotation_df)<-NULL
  clusters<-c()
  for(i in 0:(length(levels(unfiltered_seurat$seurat_clusters))-1)){
    clusters<-c(clusters,as.character(i))
  }
  annotation_df$Cluster<-clusters
  final_names<-c()
  for(i in 1:(length(levels(unfiltered_seurat$seurat_clusters)))){
    if(as.numeric(annotation_df$First_best_percent_in_cluster[i])>=50 & as.numeric(annotation_df$Average_first_best_score[i])>=0.4){
      final_names<-c(final_names, annotation_df$First_best[i])
    }
    else if(as.numeric(annotation_df$First_best_percent_in_cluster[i])>=50 & as.numeric(annotation_df$Average_first_best_score[i])<0.4){
      final_names<-c(final_names, paste(annotation_df$First_best[i], 'low_score', sep='_'))
    }
    else{
      if(annotation_df$First_best[i] == annotation_df$Second_best[i] & annotation_df$First_best[i] == annotation_df$Third_best[i]){
        final_names<-c(final_names, paste(annotation_df$First_best[i], 'low_score_and_cell_number',sep='_'))
      }
      else if(annotation_df$First_best[i] == annotation_df$Second_best[i] & annotation_df$First_best[i] != annotation_df$Third_best[i]){
        final_names<-c(final_names, paste(annotation_df$First_best[i], annotation_df$Third_best[i],sep='/'))
      }
      else if(annotation_df$First_best[i] != annotation_df$Second_best[i] & annotation_df$First_best[i] == annotation_df$Third_best[i]){
        final_names<-c(final_names, paste(annotation_df$First_best[i], annotation_df$Second_best[i],sep='/'))
      }
      else if(annotation_df$First_best[i] != annotation_df$Second_best[i] & annotation_df$Second_best[i] == annotation_df$Third_best[i]){
        final_names<-c(final_names, paste(annotation_df$First_best[i], annotation_df$Second_best[i],sep='/'))
      }
      else{
        final_names<-c(final_names, paste(annotation_df$First_best[i], annotation_df$Second_best[i],annotation_df$Third_best[i],sep='/'))
      }
    }
  }
  new_final_names<-final_names
  for(i in 1:length(unique(final_names))){
    for(j in 1:length(which(final_names==unique(final_names)[i]))){
      if(length(which(final_names==unique(final_names)[i]))>1){
        new_final_names[which(final_names==unique(final_names)[i])[j]]<-paste(final_names[which(final_names==unique(final_names)[i])[j]],as.character(j),sep='_')
      }
    }
  }
  Idents(unfiltered_seurat) <- "seurat_clusters"
  new.cluster.ids<-new_final_names
  names(new.cluster.ids)<-levels(unfiltered_seurat$seurat_clusters)
  unfiltered_seurat<-RenameIdents(unfiltered_seurat, new.cluster.ids)
  unfiltered_seurat[["cell.ID"]] <- Idents(unfiltered_seurat)
  return(unfiltered_seurat)
}
