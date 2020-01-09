######################################################################################
# name : R_otolith_microchem_elements_&_method_selection.R
# original name: R_otolith_microchemistry
# Date of creation : 2008  
#################################
#
# R Script for the calculation of the accuracy in fish habitat discrimination  
#  using otolith chemical signatures, for:
#
#  1) all possible combinations (size and content) of the N chemical elements measured 
#  2) 4 different statistical methods (LDA, QDA, Random Forest & Neural Networks)
#  -> allows selection of the statistical method & list of elements leading to optimal 
#     discrimination of the k habitats tested.
#
#################################
# To refer to this script:
#
#  Mercier et al (2011). Selecting statistical models and variable combinations for   
#  optimal classification using otolith microchemsitry. Ecol. Appl. 21(4): 1352-1364.
#  
# Original code (2008) by Leny Mercier (lenymercier@gmail.com)
# Modified (2012) by: Audrey Darnaude (audrey.darnaude@univ-montp2.fr) 
#                     Jennifer Tournois (Jennifer.Tournois@univ-montp2.fr)
#
######################################################################################
# Note before running the script...
# 
# Dataset requirements :
#   - sample size must exceed N = 10 for each of the k sites (habitats) 
#   - each of the N elements must have been measured in all the samples (so "0" values in
#     the matrix only correspond to null concentrations for an element in the corresponding 
#     sample - and not to the absence of measured data)
# 
# To run this script :
#   1) organize your data in EXCEL ("data.xls") to make a general table showing :
#        - in column 1 = the name (code) for all the samples analysed (at least 6 per site) for all sites 
#             (with the head "Sample name" on line 1)
#        - in column 2 = the sampling site (= 1 of the k habitats to discriminate) for each sample (line) 
#             (with the head "Sample origin" on line 1)
#        - in the following columns = the concentrations measured in all samples for each of the N elements
#             (with element names as column heads)
#   2) save the active worksheet that contains the table in ".csv" ("data.csv")
#          NB: this file ("data.csv") will be used to create the R objects necessary to run this script
#   3) once in R, set the "working directory" as appropriate 
#   4) enter the following script lines:
#     4.1) to download necessary files
#        > datafile <- read.csv("data.csv", header=TRUE, row.names=1, sep=";", dec=",")
#        > head(datafile) # to check that column 1 of "datafile" does indicate sample origin
#        > origin <- as.vector(datafile[,1]) 
#        > origin # to check that "origin" = column 1 of "datafile" (gives sample origin as a R vector) 
#        > test_data <- datafile[,-1]
#        > View(test_data) # to check that test_data = "datafile" without its column 1 (now in "origin")
#     4.2) to run the function "compelt" for each statistical method (here an exemple with rf) and save the results of the analysis
#        > source("R_otolith_microchem_elements_&_method_selection.R") # downloads this script
#        > res=compelt(test_data,origin,method="rf",iteration=10,ratio=3/4,ntree=50)
#          # the arguments for the "compelt" function (test_data, origin, etc) are explained below
#          # the text "Analysis in progress/Please wait" should appear 
#          # then it might take a while for the successive lists of elements to appear...
#          # the script is made to generate a namefile (eg here rf_10_50) for each analysis listing :
#             # the method tested
#             # the number of iterations chosen
#             # the number of trees chosen (see below for more details)
#          # the graph showing the min/mean/max accuracy for each combination size should appear as soon as the analysis is completed
#          # it will be automatically saved under "plot_namefile.jpeg" (with namefile = method_iteration_ntree) in your current working directory
#          # the results summary (=the average accuracy reached + the exact combinations leading to min & max accuracy for each combination size) will be saved under "Summary_Acc_namefile.rdata" 
#          # and the rescue file for this analysis will be saved under "namefile_rescue.txt"
#
###########################
#        if you ever want to retrieve the results summary for a given analysis in R
#        > load("Summary_Acc_namefile.rdata") # with "summary_Acc_namefile" being the name of the results summary file (e.g."Summary_Acc_rf_10_50.rdata")
#        > Summary_Acc_namefile # to visualize the content of the file in R (e.g.Summary_Acc_rf_10_50)
#################################################################################################################################                      
# SCRIPT

# Downloading of the required packages
  # if the packages have never been installed on your computer:
    # install.packages("mda", dependencies=TRUE) 
        # for France, choose CRAN mirror = 26 (Lyon1) or 25 (Toulouse)
    # install.packages("randomForest", dependencies=TRUE) 
    # install.packages("gtools")
    # install.packages("FactoMineR")
    # require(MASS)
    # require(class)
    # require(nnet)
   
  library(mda)
  library(randomForest)
  library(gtools)
  library(FactoMineR)
  library(MASS)
  library(class)
  library(nnet)

  
# Creation of the "compelt" function 
  # = to calculate, for each statistical method, the accuracy (minimum, maximum & average) of the discrimination 
    # between k given habitats (or sites) for all the possible combinations of 1 to N chemical elements
    # measured in > 10 otolith samples per site

  # Arguments of the "compelt" function:
    # test_data: the samples x elements data frame (with samples in rows, elements in columns) to be used for dicrimination
    # origin: a vector giving the origin (habitat or site) for each sample of test_data
    # method = the statistical method for discrimination (choice between 4)
        # "rf" = Random Forest (default value)
        # "lda" = Linear Discriminant Analysis,
        # "qda" = Quadratic Discriminant Analysis
        # "nn" = artificial Neural Networks
    # iteration: the number of replicates of the test for each element combination
        # (default is 10 -so the running of the script does not take too long- but iteration=1000 is recommended)
    # ratio: the ratio used for random selection of samples during the cross validation procedure for each iteration
        # (default is 3/4 = 75% of the samples of test_data are used for training the statistics & 25% for testing it)
    # ntree: number of trees used for the Random Forest method
        # (default is 50 -so the running of the script does not take too long- but ntree=500 is recommended)

 compelt <- function(test_data,origin,method="rf",iteration=10,ratio=3/4,ntree=50)
 {# start of "compelmt" (which includes loops for the "cacAcc" function, etc - see below)
    
  # if any, removal of the rescue file from a previous analysis with the same parameters (so the summary table for the present analysis does not include the iterations made in the previous one !!!!!)
    file_output=paste(method,iteration,ntree,"rescue.txt",sep="_")
    if (sum(file_output==list.files())!=0)
      { #start of if
      eval(parse(text =paste("file.remove ('",file_output,"')",sep="")))
      } # end of if
  
  # test for appropriate data format before starting the analysis (generation of an error message otherwise)
      if (is.data.frame(test_data)==F) stop("test_data is not a data frame \n")
      if (is.vector(origin)==F) stop("origin is not a vector \n")

  # Recall of analysis parameters (+ message "Please wait...")
    cat("Analysis parameters:\n") # n = return to the next line
    cat(" - ",method,"method used\n")
    cat(" - ", iteration,"replicates for each combination\n")
    cat(" -  Proportion of the data used for training =",ratio,"\n")
    cat("\n")
    cat("Analysis in progress\n")
    cat("   Please wait...\n")
 
  # Generation of all the possible combinations of 1 to N chemical elements (there is 2exp(N-1) possibilities)
    elmt = colnames(test_data)
    nbelmt = length(elmt)
    combelmt = lapply(seq(1:nbelmt),function(index) combinations(nbelmt,index,elmt))

  # Creation of the function for calculating the accuracy (min/mean/max) of each statistical method for all element combinations
    # + displaying the results
        # = function "calcAcc" with the same arguments than the function "compelt" (test_data,origin,method,iteration,ratio,ntree)
        # + "namecombelmt" = name of the element combination tested (autogenerated)
             
   calcAcc = function(test_data,origin,method,iteration,ratio,ntree,namecombelmt) 
   {# start of "calcAcc"
	         
    # to select & display successively all the possible combinations of the N chemical element present in test_data
	    data=as.data.frame(test_data[,namecombelmt])
	    if (length(namecombelmt)== 1) # generates the combination
       {# start of if
        colnames(test_data)= namecombelmt
       } # end of if
    	cat("\n")
      cat(namecombelmt, '\n') # displays the combination

    # Loop for the i iterations for each statistical method (lda, qda, rf or nn)
	    vectormet=vector() # loop on the total number of iterations for each method (defined below)
    	for (i in 1:iteration) 
      {# start of for i (iterations)
       
       # Random selection of the samples used for the "cross-validation" procedure in each iteration:
            # = building of the training ("dataTrain") and testing ("dataTest") datasets 
            # + the corresponding vectors indicating sample origin ("originTrain" & "originTest")
            # !!! made according to the "ratio" indicated in the "compelmt" function (defaut is 75%:25%)
              # So to ensure reliability of the accuracies calculated, the script is set to work only if 
              # >7 samples per habitat in the training dataset (which, if N total per habitat = 10,  
              # leaves 3 samples (a minimum!) per habitat for the cross-validation in each iteration
		     
        flag=FALSE
		     while(flag==FALSE) 
           {# start of while
			     v<- sample(dim(test_data)[1],round(ratio*dim(test_data)[1]),replace=FALSE)
			     dataTrain=as.data.frame(data[v,])
			     dataTest=as.data.frame(data[-v,])
			     if (length(namecombelmt)== 1) {colnames(dataTrain)= namecombelmt}
			     if (length(namecombelmt)== 1) {colnames(dataTest)= namecombelmt}
			     originTrain = origin[v,]
			     originTest = origin[-v,]
			     if(min(table(origin[v,])) > 7) # requires >7 samples per habitat in the training dataset 
			       # this condition (N>7) is valid for N total=10 per habitat to leave a minimum of N=3 for the
			       # cross-validation (if total N is greater you might want to increase it, e.g. >20 for Ntotal =30).
            {# start of if
			       cat(" sampling for cross-validation = OK for iteration",i, "\n") 
             #you might want to hide this for more than 10 iterations
			       flag = TRUE
			      }# end of if
           }# end of while
                    
        # for LDA
          if(is.na(match("lda",method))== FALSE)
          {#start of if(lda)
            try(trainlda<-lda(x=dataTrain,grouping=originTrain))
            cat("  LDA training = OK \n") #you might want to hide this for more than 10 iterations
            vectormet[i]=sum(diag(confusion(predict(trainlda,dataTest)$class,as.matrix(originTest)) ))/length(originTest)
            cat("  LDA testing = OK \n") #you might want to hide this for more than 10 iterations
          }#end of if(lda)
         
		    # for QDA
		      if(is.na(match("qda",method))== F)
		      {#start of if(qda)
		       try(trainqda<-qda(x=dataTrain,grouping=originTrain))
		       cat("  QDA training = OK \n") #you might want to hide this for more than 10 iterations
           vectormet[i]=sum(diag(confusion(predict(trainqda,dataTest)$class,as.matrix(originTest)) ))/length(originTest)
		       cat("  QDA testing = OK \n") #you might want to hide this for more than 10 iterations
          }#end of if(qda)
		     
	      # for Random Forest
      	  if(is.na(match("rf",method))== F)
	        {#start of if(rf)
    		   modelrf= randomForest(x=as.matrix(dataTrain),y=originTrain,ntree=ntree,importance=TRUE,proximity=TRUE)
    		   cat("  RF training = OK \n") #you might want to hide this for more than 10 iterations
           vectormet[i]=rF= sum(diag(confusion(predict(modelrf,dataTest),originTest)))/length(originTest)
    		   cat("  RF testing = OK \n") #you might want to hide this for more than 10 iterations
          }#end of if (rf)

	      # for Neural Network
	        if (is.na(match("nn",method))==F)
	        {#start of if(nn)
	         
            conftab = function(origin, predit) # to generate the confusion table for NN
           {#start of conftab
	          true <- max.col(origin)
	          pred <- max.col(predit)
	          a=table(true, pred)
	          b=matrix(0,dim(origin)[2],dim(origin)[2])
	          colnames(b)=rownames(b)=1:dim(origin)[2]
	          for (i in rownames(a))
		         {#start of for i
		           for (j in colnames(a))
			         {#start of for j
			          b[i,j]=a[i,j]
			         }#end of for j
		         }#end of for i
      	    colnames(b)=rownames(b)=colnames(origin)
	          return(b)
           }#end of conftab
          
           originNN=class.ind(as.matrix(origin))
	  	     modelnn= nnet(x=dataTrain,y=originNN[v,],size=25, decay=1e-4, maxit=2000)
           cat("  NN training = OK \n") #you might want to hide this for more than 10 iterations
           vectormet[i]= sum(diag(conftab(originNN[-v,],predict(modelnn, dataTest))))/sum(conftab(originNN[-v,],predict(modelnn, dataTest)))
           cat("  NN testing = OK \n") #you might want to hide this for more than 10 iterations
         }#end of if (nn)

      # building of a data frame with the results:
        resmodel=matrix(NA,1,4) # creates an empty matrix with 1 line and 4 colums 
        #names(resmodel)=c("n","Elements",sort(c(paste(method,"_mean",sep=""),paste(method,"_sd",sep=""))))
        resmodel=as.data.frame(resmodel)
        resmodel[,1]= length(namecombelmt) # gives combination size in column 1 of the matrix
        resmodel[,2]= mean(vectormet) # gives the mean value for this combination size in column 2 of the matrix
        resmodel[,3]= sd(vectormet) # gives the sd for this combination size in column 3 of the matrix
        # to give the combination corresponding to these mean & sd values in column 4 of the matrix
         nmtr="" # to give the combination corresponding to these mean & sd values in column 4 of the matrix
         for(k in 1:length(namecombelmt))
		      {#start of for k
		       nmtr<-(paste(nmtr,namecombelmt[k],sep=" "))#la liste des éléments sans les guillemets
		      }#end of for k
        resmodel[,4] = nmtr 
        for(u in 1:3)
         {#start of for u
          resmodel[,u]=as.numeric(resmodel[,u])
         }#end of for u
    
    }#end of for i (iterations)

    # to save the results
	   file_output_name=paste(method,iteration,ntree,sep="_") 
     write.table(resmodel,file=paste(file_output_name,"rescue.txt",sep="_"),append=TRUE,col.names=F,row.names=F)
	   return(resmodel)

  }#end of calcAcc

 # to get the results of "calcAcc" (min/mean/max accuracy) per combination size (irrespective of the elements included in the combination)
      # = function "applycalcAcc" with the same arguments than the function "compelt" (test_data,origin,method,iteration,ratio,ntree)
      # + "sizecombelmt" = size of element combination tested (autogenerated)
   applycalcAcc = function(sizecombelmt,origin,test_data,ratio,method,iteration,ntree) 
   {#start of applycalcAcc
    return(apply(sizecombelmt,1,calcAcc,origin=origin,test_data=test_data,ratio=ratio,method=method,iteration=iteration,ntree=ntree))
    }#end of applycalcAcc

 # to get the results of "applycalcAcc" (min/mean/max accuracy) for all possible sizes of element combinations
      # = function "applycalcAcc" with the same arguments than the function "compelt" (test_data,origin,method,iteration,ratio,ntree)
      # + "sizecombelmt" = size of element combination tested (autogenerated)
   sapplycalcAcc = function (combelmt,origin,test_data,ratio,method,iteration,ntree)
   {#start of sapplycalcAcc
    return (sapply( combelmt,applycalcAcc,origin=origin,test_data=test_data,ratio=ratio,method=method,iteration=iteration,ntree=ntree,simplify=F))
    }#end of sapplycalcAcc


 # to create a list of the results for each combination size
   reslist= sapplycalcAcc(combelmt,as.data.frame(origin),test_data,ratio,method,iteration,ntree)

 # to build the table giving a summary of the min/mean/max accuracies obtained for all combination sizes
   cat("summary table building...\n")
   file_output_name=paste(method,iteration,ntree,sep="_")
   eval(parse (text=paste("resMat = read.table('",file_output_name,"_rescue.txt')",sep=""))) 
     # saves the rescue file for this analysis
   colnames(resMat) <- c("Combination size","Mean value","SD","Corresponding combination")

    # for the maximum accuracy
      maxAcc=resMat[1,]
      for (i in 1 : dim(test_data)[2]) 
       {#start of for i
    	  bob=resMat[which(resMat[,1]==i),]
	      maxAcc[i,]=bob[which.max(bob[,2]),]
       }#end of for i

    # for the minimum accuracy
      minAcc=resMat[1,]
      for (i in 1 : dim(test_data)[2]) 
        {#start of for i
	       bob=resMat[which(resMat[,1]==i),]
	       minAcc[i,]=bob[which.min(bob[,2]),]
        }#end of for i
      
    # for the mean accuracy
      meanAcc=resMat[1,]
      meanAcc=meanAcc[,-4]
      for (i in 1:dim(test_data)[2]) 
        {#start of for i
	       meanAcc[i,1]=i
	       meanAcc[i,2]=mean(resMat[which(resMat[,1]==i),][,2])
	       meanAcc[i,3]=sd(resMat[which(resMat[,1]==i),][,2])
	      }#end of for i

 # to calculate the SD of the mean accuracy for each combination size
   moyvar<-function(x,moy,er,coul="red",form=22,siz=1) 
    {# start of moyvar
     segments(x,moy-er,x,moy+er) ; segments(x-0.05,moy-er,x+0.05,moy-er) ;
     segments(x-0.05,moy+er ,x+0.05,moy+er)
     points(x,moy,pch=form,col="black",bg=coul,cex=siz)
    } # end of moyvar

 # creation of a graph showing the summary of the results
   # = Accuracy (min, mean +-SD & max) for increasing combination size
   # + maximum accuracy obtained with the method tested (and corresponding combination of elements) 
   plot( 1:dim(test_data)[2], maxAcc[,2],type="n",xlim=c(0,dim(test_data)[2]+1),ylim=c(0,1),xlab="Combination size", ylab="Accuracy")
   polygon(rbind(maxAcc[,1:2],minAcc[dim(test_data)[2]:1,1:2]),border="red",lwd=1.5, col="#bb1f1f50")
   moyvar(meanAcc[,1],meanAcc[,2],meanAcc[,3],coul="purple")
   points(which.max(maxAcc[,2]),max(maxAcc[,2]),pch=20, bg="grey30")
   text(which.max(maxAcc[,2]),max(maxAcc[,2]),paste("max accuracy = ",round(max(maxAcc[,2]),2)), pos=3)
   text(which.max(maxAcc[,2]),0.2,paste("best combination =\n",maxAcc[which.max(maxAcc[,2]),4]),pos=3)
  # to save the graph in "jpeg"
   file_output_name=paste(method,iteration,ntree,sep="_")
   jpeg(filename = paste("Plot_",file_output_name,".jpeg",sep=""), width =20, height = 20, units = "cm",res=100)
   plot( 1:dim(test_data)[2], maxAcc[,2],type="n",xlim=c(0,dim(test_data)[2]+1),ylim=c(0,1),xlab="Combination size", ylab="Accuracy")
   polygon(rbind(maxAcc[,1:2],minAcc[dim(test_data)[2]:1,1:2]),border="red",lwd=1.5, col="#bb1f1f50")
   moyvar(meanAcc[,1],meanAcc[,2],meanAcc[,3],coul="purple")
   points(which.max(maxAcc[,2]),max(maxAcc[,2]),pch=20, bg="grey30")
   text(which.max(maxAcc[,2]),max(maxAcc[,2]),paste("maximum accuracy = ",round(max(maxAcc[,2]),2)), pos=3)
   text(which.max(maxAcc[,2]),0.2,paste("best combination =\n",maxAcc[which.max(maxAcc[,2]),4]),pos=3)
   dev.off()

 # to save the Summary table of the results in a file listing the parameters chosen for the analysis
   eval(parse(text=paste("Summary_Acc_",method,"_",iteration,"_",ntree,"= list(maxAcc=maxAcc,minAcc=minAcc,meanAcc=meanAcc)",sep="")))
   eval(parse(text=paste( "attr(Summary_Acc_",method,"_",iteration,"_",ntree,",'class') = 'compelt'",sep=""))) 
   eval(parse(text=paste("return(Summary_Acc_",method,"_",iteration,"_",ntree,")",sep="")))
   eval(parse(text= paste("save(Summary_Acc_",method,"_",iteration,"_",ntree,",file='Summary_Acc_",method,"_",iteration,"_",ntree,".rdata')",sep="")))
       

}#end of compelt
# note : if you run this script with a low number of iterations
  # the results for each method (min/mean/max Accuracy) might change from one run to the other
  # therefore, after trying the script on your dataset with 10 iterations &  50 trees (to check that it works)...
  # ... run it with 1000 iterations and 500 trees to get reliable (stable) results 