

msGARCHfit=function(eq,Data,numberMCMC=10000,numberBurn=500,GARCHmodels=c("sGARCH","sGARCH"),experiment="",estDate=Sys.Date()){
  
  
  
  
  outputData=Data
  
# Starts calculation time:  
  startCalculation=Sys.time()

  words <- strsplit(eq, "[^[:alnum:]]+") # Split the string using non-alphanumeric characters as separators
  words = words[[1]]
  words=words[-1]
  num_words <- length(words)
  
  variables=c("(Intercept)",words)
  variablesValues=rep(0,num_words)
  
  # Fit the full model 
  full.model <- lm(eq, data = Data)
  # Stepwise regression model
  step.model <- stepAIC(full.model, direction = "both", trace = FALSE)
  
  coeficientesStepWise=summary(step.model)$coefficients
  
  num_VarsStep=nrow(coeficientesStepWise)
  name_VarsStep=rownames(coeficientesStepWise)
  
  
  DBTable=data.frame(
    Date=as.character(tail(outputData$Date,1)),
    Value=rep(0,length(variables)),
    Ticker=variables,
    ModelID="factor model coefs",
    GARCHSpec="factor model coefs",
    Experiment=experiment
  )
  
  for (numVar in 1:nrow(DBTable)){
    
    idCoefRow=which(name_VarsStep[numVar]==DBTable$Ticker)
    if (length(idCoefRow)>0){
      DBTable$Value[idCoefRow]=coeficientesStepWise[numVar] 
    }

  }  
  
  DBTable=rbind(DBTable,
                data.frame(
                           Date=as.character(tail(outputData$Date,1)),
                           Value=tail(datos$Settle,1),
                           Ticker="Price at t",
                           ModelID="Price at t",
                           GARCHSpec="Price at t",
                           Experiment=experiment
                           ),
                data.frame(
                  Date=as.character(tail(outputData$Date,1)),
                  Value=tail(datos$Return,1),
                  Ticker="Return at t",
                  ModelID="Return at t",
                  GARCHSpec="Return at t",
                  Experiment=experiment
                ),
                data.frame(
                  Date=as.character(tail(outputData$Date,1)),
                  Value=tail(step.model$fitted.values,1),
                  Ticker="Return Forecast at t",
                  ModelID="Return Forecast at t",
                  GARCHSpec="Return Forecast at t",
                  Experiment=experiment
                )                   
                )
  # Extract the residuals from the model:
  residuals=step.model$residuals
  
  # Creating the DIC table:
  
  dicTable=data.frame(Date=as.character(tail(outputData$Date,1)),
                      Model=c("norm-norm",
                              "tStud-tStud",
                              "GED-GED",
                              "norm-tStud",
                              "norm-GED",
                              "tStud-GED",
                              "BestFitting"),
                      DIC=rep(Inf,7),
                      ModelID=c("c(norm,norm)",
                                  "c(std,std)",
                                  "c(ged,ged)",
                                  "c(norm,std)",
                                  "c(norm,ged)",
                                  "c(std,ged)",
                                  " "),
                      msSpec=c(paste0("MSspec", seq(1, 6)),""),
                      msMod=c(paste0("fitedMSGARCHD", seq(1, 6)),""),
                      GARCHSpec=rep(paste0("c('",paste(GARCHmodels,collapse="','"),"')"),7),
                      Experiment=experiment)
  
# Determines if the equation is a multifactor or a single one:
  MSspec1=CreateSpec(variance.spec = list(model = GARCHmodels), 
                    distribution.spec = list(distribution = c("norm","norm")),
                    switch.spec = list(do.mix = FALSE))
  
  MSspec2=CreateSpec(variance.spec = list(model = GARCHmodels), 
                     distribution.spec = list(distribution = c("std","std")),
                     switch.spec = list(do.mix = FALSE))
  
  MSspec3=CreateSpec(variance.spec = list(model = GARCHmodels), 
                     distribution.spec = list(distribution = c("ged","ged")),
                     switch.spec = list(do.mix = FALSE))  
  
  MSspec4=CreateSpec(variance.spec = list(model = GARCHmodels), 
                     distribution.spec = list(distribution = c("norm","std")),
                     switch.spec = list(do.mix = FALSE))  
  
  MSspec5=CreateSpec(variance.spec = list(model = GARCHmodels), 
                     distribution.spec = list(distribution = c("norm","ged")),
                     switch.spec = list(do.mix = FALSE))   
  
  MSspec6=CreateSpec(variance.spec = list(model = GARCHmodels), 
                     distribution.spec = list(distribution = c("std","ged")),
                     switch.spec = list(do.mix = FALSE))     

print(paste0("Estimating 2-regime ","normal","-","normal"," regime-specific LLFs MS-GARCH model"," (model 1 of 6). (",experiment,"-",estDate,"-",GARCHmodels,")"))  

  fitedMSGARCHD1 = tryCatch(FitMCMC(spec = MSspec1, data = residuals, 
                                    ctr=list(nburn=numberBurn,nmcmc=numberMCMC)) , 
                           error=function(e) NULL) 
  
  if (!is.null(fitedMSGARCHD1)){

    dicTable$DIC[1]=summary(fitedMSGARCHD1)$DIC
    meltedCoefsTable=melt(summary(fitedMSGARCHD1)$summary,id=c("Mean","SE"))
    
    DBTable=rbind(DBTable,
                  data.frame(Date=as.character(tail(outputData$Date,1)),
                       Value=meltedCoefsTable$value,
                       Ticker=paste0(meltedCoefsTable$Var1,"-",meltedCoefsTable$Var2),
                       ModelID=dicTable$Model[1],
                       GARCHSpec=dicTable$GARCHSpec[1],
                       Experiment=experiment
    )
    )
    
  }
  
  cat("\f")  
  
print(paste0("Estimating 2-regime ","tStud","-","tStud"," regime-specific LLFs MS-GARCH model"," (model 2 of 6). (",experiment,"-",estDate,"-",GARCHmodels,")"))  

  fitedMSGARCHD2 = tryCatch(FitMCMC(spec = MSspec2, data = residuals, 
                                    ctr=list(nburn=numberBurn,nmcmc=numberMCMC)) , 
                            error=function(e) NULL)   

  if (!is.null(fitedMSGARCHD2)){
    
    dicTable$DIC[2]=summary(fitedMSGARCHD2)$DIC
    meltedCoefsTable=melt(summary(fitedMSGARCHD2)$summary,id=c("Mean","SE"))
    
    DBTable=rbind(DBTable,
                  data.frame(Date=as.character(tail(outputData$Date,1)),
                  Value=meltedCoefsTable$value,
                  Ticker=paste0(meltedCoefsTable$Var1,"-",meltedCoefsTable$Var2),
                  ModelID=dicTable$Model[2],
                  GARCHSpec=dicTable$GARCHSpec[2],Experiment=experiment)
                  
                  )
    
  }  
 
  cat("\f")  
  
print(paste0("Estimating 2-regime ","GED","-","GED"," regime-specific LLFs MS-GARCH model"," (model 3 of 6). (",experiment,"-",estDate,"-",GARCHmodels,")"))  

  fitedMSGARCHD3 = tryCatch(FitMCMC(spec = MSspec3, data = residuals, 
                                    ctr=list(nburn=numberBurn,nmcmc=numberMCMC)) , 
                            error=function(e) NULL)  

  if (!is.null(fitedMSGARCHD3)){
    
    dicTable$DIC[3]=summary(fitedMSGARCHD3)$DIC
    meltedCoefsTable=melt(summary(fitedMSGARCHD3)$summary,id=c("Mean","SE"))
    
    DBTable=rbind(DBTable,
                  data.frame(Date=as.character(tail(outputData$Date,1)),
                             Value=meltedCoefsTable$value,
                             Ticker=paste0(meltedCoefsTable$Var1,"-",meltedCoefsTable$Var2),
                             ModelID=dicTable$Model[3],
                             GARCHSpec=dicTable$GARCHSpec[3],
                             Experiment=experiment)
                  
    )
    
  }
  
  cat("\f")  
  
print(paste0("Estimating 2-regime ","normal","-","tStud"," regime-specific LLFs MS-GARCH model"," (model 4 of 6). (",experiment,"-",estDate,"-",GARCHmodels,")"))  
  
  fitedMSGARCHD4 = tryCatch(FitMCMC(spec = MSspec4, data = residuals, 
                                    ctr=list(nburn=numberBurn,nmcmc=numberMCMC)) , 
                            error=function(e) NULL)    

  if (!is.null(fitedMSGARCHD4)){
    
    dicTable$DIC[4]=summary(fitedMSGARCHD4)$DIC
    meltedCoefsTable=melt(summary(fitedMSGARCHD4)$summary,id=c("Mean","SE"))
    
    DBTable=rbind(DBTable,
                  data.frame(Date=as.character(tail(outputData$Date,1)),
                             Value=meltedCoefsTable$value,
                             Ticker=paste0(meltedCoefsTable$Var1,"-",meltedCoefsTable$Var2),
                             ModelID=dicTable$Model[4],
                             GARCHSpec=dicTable$GARCHSpec[4],Experiment=experiment)
                  
    )
    
  }
  
  cat("\f")  
  
print(paste0("Estimating 2-regime ","normal","-","GED"," regime-specific LLFs MS-GARCH model"," (model 5 of 6). (",experiment,"-",estDate,"-",GARCHmodels,")"))

  fitedMSGARCHD5 = tryCatch(FitMCMC(spec = MSspec5, data = residuals, 
                                    ctr=list(nburn=numberBurn,nmcmc=numberMCMC)) , 
                            error=function(e) NULL)  

  if (!is.null(fitedMSGARCHD5)){
    
    dicTable$DIC[5]=summary(fitedMSGARCHD5)$DIC
    meltedCoefsTable=melt(summary(fitedMSGARCHD5)$summary,id=c("Mean","SE"))
    
    DBTable=rbind(DBTable,
                  data.frame(Date=as.character(tail(outputData$Date,1)),
                             Value=meltedCoefsTable$value,
                             Ticker=paste0(meltedCoefsTable$Var1,"-",meltedCoefsTable$Var2),
                             ModelID=dicTable$Model[5],
                             GARCHSpec=dicTable$GARCHSpec[5],Experiment=experiment)
    )
    
  }
  
  cat("\f")  

print(paste0("Estimating 2-regime ","tStud","-","GED"," regime-specific LLFs MS-GARCH model"," (model 6 of 6). (",experiment,"-",estDate,"-",GARCHmodels,")"))      

  fitedMSGARCHD6 = tryCatch(FitMCMC(spec = MSspec6, data = residuals, 
                                    ctr=list(nburn=numberBurn,nmcmc=numberMCMC)) , 
                            error=function(e) NULL)    
  
  if (!is.null(fitedMSGARCHD6)){
    
    dicTable$DIC[6]=summary(fitedMSGARCHD6)$DIC
    meltedCoefsTable=melt(summary(fitedMSGARCHD6)$summary,id=c("Mean","SE"))
    
    DBTable=rbind(DBTable,
                  data.frame(Date=as.character(tail(outputData$Date,1)),
                             Value=meltedCoefsTable$value,
                             Ticker=paste0(meltedCoefsTable$Var1,"-",meltedCoefsTable$Var2),
                             ModelID=dicTable$Model[6],
                             GARCHSpec=dicTable$GARCHSpec[6],
                             Experiment=experiment)
    )
    
  }  
  
  cat("\f")  
  
print("Determining the best fitting model...") 

bestModelRow=which(dicTable$DIC==min(dicTable$DIC))  
  
dicTable$DIC[7]=dicTable$DIC[bestModelRow]
  
dicTable$ModelID[7]=dicTable$ModelID[bestModelRow]  
  
dicTable$msSpec[7]=dicTable$msSpec[bestModelRow]
  
dicTable$msMod[7]=dicTable$msMod[bestModelRow]

#dicTable$Model[7]=paste0(dicTable$Model[7]," (",dicTable$Model[bestModelRow],")")

dicTableDB=cbind(data.frame(Date=tail(outputData$Date,1)),dicTable)

# Historiacl best fitting model DB table:

eval(parse(text=paste0("fitedMSGARCHDBest=fitedMSGARCHD",bestModelRow)))

if (!is.null(fitedMSGARCHDBest)){
  bestFittingModelId=paste0("fitedMSGARCHD",bestModelRow)
  
  eval(parse(text=paste0("meltedCoefsTable=melt(summary(",
  bestFittingModelId,
  ")$summary,id=c('Mean','SE'))")))
  
  DBTable=rbind(DBTable,
                data.frame(Date=as.character(tail(outputData$Date,1)),
                           Value=meltedCoefsTable$value,
                           Ticker=paste0(meltedCoefsTable$Var1,"-",meltedCoefsTable$Var2),
                           ModelID=dicTable$Model[bestModelRow],
                           GARCHSpec=dicTable$GARCHSpec[bestModelRow],
                           Experiment=experiment)
  )
  
}  

# Generating the Smoothed, transtition and forecasted probabilities for model 1:

cat("\f")  

print(paste0("Estimating regime-specific smoothed probs. ","norm","-","norm"," (model 1 of 6). (",experiment,"-",estDate,"-",GARCHmodels,")")) 

# Smoothed probabilites:
Smooth.probs1 = State(fitedMSGARCHD1)$SmoothProb[,1, 1:MSspec1$K, drop = TRUE]

# Transition probability matrix:
transprob1=summary(fitedMSGARCHD1)$post.trans.mat

# Forecasted prpbabilities:

Predprobs1 = rbind(Smooth.probs1[nrow(Smooth.probs1),]%*%transprob1,
                  Smooth.probs1[nrow(Smooth.probs1),]%*%transprob1^2,
                  Smooth.probs1[nrow(Smooth.probs1),]%*%transprob1^3,
                  Smooth.probs1[nrow(Smooth.probs1),]%*%transprob1^4,
                  Smooth.probs1[nrow(Smooth.probs1),]%*%transprob1^5)

# smooth probs Data Base table:
DBTable=rbind(DBTable,
            data.frame(Date=as.character(tail(outputData$Date,1)),
           Value=c(Predprobs1[,1],Predprobs1[,2]),
           Ticker=c(paste0("Calm prob. forecast t+",seq(1,5)),
                    paste0("Crisis prob. forecast t+",seq(1,5))),
           ModelID=dicTable$Model[1],
           GARCHSpec=dicTable$GARCHSpec[1],
           Experiment=experiment
           )
)

# Adding smooth probs to output data:

outputData$NormNormSmoothProbR1=Smooth.probs1[2:nrow(Smooth.probs1),1]
outputData$NormNormSmoothProbR2=Smooth.probs1[2:nrow(Smooth.probs1),2]

# Creating the regime-specific chart data:

#regime1ChartDataPrice=data.frame(Date=Data$Date,
          #                       Value=Data$Settle,
           #                      Ticker="Historical price")

regime1ChartDataPrice=data.frame(Date=Data$Date,
                                 Value=Smooth.probs1[2:nrow(Smooth.probs1),1],
                                 Ticker=paste0(dicTable$Model[1]," LLFs"))

regime2ChartDataPrice=data.frame(Date=Data$Date,
                                       Value=Smooth.probs1[2:nrow(Smooth.probs1),2],
                                       Ticker=paste0(dicTable$Model[1]," LLFs"))


cat("\f")  

print(paste0("Estimating regime-specific smoothed probs. ","tStud","-","tSud"," (model 2 of 6). (",experiment,"-",estDate,"-",GARCHmodels,")")) 

# Smoothed, transtition and forecasted probabilities for model 2:

# Smoothed probabilites:
Smooth.probs2 = State(fitedMSGARCHD2)$SmoothProb[,1, 1:MSspec2$K, drop = TRUE]

# Transition probability matrix:
transprob2 = summary(fitedMSGARCHD2)$post.trans.mat

# Forecasted prpbabilities:

Predprobs2 = rbind(Smooth.probs2[nrow(Smooth.probs2),]%*%transprob2,
                   Smooth.probs2[nrow(Smooth.probs2),]%*%transprob2^2,
                   Smooth.probs2[nrow(Smooth.probs2),]%*%transprob2^3,
                   Smooth.probs2[nrow(Smooth.probs2),]%*%transprob2^4,
                   Smooth.probs2[nrow(Smooth.probs2),]%*%transprob2^5)

DBTable=rbind(DBTable,
              data.frame(Date=as.character(tail(outputData$Date,1)),
                   Value=c(Predprobs2[,1],Predprobs2[,2]),
                   Ticker=c(paste0("Calm prob. forecast t+",seq(1,5)),
                            paste0("Crisis prob. forecast t+",seq(1,5))),
                   ModelID=dicTable$Model[2],
                   GARCHSpec=dicTable$GARCHSpec[2],
                   Experiment=experiment
)
)
# Adding smooth probs to output data:

outputData$tStudtStudSmoothProbR1=Smooth.probs2[2:nrow(Smooth.probs2),1]
outputData$tStudtStudSmoothProbR2=Smooth.probs2[2:nrow(Smooth.probs2),2]

# Creating the regime-specific chart data:

regime1ChartDataPrice=rbind(regime1ChartDataPrice,
                            data.frame(Date=Data$Date,
                                       Value=Smooth.probs2[2:nrow(Smooth.probs2),1],
                                       Ticker=paste0(dicTable$Model[2]," LLFs"))
)

regime2ChartDataPrice=rbind(regime2ChartDataPrice,
                            data.frame(Date=Data$Date,
                                       Value=Smooth.probs2[2:nrow(Smooth.probs2),2],
                                       Ticker=paste0(dicTable$Model[2]," LLFs"))
)

cat("\f")  

print(paste0("Estimating regime-specific smoothed probs. ","GED","-","GED"," (model 3 of 6). (",experiment,"-",estDate,"-",GARCHmodels,")")) 

# Smoothed, transtition and forecasted probabilities for model 3:

# Smoothed probabilites:
Smooth.probs3 = State(fitedMSGARCHD3)$SmoothProb[,1, 1:MSspec3$K, drop = TRUE]

# Transition probability matrix:
transprob3 = summary(fitedMSGARCHD3)$post.trans.mat

# Forecasted prpbabilities:

Predprobs3 = rbind(Smooth.probs3[nrow(Smooth.probs3),]%*%transprob3,
                   Smooth.probs3[nrow(Smooth.probs3),]%*%transprob3^2,
                   Smooth.probs3[nrow(Smooth.probs3),]%*%transprob3^3,
                   Smooth.probs3[nrow(Smooth.probs3),]%*%transprob3^4,
                   Smooth.probs3[nrow(Smooth.probs3),]%*%transprob3^5)

DBTable=rbind(DBTable,
              data.frame(Date=as.character(tail(outputData$Date,1)),
                         Value=c(Predprobs3[,1],Predprobs3[,2]),
                         Ticker=c(paste0("Calm prob. forecast t+",seq(1,5)),
                                  paste0("Crisis prob. forecast t+",seq(1,5))),
                         ModelID=dicTable$Model[3],
                         GARCHSpec=dicTable$GARCHSpec[3],
                         Experiment=experiment
              )
)

# Adding smooth probs to output data:

outputData$GEDGEDSmoothProbR1=Smooth.probs3[2:nrow(Smooth.probs3),1]
outputData$GEDGEDSmoothProbR2=Smooth.probs3[2:nrow(Smooth.probs3),2]

# chart data tables update:
regime1ChartDataPrice=rbind(regime1ChartDataPrice,
                            data.frame(Date=Data$Date,
                                       Value=Smooth.probs3[2:nrow(Smooth.probs3),1],
                                       Ticker=paste0(dicTable$Model[3]," LLFs"))
)

regime2ChartDataPrice=rbind(regime2ChartDataPrice,
                            data.frame(Date=Data$Date,
                                       Value=Smooth.probs3[2:nrow(Smooth.probs3),2],
                                       Ticker=paste0(dicTable$Model[3]," LLFs"))
)

cat("\f")  

print(paste0("Estimating regime-specific smoothed probs. ","norm","-","tStud"," (model 4 of 6). (",experiment,"-",estDate,"-",GARCHmodels,")")) 


# Smoothed, transtition and forecasted probabilities for model 4:

# Smoothed probabilites:
Smooth.probs4 = State(fitedMSGARCHD4)$SmoothProb[,1, 1:MSspec4$K, drop = TRUE]

# Transition probability matrix:
transprob4 = summary(fitedMSGARCHD4)$post.trans.mat

# Forecasted prpbabilities:

Predprobs4 = rbind(Smooth.probs4[nrow(Smooth.probs4),]%*%transprob4,
                   Smooth.probs4[nrow(Smooth.probs4),]%*%transprob4^2,
                   Smooth.probs4[nrow(Smooth.probs4),]%*%transprob4^3,
                   Smooth.probs4[nrow(Smooth.probs4),]%*%transprob4^4,
                   Smooth.probs4[nrow(Smooth.probs4),]%*%transprob4^5)

DBTable=rbind(DBTable,
              data.frame(Date=as.character(tail(outputData$Date,1)),
                         Value=c(Predprobs4[,1],Predprobs4[,2]),
                         Ticker=c(paste0("Calm prob. forecast t+",seq(1,5)),
                                  paste0("Crisis prob. forecast t+",seq(1,5))),
                         ModelID=dicTable$Model[4],
                         GARCHSpec=dicTable$GARCHSpec[4],
                         Experiment=experiment
              )
)
# Adding smooth probs to output data:

outputData$NormtStudSmoothProbR1=Smooth.probs4[2:nrow(Smooth.probs4),1]
outputData$NormtStudSmoothProbR2=Smooth.probs4[2:nrow(Smooth.probs4),2]

# chart data tables update:
regime1ChartDataPrice=rbind(regime1ChartDataPrice,
                            data.frame(Date=Data$Date,
                                       Value=Smooth.probs4[2:nrow(Smooth.probs4),1],
                                       Ticker=paste0(dicTable$Model[4]," LLFs"))
)

regime2ChartDataPrice=rbind(regime2ChartDataPrice,
                            data.frame(Date=Data$Date,
                                       Value=Smooth.probs4[2:nrow(Smooth.probs4),2],
                                       Ticker=paste0(dicTable$Model[4]," LLFs"))
)

cat("\f")  

print(paste0("Estimating regime-specific smoothed probs. ","norm","-","GED"," (model 5 of 6). (",experiment,"-",estDate,"-",GARCHmodels,")")) 


# Smoothed, transtition and forecasted probabilities for model 5:

# Smoothed probabilites:
Smooth.probs5 = State(fitedMSGARCHD5)$SmoothProb[,1, 1:MSspec5$K, drop = TRUE]

# Transition probability matrix:
transprob5 = summary(fitedMSGARCHD5)$post.trans.mat

# Forecasted prpbabilities:

Predprobs5 = rbind(Smooth.probs5[nrow(Smooth.probs5),]%*%transprob5,
                   Smooth.probs5[nrow(Smooth.probs5),]%*%transprob5^2,
                   Smooth.probs5[nrow(Smooth.probs5),]%*%transprob5^3,
                   Smooth.probs5[nrow(Smooth.probs5),]%*%transprob5^4,
                   Smooth.probs5[nrow(Smooth.probs5),]%*%transprob5^5)

DBTable=rbind(DBTable,
              data.frame(Date=as.character(tail(outputData$Date,1)),
                         Value=c(Predprobs5[,1],Predprobs5[,2]),
                         Ticker=c(paste0("Calm prob. forecast t+",seq(1,5)),
                                  paste0("Crisis prob. forecast t+",seq(1,5))),
                         ModelID=dicTable$Model[5],
                         GARCHSpec=dicTable$GARCHSpec[5],
                         Experiment=experiment
              )
)

# Adding smooth probs to output data:

outputData$NormGEDSmoothProbR1=Smooth.probs5[2:nrow(Smooth.probs5),1]
outputData$NormGEDSmoothProbR2=Smooth.probs5[2:nrow(Smooth.probs5),2]

# chart data tables update:
regime1ChartDataPrice=rbind(regime1ChartDataPrice,
                            data.frame(Date=Data$Date,
                                       Value=Smooth.probs5[2:nrow(Smooth.probs5),1],
                                       Ticker=paste0(dicTable$Model[5]," LLFs"))
)

regime2ChartDataPrice=rbind(regime2ChartDataPrice,
                            data.frame(Date=Data$Date,
                                       Value=Smooth.probs5[2:nrow(Smooth.probs5),2],
                                       Ticker=paste0(dicTable$Model[5]," LLFs"))
)

cat("\f")  

print(paste0("Estimating regime-specific smoothed probs. ","tStud","-","GED"," (model 6 of 6). (",experiment,"-",estDate,"-",GARCHmodels,")")) 

# Smoothed, transtition and forecasted probabilities for model 6:

# Smoothed probabilites:
Smooth.probs6 = State(fitedMSGARCHD6)$SmoothProb[,1, 1:MSspec6$K, drop = TRUE]

# Transition probability matrix:
transprob6 = summary(fitedMSGARCHD6)$post.trans.mat

# Forecasted prpbabilities:

Predprobs6 = rbind(Smooth.probs6[nrow(Smooth.probs6),]%*%transprob6,
                   Smooth.probs6[nrow(Smooth.probs6),]%*%transprob6^2,
                   Smooth.probs6[nrow(Smooth.probs6),]%*%transprob6^3,
                   Smooth.probs6[nrow(Smooth.probs6),]%*%transprob6^4,
                   Smooth.probs6[nrow(Smooth.probs6),]%*%transprob6^5)

DBTable=rbind(DBTable,
              data.frame(Date=as.character(tail(outputData$Date,1)),
                         Value=c(Predprobs6[,1],Predprobs6[,2]),
                         Ticker=c(paste0("Calm prob. forecast t+",seq(1,5)),
                                  paste0("Crisis prob. forecast t+",seq(1,5))),
                         ModelID=dicTable$Model[6],
                         GARCHSpec=dicTable$GARCHSpec[6],
                         Experiment=experiment
              )
)

# Adding smooth probs to output data:

outputData$tStudGEDSmoothProbR1=Smooth.probs6[2:nrow(Smooth.probs6),1]
outputData$tStudGEDSmoothProbR2=Smooth.probs6[2:nrow(Smooth.probs6),2]

# chart data tables update:
regime1ChartDataPrice=rbind(regime1ChartDataPrice,
                            data.frame(Date=Data$Date,
                                       Value=Smooth.probs6[2:nrow(Smooth.probs6),1],
                                       Ticker=paste0(dicTable$Model[6]," LLFs"))
)

regime2ChartDataPrice=rbind(regime2ChartDataPrice,
                            data.frame(Date=Data$Date,
                                       Value=Smooth.probs6[2:nrow(Smooth.probs6),2],
                                       Ticker=paste0(dicTable$Model[6]," LLFs"))
)

cat("\f")  

print(paste0("Estimating regime-specific smoothed probs. ","(best fiting model)")) 

msBestmodel=dicTable$msMod[7]

# Smoothed, transtition and forecasted probabilities for best model:

# Smoothed probabilites:
eval(parse(text=paste0("Smooth.probsBest = State(",msBestmodel,")$SmoothProb[,1, 1:MSspec5$K, drop = TRUE]")))

# Transition probability matrix:
eval(parse(text=paste0("transprobBest = summary(",msBestmodel,")$post.trans.mat")))

# Forecasted prpbabilities:

PredprobsBest = rbind(Smooth.probsBest[nrow(Smooth.probsBest),]%*%transprobBest,
                      Smooth.probsBest[nrow(Smooth.probsBest),]%*%transprobBest^2,
                      Smooth.probsBest[nrow(Smooth.probsBest),]%*%transprobBest^3,
                      Smooth.probsBest[nrow(Smooth.probsBest),]%*%transprobBest^4,
                      Smooth.probsBest[nrow(Smooth.probsBest),]%*%transprobBest^5)

# Adding smooth probs to output data:

outputData$BestModelSmoothProbR1=Smooth.probsBest[2:nrow(Smooth.probsBest),1]
outputData$BestModelSmoothProbR2=Smooth.probsBest[2:nrow(Smooth.probsBest),2]

# chart data tables update:
regime1ChartDataPrice=rbind(regime1ChartDataPrice,
                            data.frame(Date=Data$Date,
                                       Value=Smooth.probsBest[2:nrow(Smooth.probsBest),1],
                                       Ticker=paste0(dicTable$Model[7]," (",
                                                     dicTable$Model[bestModelRow],")"))
)

regime2ChartDataPrice=rbind(regime2ChartDataPrice,
                            data.frame(Date=Data$Date,
                                       Value=Smooth.probsBest[2:nrow(Smooth.probsBest),2],
                                       Ticker=paste0(dicTable$Model[7]," (",
                                                     dicTable$Model[bestModelRow],")"))
)

# Estimates forecasted volatility and VaR at t:

cat("\f")
print(paste0("Estimating expected volatility at t (model 1). (",experiment,"-",estDate,")"))

pred1 <- predict(fitedMSGARCHD1, nahead = 5, do.return.draws = FALSE)


DBTable=rbind(DBTable,
              data.frame(
                Date=as.character(tail(outputData$Date,1)),
                Value=pred1$vol,
                Ticker=paste0(paste0("Volatility at t+",seq(1:5))," model 1"),
                ModelID=dicTable$Model[1],
                GARCHSpec=dicTable$GARCHSpec[1],
                Experiment=experiment
              )           
)

print(paste0("Estimating expected volatility at t (model 2). (",experiment,"-",estDate,"-",GARCHmodels,")"))

pred2 <- predict(fitedMSGARCHD2, nahead = 5, do.return.draws = FALSE)


DBTable=rbind(DBTable,
              data.frame(
                Date=as.character(tail(outputData$Date,1)),
                Value=pred2$vol,
                Ticker=paste0(paste0("Volatility at t+",seq(1:5))," model 2"),
                ModelID=dicTable$Model[2],
                GARCHSpec=dicTable$GARCHSpec[2],
                Experiment=experiment
              )           
)

print(paste0("Estimating expected volatility at t (model 3). (",experiment,"-",estDate,"-",GARCHmodels,")"))

pred3 <- predict(fitedMSGARCHD3, nahead = 5, do.return.draws = FALSE)


DBTable=rbind(DBTable,
              data.frame(
                Date=as.character(tail(outputData$Date,1)),
                Value=pred3$vol,
                Ticker=paste0(paste0("Volatility at t+",seq(1:5))," model 3"),
                ModelID=dicTable$Model[3],
                GARCHSpec=dicTable$GARCHSpec[3],
                Experiment=experiment
              )           
)

print(paste0("Estimating expected volatility at t (model 4). (",experiment,"-",estDate,"-",GARCHmodels,")"))

pred4 <- predict(fitedMSGARCHD4, nahead = 5, do.return.draws = FALSE)


DBTable=rbind(DBTable,
              data.frame(
                Date=as.character(tail(outputData$Date,1)),
                Value=pred4$vol,
                Ticker=paste0(paste0("Volatility at t+",seq(1:5))," model 4"),
                ModelID=dicTable$Model[4],
                GARCHSpec=dicTable$GARCHSpec[4],
                Experiment=experiment
              )           
)

print(paste0("Estimating expected volatility at t (model 5). (",experiment,"-",estDate,"-",GARCHmodels,")"))

pred5 <- predict(fitedMSGARCHD5, nahead = 5, do.return.draws = FALSE)


DBTable=rbind(DBTable,
              data.frame(
                Date=as.character(tail(outputData$Date,1)),
                Value=pred5$vol,
                Ticker=paste0(paste0("Volatility at t+",seq(1:5))," model 5"),
                ModelID=dicTable$Model[5],
                GARCHSpec=dicTable$GARCHSpec[5],
                Experiment=experiment
              )           
)

print(paste0("Estimating expected volatility at t (model 6). (",experiment,"-",estDate,"-",GARCHmodels,")"))

pred6 <- predict(fitedMSGARCHD6, nahead = 5, do.return.draws = FALSE)


DBTable=rbind(DBTable,
              data.frame(
                Date=as.character(tail(outputData$Date,1)),
                Value=pred1$vol,
                Ticker=paste0(paste0("Volatility at t+",seq(1:5))," model 6"),
                ModelID=dicTable$Model[6],
                GARCHSpec=dicTable$GARCHSpec[6],
                Experiment=experiment
              )           
)

print(paste0("Estimating expected volatility at t (Best model). (",experiment,"-",estDate,"-",GARCHmodels,")"))

predBest <- predict(fitedMSGARCHDBest, nahead = 5, do.return.draws = FALSE)


DBTable=rbind(DBTable,
              data.frame(
                Date=as.character(tail(outputData$Date,1)),
                Value=predBest$vol,
                Ticker=paste0(paste0("Volatility at t+",seq(1:5))," Best model"),
                ModelID=dicTable$Model[7],
                GARCHSpec=dicTable$GARCHSpec[7],
                Experiment=experiment
              )           
)

# Closes the code:  
  tiempoPasado=as.numeric(Sys.time()-startCalculation)
  
  elapsTimeMsg=paste0("Estimation elapsed time: ",
                      round(tiempoPasado/360,0),
                      ":",
                      round(tiempoPasado/60,0),":",
                      round(tiempoPasado,3))
  
  MSGARCHResults=list(
    outputData=outputData,
    dicDBTable=dicTable,
    dicTableDB=dicTableDB,
    regime1ChartSmoothProbs=regime1ChartDataPrice,
    regime2ChartSmoothProbs=regime2ChartDataPrice,
    probsDBtable=DBTable,
    elapsTimeMsg=elapsTimeMsg,
    elapsTime=tiempoPasado
  )
  
  cat("\f")
  print(elapsTimeMsg)
  
  return(MSGARCHResults)

}
