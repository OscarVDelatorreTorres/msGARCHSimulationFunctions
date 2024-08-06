if (!require(fGarch)){install.packages("fGarch") 
  require(fGarch)} else {require(fGarch)}

if (!require(DescTools)){install.packages("DescTools") 
  require(DescTools)} else {require(DescTools)}

if (!require(MASS)){install.packages("MASS") 
  require(MASS)} else {require(MASS)}

if (!require(data.table)){install.packages("data.table") 
  require(data.table)} else {require(data.table)}

if (!require(openxlsx)){install.packages("openxlsx") 
  require(openxlsx)} else {require(openxlsx)}

if (!require(MSGARCH)){install.packages("MSGARCH") 
  require(MSGARCH)} else {require(MSGARCH)}

if (!require(zoo)){install.packages("zoo") 
  require(zoo)} else {require(zoo)}

emMCMCMSGARCHFit=function(eq,data,numberMCMC=10000,numberBurn=500,GARCHmodels=c("sGARCH","sGARCH"),pdfFunct=c("norm","norm"),experiment="",timeFixed=FALSE){
  # Starts calculation time:
  startCalculation=Sys.time()
  
  words <- strsplit(eq, "[^[:alnum:]]+") # Split the string using non-alphanumeric characters as separators
  words = words[[1]]
  words=words[-1]
  num_words <- length(words)
  
  # Column with the dependent variable in the data.frame:
  depVar=strsplit(eq,"~")[[1]][1]
  depVarPos=which(colnames(data)==depVar)
  
  variables=c("(Intercept)",words)
  variablesValues=rep(0,num_words)
  
  # Fit the full model
  full.model <- lm(eq, data = data)
  meanEqCoefsTable=summary(full.model)$coefficients
  meanEqR2=summary(full.model)$adj.r.squared
  maenEqSigma=summary(full.model)$sigma
  # Stepwise regression model
  #step.model <- stepAIC(full.model, direction = "forward", trace = F)
  step.model=full.model
  coeficientesStepWise=summary(step.model)$coefficients
  
  num_VarsStep=nrow(coeficientesStepWise)
  name_VarsStep=rownames(coeficientesStepWise)
  
  # Coefficients tbale in DBTable
  DBTable=data.frame(
    Date=as.character(tail(data$Date,1)),
    Value=step.model$coefficients,
    Ticker="factor model coefs",
    Experiment=experiment
  )
  
  # Forecast in DBTable:
  
  meanForecast=tail(step.model$fitted.values,1)
  
  for (numVar in 1:nrow(DBTable)){
    
    idCoefRow=which(name_VarsStep[numVar]==DBTable$Ticker)
    if (length(idCoefRow)>0){
      DBTable$Value[idCoefRow]=coeficientesStepWise[numVar]
    }
    
  }
  
  if (tail(data[depVarPos],1)>0){
    returnScenarioText="Bullish"
    returnScenario=1
  } else {
    returnScenarioText="Bearish"
    returnScenario=-1
  }
  
  if (meanForecast>0){
    expectedReturnScenarioText="Bullish"
    expectedReturnScenario=1
  } else {
    expectedReturnScenarioText="Bearish"
    expectedReturnScenario=-1
  }
  
  DBTable=rbind(DBTable,
                data.frame(
                  Date=as.character(tail(data$Date,1)),
                  Value=as.numeric(tail(data[depVarPos],1)),
                  Ticker="Return at t",
                  Experiment=experiment
                ),
                data.frame(
                  Date=as.character(tail(data$Date,1)),
                  Value=meanForecast,
                  Ticker="Return Forecast at t",
                  Experiment=experiment
                ),
                data.frame(
                  Date=as.character(tail(data$Date,1)),
                  Value=returnScenario,
                  Ticker="Return scenario",
                  Experiment=experiment
                ),
                data.frame(
                  Date=as.character(tail(data$Date,1)),
                  Value=returnScenarioText,
                  Ticker="Return scenario text",
                  Experiment=experiment
                ),
                data.frame(
                  Date=as.character(tail(data$Date,1)),
                  Value=expectedReturnScenario,
                  Ticker="Expected return scenario",
                  Experiment=experiment
                ),
                data.frame(
                  Date=as.character(tail(data$Date,1)),
                  Value=expectedReturnScenarioText,
                  Ticker="Expected return scenario text",
                  Experiment=experiment
                )
  )
  
# Estimates the E-M MS-GARCH model as a first step:====
  
  # Extract the residuals from the model:
  residuals=step.model$residuals
  
  # Starts calculation time:
  startCalculation=Sys.time()
  
  # Determines if the equation is a multifactor or a single one:
  
  if (isTRUE(timeFixed)){
    
    MSspec=CreateSpec(variance.spec = list(model = c("sARCH","sARCH")),
                      distribution.spec = list(distribution = pdfFunct),
                      switch.spec = list(do.mix = FALSE),
                      constraint.spec = list(fixed = list(alpha1_1 = 1e-06,
                                                          alpha1_2 = 1e-06)))
  } else {
    
    MSspec=CreateSpec(variance.spec = list(model = GARCHmodels),
                      distribution.spec = list(distribution = pdfFunct),
                      switch.spec = list(do.mix = FALSE))
    
  }
  
  # Estimates the MSGARH E-M model:
  cat("\f")
print(paste0("Estimating MS-GARCH (EM method) model for experiment ", experiment,", date: ",tail(data$Date,1)))
  
# ML estimation:
  simmethod="E-M"
  fittedMSGARCHD = tryCatch(FitML(spec = MSspec, data = residuals) ,
                            error=function(e) NULL)
  
  if (!is.null(fittedMSGARCHD)){
    
    meltedCoefsTable=melt(fittedMSGARCHD$Inference$MatCoef)
    msGARCHOutPutCoefsTable=fittedMSGARCHD$Inference$MatCoef
    # Parameters coefficients:
    DBTable=rbind(DBTable,
                  data.frame(Date=as.character(tail(data$Date,1)),
                             Value=meltedCoefsTable$value,
                             Ticker=paste0(meltedCoefsTable$Var1,"-",meltedCoefsTable$Var2),
                             Experiment=experiment
                  )
    )
    # AIC values:
    DBTable=rbind(DBTable,
                  data.frame(Date=as.character(tail(data$Date,1)),
                             Value=summary(fittedMSGARCHD)$AIC,
                             Ticker="AIC-ML estimation",
                             Experiment=experiment
                  )
    )
    # estimation method
    DBTable=rbind(DBTable,
                  data.frame(Date=as.character(tail(data$Date,1)),
                             Value="EM",
                             Ticker="MSGARCH estimation method",
                             Experiment=experiment
                  )
    )
    
#Smoothed probabilities and forecasts:
    # Smoothed probabilites:
    stateMS=State(fittedMSGARCHD)
    
    Smooth.probs1 = stateMS$SmoothProb[,1, 1:MSspec$K, drop = TRUE]
    
    # Transition probability matrix:
    
    if (!is.null(fittedMSGARCHD)){
      transprob1=summary(fittedMSGARCHD)$trans.mat
    } else {
      transprob1=summary(fittedMSGARCHD)$post.trans.mat
    }
    
    
# Forecasted probabilities:
    
    Predprobs1 = rbind(Smooth.probs1[nrow(Smooth.probs1),]%*%transprob1,
                       Smooth.probs1[nrow(Smooth.probs1),]%*%transprob1^2,
                       Smooth.probs1[nrow(Smooth.probs1),]%*%transprob1^3,
                       Smooth.probs1[nrow(Smooth.probs1),]%*%transprob1^4,
                       Smooth.probs1[nrow(Smooth.probs1),]%*%transprob1^5)
    
# smooth probs data Base table:
    DBTable=rbind(DBTable,
                  data.frame(Date=as.character(tail(data$Date,1)),
                             Value=c(Predprobs1[,1],Predprobs1[,2]),
                             Ticker=c(paste0("Calm prob. forecast t+",seq(1,5)),
                                      paste0("Crisis prob. forecast t+",seq(1,5))),
                             Experiment=experiment
                  )
    )
    
#  Actual return and actual smooth probability regime:
    
    if (tail(data[,depVarPos],1)>0){
      if (tail(Smooth.probs1[,2],1)<=0.5){
        msTradingScenarioActual=1
        msTradingScenarioTextActual="normal bullish"
        msVolatilityScenarioTextActual="normal"
        msVolatilityScenarioActual=1
      } else {
        msTradingScenarioActual=2
        msTradingScenarioTextActual="crisis bullish"
        msVolatilityScenarioTextActual="crisis"
        msVolatilityScenarioActual=2
      }
    } else {
      if (tail(Smooth.probs1[,2],1)<=0.5){
        msTradingScenarioActual=3
        msTradingScenarioTextActual="normal bearish"
        msVolatilityScenarioTextActual="normal"
        msVolatilityScenarioActual=1
      } else {
        msTradingScenarioActual=4
        msTradingScenarioTextActual="crisis bearish"
        msVolatilityScenarioTextActual="crisis"
        msVolatilityScenarioActual=2
      }
    }
    
# Expected return and actual smooth probability forecasted regimes:
    
    if (expectedReturnScenario>0){
      if (tail(Smooth.probs1[,2],1)<=0.5){
        msTradingScenarioF1=1
        msTradingScenarioTextF1="normal bullish"
        msVolatilityScenarioTextF1="normal"
        msVolatilityScenarioF1=1
      } else {
        msTradingScenarioF1=2
        msTradingScenarioTextF1="crisis bullish"
        msVolatilityScenarioTextF1="crisis"
        msVolatilityScenarioF1=2
      }
    } else {
      if (tail(Smooth.probs1[,2],1)<=0.5){
        msTradingScenarioF1=3
        msTradingScenarioTextF1="normal bearish"
        msVolatilityScenarioTextF1="normal"
        msVolatilityScenarioF1=1
      } else {
        msTradingScenarioF1=4
        msTradingScenarioTextF1="crisis bearish"
        msVolatilityScenarioTextF1="crisis"
        msVolatilityScenarioF1=2
      }
    }
    
# Expected return  and t+1 forecasted smooth probability forecasted regimes:
    
    if (expectedReturnScenario>0){
      if (Predprobs1[1,2]<=0.5){
        msTradingScenarioF2=1
        msTradingScenarioTextF2="normal bullish"
        msVolatilityScenarioTextF2="normal"
        msVolatilityScenarioF2=1
      } else {
        msTradingScenarioF2=2
        msTradingScenarioTextF2="crisis bullish"
        msVolatilityScenarioTextF2="crisis"
        msVolatilityScenarioF2=2
      }
    } else {
      if (Predprobs1[1,2]<=0.5){
        msTradingScenarioF2=3
        msTradingScenarioTextF2="normal bearish"
        msVolatilityScenarioTextF2="normal"
        msVolatilityScenarioF2=1
      } else {
        msTradingScenarioF2=4
        msTradingScenarioTextF2="crisis bearish"
        msVolatilityScenarioTextF2="crisis"
        msVolatilityScenarioF2=2
      }
    }
    
    DBTable=rbind(DBTable,
                  data.frame(
                    Date=as.character(tail(data$Date,1)),
                    Value=c(msTradingScenarioActual,
                            msTradingScenarioTextActual,
                            msVolatilityScenarioActual,
                            msVolatilityScenarioTextActual),
                    Ticker=c("MS trading scenario Actual",
                             "MS trading scenario text Actual",
                             "MS volatility scenario Actual",
                             "MS volatility scenario text Actual"),
                    Experiment=experiment
                  ),
                  
                  data.frame(
                    Date=as.character(tail(data$Date,1)),
                    Value=c(msTradingScenarioF1,
                            msTradingScenarioTextF1,
                            msVolatilityScenarioF1,
                            msVolatilityScenarioTextF1),
                    Ticker=c("MS trading scenario Actual vol forecast",
                             "MS trading scenario text Actual vol forecast",
                             "MS volatility scenario Actual vol forecast",
                             "MS volatility scenario text Actual vol forecast"),
                    Experiment=experiment
                  ),
                  
                  data.frame(
                    Date=as.character(tail(data$Date,1)),
                    Value=c(msTradingScenarioF2,
                            msTradingScenarioTextF2,
                            msVolatilityScenarioF2,
                            msVolatilityScenarioTextF2),
                    Ticker=c("MS trading scenario all Forecast",
                             "MS trading scenario text all Forecast",
                             "MS volatility scenario all Forecast",
                             "MS volatility scenario text all Forecast"),
                    Experiment=experiment
                  )
    )
    
# Estimates forecasted volatility and VaR at t:
    cat("\f")   
print(paste0("Estimating expected volatility at t. (",experiment,"-",GARCHmodels,")"))
    
    pred1 <- predict(fittedMSGARCHD, nahead = 5, do.return.draws = FALSE)
    checkNARisk=which(is.na(pred1$vol))
    
    
    if (length(checkNARisk)>0){
      pred1$vol=na.locf(pred1$vol)
    }
    
    DBTable=rbind(DBTable,
                  data.frame(
                    Date=as.character(tail(data$Date,1)),
                    Value=pred1$vol,
                    Ticker=paste0(paste0("Volatility at t+",seq(1:5))," model 1"),
                    Experiment=experiment
                  )
    )
    
# general model LLF and AIC estimation from data:
    
    sigmaLLF=pred1$vol[1]
    muLLF=meanForecast
    
    llfFunct=pdfFunct[1]
    
    switch(llfFunct,
           "norm"={
             kMS=sum(fittedMSGARCHD$spec$n.params)
             
             logLikelihoodF=sum(dnorm(residuals,mean=0,sd=sigmaLLF,log=TRUE))
             aicCrit=(2*kMS)-(2*logLikelihoodF)
           },
           "std"={
             
             if (!is.null(fittedMSGARCHD)){
               # for ML models:
               coefsMS=as.data.frame(summary(fittedMSGARCHD)$estimate)
               stableProbs=summary(fittedMSGARCHD)$stable.prob
             } else {
               # for MCMC models:
               coefsMS=as.data.frame(summary(fittedMSGARCHD)$summary )
               stableProbs=summary(fittedMSGARCHD)$post.stable.prob
             }
             
             
             nuMS=c(coefsMS[which(rownames(coefsMS)=="nu_1"),1],
                    coefsMS[which(rownames(coefsMS)=="nu_2"),2])
             nuMS=as.numeric(nuMS%*%stableProbs)
             
             kMS=sum(fittedMSGARCHD$spec$n.params)
             
             logLikelihoodF=sum(dstd((residuals/1),mean=0,sd=sigmaLLF,nu=nuMS,log=TRUE))
             aicCrit=(2*kMS)-(2*logLikelihoodF)
           },
           "ged"={
             if (!is.null(fittedMSGARCHD)){
               # for ML models:
               coefsMS=as.data.frame(summary(fittedMSGARCHD)$estimate)
               stableProbs=summary(fittedMSGARCHD)$stable.prob
             } else {
               # for MCMC models:
               coefsMS=as.data.frame(summary(fittedMSGARCHD)$summary )
               stableProbs=summary(fittedMSGARCHD)$post.stable.prob
             }
             
             nuMS=c(coefsMS[which(rownames(coefsMS)=="nu_1"),1],
                    coefsMS[which(rownames(coefsMS)=="nu_2"),2])
             nuMS=as.numeric(nuMS%*%stableProbs)
             
             kMS=sum(fittedMSGARCHD$spec$n.params)
             
             logLikelihoodF=sum(dged(residuals,mean=0,sd=sigmaLLF,nu=nuMS,log=TRUE))
             aicCrit=(2*kMS)-(2*logLikelihoodF)
           }
    )
    
    DBTable=rbind(DBTable,
                  data.frame(Date=as.character(tail(data$Date,1)),
                             Value=c(logLikelihoodF,
                                     aicCrit,
                                     kMS),
                             Ticker=c("general LLF",
                                      "general Akaike",
                                      "number of parameters"),
                             Experiment=experiment)
    )      
    
# Closes the if statement for the ML estimation in the fittedMSGARCHD object:====
  } 

# The case in which is necessary to use MCMC estimation:====
  
if (is.null(fittedMSGARCHD)){
  cat("\f")
print("E-M estimation not feasible. Estimating with MCMC method instead...")
  
  simmethod="MCMC"
  fittedMSGARCHD = tryCatch(FitMCMC(spec = MSspec, data = residuals,
                                    ctr=list(nburn=numberBurn,nmcmc=numberMCMC)) ,
                            error=function(e) NULL)
    
  meltedCoefsTable=melt(summary(fittedMSGARCHD)$summary,id=c("Mean","SE"))
  msGARCHOutPutCoefsTable=summary(fittedMSGARCHD)$summary
    # Parameters coefficients:
    DBTable=rbind(DBTable,
                  data.frame(Date=as.character(tail(data$Date,1)),
                             Value=meltedCoefsTable$value,
                             Ticker=paste0(meltedCoefsTable$Var1,"-",meltedCoefsTable$Var2),
                             Experiment=experiment
                  )
    )
    # DIC values:
    DBTable=rbind(DBTable,
                  data.frame(Date=as.character(tail(data$Date,1)),
                             Value=summary(fittedMSGARCHD)$DIC,
                             Ticker="DIC-MCMC estimation",
                             Experiment=experiment
                  )
    )
    # estimation method
    DBTable=rbind(DBTable,
                  data.frame(Date=as.character(tail(data$Date,1)),
                             Value="MCMC",
                             Ticker="MSGARCH estimation method",
                             Experiment=experiment
                  )
    )
    
    #Smoothed probabilities and forecasts:
    # Smoothed probabilites:
    stateMS=State(fittedMSGARCHD)
    
    Smooth.probs1 = stateMS$SmoothProb[,1, 1:MSspec$K, drop = TRUE]
    
    # Transition probability matrix:
    
    transprob1=summary(fittedMSGARCHD)$post.trans.mat
    
    # Forecasted probabilities:
    
    Predprobs1 = rbind(Smooth.probs1[nrow(Smooth.probs1),]%*%transprob1,
                       Smooth.probs1[nrow(Smooth.probs1),]%*%transprob1^2,
                       Smooth.probs1[nrow(Smooth.probs1),]%*%transprob1^3,
                       Smooth.probs1[nrow(Smooth.probs1),]%*%transprob1^4,
                       Smooth.probs1[nrow(Smooth.probs1),]%*%transprob1^5)
    
    # smooth probs data Base table:
    DBTable=rbind(DBTable,
                  data.frame(Date=as.character(tail(data$Date,1)),
                             Value=c(Predprobs1[,1],Predprobs1[,2]),
                             Ticker=c(paste0("Calm prob. forecast t+",seq(1,5)),
                                      paste0("Crisis prob. forecast t+",seq(1,5))),
                             Experiment=experiment
                  )
    )
    
    #  Actual return and actual smooth probability regime:
    
    if (tail(data[,depVarPos],1)>0){
      if (tail(Smooth.probs1[,2],1)<=0.5){
        msTradingScenarioActual=1
        msTradingScenarioTextActual="normal bullish"
        msVolatilityScenarioTextActual="normal"
        msVolatilityScenarioActual=1
      } else {
        msTradingScenarioActual=2
        msTradingScenarioTextActual="crisis bullish"
        msVolatilityScenarioTextActual="crisis"
        msVolatilityScenarioActual=2
      }
    } else {
      if (tail(Smooth.probs1[,2],1)<=0.5){
        msTradingScenarioActual=3
        msTradingScenarioTextActual="normal bearish"
        msVolatilityScenarioTextActual="normal"
        msVolatilityScenarioActual=1
      } else {
        msTradingScenarioActual=4
        msTradingScenarioTextActual="crisis bearish"
        msVolatilityScenarioTextActual="crisis"
        msVolatilityScenarioActual=2
      }
    }
    
    # Expected return and actual smooth probability forecasted regimes:
    
    if (expectedReturnScenario>0){
      if (tail(Smooth.probs1[,2],1)<=0.5){
        msTradingScenarioF1=1
        msTradingScenarioTextF1="normal bullish"
        msVolatilityScenarioTextF1="normal"
        msVolatilityScenarioF1=1
      } else {
        msTradingScenarioF1=2
        msTradingScenarioTextF1="crisis bullish"
        msVolatilityScenarioTextF1="crisis"
        msVolatilityScenarioF1=2
      }
    } else {
      if (tail(Smooth.probs1[,2],1)<=0.5){
        msTradingScenarioF1=3
        msTradingScenarioTextF1="normal bearish"
        msVolatilityScenarioTextF1="normal"
        msVolatilityScenarioF1=1
      } else {
        msTradingScenarioF1=4
        msTradingScenarioTextF1="crisis bearish"
        msVolatilityScenarioTextF1="crisis"
        msVolatilityScenarioF1=2
      }
    }
    
    # Expected return  and t+1 forecasted smooth probability forecasted regimes:
    
    if (expectedReturnScenario>0){
      if (Predprobs1[1,2]<=0.5){
        msTradingScenarioF2=1
        msTradingScenarioTextF2="normal bullish"
        msVolatilityScenarioTextF2="normal"
        msVolatilityScenarioF2=1
      } else {
        msTradingScenarioF2=2
        msTradingScenarioTextF2="crisis bullish"
        msVolatilityScenarioTextF2="crisis"
        msVolatilityScenarioF2=2
      }
    } else {
      if (Predprobs1[1,2]<=0.5){
        msTradingScenarioF2=3
        msTradingScenarioTextF2="normal bearish"
        msVolatilityScenarioTextF2="normal"
        msVolatilityScenarioF2=1
      } else {
        msTradingScenarioF2=4
        msTradingScenarioTextF2="crisis bearish"
        msVolatilityScenarioTextF2="crisis"
        msVolatilityScenarioF2=2
      }
    }
    
    DBTable=rbind(DBTable,
                  data.frame(
                    Date=as.character(tail(data$Date,1)),
                    Value=c(msTradingScenarioActual,
                            msTradingScenarioTextActual,
                            msVolatilityScenarioActual,
                            msVolatilityScenarioTextActual),
                    Ticker=c("MS trading scenario Actual",
                             "MS trading scenario text Actual",
                             "MS volatility scenario Actual",
                             "MS volatility scenario text Actual"),
                    Experiment=experiment
                  ),
                  
                  data.frame(
                    Date=as.character(tail(data$Date,1)),
                    Value=c(msTradingScenarioF1,
                            msTradingScenarioTextF1,
                            msVolatilityScenarioF1,
                            msVolatilityScenarioTextF1),
                    Ticker=c("MS trading scenario Actual vol forecast",
                             "MS trading scenario text Actual vol forecast",
                             "MS volatility scenario Actual vol forecast",
                             "MS volatility scenario text Actual vol forecast"),
                    Experiment=experiment
                  ),
                  
                  data.frame(
                    Date=as.character(tail(data$Date,1)),
                    Value=c(msTradingScenarioF2,
                            msTradingScenarioTextF2,
                            msVolatilityScenarioF2,
                            msVolatilityScenarioTextF2),
                    Ticker=c("MS trading scenario all Forecast",
                             "MS trading scenario text all Forecast",
                             "MS volatility scenario all Forecast",
                             "MS volatility scenario text all Forecast"),
                    Experiment=experiment
                  )
    )
    
# Estimates forecasted volatility and VaR at t:
    cat("\f")
print(paste0("Estimating expected volatility at t. (",experiment,"-",GARCHmodels,")"))
    
    pred1 <- predict(fittedMSGARCHD, nahead = 5, do.return.draws = FALSE)
    checkNARisk=which(is.na(pred1$vol))
    
    
    if (length(checkNARisk)>0){
      pred1$vol=na.locf(pred1$vol)
    }
    
    DBTable=rbind(DBTable,
                  data.frame(
                    Date=as.character(tail(data$Date,1)),
                    Value=pred1$vol,
                    Ticker=paste0(paste0("Volatility at t+",seq(1:5))," model 1"),
                    Experiment=experiment
                  )
    )
    
    # general model LLF and AIC estimation from data:
    
    sigmaLLF=pred1$vol[1]
    muLLF=meanForecast
    
    llfFunct=pdfFunct[1]
    
    switch(llfFunct,
           "norm"={
             kMS=sum(fittedMSGARCHD$spec$n.params)
             
             logLikelihoodF=sum(dnorm(residuals,mean=0,sd=sigmaLLF,log=TRUE))
             aicCrit=(2*kMS)-(2*logLikelihoodF)
           },
           "std"={
             
             # for MCMC models:
             coefsMS=as.data.frame(summary(fittedMSGARCHD)$summary )
             stableProbs=summary(fittedMSGARCHD)$post.stable.prob
             nuMS=c(coefsMS[which(rownames(coefsMS)=="nu_1"),1],
                    coefsMS[which(rownames(coefsMS)=="nu_2"),2])
             nuMS=as.numeric(nuMS%*%stableProbs)
             
             kMS=sum(fittedMSGARCHD$spec$n.params)
             
             logLikelihoodF=sum(dstd((residuals/1),mean=0,sd=sigmaLLF,nu=nuMS,log=TRUE))
             aicCrit=(2*kMS)-(2*logLikelihoodF)
           },
           "ged"={
             # for MCMC models:
             coefsMS=as.data.frame(summary(fittedMSGARCHD)$summary )
             stableProbs=summary(fittedMSGARCHD)$post.stable.prob
             
             nuMS=c(coefsMS[which(rownames(coefsMS)=="nu_1"),1],
                    coefsMS[which(rownames(coefsMS)=="nu_2"),2])
             nuMS=as.numeric(nuMS%*%stableProbs)
             
             kMS=sum(fittedMSGARCHD$spec$n.params)
             
             logLikelihoodF=sum(dged(residuals,mean=0,sd=sigmaLLF,nu=nuMS,log=TRUE))
             aicCrit=(2*kMS)-(2*logLikelihoodF)
           }
    )
    
    DBTable=rbind(DBTable,
                  data.frame(Date=as.character(tail(data$Date,1)),
                             Value=c(logLikelihoodF,
                                     aicCrit,
                                     kMS),
                             Ticker=c("general LLF",
                                      "general Akaike",
                                      "number of parameters"),
                             Experiment=experiment)
    )      
    
  }  
  

# Closes the simulation records:
tiempoPasado=as.numeric(Sys.time()-startCalculation)

elapsTimeMsg=paste0("Estimation elapsed time: ",
                    round(tiempoPasado/360,0),
                    ":",
                    round(tiempoPasado/60,0),":",
                    round(tiempoPasado,3))
DBTable=rbind(DBTable,
              data.frame(Date=as.character(tail(data$Date,1)),
                         Value=c(tiempoPasado,
                                 elapsTimeMsg),
                         Ticker=c("elapsed time",
                                  "elapsed time message"),
                         Experiment=experiment))

# Output list object:
output=list(msGarchFittedModel=fittedMSGARCHD,
            dbInputTable=DBTable,
            smoothProbs=Smooth.probs1,
            predProbs=Predprobs1,
            estimationTime=elapsTimeMsg,
            msGARCHOutPutCoefsTable=msGARCHOutPutCoefsTable,
            meanEqCoefsTable=meanEqCoefsTable,
            meanEqR2=meanEqR2,
            maenEqSigma=maenEqSigma
            )
cat("\f")
print(paste0("MS-GARCH estimation concluded with ",simmethod," method. ",elapsTimeMsg))
# Returns the output object:
return(output)
# Function ends here:   
}