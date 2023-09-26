
msMCMCGARCHFit=function(eq=EQ,data=Data,numberMCMC=10000,numberBurn=500,GARCHmodels=c("sGARCH","sGARCH"),pdfFunct=c("norm","norm"),experiment="",timeFixed=FALSE){

  #outputData=Data

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

# Coefficients tbale in DBTable
  DBTable=data.frame(
    Date=as.character(tail(Data$Date,1)),
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

  if (tail(datos$Return,1)>0){
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
                  Date=as.character(tail(Data$Date,1)),
                  Value=tail(datos$Settle,1),
                  Ticker="Price at t",
                  Experiment=experiment
                ),
                data.frame(
                  Date=as.character(tail(Data$Date,1)),
                  Value=tail(datos$Return,1),
                  Ticker="Return at t",
                  Experiment=experiment
                ),
                data.frame(
                  Date=as.character(tail(Data$Date,1)),
                  Value=meanForecast,
                  Ticker="Return Forecast at t",
                  Experiment=experiment
                ),
                data.frame(
                  Date=as.character(tail(Data$Date,1)),
                  Value=returnScenario,
                  Ticker="Return scenario",
                  Experiment=experiment
                ),
                data.frame(
                  Date=as.character(tail(Data$Date,1)),
                  Value=returnScenarioText,
                  Ticker="Return scenario text",
                  Experiment=experiment
                ),
                data.frame(
                  Date=as.character(tail(Data$Date,1)),
                  Value=expectedReturnScenario,
                  Ticker="Expected return scenario",
                  Experiment=experiment
                ),
                data.frame(
                  Date=as.character(tail(Data$Date,1)),
                  Value=expectedReturnScenarioText,
                  Ticker="Expected return scenario text",
                  Experiment=experiment
                )
  )
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

# Estimates the MSGARH model:

print(paste0("Estimating MS-GARCH (EM method) model for experiment ", experiment,", date: ",tail(Data$Date,1)))

# ML estimation:
simmethod="E-M"
fittedMSGARCHD = tryCatch(FitML(spec = MSspec, data = residuals) ,
                          error=function(e) NULL)

if (!is.null(fittedMSGARCHD)){

  meltedCoefsTable=melt(fittedMSGARCHD$Inference$MatCoef)
  # Parameters coefficients:
  DBTable=rbind(DBTable,
                data.frame(Date=as.character(tail(Data$Date,1)),
                           Value=meltedCoefsTable$value,
                           Ticker=paste0(meltedCoefsTable$Var1,"-",meltedCoefsTable$Var2),
                           Experiment=experiment
                )
  )
  # AIC values:
  DBTable=rbind(DBTable,
                data.frame(Date=as.character(tail(Data$Date,1)),
                           Value=summary(fittedMSGARCHD)$AIC,
                           Ticker="AIC-ML estimation",
                           Experiment=experiment
                            )
                )
  # estimation method
  DBTable=rbind(DBTable,
                data.frame(Date=as.character(tail(Data$Date,1)),
                           Value="EM",
                           Ticker="MSGARCH estimation method",
                           Experiment=experiment
                            )
                )

MLestimation=TRUE
} else {

MLestimation=FALSE
}


# MCMC estimation:

if (!isTRUE(MLestimation)){

print(paste0("Estimating MS-GARCH (MCMC method) model for", experiment,", date: ",tail(Data$Date,1)))
  simmethod="MCMC"
  fittedMSGARCHD = tryCatch(FitMCMC(spec = MSspec, data = residuals,
                                    ctr=list(nburn=numberBurn,nmcmc=numberMCMC)) ,
                            error=function(e) NULL)

  if (!is.null(fittedMSGARCHD)){

    meltedCoefsTable=melt(summary(fittedMSGARCHD)$summary,id=c("Mean","SE"))
    # Parameters coefficients:
    DBTable=rbind(DBTable,
                  data.frame(Date=as.character(tail(Data$Date,1)),
                             Value=meltedCoefsTable$value,
                             Ticker=paste0(meltedCoefsTable$Var1,"-",meltedCoefsTable$Var2),
                             Experiment=experiment
                  )
    )
    # DIC values:
    DBTable=rbind(DBTable,
                  data.frame(Date=as.character(tail(Data$Date,1)),
                             Value=summary(fittedMSGARCHD)$DIC,
                             Ticker="DIC-MCMC estimation",
                             Experiment=experiment
                  )
    )

    # estimation method
    DBTable=rbind(DBTable,
                  data.frame(Date=as.character(tail(Data$Date,1)),
                             Value="MCMC",
                             Ticker="MSGARCH estimation method",
                             Experiment=experiment
                  )
    )

  }

}


# Generating the Smoothed, transtition and forecasted probabilities for model 1:

print(paste0("Estimating regime-specific smoothed probs. ","norm","-","norm"," (model ",b,"). (",experiment,"-",GARCHmodels,")"))

# Smoothed probabilites:
stateMS=State(fittedMSGARCHD)

Smooth.probs1 = stateMS$SmoothProb[,1, 1:MSspec$K, drop = TRUE]

# Transition probability matrix:

if (isTRUE(MLestimation)){
  transprob1=summary(fittedMSGARCHD)$trans.mat
} else {
  transprob1=summary(fittedMSGARCHD)$post.trans.mat
}


# Forecasted prpbabilities:

Predprobs1 = rbind(Smooth.probs1[nrow(Smooth.probs1),]%*%transprob1,
                  Smooth.probs1[nrow(Smooth.probs1),]%*%transprob1^2,
                  Smooth.probs1[nrow(Smooth.probs1),]%*%transprob1^3,
                  Smooth.probs1[nrow(Smooth.probs1),]%*%transprob1^4,
                  Smooth.probs1[nrow(Smooth.probs1),]%*%transprob1^5)

# smooth probs Data Base table:
DBTable=rbind(DBTable,
            data.frame(Date=as.character(tail(Data$Date,1)),
           Value=c(Predprobs1[,1],Predprobs1[,2]),
           Ticker=c(paste0("Calm prob. forecast t+",seq(1,5)),
                    paste0("Crisis prob. forecast t+",seq(1,5))),
           Experiment=experiment
           )
)

#  Actual return and actual smooth probability regime:

if (tail(Data$Return,1)>0){
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
                Date=as.character(tail(Data$Date,1)),
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
                Date=as.character(tail(Data$Date,1)),
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
                Date=as.character(tail(Data$Date,1)),
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

print(paste0("Estimating expected volatility at t. (",experiment,"-",GARCHmodels,")"))

pred1 <- predict(fittedMSGARCHD, nahead = 5, do.return.draws = FALSE)


DBTable=rbind(DBTable,
              data.frame(
                Date=as.character(tail(Data$Date,1)),
                Value=pred1$vol,
                Ticker=paste0(paste0("Volatility at t+",seq(1:5))," model 1"),
                Experiment=experiment
              )
)

# Estimates CVaR

RiskPred=Risk(fittedMSGARCHD,alpha=c(0.02,0.05),nahead=5)

# Risk forecast table:

DBTable=rbind(DBTable,
              data.frame(Date=as.character(tail(Data$Date,1)),
                         Value=c(RiskPred$VaR[,1],
                                 RiskPred$VaR[,2]),
                         Ticker=c(paste0("VaR 98% t+",seq(1,5)),
                                  paste0("VaR 95% t+",seq(1,5))),
                         Experiment=experiment
                          ),
              data.frame(Date=as.character(tail(Data$Date,1)),
                         Value=c(RiskPred$ES[,1],
                                 RiskPred$ES[,2]),
                         Ticker=c(paste0("CVaR 98% t+",seq(1,5)),
                                  paste0("CVaR 95% t+",seq(1,5))),
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

 if (isTRUE(MLestimation)){
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
           if (isTRUE(MLestimation)){
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
                data.frame(Date=as.character(tail(Data$Date,1)),
                           Value=c(logLikelihoodF,
                                   aicCrit,
                                   kMS),
                           Ticker=c("general LLF",
                                    "general Akaike",
                                    "number of parameters"),
                           Experiment=experiment)
  )

  # Closes the code:
  tiempoPasado=as.numeric(Sys.time()-startCalculation)

  elapsTimeMsg=paste0("Estimation elapsed time: ",
                      round(tiempoPasado/360,0),
                      ":",
                      round(tiempoPasado/60,0),":",
                      round(tiempoPasado,3))
  DBTable=rbind(DBTable,
                data.frame(Date=as.character(tail(Data$Date,1)),
                           Value=c(tiempoPasado,
                                   elapsTimeMsg),
                           Ticker=c("elapsed time",
                                    "elapsed time message"),
                           Experiment=experiment)
  )

# Creating the output list object:
  MSGARCHResults=list(
    outputData=Data,
    outPutDBtable=DBTable,
    elapsTimeMsg=elapsTimeMsg,
    elapsTime=tiempoPasado,
    simMethod=simmethod
  )

  cat("\f")
  print(elapsTimeMsg)

# Return object:
  return(MSGARCHResults)

}
