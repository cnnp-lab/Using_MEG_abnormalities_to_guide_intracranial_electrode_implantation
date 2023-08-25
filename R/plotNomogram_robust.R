#Section Start: Load in the packages ----

if(!('pacman' %in% rownames(installed.packages()))){
  install.packages('pacman')
}

pacman::p_load('rms',
               'dplyr',
               'pROC')


#Section Start: Load in the data ----
patient_metrics=read.csv('../data/nomogram_data.csv',stringsAsFactors = FALSE)

#convert binary outcome to ILAE scores
patient_metrics=patient_metrics%>%
  mutate(outcome=ifelse(outcome==0,'ILAE2+','ILAE1'))

#convert to a factor
patient_metrics$outcome=factor(patient_metrics$outcome,levels=c("ILAE2+","ILAE1"))


patient_metrics_subset=patient_metrics[,c('meg_drs_electrode_coverage','electrode_concordance','eeg_drs','outcome')]

colnames(patient_metrics_subset)=c('MEG_DRS','EPC','EEG_DRS','outcome')



t.data=datadist(patient_metrics_subset)

options(datadist='t.data')


var.labels <- c(MEG_DRS = "MEG DRS",
                EPC = "Electrode Concordance",
                EEG_DRS="iEEG DRS",
                outcome = "Surgical Outcome")




label(patient_metrics_subset) = lapply(names(var.labels),
                   function(x) label(patient_metrics_subset[,x]) <- var.labels[x])


#################################################
#Functions
#################################################

#####################################################################

formula_lp <- function(nomogram,power,digits=6){
  #total.points always appears in names of nomogram
  #total.points can only be changed in plot using points.label
  variable_part=nomogram["lp" == names(nomogram)]
  if (missing(power)){
    #missing power : choose power automatically
    power = 0
    test=data.frame(R2=0.5)
    while (test$R2<1) {
      power=power+1
      ######get 2 variables
      #1
      points=as.numeric(unlist(variable_part$lp[1]))
      #2
      lp=as.numeric(unlist(variable_part$lp[2]))
      ######calculate
      formu=as.formula(paste0('points~',
                              inner_Add_Symbol(paste0("I(lp^",1:power,")"))))
      #regressiong
      reg=lm(formu)
      #formula
      lm.result=data.frame(t(reg$coefficients))
      rownames(lm.result)="linear predictor"
      colnames(lm.result)=c("b0",paste0("x^",1:power))
      #real,fit,diff
      fit=predict(reg)
      diff=points-fit
      real_fit=t(data.frame(nomogram=points,fit,diff))
      colnames(real_fit)=lp
      # test
      R2=suppressWarnings(summary(reg)$r.squared)
      RMSE=(mean(predict(reg)-points)^2)^(1/2)
      test=data.frame(R2,RMSE)
    }
  }else{
    #exist power
    if (power<1) stop("power must not be less 1")
    ######get 2 variables
    #1
    points=as.numeric(unlist(variable_part$lp[1]))
    #2
    lp=as.numeric(unlist(variable_part$lp[2]))
    ######calculate
    formu=as.formula(paste0('points~',
                            inner_Add_Symbol(paste0("I(lp^",1:power,")"))))
    #regressiong
    reg=lm(formu)
    #formula
    lm.result=data.frame(t(reg$coefficients))
    rownames(lm.result)="linear predictor"
    colnames(lm.result)=c("b0",paste0("x^",1:power))
    #real,fit,diff
    fit=predict(reg)
    diff=points-fit
    real_fit=t(data.frame(nomogram=points,fit,diff))
    colnames(real_fit)=lp
    # test
    R2=suppressWarnings(summary(reg)$r.squared)
    RMSE=(mean(predict(reg)-points)^2)^(1/2)
    test=data.frame(R2,RMSE)
  }
  rownames(test)="linear predictor"
  result=list(formula=round(lm.result,digits),
              test=round(test,digits),
              diff=round(real_fit,digits))
  return(result)
}


points_cal <- function(formula,rd,lp,digits=6){
  #data
  if (!missing(rd)){
    if (any(is.na(rd))) stop("data must not have NAs.")
    nomoF.matrix=as.matrix(formula)
    data1=rd[,rownames(nomoF.matrix)]
    for (i in 1:ncol(data1)) {
      if (i==1) score=rep(0,nrow(data1))
      beta_for_x=nomoF.matrix[i,]
      x_for_beta=data1[,i]
      if (is.factor(x_for_beta)){
        x_for_beta=as.numeric(as.character(x_for_beta))
      }
      for (j in 1:length(nomoF.matrix[i,])) {
        if (is.na(beta_for_x[j])) next(j)
        if (j==1){
          each.var.score=(beta_for_x[j])*(x_for_beta^(j-1))
        }else{
          each.var.score=each.var.score+(beta_for_x[j])*(x_for_beta^(j-1))
        }
      }
      score=score+each.var.score
    }
    
    return(score)
    #lp
  }else if (!missing(lp)){
    for (j in 1:ncol(formula)) {
      if (is.na(formula[1,j])) next(j)
      if (j==1){
        score=(formula[1,j])*(lp^(j-1))
      }else{
        score=score+(formula[1,j])*(lp^(j-1))
      }
    }
    score=round(score,digits)
    return(score)
  }
}


inner_Add_Symbol <- function(character,symbol="+"){
  if (length(character)>=2){
    for (character.i in 1:length(character)) {
      if (character.i==1){
        adj=character[1]
      }else{
        adj=paste0(adj,symbol,character[character.i])
      }
    }
  }else{
    adj=character
  }
  adj
}


leave_one_out=function(training, test){
  
  n=dim(training)[1]
  
  weights=case_when(training$outcome=='ILAE2+'~freq['ILAE1']/n,TRUE~freq['ILAE2+']/n)
  
  fit=lrm(outcome~MEG_DRS+EPC+EEG_DRS,data=training,scale=FALSE,weights=weights)
  
  nomo=nomogram(fit)
  
  results <- formula_lp(nomo)
  
  
  
  lin_pred=fit$coefficients[1]+
    sum(fit$coefficients[-1]*test[1,c(1,2,3)])
  
  
  pts=points_cal(formula = results$formula,lp=lin_pred)
  
  pts_training=points_cal(formula = results$formula,lp=fit$linear.predictors)
  
  
  auc=roc(training$outcome,pts_training)
  
  threshold=points_cal(formula = results$formula,lp=0)
  
  #return(threshold)
  return(auc$auc)
  #return(pts)
}


#Section Start: Create a nomogram ----

# Run the nomogram normally
freq=table(patient_metrics_subset$outcome)

weights=case_when(patient_metrics_subset$outcome=='ILAE2+'~freq['ILAE1']/nrow(patient_metrics),TRUE~freq['ILAE2+']/nrow(patient_metrics))


#Fit the model
fit=lrm(outcome~EPC+MEG_DRS+EEG_DRS+0,data=patient_metrics_subset,scale=FALSE,weights = weights)

#compute the nomogram
nomo=nomogram(fit)

#create dir if doesn't exist
if (!dir.exists("../figures/figure4")){
  dir.create("../figures/figure4")
}

#uncomment these lines if you want to plot and save the nomogram
pdf(file = "../figures//figure4/drs_nomogram.pdf")
plot(nomo,cex.axis = 2,cex.var=2,tcl=-0.5)
dev.off()


results <- formula_lp(nomo)

pts=points_cal(formula = results$formula,lp=fit$linear.predictors)


patient_metrics['pts']=pts

#now run it using leave one out
patient_metrics['auc_robust']=sapply(1:32,function(i) leave_one_out(patient_metrics_subset[-i,],
                                                                    patient_metrics_subset[i,]))



#save the results to a new csv file
write.csv(patient_metrics,file='../data/nomogram_predict_full.csv',row.names =FALSE)

