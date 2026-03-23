#####################subgroup analysis
rm(list = ls())
load("primary analysis data (income level).RDATA")
###########
analysis<-function(data,data_name){
  data$age_c2<-ifelse(data$age_wave1<70,0,ifelse(
    data$age_wave1>=70,1,NA
  ))
  data$age_c2<-as.factor(data$age_c2)
  data$education_c2<-ifelse(data$education_wave1==0,0,ifelse(
    data$education_wave1 %in% c(1,2),1,NA
  ))
  data$education_c2<-as.factor(data$education_c2)
  var<-c("age_c2","sex_wave1","education_c2","marital_status_wave1","smoke_wave1",
         "drink_wave1","pa_wave1","diabetes_wave1","hypertension_wave1","cvd_wave1")
  result<-NULL
  for (i in var) {
    a<-unique(data[[i]])
    for (j in a) {
      data1<-subset(data,data[[i]]==j)
      fit<-coxme(Surv(follow_year,dementia_final)~hearing_aids_c2+(1|cohort),data=data1,weights = w)
      a<-data.frame(hr=exp(coef(fit)),
                    cilower=exp(confint(fit))[1],
                    ciupper=exp(confint(fit))[2],
                    p=summary(fit)$coefficients[5],
                    variable="hearing aids")
      a$hrci<-sprintf("%0.2f (%0.2f, %0.2f)",
                      a$hr,a$cilower,a$ciupper)
      
      ####hearing aids wearing effect
      fit<-coxme(Surv(follow_year,dementia_final)~hearing_loss_aids_c5+(1|cohort),data=data1,weights = w1)
      b<-data.frame(hr=exp(coef(fit)),
                    cilower=exp(confint(fit))[,1],
                    ciupper=exp(confint(fit))[,2],
                    p=summary(fit)$coefficients[,5],
                    variable="hearing aids effect")
      b$hrci<-sprintf("%0.2f (%0.2f, %0.2f)",
                      b$hr,b$cilower,b$ciupper)
      ab<-rbind(a,b)
      ab$cohort<-data_name
      ab$variable<-i
      ab$subgroup<-j
      result<-rbind(result,ab)
    }
  }
  return(result)
}
result_list<-list()
for (i in seq_along(data_list1)) {
  data_name<-names(data_list1)[[i]]
  result_list[[data_name]]<-analysis(data_list1[[i]],data_name)
}
result<-do.call(bind_rows,result_list)


####################p for interaction
analysis<-function(data,data_name){
  data$age_c2<-ifelse(data$age_wave1<70,0,ifelse(
    data$age_wave1>=70,1,NA
  ))
  data$age_c2<-as.factor(data$age_c2)
  data$education_c2<-ifelse(data$education_wave1==0,0,ifelse(
    data$education_wave1 %in% c(1,2),1,NA
  ))
  data$education_c2<-as.factor(data$education_c2)
  var<-c("age_c2","sex_wave1","education_c2","marital_status_wave1","smoke_wave1",
         "drink_wave1","pa_wave1","diabetes_wave1","hypertension_wave1","cvd_wave1")
  result<-NULL
  for (i in var) {
    fit<-coxme(Surv(follow_year,dementia_final)~hearing_aids_c2+data[[i]]:hearing_aids_c2+(1|cohort),data=data,weights = w)
    a<-summary(fit)$coefficients
    a<-as.data.frame(a)
    a<-a[-1,]
    a$variable<-i
    a$exp<-"hearing aids"
    a<-a[,c(5,6,7)]
    ####hearing aids wearing effect
    fit<-coxme(Surv(follow_year,dementia_final)~hearing_loss_aids_c5+hearing_loss_aids_c5:data[[i]]+(1|cohort),data=data,weights = w1)
    b<-summary(fit)$coefficients
    b<-as.data.frame(b)
    b<-b[-c(1,2),]
    b$variable<-i
    b$exp<-"hearing aids effect"
    b<-b[,c(5,6,7)]
    ab<-rbind(a,b)
    ab$cohort<-data_name
    result<-rbind(result,ab)
  }
  return(result)
}
result_list<-list()
for (i in seq_along(data_list1)) {
  data_name<-names(data_list1)[[i]]
  result_list[[data_name]]<-analysis(data_list1[[i]],data_name)
}
result<-do.call(bind_rows,result_list)
