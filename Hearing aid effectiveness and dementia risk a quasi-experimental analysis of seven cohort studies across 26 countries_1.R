rm(list = ls())
load("hrs_final_data.RDATA")
load("MHAS_final_data.RDATA")
load("charls_final_data.RDATA")
load("klosa_final_data.RDATA")
load("share_final_data.RDATA")
load("elsa_final_data.RDATA")
load("TILDA_final_data.RDATA")

########weighted data using IPTW
library(ipw)
data_list<-list(elsa,share,TILDA,charls,klosa,hrs,MHAS)
all<-do.call(bind_rows,data_list)

weight<-function(data){
  data<-subset(data,follow_year>0)
  data<-subset(data,age_wave1>=55)
  data<-subset(data,!is.na(hearing_loss_aids_c5))
  data$hearing_loss_aids_c5<-ifelse(data$hearing_loss_aids_c5==3,0,ifelse(
    data$hearing_loss_aids_c5==1,2,ifelse(
      data$hearing_loss_aids_c5==2,1,NA
    )
  ))
  data$hearing_loss_aids_c5<-as.factor(data$hearing_loss_aids_c5)
  data$hearing_aids_c2<-as.factor(data$hearing_aids_c2)
  fit<-glm(hearing_aids_c2~age_wave1+sex_wave1+education_wave1+marital_status_wave1+smoke_wave1+
             drink_wave1+pa_wave1+diabetes_wave1+hypertension_wave1+cvd_wave1,data=data, 
           family = binomial(link = "logit"))
  
  data$ps<-predict(fit,type = "response")
  data$w <- ifelse(data$hearing_aids_c2==1,1/data$ps,1/(1-data$ps))
  
  w1<-ipwpoint(exposure = hearing_loss_aids_c5,
               family = "multinomial",
               numerator = ~1,
               denominator = ~age_wave1+sex_wave1+education_wave1+marital_status_wave1+smoke_wave1+
                 drink_wave1+pa_wave1+diabetes_wave1+hypertension_wave1+cvd_wave1,
               data = data)
  data$w1<-w1$ipw.weights
  return(data)
}
all<-lapply(all, weight)


#########eastern asia
asia<-subset(all, cohort %in% c(4,5))
#########europe
europe<-subset(all, cohort %in% c(1,2,3))
#########north america
america<-subset(all, cohort %in% c(6,7))
########high-income countries
high_income<-subset(all, cohort %in% c(1,2,3,5,6))
########middle-income countries
middle_income<-subset(all, cohort %in% c(4,7))

data_list<-setNames(list(all,asia,europe,america,high_income,middle_income),
                    c("all","asia","europe","america","high_income","middle_income"))

###########primary analysis
library(coxme)
analysis1<-function(data,data_name){
  fit<-coxme(Surv(follow_year,dementia_final)~hearing_aids_c2+(1|cohort),data=data,weights = w)
  a<-data.frame(hr=exp(coef(fit)),
                cilower=exp(confint(fit))[1],
                ciupper=exp(confint(fit))[2],
                p=summary(fit)$coefficients[5],
                variable="hearing aids")
  a$hrci<-sprintf("%0.2f (%0.2f, %0.2f)",
                  a$hr,a$cilower,a$ciupper)
  
  ####hearing aids wearing effect
  fit<-coxme(Surv(follow_year,dementia_final)~hearing_loss_aids_c5+(1|cohort),data=data,weights = w1)
  b<-data.frame(hr=exp(coef(fit)),
                cilower=exp(confint(fit))[,1],
                ciupper=exp(confint(fit))[,2],
                p=summary(fit)$coefficients[,5],
                variable="hearing aids effect")
  b$hrci<-sprintf("%0.2f (%0.2f, %0.2f)",
                  b$hr,b$cilower,b$ciupper)
  ab<-rbind(a,b)
  ab$cohort<-data_name
  return(ab)
}
result_list<-list()
for (i in seq_along(data_list)) {
  data_name<-names(data_list)[[i]]
  result_list[[data_name]]<-analysis1(data_list[[i]],data_name)
}
result<-do.call(bind_rows,result_list)


###########sensitivity analysis (aged over 60 years old)
analysis<-function(data,data_name){
  data<-subset(data,age_wave1>=60)
  fit<-coxme(Surv(follow_year,dementia_final)~hearing_aids_c2+(1|cohort),data=data,weights = w)
  a<-data.frame(hr=exp(coef(fit)),
                cilower=exp(confint(fit))[1],
                ciupper=exp(confint(fit))[2],
                p=summary(fit)$coefficients[5],
                variable="hearing aids")
  a$hrci<-sprintf("%0.2f (%0.2f, %0.2f)",
                  a$hr,a$cilower,a$ciupper)
  
  ####hearing aids wearing effect
  fit<-coxme(Surv(follow_year,dementia_final)~hearing_loss_aids_c5+(1|cohort),data=data,weights = w1)
  b<-data.frame(hr=exp(coef(fit)),
                cilower=exp(confint(fit))[,1],
                ciupper=exp(confint(fit))[,2],
                p=summary(fit)$coefficients[,5],
                variable="hearing aids effect")
  b$hrci<-sprintf("%0.2f (%0.2f, %0.2f)",
                  b$hr,b$cilower,b$ciupper)
  ab<-rbind(a,b)
  ab$cohort<-data_name
  return(ab)
}
result_list<-list()
for (i in seq_along(data_list)) {
  data_name<-names(data_list)[[i]]
  result_list[[data_name]]<-analysis(data_list[[i]],data_name)
}
result<-do.call(bind_rows,result_list)

###########sensitivity analysis (five years after baseline)
analysis<-function(data,data_name){
  data<-subset(data,follow_year>=5)
  fit<-coxme(Surv(follow_year,dementia_final)~hearing_aids_c2+(1|cohort),data=data,weights = w)
  a<-data.frame(hr=exp(coef(fit)),
                cilower=exp(confint(fit))[1],
                ciupper=exp(confint(fit))[2],
                p=summary(fit)$coefficients[5],
                variable="hearing aids")
  a$hrci<-sprintf("%0.2f (%0.2f, %0.2f)",
                  a$hr,a$cilower,a$ciupper)
  
  ####hearing aids wearing effect
  fit<-coxme(Surv(follow_year,dementia_final)~hearing_loss_aids_c5+(1|cohort),data=data,weights = w1)
  b<-data.frame(hr=exp(coef(fit)),
                cilower=exp(confint(fit))[,1],
                ciupper=exp(confint(fit))[,2],
                p=summary(fit)$coefficients[,5],
                variable="hearing aids effect")
  b$hrci<-sprintf("%0.2f (%0.2f, %0.2f)",
                  b$hr,b$cilower,b$ciupper)
  ab<-rbind(a,b)
  ab$cohort<-data_name
  return(ab)
}
result_list<-list()
for (i in seq_along(data_list)) {
  data_name<-names(data_list)[[i]]
  result_list[[data_name]]<-analysis(data_list[[i]],data_name)
}
result<-do.call(bind_rows,result_list)

#########sensitivity analysis4 (dementia defined in a more rigorous way)
rm(list = ls())
load("hrs_final_data(transient dementia).RDATA")
load("MHAS_final_data(transient dementia).RDATA")
load("charls_final_data(transient dementia).RDATA")
load("klosa_final_data(transient dementia).RDATA")
load("share_final_data(transient dementia).RDATA")
load("elsa_final_data(transient dementia).RDATA")
load("TILDA_final_data(transient dementia).RDATA")


######IPTW
data_list<-list(elsa,share,TILDA,charls,klosa,hrs,MHAS)
all<-do.call(bind_rows,data_list)
all<-lapply(all, weight)
#########eastern asia
asia<-subset(all, cohort %in% c(4,5))
#########europe
europe<-subset(all, cohort %in% c(1,2,3))
#########north america
america<-subset(all, cohort %in% c(6,7))
########high-income countries
high_income<-subset(all, cohort %in% c(1,2,3,5,6))
########middle-income countries
middle_income<-subset(all, cohort %in% c(4,7))

data_list<-setNames(list(all,asia,europe,america,high_income,middle_income),
                    c("all","asia","europe","america","high_income","middle_income"))

result_list<-list()
for (i in seq_along(data_list)) {
  data_name<-names(data_list)[[i]]
  result_list[[data_name]]<-analysis1(data_list[[i]],data_name)
}
result<-do.call(bind_rows,result_list)

#######################sensitivity analysis 3 (at least 3 waves)
rm(list = ls())
load("hrs_final_data(3 waves).RDATA")
load("MHAS_final_data(3 waves).RDATA")
load("charls_final_data(3 waves).RDATA")
load("klosa_final_data(3 waves).RDATA")
load("share_final_data(3 waves).RDATA")
load("elsa_final_data(3 waves).RDATA")
load("TILDA_final_data(3 waves).RDATA")

######IPTW
data_list<-list(elsa,share,TILDA,charls,klosa,hrs,MHAS)
all<-do.call(bind_rows,data_list)
all<-lapply(all, weight)
#########eastern asia
asia<-subset(all, cohort %in% c(4,5))
#########europe
europe<-subset(all, cohort %in% c(1,2,3))
#########north america
america<-subset(all, cohort %in% c(6,7))
########high-income countries
high_income<-subset(all, cohort %in% c(1,2,3,5,6))
########middle-income countries
middle_income<-subset(all, cohort %in% c(4,7))

data_list<-setNames(list(all,asia,europe,america,high_income,middle_income),
                    c("all","asia","europe","america","high_income","middle_income"))

###########
result_list<-list()
for (i in seq_along(data_list)) {
  data_name<-names(data_list)[[i]]
  result_list[[data_name]]<-analysis1(data_list[[i]],data_name)
}
result<-do.call(bind_rows,result_list)



#######################sensitivity analysis 5 (excluding hearing aid in follow up)
rm(list = ls())
load("hrs_final_data (hearing aid sustain).RDATA")
load("MHAS_final_data (hearing aid sustain).RDATA")
load("charls_final_data (hearing aid sustain).RDATA")
load("klosa_final_data (hearing aid sustain).RDATA")
load("share_final_data (hearing aid sustain).RDATA")


######IPTW
data_list<-list(share,charls,klosa,hrs,MHAS)
all<-do.call(bind_rows,data_list)
all<-lapply(all, weight)
#########eastern asia
asia<-subset(all, cohort %in% c(4,5))
#########europe
europe<-subset(all, cohort %in% c(2))
#########north america
america<-subset(all, cohort %in% c(6,7))
########high-income countries
high_income<-subset(all, cohort %in% c(2,5,6))
########middle-income countries
middle_income<-subset(all, cohort %in% c(4,7))

data_list<-setNames(list(all,asia,europe,america,high_income,middle_income),
                    c("all","asia","europe","america","high_income","middle_income"))


result_list<-list()
for (i in seq_along(data_list)) {
  data_name<-names(data_list)[[i]]
  result_list[[data_name]]<-analysis1(data_list[[i]],data_name)
}
result<-do.call(bind_rows,result_list)

#########sensitivity analysis6 (dementia defined without iqcode)
rm(list = ls())
load("hrs_final_data_noiqcode.RDATA")
load("MHAS_final_data_noiqcode.RDATA")
load("charls_final_data.RDATA")
load("klosa_final_data.RDATA")
load("share_final_data.RDATA")
load("elsa_final_data_noiqcode.RDATA")
load("TILDA_final_data.RDATA")

######IPTW
data_list<-list(elsa,share,TILDA,charls,klosa,hrs,MHAS)
all<-do.call(bind_rows,data_list)
all<-lapply(all, weight)
#########eastern asia
asia<-subset(all, cohort %in% c(4,5))
#########europe
europe<-subset(all, cohort %in% c(1,2,3))
#########north america
america<-subset(all, cohort %in% c(6,7))
########high-income countries
high_income<-subset(all, cohort %in% c(1,2,3,5,6))
########middle-income countries
middle_income<-subset(all, cohort %in% c(4,7))

data_list<-setNames(list(all,asia,europe,america,high_income,middle_income),
                    c("all","asia","europe","america","high_income","middle_income"))


result_list<-list()
for (i in seq_along(data_list)) {
  data_name<-names(data_list)[[i]]
  result_list[[data_name]]<-analysis1(data_list[[i]],data_name)
}
result<-do.call(bind_rows,result_list)

#########sensitivity analysis7 (dementia defined without imputation)
rm(list = ls())
load("hrs_final_data_noimp.RDATA")
load("MHAS_final_data_noimp.RDATA")
load("charls_final_data_noimp.RDATA")
load("klosa_final_data_noimp.RDATA")
load("share_final_data_noimp.RDATA")
load("elsa_final_data_noimp.RDATA")
load("TILDA_final_data_noimp.RDATA")


######IPTW
data_list<-list(elsa,share,TILDA,charls,klosa,hrs,MHAS)
all<-do.call(bind_rows,data_list)
all<-lapply(all, weight)
#########eastern asia
asia<-subset(all, cohort %in% c(4,5))
#########europe
europe<-subset(all, cohort %in% c(1,2,3))
#########north america
america<-subset(all, cohort %in% c(6,7))
########high-income countries
high_income<-subset(all, cohort %in% c(1,2,3,5,6))
########middle-income countries
middle_income<-subset(all, cohort %in% c(4,7))

data_list<-setNames(list(all,asia,europe,america,high_income,middle_income),
                    c("all","asia","europe","america","high_income","middle_income"))


result_list<-list()
for (i in seq_along(data_list)) {
  data_name<-names(data_list)[[i]]
  result_list[[data_name]]<-analysis1(data_list[[i]],data_name)
}
result<-do.call(bind_rows,result_list)


