library(foreach)
library(doParallel)
library(hrbrthemes)
library(ggplot2)

source("C:/Users/lizha/Desktop/Content.R")
con_critical_cold<-10^seq(from=3,to=10,by=0.02)


cl<-makeCluster(6)
registerDoParallel(cl)

Result_cold_critical <- foreach(j=1:351,.combine = 'rbind',.packages = 'deSolve') %dopar% {
  Dose_cold_conv_overall<-c()
  for (i in 1:1000) {
    state <- c(C_aer_1 = 0,C_aer_2=0,C_aer_3=0,C_aer_4=0,C_aer_5=0,Dose=0)
    differetial <- function(t, state, parms){
      C_aer_1<-state[1]
      C_aer_2<-state[2]
      C_aer_3<-state[3]
      C_aer_4<-state[4]
      C_aer_5<-state[5]
      Dose<-state[6]
      Q_air<-if (t<10) {Ventilation_cold_conv_1[i]} else {Ventilation_cold_conv_2[i]}
      G_1<-if (t<10) {Generation_cold_conv_1[i]} else {0}
      numda_1<-if (t<10) {removal_rate_other_1a_cold_conv[i]} else {removal_rate_other_1b_cold_conv[i]}
      dC_aer_1<-G_1/Volume-Flow/Volume*C_aer_1-Q_air/Volume*C_aer_1-v_floor[1]*Aera/Volume*C_aer_1-v_vertical[1]*(Length*Height*2+Width*Height*2)/Volume*C_aer_1-numda_1*C_aer_1
      G_2<-if (t<10) {Generation_cold_conv_2[i]} else {0}
      numda_2<-if (t<10) {removal_rate_other_2a_cold_conv[i]} else {removal_rate_other_2b_cold_conv[i]}
      dC_aer_2<-G_2/Volume-Flow/Volume*C_aer_2-Q_air/Volume*C_aer_2-v_floor[2]*Aera/Volume*C_aer_2-v_vertical[2]*(Length*Height*2+Width*Height*2)/Volume*C_aer_2-numda_2*C_aer_2
      G_3<-if (t<10) {Generation_cold_conv_3[i]} else {0}
      numda_3<-if (t<10) {removal_rate_other_3a_cold_conv[i]} else {removal_rate_other_3b_cold_conv[i]}
      dC_aer_3<-G_3/Volume-Flow/Volume*C_aer_3-Q_air/Volume*C_aer_3-v_floor[3]*Aera/Volume*C_aer_3-v_vertical[3]*(Length*Height*2+Width*Height*2)/Volume*C_aer_3-numda_3*C_aer_3
      G_4<-if (t<10) {Generation_cold_conv_4[i]} else {0}
      numda_4<-if (t<10) {removal_rate_other_4a_cold_conv[i]} else {removal_rate_other_4b_cold_conv[i]}
      dC_aer_4<-G_4/Volume-Flow/Volume*C_aer_4-Q_air/Volume*C_aer_4-v_floor[4]*Aera/Volume*C_aer_3-v_vertical[4]*(Length*Height*2+Width*Height*2)/Volume*C_aer_4-numda_4*C_aer_4
      G_5<-if (t<10) {Generation_cold_conv_5[i]} else {0}
      numda_5<-if (t<10) {removal_rate_other_5a_cold_conv[i]} else {removal_rate_other_5b_cold_conv[i]}
      dC_aer_5<-G_5/Volume-Flow/Volume*C_aer_5-Q_air/Volume*C_aer_5-v_floor[5]*Aera/Volume*C_aer_5-v_vertical[5]*(Length*Height*2+Width*Height*2)/Volume*C_aer_5-numda_5*C_aer_5
      dDose<-(C_aer_1*Fraction[1]*DE1[i]+C_aer_2*Fraction[2]*DE2[i]+C_aer_3*Fraction[3]*DE3[i]+C_aer_4*Fraction[4]*DE4[i]+C_aer_5*Fraction[5]*DE5[i])*con_critical_cold[j]*Inhalation[i]
      return(list(c(dC_aer_1,dC_aer_2,dC_aer_3,dC_aer_4,dC_aer_5,dDose)))
    }
    
    
    out_cold_conv <- ode.1D(y = state, times = time, func = differetial, parms = NULL,nspec=1)
    Dose_cold_conv<-out_cold_conv[duration/time_step+1,"Dose"]
    Dose_cold_conv_overall<-c(Dose_cold_conv_overall,Dose_cold_conv)
  }
  Risk_cold_conv<-1-exp(-r_inf*Dose_cold_conv_overall)
  Risk_cold_conv_annual<-1-(1-Risk_cold_conv)^365
  return(Risk_cold_conv_annual)
}

Cold_critical<-as.data.frame(Result_cold_critical)
Cold_critical<-log10(Cold_critical)

critical<-foreach(i=1:351,.combine = 'c') %dopar% {
  number<-sum(Cold_critical[i,]>(-4.1) & Cold_critical[i,]<(-3.9))
  return(number)
}

data_critical_cold<-c()
for (i in 1:351) {
  data_critical_cold<-c(data_critical_cold,rep(con_critical_cold[i],critical[i]))
}
data_critical_cold<-as.data.frame(data_critical_cold)
ggplot(data_critical_cold, aes(x=data_critical_cold))+
  geom_density(color="darkblue", fill="lightblue")+scale_x_continuous(transform = "log10")


stopCluster(cl)
