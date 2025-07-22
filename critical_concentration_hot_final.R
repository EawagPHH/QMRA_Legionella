library(foreach)
library(doParallel)
library(hrbrthemes)

source("C:/Users/lizha/Desktop/Content.R")
con_critical_hot<-10^seq(from=3,to=10,by=0.02)

cl<-makeCluster(6)
registerDoParallel(cl)

Result_hot_critical <- foreach(j=1:351,.combine = 'rbind',.packages = 'deSolve') %dopar% {
  Dose_hot_conv_overall<-c()
  for (i in 1:1000) {
    state <- c(C_aer_1 = 0,C_aer_2=0,C_aer_3=0,C_aer_4=0,C_aer_5=0,C_aer_6=0,C_aer_7=0,C_aer_8=0,C_aer_9=0,Dose=0)
    differetial <- function(t, state, parms,G,Q_air,numda){
      C_aer_1<-state[1]
      C_aer_2<-state[2]
      C_aer_3<-state[3]
      C_aer_4<-state[4]
      C_aer_5<-state[5]
      C_aer_6<-state[6]
      C_aer_7<-state[7]
      C_aer_8<-state[8]
      C_aer_9<-state[9]
      Dose<-state[10]
      Q_air<-if (t<10) {Ventilation_hot_conv_1[i]} else {Ventilation_hot_conv_2[i]}
      G_1<-if (t<1) {Generation_hot_conv_1a[i]} else if (t>1 && t<10) {Generation_hot_conv_1b[i]} else {0}
      numda_1<-if (t<1) {0} else if (t>1 && t<10) {removal_rate_other_1a_hot_conv[i]} else {removal_rate_other_1b_hot_conv[i]}
      dC_aer_1<-G_1/Volume-Flow/Volume*C_aer_1-Q_air/Volume*C_aer_1-v_floor[1]*Aera/Volume*C_aer_1-v_vertical[1]*(Length*Height*2+Width*Height*2)/Volume*C_aer_1-numda_1*C_aer_1
      G_2<-if (t<1) {Generation_hot_conv_2a[i]} else if (t>1 && t<10) {Generation_hot_conv_2b[i]} else {0}
      numda_2<-if (t<1) {0} else if (t>1 && t<10) {removal_rate_other_2a_hot_conv[i]} else {removal_rate_other_2b_hot_conv[i]}
      dC_aer_2<-G_2/Volume-Flow/Volume*C_aer_2-Q_air/Volume*C_aer_2-v_floor[2]*Aera/Volume*C_aer_2-v_vertical[2]*(Length*Height*2+Width*Height*2)/Volume*C_aer_2-numda_2*C_aer_2
      G_3<-if (t<1) {Generation_hot_conv_3a[i]} else if (t>1 && t<10) {Generation_hot_conv_3b[i]} else {0}
      numda_3<-if (t<1) {0} else if (t>1 && t<10) {removal_rate_other_3a_hot_conv[i]} else {removal_rate_other_3b_hot_conv[i]}
      dC_aer_3<-G_3/Volume-Flow/Volume*C_aer_3-Q_air/Volume*C_aer_3-v_floor[3]*Aera/Volume*C_aer_3-v_vertical[3]*(Length*Height*2+Width*Height*2)/Volume*C_aer_3-numda_3*C_aer_3
      G_4<-if (t<1) {Generation_hot_conv_4a[i]} else if (t>1 && t<10) {Generation_hot_conv_4b[i]} else {0}
      numda_4<-if (t<1) {0} else if (t>1 && t<10) {removal_rate_other_4a_hot_conv[i]} else {removal_rate_other_4b_hot_conv[i]}
      dC_aer_4<-G_4/Volume-Flow/Volume*C_aer_4-Q_air/Volume*C_aer_4-v_floor[4]*Aera/Volume*C_aer_3-v_vertical[4]*(Length*Height*2+Width*Height*2)/Volume*C_aer_4-numda_4*C_aer_4
      G_5<-if (t<1.5) {Generation_hot_conv_5a[i]} else if (t>1.5 && t<10) {Generation_hot_conv_5b[i]} else {0}
      numda_5<-if (t<1) {0} else if (t>1 && t<10) {removal_rate_other_5a_hot_conv[i]} else {removal_rate_other_5b_hot_conv[i]}
      dC_aer_5<-G_5/Volume-Flow/Volume*C_aer_5-Q_air/Volume*C_aer_5-v_floor[5]*Aera/Volume*C_aer_5-v_vertical[5]*(Length*Height*2+Width*Height*2)/Volume*C_aer_5-numda_5*C_aer_5
      G_6<-if (t<1) {Generation_hot_conv_6a[i]} else if (t>1 && t<10) {Generation_hot_conv_6b[i]} else {0}
      numda_6<-if (t<1) {0} else if (t>1 && t<10) {removal_rate_other_6a_hot_conv[i]} else {removal_rate_other_6b_hot_conv[i]}
      dC_aer_6<-G_6/Volume-Flow/Volume*C_aer_6-Q_air/Volume*C_aer_6-v_floor[6]*Aera/Volume*C_aer_6-v_vertical[6]*(Length*Height*2+Width*Height*2)/Volume*C_aer_6-numda_6*C_aer_6
      G_7<-if (t<1) {Generation_hot_conv_7a[i]} else if (t>1 && t<10) {Generation_hot_conv_7b[i]} else {0}
      numda_7<-if (t<1) {0} else if (t>1 && t<10) {removal_rate_other_7a_hot_conv[i]} else {removal_rate_other_7b_hot_conv[i]}
      dC_aer_7<-G_7/Volume-Flow/Volume*C_aer_7-Q_air/Volume*C_aer_7-v_floor[7]*Aera/Volume*C_aer_7-v_vertical[7]*(Length*Height*2+Width*Height*2)/Volume*C_aer_7-numda_7*C_aer_7
      G_8<-if (t<1) {Generation_hot_conv_8a[i]} else if (t>1 && t<10) {Generation_hot_conv_8b[i]} else {0}
      numda_8<-if (t<1) {0} else if (t>1 && t<10) {removal_rate_other_8a_hot_conv[i]} else {removal_rate_other_8b_hot_conv[i]}
      dC_aer_8<-G_8/Volume-Flow/Volume*C_aer_8-Q_air/Volume*C_aer_8-v_floor[8]*Aera/Volume*C_aer_8-v_vertical[8]*(Length*Height*2+Width*Height*2)/Volume*C_aer_8-numda_8*C_aer_8
      G_9<-if (t<1) {Generation_hot_conv_9a[i]} else if (t>1 && t<10) {Generation_hot_conv_9b[i]} else {0}
      numda_9<-if (t<1) {0} else if (t>1 && t<10) {removal_rate_other_9a_hot_conv[i]} else {removal_rate_other_9b_hot_conv[i]}
      dC_aer_9<-G_9/Volume-Flow/Volume*C_aer_9-Q_air/Volume*C_aer_9-v_floor[9]*Aera/Volume*C_aer_9-v_vertical[9]*(Length*Height*2+Width*Height*2)/Volume*C_aer_9-numda_9*C_aer_9
      
      dDose<-(C_aer_1*Fraction[1]*DE1[i]+C_aer_2*Fraction[2]*DE2[i]+C_aer_3*Fraction[3]*DE3[i]+C_aer_4*Fraction[4]*DE4[i]+C_aer_5*Fraction[5]*DE5[i]+C_aer_6*Fraction[6]*DE6[i]+C_aer_7*Fraction[7]*DE7[i]+C_aer_8*Fraction[8]*DE8[i]+C_aer_9*Fraction[9]*DE9[i])*con_critical_hot[j]*Inhalation[i]
      return(list(c(dC_aer_1,dC_aer_2,dC_aer_3,dC_aer_4,dC_aer_5,dC_aer_6,dC_aer_7,dC_aer_8,dC_aer_9,dDose)))
    }
    
    out_hot_conv <- ode.1D(y = state, times = time, func = differetial, parms = NULL,nspec=1)
    Dose_hot_conv<-out_hot_conv[duration/time_step+1,"Dose"]
    Dose_hot_conv_overall<-c(Dose_hot_conv_overall,Dose_hot_conv)
    
  }
  Risk_hot_conv<-1-exp(-r_inf*Dose_hot_conv_overall)
  Risk_hot_conv_annual<-1-(1-Risk_hot_conv)^365
  return(Risk_hot_conv_annual)
}

Hot_critical<-as.data.frame(Result_hot_critical)
Hot_critical<-log10(Hot_critical)

critical_hot<-foreach(i=1:351,.combine = 'c') %dopar% {
  number_hot<-sum(Hot_critical[i,]>(-4.1) & Hot_critical[i,]<(-3.9))
  return(number_hot)
}

data_critical_hot<-c()
for (i in 1:351) {
  data_critical_hot<-c(data_critical_hot,rep(con_critical_hot[i],critical_hot[i]))
}
data_critical_hot<-as.data.frame(data_critical_hot)
ggplot(data_critical_hot, aes(x=data_critical_hot))+geom_histogram(aes(y=..density..), colour="black", fill="pink")+geom_density(alpha=0.2, fill="pink")+scale_x_continuous(trans="log")




data_critical_cold[,2]<-rep("Cold water",9999)
colnames(data_critical_cold)<-c("Concentration","Temperature")
data_critical_hot[,2]<-rep("Hot water",9998)
colnames(data_critical_hot)<-c("Concentration","Temperature")
data_critical_overall<-rbind(data_critical_cold,data_critical_hot)
data_critical_overall[,1]<-data_critical_overall[,1]/1000
ggplot(data_critical_overall, aes(x=Concentration, fill=Temperature))+geom_histogram(aes(y=..density..), alpha=0.5, position="identity",colour="black",binwidth = .1)+geom_density(alpha=0.4)+xlab("Critical concentration (CFU/L)")+scale_x_continuous(transform = "log10")+scale_fill_manual(values=c("lightblue", "pink"))+theme(panel.grid.major = element_blank(), 
                                                                                                                                                                                                                                                         panel.grid.minor = element_blank(),panel.background = element_blank(),
                                                                                                                                                                                                                                                         panel.border = element_rect(colour = "black", fill=NA, size=2),text = element_text(size = 18))


stopCluster(cl)

