library(foreach)
library(doParallel)
library(hrbrthemes)
library(ggplot2)

source("C:/Users/lizha/Desktop/Content.R")
Inlet_leg<-10^seq(from=3,to=10,by=0.2)
Init_leg<-10^seq(from=3,to=10,by=0.2)
combinations_cold<- expand.grid(Inlet_leg = Inlet_leg, Init_leg = Init_leg)

cl<-makeCluster(6)
registerDoParallel(cl)

heatmap_result_cold <- foreach(j=1:nrow(combinations_cold),.combine = 'c',.packages = 'deSolve') %dopar% {
  Dose_cold_conv_overall<-c()
  for (i in 1:100) {
    state <- c(C_aer_1 = 0,C_aer_2=0,C_aer_3=0,C_aer_4=0,C_aer_5=0,C_leg=combinations_cold$Init_leg[j],Dose=0)
    differetial <- function(t, state, parms){
      C_aer_1<-state[1]
      C_aer_2<-state[2]
      C_aer_3<-state[3]
      C_aer_4<-state[4]
      C_aer_5<-state[5]
      C_leg<-state[6]
      Dose<-state[7]
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
      dC_leg<-release_constant*(combinations_cold$Inlet_leg[j]-C_leg)
      dDose<-(C_aer_1*Fraction[1]*DE1[i]+C_aer_2*Fraction[2]*DE2[i]+C_aer_3*Fraction[3]*DE3[i]+C_aer_4*Fraction[4]*DE4[i]+C_aer_5*Fraction[5]*DE5[i])*C_leg*Inhalation[i]
      return(list(c(dC_aer_1,dC_aer_2,dC_aer_3,dC_aer_4,dC_aer_5,dC_leg,dDose)))
    }
    
    
    out_cold_conv <- ode.1D(y = state, times = time, func = differetial, parms = NULL,nspec=1)
    Dose_cold_conv<-out_cold_conv[duration/time_step+1,"Dose"]
    Dose_cold_conv_overall<-c(Dose_cold_conv_overall,Dose_cold_conv)
  }
  Risk_cold_conv<-1-exp(-r_inf*Dose_cold_conv_overall)
  Risk_cold_conv_annual<-1-(1-Risk_cold_conv)^365
  Risk_cold_conv_annual_median<-median(Risk_cold_conv_annual)
  return(Risk_cold_conv_annual_median)
}

stopCluster(cl)

combinations_cold$Risk<-heatmap_result_cold
Legionella_heatmap <- read_excel("Legionella_heatmap") #import dataset named Legionella_heatmap
Legionella_heatmap_cold<- dplyr::filter(Legionella_heatmap, Temperature=="Cold")
Legionella_heatmap_cold$Risk<-rep(1e-4,42)
combinations_cold$Risk_log<-log10(combinations_cold$Risk)
combinations_4_cold<-combinations_cold[combinations_cold$Risk_log>(-4.09) & combinations_cold$Risk_log<(-3.91),]
ggplot(combinations_cold,aes(x=log10(Init_leg)-3,y=log10(Inlet_leg)-3,fill=log10(Risk)))+geom_tile()+scale_fill_distiller(palette = "Blues",direction = 1)+theme_ipsum()+geom_point(Legionella_heatmap_cold,mapping=aes(log10(Concentration_init),log10(Concentration_steady),shape  = Location),size=3)+geom_line(combinations_4_cold,mapping=aes(x=log10(Init_leg)-3,y=log10(Inlet_leg)-3))