library(foreach)
library(doParallel)
library(hrbrthemes)

source("C:/Users/lizha/Desktop/Content.R")
Inlet_leg<-10^seq(from=3,to=10,by=0.2)
Init_leg<-10^seq(from=3,to=10,by=0.2)
combinations<- expand.grid(Inlet_leg = Inlet_leg, Init_leg = Init_leg)

cl<-makeCluster(6)
registerDoParallel(cl)

heatmap_result_hot <- foreach(j=1:nrow(combinations),.combine = 'c',.packages = 'deSolve') %dopar% {
  Dose_hot_conv_overall<-c()
  for (i in 1:100) {
    state <- c(C_aer_1 = 0,C_aer_2=0,C_aer_3=0,C_aer_4=0,C_aer_5=0,C_aer_6=0,C_aer_7=0,C_aer_8=0,C_aer_9=0,C_leg=combinations$Init_leg[j],Dose=0)
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
      C_leg<-state[10]
      Dose<-state[11]
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
      
      dC_leg<-release_constant*(combinations$Inlet_leg[j]-C_leg)
      dDose<-(C_aer_1*Fraction[1]*DE1[i]+C_aer_2*Fraction[2]*DE2[i]+C_aer_3*Fraction[3]*DE3[i]+C_aer_4*Fraction[4]*DE4[i]+C_aer_5*Fraction[5]*DE5[i]+C_aer_6*Fraction[6]*DE6[i]+C_aer_7*Fraction[7]*DE7[i]+C_aer_8*Fraction[8]*DE8[i]+C_aer_9*Fraction[9]*DE9[i])*C_leg*Inhalation[i]
      return(list(c(dC_aer_1,dC_aer_2,dC_aer_3,dC_aer_4,dC_aer_5,dC_aer_6,dC_aer_7,dC_aer_8,dC_aer_9,dC_leg,dDose)))
    }
    
    out_hot_conv <- ode.1D(y = state, times = time, func = differetial, parms = NULL,nspec=1)
    Dose_hot_conv<-out_hot_conv[duration/time_step+1,"Dose"]
    Dose_hot_conv_overall<-c(Dose_hot_conv_overall,Dose_hot_conv)
    
  }
  Risk_hot_conv<-1-exp(-r_inf*Dose_hot_conv_overall)
  Risk_hot_conv_annual<-1-(1-Risk_hot_conv)^365
  Risk_hot_conv_annual_median<-median(Risk_hot_conv_annual)
  return(Risk_hot_conv_annual_median)
}

stopCluster(cl)

Legionella_heatmap <- Legionella_heatmap #imported dataset named Legionella_heatmap
Legionella_heatmap_hot<- dplyr::filter(Legionella_heatmap, Temperature=="Hot")
Legionella_heatmap_hot$Risk<-rep(1e-4,41)
combinations$Risk<-heatmap_result_hot
combinations$Risk_log<-log10(combinations$Risk)
combinations_4<-combinations[combinations$Risk_log>(-4.06) & combinations$Risk_log<(-3.94),]
ggplot(combinations,mapping=aes(x=log10(Init_leg)-3,y=log10(Inlet_leg)-3,fill=log10(Risk)))+geom_tile()+scale_fill_distiller(palette = "OrRd",direction = 1)+theme_ipsum()+geom_point(Legionella_heatmap_hot,mapping=aes(log10(Concentration_init+1),log10(Concentration_steady+1),shape = Location),size=3)+geom_line(combinations_4,mapping=aes(x=log10(Init_leg)-3,y=log10(Inlet_leg)-3))+theme(text=element_text(size=20))
