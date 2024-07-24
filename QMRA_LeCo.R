library(deSolve)
library(ggplot2)
library(ReacTran)
library(reshape2)
library(mc2d)
scientific_10 <- function(x) {
  parse(text=gsub("e", " %*% 10^", scales::scientific_format()(x)))
}
RH<-0.5 #Relative humitidy
Tem<-298 #Temperature
Length<-1.1 #Length of shower stall(m)
Width<-1.2 #Width of shower stall(m)
Height<-2.9 #Height of shower stall(m)
Aera<-Length*Width #Cross sectional aera of shower stall(m)
Volume<-Aera*Height #Volume of shower stall(m3)
Flow<-0.015 # flow rate of APS (m3/min)
Ventilation<-1.56155 #ventiltation rate (m3/min)
path<-0.066e-6 #mean free path of air (m)
vis<-1.57e-5 #Air viscosity
v_air<-(Ventilation+Flow)/Aera #Mean air velocity (m/min)
Dh<-2*Aera/(Length+Width) #Hydraulic diameter (m)
Re_air<-v_air/60*Dh/vis #Reynolds number
Cf<-0.027/(Re_air^(1/7)) #skin friction coefficient
u_friction<-(Cf/2)^0.5*v_air #friction velocity (m/min)


Dp<-c(1.5e-6,2.5e-6,3.5e-6,4.5e-6,5.5e-6,6.5e-6,7.5e-6,8.5e-6,9.5e-6) #aerosol diameter
Cc<-1+2*path/Dp*(1.257+0.4*exp(-1.1*Dp/2/path))
D<-(1.38e-23)*Tem*Cc/(3*3.14*1.81e-5*Dp)
Sc<-vis/D
r<-Dp*u_friction/60/2/vis
a<-0.5*log((10.92*Sc^(-1/3)+4.3)^3/(Sc^(-1)+0.0609))+3^0.5*atan((8.6-10.92*Sc^(-1/3))/(3^0.5*10.92*Sc^(-1/3)))
b<-0.5*log((10.92*Sc^(-1/3)+r)^3/(Sc^(-1)+7.669e-4*r^3))+3^0.5*atan((2*r-10.92*Sc^(-1/3))/(3^0.5*10.92*Sc^(-1/3)))
I<-3.64*Sc*(2/3)*(a-b)+39
v_vertical<-u_friction/I #deposition velocity on vertical surface
vs<-Dp^2*1000*Cc*9.81/18/1.81e-5*60 #settling velocity
v_floor=vs/(1-exp(-vs*I/u_friction)) #deposition velocity on floor

set.seed(100)
# Parameter bacteria
release_constant<-1.3
Inlet_con<-1000*rlnorm(1000,meanlog=4.414036,sdlog = 3.931661)
Init_con<-1000*rlnorm(1000,meanlog=3.913390,sdlog=3.312003)

Inhalation<-runif(1000,min=0.013,max=0.017)
Fraction<-c(0.175,0.1639,0.1556,0.0667,0.0389,0.025,0.0278,0.05,0.0528)
DE1<-runif(1000,min=0.23,max=0.25)
DE2<-runif(1000,min=0.4,max=0.53)
DE3<-runif(1000,min=0.36,max=0.62)
DE4<-runif(1000,min=0.29,max=0.61)
DE5<-runif(1000,min=0.19,max=0.52)
DE6<-runif(1000,min=0.1,max=0.4)
DE7<-runif(1000,min=0.06,max=0.29)
DE8<-runif(1000,min=0.03,max=0.19)
DE9<-runif(1000,min=0.01,max=0.12)
DE<-cbind(DE1=DE1,DE2=DE2,DE3=DE3,DE4=DE4,DE5=DE5,DE6=DE6,DE7=DE7,DE8=DE8,DE9=DE9)
r_inf<-rlnorm(1000,meanlog=-2.93,sdlog=0.49)

duration<-15
time_step<-0.01
time <- seq(from=0, to=duration, by = time_step)

#demographic parameter
Morbidity_child<-0.0077
Morbidity_adult<-0.097
Morbidity_elderly<-0.75


#cold water conventional showerhead
Generation_cold_conv_1<-runif(1000,min=2.57e-11,max=3.91e-11)
Generation_cold_conv_2<-runif(1000,min=3.4e-11,max=7.56e-11)
Generation_cold_conv_3<-runif(1000,min=3.95e-11,max=1.5e-10)
Generation_cold_conv_4<-runif(1000,min=3.01e-11,max=1.21e-10)
Generation_cold_conv_5<-runif(1000,min=2.54e-11,max=2.07e-10)

removal_rate_other_1a_cold_conv<-runif(1000,min=-0.07,max=0.216)
removal_rate_other_2a_cold_conv<-runif(1000,min=-0.124,max=0.058)
removal_rate_other_3a_cold_conv<-runif(1000,min=-0.11,max=0.1)
removal_rate_other_4a_cold_conv<-runif(1000,min=-0.16672,max=0.013)
removal_rate_other_5a_cold_conv<-runif(1000,min=-0.13,max=0.106)

removal_rate_other_1b_cold_conv<-runif(1000,min=0.03,max=0.45)
removal_rate_other_2b_cold_conv<-runif(1000,min=0.11,max=0.654)
removal_rate_other_3b_cold_conv<-runif(1000,min=0.196,max=0.93)
removal_rate_other_4b_cold_conv<-runif(1000,min=0.22,max=1.388)
removal_rate_other_5b_cold_conv<-runif(1000,min=0.318,max=2.165)

Ventilation_cold_conv_1<-runif(1000,min=0.38187,max=0.885)
Ventilation_cold_conv_2<-runif(1000,min=0.291,max=0.95)

Risk_cold_conv_overall<-numeric()
Risk_cold_conv_child_ill<-c()
Risk_cold_conv_adult_ill<-c()
Risk_cold_conv_elderly_ill<-c()
for (i in 1:1000) {
state <- c(C_aer_1 = 0,C_aer_2=0,C_aer_3=0,C_aer_4=0,C_aer_5=0,C_leg=Init_con[i],Dose=0)
differetial <- function(t, state, parms,G,Q_air,numda){
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
  dC_leg<-release_constant*(Inlet_con[i]-C_leg)
  dDose<-(C_aer_1*Fraction[1]*DE1[i]+C_aer_2*Fraction[2]*DE2[i]+C_aer_3*Fraction[3]*DE3[i]+C_aer_4*Fraction[4]*DE4[i]+C_aer_5*Fraction[5]*DE5[i])*C_leg*Inhalation[i]
  return(list(c(dC_aer_1,dC_aer_2,dC_aer_3,dC_aer_4,dC_aer_5,dC_leg,dDose)))
}

out_cold_conv <- ode.1D(y = state, times = time, func = differetial, parms = NULL,nspec=1)
Dose_cold_conv<-out_cold_conv[,"Dose"]
Risk_cold_conv<-1-exp(-r_inf[i]*Dose_cold_conv)
Risk_cold_conv_child<-Morbidity_child*Risk_cold_conv[duration/time_step+1]
Risk_cold_conv_adult<-Morbidity_adult*Risk_cold_conv[duration/time_step+1]
Risk_cold_conv_elderly<-Morbidity_elderly*Risk_cold_conv[duration/time_step+1]
Risk_cold_conv_annual<-1-(1-Risk_cold_conv)^365
Risk_cold_conv_child_annual<-1-(1-Risk_cold_conv_child)^365
Risk_cold_conv_adult_annual<-1-(1-Risk_cold_conv_adult)^365
Risk_cold_conv_elderly_annual<-1-(1-Risk_cold_conv_elderly)^365
Risk_cold_conv_child_ill<-c(Risk_cold_conv_child_ill,Risk_cold_conv_child_annual)
Risk_cold_conv_adult_ill<-c(Risk_cold_conv_adult_ill,Risk_cold_conv_adult_annual)
Risk_cold_conv_elderly_ill<-c(Risk_cold_conv_elderly_ill,Risk_cold_conv_elderly_annual)
Risk_cold_conv_overall<-cbind(Risk_cold_conv_overall,Risk_cold_conv_annual)
}
mean_risk_cold_conv<-apply(Risk_cold_conv_overall,1,median)
percentile_risk_low_cold_conv<-apply(Risk_cold_conv_overall,1,quantile,probs=c(.25))
percentile_risk_high_cold_conv<-apply(Risk_cold_conv_overall,1,quantile,probs=c(.75))
risk_overall_cold_conv<-data.frame(Time=seq(from=0,to=duration,by=time_step),risk_mean=mean_risk_cold_conv,risk_lower=percentile_risk_low_cold_conv,risk_higher=percentile_risk_high_cold_conv,Condition=rep("Cold water",1501),Type=rep("Conventional showerhead",1501))
Final_risk_cold_conv_child<-data.frame(Risk=Risk_cold_conv_child_ill,Group=rep("Child",1000),Condition=rep("Cold water_conventional showerhead",1000))
Final_risk_cold_conv_adult<-data.frame(Risk=Risk_cold_conv_adult_ill,Group=rep("Aldult",1000),Condition=rep("Cold water_conventional showerhead",1000))
Final_risk_cold_conv_elderly<-data.frame(Risk=Risk_cold_conv_elderly_ill,Group=rep("Elderly",1000),Condition=rep("Cold water_conventional showerhead",1000))


#cold water rain showerhead
Generation_cold_rain_1<-runif(1000,min=2.65e-11,max=3.16e-11)
Generation_cold_rain_2<-runif(1000,min=3.5e-11,max=4.63e-11)
Generation_cold_rain_3<-runif(1000,min=4.12e-11,max=5.43e-11)
Generation_cold_rain_4<-runif(1000,min=2.49e-11,max=4.24e-11)
Generation_cold_rain_5<-runif(1000,min=1.58e-11,max=4.35e-11)

removal_rate_other_1a_cold_rain<-runif(1000,min=-0.1648,max=0.10572)
removal_rate_other_2a_cold_rain<-runif(1000,min=-0.17824,max=0.09559)
removal_rate_other_3a_cold_rain<-runif(1000,min=-0.04965,max=0.11029)
removal_rate_other_4a_cold_rain<-runif(1000,min=-0.04413,max=0.28227)
removal_rate_other_5a_cold_rain<-runif(1000,min=-0.05096,max=0.90755)

removal_rate_other_1b_cold_rain<-runif(1000,min=-0.08575,max=0.20062)
removal_rate_other_2b_cold_rain<-runif(1000,min=-0.01866,max=0.2497)
removal_rate_other_3b_cold_rain<-runif(1000,min=-0.04,max=0.232)
removal_rate_other_4b_cold_rain<-runif(1000,min=-0.02627,max=0.31202)
removal_rate_other_5b_cold_rain<-runif(1000,min=-0.02328,max=0.22278)

Ventilation_cold_rain_1<-runif(1000,min=1.06921,max=2.39117)
Ventilation_cold_rain_2<-runif(1000,min=0,max=1)

Risk_cold_rain_overall<-numeric()
Risk_cold_rain_child_ill<-c()
Risk_cold_rain_adult_ill<-c()
Risk_cold_rain_elderly_ill<-c()
for (i in 1:1000) {
  state <- c(C_aer_1 = 0,C_aer_2=0,C_aer_3=0,C_aer_4=0,C_aer_5=0,C_leg=Init_con[i],Dose=0)
  differetial <- function(t, state, parms,G,Q_air,numda){
    C_aer_1<-state[1]
    C_aer_2<-state[2]
    C_aer_3<-state[3]
    C_aer_4<-state[4]
    C_aer_5<-state[5]
    C_leg<-state[6]
    Dose<-state[7]
    Q_air<-if (t<10) {Ventilation_cold_rain_1[i]} else {Ventilation_cold_rain_2[i]}
    G_1<-if (t<10) {Generation_cold_rain_1[i]} else {0}
    numda_1<-if (t<10) {removal_rate_other_1a_cold_rain[i]} else {removal_rate_other_1b_cold_rain[i]}
    dC_aer_1<-G_1/Volume-Flow/Volume*C_aer_1-Q_air/Volume*C_aer_1-v_floor[1]*Aera/Volume*C_aer_1-v_vertical[1]*(Length*Height*2+Width*Height*2)/Volume*C_aer_1-numda_1*C_aer_1
    G_2<-if (t<10) {Generation_cold_rain_2[i]} else {0}
    numda_2<-if (t<10) {removal_rate_other_2a_cold_rain[i]} else {removal_rate_other_2b_cold_rain[i]}
    dC_aer_2<-G_2/Volume-Flow/Volume*C_aer_2-Q_air/Volume*C_aer_2-v_floor[2]*Aera/Volume*C_aer_2-v_vertical[2]*(Length*Height*2+Width*Height*2)/Volume*C_aer_2-numda_2*C_aer_2
    G_3<-if (t<10) {Generation_cold_rain_3[i]} else {0}
    numda_3<-if (t<10) {removal_rate_other_3a_cold_rain[i]} else {removal_rate_other_3b_cold_rain[i]}
    dC_aer_3<-G_3/Volume-Flow/Volume*C_aer_3-Q_air/Volume*C_aer_3-v_floor[3]*Aera/Volume*C_aer_3-v_vertical[3]*(Length*Height*2+Width*Height*2)/Volume*C_aer_3-numda_3*C_aer_3
    G_4<-if (t<10) {Generation_cold_rain_4[i]} else {0}
    numda_4<-if (t<10) {removal_rate_other_4a_cold_rain[i]} else {removal_rate_other_4b_cold_rain[i]}
    dC_aer_4<-G_4/Volume-Flow/Volume*C_aer_4-Q_air/Volume*C_aer_4-v_floor[4]*Aera/Volume*C_aer_3-v_vertical[4]*(Length*Height*2+Width*Height*2)/Volume*C_aer_4-numda_4*C_aer_4
    G_5<-if (t<10) {Generation_cold_rain_5[i]} else {0}
    numda_5<-if (t<10) {removal_rate_other_5a_cold_rain[i]} else {removal_rate_other_5b_cold_rain[i]}
    dC_aer_5<-G_5/Volume-Flow/Volume*C_aer_5-Q_air/Volume*C_aer_5-v_floor[5]*Aera/Volume*C_aer_5-v_vertical[5]*(Length*Height*2+Width*Height*2)/Volume*C_aer_5-numda_5*C_aer_5
    dC_leg<-release_constant*(Inlet_con[i]-C_leg)
    dDose<-(C_aer_1*Fraction[1]*DE1[i]+C_aer_2*Fraction[2]*DE2[i]+C_aer_3*Fraction[3]*DE3[i]+C_aer_4*Fraction[4]*DE4[i]+C_aer_5*Fraction[5]*DE5[i])*C_leg*Inhalation[i]
    return(list(c(dC_aer_1,dC_aer_2,dC_aer_3,dC_aer_4,dC_aer_5,dC_leg,dDose)))
  }
  
  out_cold_rain <- ode.1D(y = state, times = time, func = differetial, parms = NULL,nspec=1)
  Dose_cold_rain<-out_cold_rain[,"Dose"]
  Risk_cold_rain<-1-exp(-r_inf[i]*Dose_cold_rain)
  Risk_cold_rain_child<-Morbidity_child*Risk_cold_rain[duration/time_step+1]
  Risk_cold_rain_adult<-Morbidity_adult*Risk_cold_rain[duration/time_step+1]
  Risk_cold_rain_elderly<-Morbidity_elderly*Risk_cold_rain[duration/time_step+1]
  Risk_cold_rain_annual<-1-(1-Risk_cold_rain)^365
  Risk_cold_rain_child_annual<-1-(1-Risk_cold_rain_child)^365
  Risk_cold_rain_adult_annual<-1-(1-Risk_cold_rain_adult)^365
  Risk_cold_rain_elderly_annual<-1-(1-Risk_cold_rain_elderly)^365
  Risk_cold_rain_child_ill<-c(Risk_cold_rain_child_ill,Risk_cold_rain_child_annual)
  Risk_cold_rain_adult_ill<-c(Risk_cold_rain_adult_ill,Risk_cold_rain_adult_annual)
  Risk_cold_rain_elderly_ill<-c(Risk_cold_rain_elderly_ill,Risk_cold_rain_elderly_annual)
  Risk_cold_rain_overall<-cbind(Risk_cold_rain_overall,Risk_cold_rain_annual)
}
mean_risk_cold_rain<-apply(Risk_cold_rain_overall,1,median)
percentile_risk_low_cold_rain<-apply(Risk_cold_rain_overall,1,quantile,probs=c(.25))
percentile_risk_high_cold_rain<-apply(Risk_cold_rain_overall,1,quantile,probs=c(.75))
risk_overall_cold_rain<-data.frame(Time=seq(from=0,to=duration,by=time_step),risk_mean=mean_risk_cold_rain,risk_lower=percentile_risk_low_cold_rain,risk_higher=percentile_risk_high_cold_rain,Condition=rep("Cold water",1501),Type=rep("Rain showerhead",1501))
Final_risk_cold_rain_child<-data.frame(Risk=Risk_cold_rain_child_ill,Group=rep("Child",1000),Condition=rep("Cold water_rain showerhead",1000))
Final_risk_cold_rain_adult<-data.frame(Risk=Risk_cold_rain_adult_ill,Group=rep("Aldult",1000),Condition=rep("Cold water_rain showerhead",1000))
Final_risk_cold_rain_elderly<-data.frame(Risk=Risk_cold_rain_elderly_ill,Group=rep("Elderly",1000),Condition=rep("Cold water_rain showerhead",1000))


#hot water conventional showerhead
Generation_hot_conv_1a<-runif(1000,min=0.000000000195,max=0.00000000145)
Generation_hot_conv_2a<-runif(1000,min=0.00000000027098,max=0.0000000096)
Generation_hot_conv_3a<-runif(1000,min=0.00000000040459,max=0.000000045)
Generation_hot_conv_4a<-runif(1000,min=0.00000000047,max=0.000000115)
Generation_hot_conv_5a<-runif(1000,min=0.00000000057,max=0.000000323)
Generation_hot_conv_6a<-runif(1000,min=0.0000000002277,max=0.00000024)
Generation_hot_conv_7a<-runif(1000,min=0.00000000017309,max=0.000000215)
Generation_hot_conv_8a<-runif(1000,min=0.00000000014,max=0.00000022)
Generation_hot_conv_9a<-runif(1000,min=0.000000000075043,max=0.000000175)

Generation_hot_conv_1b<-runif(1000,min=0.000000000051779,max=0.00000000024213)
Generation_hot_conv_2b<-runif(1000,min=0.00000000027098,max=0.0000000096)
Generation_hot_conv_3b<-runif(1000,min=0.0000000001176,max=0.00000001434)
Generation_hot_conv_4b<-runif(1000,min=0.00000000016276,max=0.000000026851)
Generation_hot_conv_5b<-runif(1000,min=0.00000000018086,max=0.00000012032)
Generation_hot_conv_6b<-runif(1000,min=0.000000000045372,max=0.000000042755)
Generation_hot_conv_7b<-runif(1000,min=0.000000000025083,max=0.000000072289)
Generation_hot_conv_8b<-runif(1000,min=0.000000000011731,max=0.000000073568)
Generation_hot_conv_9b<-runif(1000,min=0,max=0.000000080369)

removal_rate_other_1a_hot_conv<-runif(1000,min=0.26589,max=1.10439)
removal_rate_other_2a_hot_conv<-runif(1000,min=0.5,max=0.89)
removal_rate_other_3a_hot_conv<-runif(1000,min=0.41,max=0.63656)
removal_rate_other_4a_hot_conv<-runif(1000,min=0.29,max=0.5)
removal_rate_other_5a_hot_conv<-runif(1000,min=0.06,max=1)
removal_rate_other_6a_hot_conv<-runif(1000,min=0.42,max=0.74)
removal_rate_other_7a_hot_conv<-runif(1000,min=0.16,max=0.6)
removal_rate_other_8a_hot_conv<-runif(1000,min=0.048,max=0.34)
removal_rate_other_9a_hot_conv<-runif(1000,min=-0.096,max=0.46235)


removal_rate_other_1b_hot_conv<-runif(1000,min=1.12,max=1.5)
removal_rate_other_2b_hot_conv<-runif(1000,min=1.22,max=2.85)
removal_rate_other_3b_hot_conv<-runif(1000,min=1.67,max=2.5)
removal_rate_other_4b_hot_conv<-runif(1000,min=1.57,max=3.3)
removal_rate_other_5b_hot_conv<-runif(1000,min=1.5,max=3.1)
removal_rate_other_6b_hot_conv<-runif(1000,min=1.2,max=3.2)
removal_rate_other_7b_hot_conv<-runif(1000,min=1.2,max=3.1)
removal_rate_other_8b_hot_conv<-runif(1000,min=1.2,max=2.84)
removal_rate_other_9b_hot_conv<-runif(1000,min=1.2,max=2.5)


Ventilation_hot_conv_1<-runif(1000,min=0.49,max=1.5)
Ventilation_hot_conv_2<-runif(1000,min=0.95,max=1.81)

Risk_hot_conv_overall<-numeric()
Risk_hot_conv_child_ill<-c()
Risk_hot_conv_adult_ill<-c()
Risk_hot_conv_elderly_ill<-c()
for (i in 1:1000) {
  state <- c(C_aer_1 = 0,C_aer_2=0,C_aer_3=0,C_aer_4=0,C_aer_5=0,C_aer_6=0,C_aer_7=0,C_aer_8=0,C_aer_9=0,C_leg=Init_con[i],Dose=0)
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
    G_1<-if (t<1.5) {Generation_hot_conv_1a[i]} else if (t>1.5 && t<10) {Generation_hot_conv_1b[i]} else {0}
    numda_1<-if (t<1.5) {0} else if (t>1.5 && t<10) {removal_rate_other_1a_hot_conv[i]} else {removal_rate_other_1b_hot_conv[i]}
    dC_aer_1<-G_1/Volume-Flow/Volume*C_aer_1-Q_air/Volume*C_aer_1-v_floor[1]*Aera/Volume*C_aer_1-v_vertical[1]*(Length*Height*2+Width*Height*2)/Volume*C_aer_1-numda_1*C_aer_1
    G_2<-if (t<1.5) {Generation_hot_conv_2a[i]} else if (t>1.5 && t<10) {Generation_hot_conv_2b[i]} else {0}
    numda_2<-if (t<1.5) {0} else if (t>1.5 && t<10) {removal_rate_other_2a_hot_conv[i]} else {removal_rate_other_2b_hot_conv[i]}
    dC_aer_2<-G_2/Volume-Flow/Volume*C_aer_2-Q_air/Volume*C_aer_2-v_floor[2]*Aera/Volume*C_aer_2-v_vertical[2]*(Length*Height*2+Width*Height*2)/Volume*C_aer_2-numda_2*C_aer_2
    G_3<-if (t<1.5) {Generation_hot_conv_3a[i]} else if (t>1.5 && t<10) {Generation_hot_conv_3b[i]} else {0}
    numda_3<-if (t<1.5) {0} else if (t>1.5 && t<10) {removal_rate_other_3a_hot_conv[i]} else {removal_rate_other_3b_hot_conv[i]}
    dC_aer_3<-G_3/Volume-Flow/Volume*C_aer_3-Q_air/Volume*C_aer_3-v_floor[3]*Aera/Volume*C_aer_3-v_vertical[3]*(Length*Height*2+Width*Height*2)/Volume*C_aer_3-numda_3*C_aer_3
    G_4<-if (t<1.5) {Generation_hot_conv_4a[i]} else if (t>1.5 && t<10) {Generation_hot_conv_4b[i]} else {0}
    numda_4<-if (t<1.5) {0} else if (t>1.5 && t<10) {removal_rate_other_4a_hot_conv[i]} else {removal_rate_other_4b_hot_conv[i]}
    dC_aer_4<-G_4/Volume-Flow/Volume*C_aer_4-Q_air/Volume*C_aer_4-v_floor[4]*Aera/Volume*C_aer_3-v_vertical[4]*(Length*Height*2+Width*Height*2)/Volume*C_aer_4-numda_4*C_aer_4
    G_5<-if (t<1.5) {Generation_hot_conv_5a[i]} else if (t>1.5 && t<10) {Generation_hot_conv_5b[i]} else {0}
    numda_5<-if (t<1.5) {0} else if (t>1.5 && t<10) {removal_rate_other_5a_hot_conv[i]} else {removal_rate_other_5b_hot_conv[i]}
    dC_aer_5<-G_5/Volume-Flow/Volume*C_aer_5-Q_air/Volume*C_aer_5-v_floor[5]*Aera/Volume*C_aer_5-v_vertical[5]*(Length*Height*2+Width*Height*2)/Volume*C_aer_5-numda_5*C_aer_5
    G_6<-if (t<1.5) {Generation_hot_conv_6a[i]} else if (t>1.5 && t<10) {Generation_hot_conv_6b[i]} else {0}
    numda_6<-if (t<1.5) {0} else if (t>1.5 && t<10) {removal_rate_other_6a_hot_conv[i]} else {removal_rate_other_6b_hot_conv[i]}
    dC_aer_6<-G_6/Volume-Flow/Volume*C_aer_6-Q_air/Volume*C_aer_6-v_floor[6]*Aera/Volume*C_aer_6-v_vertical[6]*(Length*Height*2+Width*Height*2)/Volume*C_aer_6-numda_6*C_aer_6
    G_7<-if (t<1.5) {Generation_hot_conv_7a[i]} else if (t>1.5 && t<10) {Generation_hot_conv_7b[i]} else {0}
    numda_7<-if (t<1.5) {0} else if (t>1.5 && t<10) {removal_rate_other_7a_hot_conv[i]} else {removal_rate_other_7b_hot_conv[i]}
    dC_aer_7<-G_7/Volume-Flow/Volume*C_aer_7-Q_air/Volume*C_aer_7-v_floor[7]*Aera/Volume*C_aer_7-v_vertical[7]*(Length*Height*2+Width*Height*2)/Volume*C_aer_7-numda_7*C_aer_7
    G_8<-if (t<1.5) {Generation_hot_conv_8a[i]} else if (t>1.5 && t<10) {Generation_hot_conv_8b[i]} else {0}
    numda_8<-if (t<1.5) {0} else if (t>1.5 && t<10) {removal_rate_other_8a_hot_conv[i]} else {removal_rate_other_8b_hot_conv[i]}
    dC_aer_8<-G_8/Volume-Flow/Volume*C_aer_8-Q_air/Volume*C_aer_8-v_floor[8]*Aera/Volume*C_aer_8-v_vertical[8]*(Length*Height*2+Width*Height*2)/Volume*C_aer_8-numda_8*C_aer_8
    G_9<-if (t<1.5) {Generation_hot_conv_9a[i]} else if (t>1.5 && t<10) {Generation_hot_conv_9b[i]} else {0}
    numda_9<-if (t<1.5) {0} else if (t>1.5 && t<10) {removal_rate_other_9a_hot_conv[i]} else {removal_rate_other_9b_hot_conv[i]}
    dC_aer_9<-G_9/Volume-Flow/Volume*C_aer_9-Q_air/Volume*C_aer_9-v_floor[9]*Aera/Volume*C_aer_9-v_vertical[9]*(Length*Height*2+Width*Height*2)/Volume*C_aer_9-numda_9*C_aer_9
    
    dC_leg<-release_constant*(Inlet_con[i]-C_leg)
    dDose<-(C_aer_1*Fraction[1]*DE1[i]+C_aer_2*Fraction[2]*DE2[i]+C_aer_3*Fraction[3]*DE3[i]+C_aer_4*Fraction[4]*DE4[i]+C_aer_5*Fraction[5]*DE5[i]+C_aer_6*Fraction[6]*DE6[i]+C_aer_7*Fraction[7]*DE7[i]+C_aer_8*Fraction[8]*DE8[i]+C_aer_9*Fraction[9]*DE9[i])*C_leg*Inhalation[i]
    return(list(c(dC_aer_1,dC_aer_2,dC_aer_3,dC_aer_4,dC_aer_5,dC_aer_6,dC_aer_7,dC_aer_8,dC_aer_9,dC_leg,dDose)))
  }
  
  out_hot_conv <- ode.1D(y = state, times = time, func = differetial, parms = NULL,nspec=1)
  Dose_hot_conv<-out_hot_conv[,"Dose"]
  Risk_hot_conv<-1-exp(-r_inf[i]*Dose_hot_conv)
  Risk_hot_conv_child<-Morbidity_child*Risk_hot_conv[duration/time_step+1]
  Risk_hot_conv_adult<-Morbidity_adult*Risk_hot_conv[duration/time_step+1]
  Risk_hot_conv_elderly<-Morbidity_elderly*Risk_hot_conv[duration/time_step+1]
  Risk_hot_conv_annual<-1-(1-Risk_hot_conv)^365
  Risk_hot_conv_child_annual<-1-(1-Risk_hot_conv_child)^365
  Risk_hot_conv_adult_annual<-1-(1-Risk_hot_conv_adult)^365
  Risk_hot_conv_elderly_annual<-1-(1-Risk_hot_conv_elderly)^365
  Risk_hot_conv_child_ill<-c(Risk_hot_conv_child_ill,Risk_hot_conv_child_annual)
  Risk_hot_conv_adult_ill<-c(Risk_hot_conv_adult_ill,Risk_hot_conv_adult_annual)
  Risk_hot_conv_elderly_ill<-c(Risk_hot_conv_elderly_ill,Risk_hot_conv_elderly_annual)
  Risk_hot_conv_overall<-cbind(Risk_hot_conv_overall,Risk_hot_conv_annual)
}
mean_risk_hot_conv<-apply(Risk_hot_conv_overall,1,median)
percentile_risk_low_hot_conv<-apply(Risk_hot_conv_overall,1,quantile,probs=c(.25))
percentile_risk_high_hot_conv<-apply(Risk_hot_conv_overall,1,quantile,probs=c(.75))
risk_overall_hot_conv<-data.frame(Time=seq(from=0,to=duration,by=time_step),risk_mean=mean_risk_hot_conv,risk_lower=percentile_risk_low_hot_conv,risk_higher=percentile_risk_high_hot_conv,Condition=rep("Hot water",1501),Type=rep("Conventional showerhead",1501))
Final_risk_hot_conv_child<-data.frame(Risk=Risk_hot_conv_child_ill,Group=rep("Child",1000),Condition=rep("Hot water_conventional showerhead",1000))
Final_risk_hot_conv_adult<-data.frame(Risk=Risk_hot_conv_adult_ill,Group=rep("Aldult",1000),Condition=rep("Hot water_conventional showerhead",1000))
Final_risk_hot_conv_elderly<-data.frame(Risk=Risk_hot_conv_elderly_ill,Group=rep("Elderly",1000),Condition=rep("Hot water_conventional showerhead",1000))

#hot water rain showerhead
Generation_hot_rain_1a<-runif(1000,min=1.87e-10,max=1.83e-9)
Generation_hot_rain_2a<-runif(1000,min=0.00000000119,max=0.0000000128)
Generation_hot_rain_3a<-runif(1000,min=0.0000000046,max=6.1e-8)
Generation_hot_rain_4a<-runif(1000,min=8.9e-9,max=1.55e-7)
Generation_hot_rain_5a<-runif(1000,min=2.4e-8,max=3.7e-7)
Generation_hot_rain_6a<-runif(1000,min=1.69e-8,max=0.00000029)
Generation_hot_rain_7a<-runif(1000,min=2.65e-8,max=2.5e-7)
Generation_hot_rain_8a<-runif(1000,min=4.4e-8,max=2.7e-7)
Generation_hot_rain_9a<-runif(1000,min=4.8e-8,max=2.4e-7)

Generation_hot_rain_1b<-runif(1000,min=3.73e-11,max=8.66e-10)
Generation_hot_rain_2b<-runif(1000,min=2.2e-10,max=6.5e-9)
Generation_hot_rain_3b<-runif(1000,min=1.12e-9,max=3.4e-8)
Generation_hot_rain_4b<-runif(1000,min=4.17e-9,max=5.99e-8)
Generation_hot_rain_5b<-runif(1000,min=1.66e-8,max=3.33e-7)
Generation_hot_rain_6b<-runif(1000,min=1.25e-8,max=3.52e-8)
Generation_hot_rain_7b<-runif(1000,min=2.3e-8,max=7.4248e-8)
Generation_hot_rain_8b<-runif(1000,min=3.1581e-8,max=1.3856e-7)
Generation_hot_rain_9b<-runif(1000,min=5.9967e-8,max=0.00000011)

removal_rate_other_1a_hot_rain<-runif(1000,min=1.13956,max=4.15509)
removal_rate_other_2a_hot_rain<-runif(1000,min=1.03927,max=4.02613)
removal_rate_other_3a_hot_rain<-runif(1000,min=0.95776,max=2.95639)
removal_rate_other_4a_hot_rain<-runif(1000,min=0.97491,max=1.84009)
removal_rate_other_5a_hot_rain<-runif(1000,min=0.99042,max=4.19891)
removal_rate_other_6a_hot_rain<-runif(1000,min=0.20946,max=2.28459)
removal_rate_other_7a_hot_rain<-runif(1000,min=0.40944,max=2.17553)
removal_rate_other_8a_hot_rain<-runif(1000,min=0.61541,max=1.91014)
removal_rate_other_9a_hot_rain<-runif(1000,min=0.41268,max=1.8856)


removal_rate_other_1b_hot_rain<-runif(1000,min=0.21243,max=1.45867)
removal_rate_other_2b_hot_rain<-runif(1000,min=1.19236,max=2.75681)
removal_rate_other_3b_hot_rain<-runif(1000,min=1.8,max=2.9)
removal_rate_other_4b_hot_rain<-runif(1000,min=1.8,max=4.1)
removal_rate_other_5b_hot_rain<-runif(1000,min=1.8,max=5.2)
removal_rate_other_6b_hot_rain<-runif(1000,min=1.8,max=5)
removal_rate_other_7b_hot_rain<-runif(1000,min=1.8,max=5.2)
removal_rate_other_8b_hot_rain<-runif(1000,min=1.8,max=5.1)
removal_rate_other_9b_hot_rain<-runif(1000,min=1.8,max=5.8)


Ventilation_hot_rain_1<-runif(1000,min=0.63306,max=0.90507)
Ventilation_hot_rain_2<-runif(1000,min=0.01853,max=0.0624)

Risk_hot_rain_overall<-numeric()
Risk_hot_rain_child_ill<-c()
Risk_hot_rain_adult_ill<-c()
Risk_hot_rain_elderly_ill<-c()
for (i in 1:1000) {
  state <- c(C_aer_1 = 0,C_aer_2=0,C_aer_3=0,C_aer_4=0,C_aer_5=0,C_aer_6=0,C_aer_7=0,C_aer_8=0,C_aer_9=0,C_leg=Init_con[i],Dose=0)
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
    Q_air<-if (t<10) {Ventilation_hot_rain_1[i]} else {Ventilation_hot_rain_2[i]}
    G_1<-if (t<1.5) {Generation_hot_rain_1a[i]} else if (t>1.5 && t<10) {Generation_hot_rain_1b[i]} else {0}
    numda_1<-if (t<1.5) {0} else if (t>1.5 && t<10) {removal_rate_other_1a_hot_rain[i]} else {removal_rate_other_1b_hot_rain[i]}
    dC_aer_1<-G_1/Volume-Flow/Volume*C_aer_1-Q_air/Volume*C_aer_1-v_floor[1]*Aera/Volume*C_aer_1-v_vertical[1]*(Length*Height*2+Width*Height*2)/Volume*C_aer_1-numda_1*C_aer_1
    G_2<-if (t<1.5) {Generation_hot_rain_2a[i]} else if (t>1.5 && t<10) {Generation_hot_rain_2b[i]} else {0}
    numda_2<-if (t<1.5) {0} else if (t>1.5 && t<10) {removal_rate_other_2a_hot_rain[i]} else {removal_rate_other_2b_hot_rain[i]}
    dC_aer_2<-G_2/Volume-Flow/Volume*C_aer_2-Q_air/Volume*C_aer_2-v_floor[2]*Aera/Volume*C_aer_2-v_vertical[2]*(Length*Height*2+Width*Height*2)/Volume*C_aer_2-numda_2*C_aer_2
    G_3<-if (t<1.5) {Generation_hot_rain_3a[i]} else if (t>1.5 && t<10) {Generation_hot_rain_3b[i]} else {0}
    numda_3<-if (t<1.5) {0} else if (t>1.5 && t<10) {removal_rate_other_3a_hot_rain[i]} else {removal_rate_other_3b_hot_rain[i]}
    dC_aer_3<-G_3/Volume-Flow/Volume*C_aer_3-Q_air/Volume*C_aer_3-v_floor[3]*Aera/Volume*C_aer_3-v_vertical[3]*(Length*Height*2+Width*Height*2)/Volume*C_aer_3-numda_3*C_aer_3
    G_4<-if (t<1.5) {Generation_hot_rain_4a[i]} else if (t>1.5 && t<10) {Generation_hot_rain_4b[i]} else {0}
    numda_4<-if (t<1.5) {0} else if (t>1.5 && t<10) {removal_rate_other_4a_hot_rain[i]} else {removal_rate_other_4b_hot_rain[i]}
    dC_aer_4<-G_4/Volume-Flow/Volume*C_aer_4-Q_air/Volume*C_aer_4-v_floor[4]*Aera/Volume*C_aer_3-v_vertical[4]*(Length*Height*2+Width*Height*2)/Volume*C_aer_4-numda_4*C_aer_4
    G_5<-if (t<1.5) {Generation_hot_rain_5a[i]} else if (t>1.5 && t<10) {Generation_hot_rain_5b[i]} else {0}
    numda_5<-if (t<1.5) {0} else if (t>1.5 && t<10) {removal_rate_other_5a_hot_rain[i]} else {removal_rate_other_5b_hot_rain[i]}
    dC_aer_5<-G_5/Volume-Flow/Volume*C_aer_5-Q_air/Volume*C_aer_5-v_floor[5]*Aera/Volume*C_aer_5-v_vertical[5]*(Length*Height*2+Width*Height*2)/Volume*C_aer_5-numda_5*C_aer_5
    G_6<-if (t<1.5) {Generation_hot_rain_6a[i]} else if (t>1.5 && t<10) {Generation_hot_rain_6b[i]} else {0}
    numda_6<-if (t<1.5) {0} else if (t>1.5 && t<10) {removal_rate_other_6a_hot_rain[i]} else {removal_rate_other_6b_hot_rain[i]}
    dC_aer_6<-G_6/Volume-Flow/Volume*C_aer_6-Q_air/Volume*C_aer_6-v_floor[6]*Aera/Volume*C_aer_6-v_vertical[6]*(Length*Height*2+Width*Height*2)/Volume*C_aer_6-numda_6*C_aer_6
    G_7<-if (t<1.5) {Generation_hot_rain_7a[i]} else if (t>1.5 && t<10) {Generation_hot_rain_7b[i]} else {0}
    numda_7<-if (t<1.5) {0} else if (t>1.5 && t<10) {removal_rate_other_7a_hot_rain[i]} else {removal_rate_other_7b_hot_rain[i]}
    dC_aer_7<-G_7/Volume-Flow/Volume*C_aer_7-Q_air/Volume*C_aer_7-v_floor[7]*Aera/Volume*C_aer_7-v_vertical[7]*(Length*Height*2+Width*Height*2)/Volume*C_aer_7-numda_7*C_aer_7
    G_8<-if (t<1.5) {Generation_hot_rain_8a[i]} else if (t>1.5 && t<10) {Generation_hot_rain_8b[i]} else {0}
    numda_8<-if (t<1.5) {0} else if (t>1.5 && t<10) {removal_rate_other_8a_hot_rain[i]} else {removal_rate_other_8b_hot_rain[i]}
    dC_aer_8<-G_8/Volume-Flow/Volume*C_aer_8-Q_air/Volume*C_aer_8-v_floor[8]*Aera/Volume*C_aer_8-v_vertical[8]*(Length*Height*2+Width*Height*2)/Volume*C_aer_8-numda_8*C_aer_8
    G_9<-if (t<1.5) {Generation_hot_rain_9a[i]} else if (t>1.5 && t<10) {Generation_hot_rain_9b[i]} else {0}
    numda_9<-if (t<1.5) {0} else if (t>1.5 && t<10) {removal_rate_other_9a_hot_rain[i]} else {removal_rate_other_9b_hot_rain[i]}
    dC_aer_9<-G_9/Volume-Flow/Volume*C_aer_9-Q_air/Volume*C_aer_9-v_floor[9]*Aera/Volume*C_aer_9-v_vertical[9]*(Length*Height*2+Width*Height*2)/Volume*C_aer_9-numda_9*C_aer_9
    
    dC_leg<-release_constant*(Inlet_con[i]-C_leg)
    dDose<-(C_aer_1*Fraction[1]*DE1[i]+C_aer_2*Fraction[2]*DE2[i]+C_aer_3*Fraction[3]*DE3[i]+C_aer_4*Fraction[4]*DE4[i]+C_aer_5*Fraction[5]*DE5[i]+C_aer_6*Fraction[6]*DE6[i]+C_aer_7*Fraction[7]*DE7[i]+C_aer_8*Fraction[8]*DE8[i]+C_aer_9*Fraction[9]*DE9[i])*C_leg*Inhalation[i]
    return(list(c(dC_aer_1,dC_aer_2,dC_aer_3,dC_aer_4,dC_aer_5,dC_aer_6,dC_aer_7,dC_aer_8,dC_aer_9,dC_leg,dDose)))
  }
  
  out_hot_rain <- ode.1D(y = state, times = time, func = differetial, parms = NULL,nspec=1)
  Dose_hot_rain<-out_hot_rain[,"Dose"]
  Risk_hot_rain<-1-exp(-r_inf[i]*Dose_hot_rain)
  Risk_hot_rain_child<-Morbidity_child*Risk_hot_rain[duration/time_step+1]
  Risk_hot_rain_adult<-Morbidity_adult*Risk_hot_rain[duration/time_step+1]
  Risk_hot_rain_elderly<-Morbidity_elderly*Risk_hot_rain[duration/time_step+1]
  Risk_hot_rain_annual<-1-(1-Risk_hot_rain)^365
  Risk_hot_rain_child_annual<-1-(1-Risk_hot_rain_child)^365
  Risk_hot_rain_adult_annual<-1-(1-Risk_hot_rain_adult)^365
  Risk_hot_rain_elderly_annual<-1-(1-Risk_hot_rain_elderly)^365
  Risk_hot_rain_child_ill<-c(Risk_hot_rain_child_ill,Risk_hot_rain_child_annual)
  Risk_hot_rain_adult_ill<-c(Risk_hot_rain_adult_ill,Risk_hot_rain_adult_annual)
  Risk_hot_rain_elderly_ill<-c(Risk_hot_rain_elderly_ill,Risk_hot_rain_elderly_annual)
  Risk_hot_rain_overall<-cbind(Risk_hot_rain_overall,Risk_hot_rain_annual)
}
mean_risk_hot_rain<-apply(Risk_hot_rain_overall,1,median)
percentile_risk_low_hot_rain<-apply(Risk_hot_rain_overall,1,quantile,probs=c(.25))
percentile_risk_high_hot_rain<-apply(Risk_hot_rain_overall,1,quantile,probs=c(.75))
risk_overall_hot_rain<-data.frame(Time=seq(from=0,to=duration,by=time_step),risk_mean=mean_risk_hot_rain,risk_lower=percentile_risk_low_hot_rain,risk_higher=percentile_risk_high_hot_rain,Condition=rep("Hot water",1501),Type=rep("Rain showerhead",1501))
Final_risk_hot_rain_child<-data.frame(Risk=Risk_hot_rain_child_ill,Group=rep("Child",1000),Condition=rep("Hot water_rain showerhead",1000))
Final_risk_hot_rain_adult<-data.frame(Risk=Risk_hot_rain_adult_ill,Group=rep("Aldult",1000),Condition=rep("Hot water_rain showerhead",1000))
Final_risk_hot_rain_elderly<-data.frame(Risk=Risk_hot_rain_elderly_ill,Group=rep("Elderly",1000),Condition=rep("Hot water_rain showerhead",1000))




Risk_over_time<-rbind(risk_overall_cold_conv,risk_overall_hot_conv,risk_overall_cold_rain,risk_overall_hot_rain)

         

Final_risk_overall<-rbind(Final_risk_cold_rain_child,Final_risk_cold_rain_adult,Final_risk_cold_rain_elderly,Final_risk_hot_rain_child,Final_risk_hot_rain_adult,Final_risk_hot_rain_elderly,Final_risk_cold_conv_child,Final_risk_cold_conv_adult,Final_risk_cold_conv_elderly,Final_risk_hot_conv_child,Final_risk_hot_conv_adult,Final_risk_hot_conv_elderly)

ggplot(Risk_over_time,aes(x=Time,y=risk_mean))+geom_line()+scale_y_continuous(trans='log10',breaks=c(1e-12,1e-10,1e-8,1e-6,1e-4,1e-2,1),label=scientific_10)+scale_x_continuous(breaks = seq(0,duration,b=1))+xlab("Time (min)")+ylab("Cumulative risk of infection")+geom_ribbon(aes(x=Time,ymax=risk_higher,ymin=risk_lower,fill=Condition),linetype="dashed",alpha=0.5,col="black")+scale_fill_brewer(palette="Pastel1",direction=-1)+facet_grid(Condition~Type)+theme(panel.grid.major = element_blank(), 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        panel.grid.minor = element_blank(),panel.background = element_blank(),
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        panel.border = element_rect(colour = "black", fill=NA, size=2),,text=element_text(size=15))

ggplot(Final_risk_overall,aes(Group,Risk,fill=Condition))+stat_boxplot(geom = "errorbar", width = 0.4,position=position_dodge(0.8))+geom_boxplot(outlier.shape = NA,width=0.6,position =position_dodge(0.8))+ylab("Risk of illness")+scale_y_continuous(trans='log10',breaks = c(1e-12,1e-10,1e-8,1e-6,1e-4,1e-2,1),label=scientific_10)+theme(panel.grid.major = element_blank(), 
                                                                                                                                                                                                                                                                                                                                                             panel.grid.minor = element_blank(),panel.background = element_blank(),
                                                                                                                                                                                                                                                                                                                                                             panel.border = element_rect(colour = "black", fill=NA, size=2),,text=element_text(size=15))+scale_fill_brewer(palette="RdBu",direction=-1)
# sensitivity analysis
Spearman_cold_conv<-data.frame(DE1=DE1,DE2=DE2,DE3=DE3,DE4=DE4,DE5=DE5,Inhalation_rate=Inhalation,Dose_response=r_inf,Generation_rate_1=Generation_cold_conv_1,Generation_rate_2=Generation_cold_conv_2,Generation_rate_3=Generation_cold_conv_3,Generation_rate_4=Generation_cold_conv_4,Generation_rate_5=Generation_cold_conv_5,Removal_rate_1a=removal_rate_other_1a_cold_conv,Removal_rate_2a=removal_rate_other_2a_cold_conv,Removal_rate_3a=removal_rate_other_3a_cold_conv,Removal_rate_4a=removal_rate_other_4a_cold_conv,Removal_rate_5a=removal_rate_other_5a_cold_conv,Removal_rate_1b=removal_rate_other_1b_cold_conv,Removal_rate_2b=removal_rate_other_2b_cold_conv,Removal_rate_3b=removal_rate_other_3b_cold_conv,Removal_rate_4b=removal_rate_other_4b_cold_conv,Removal_rate_5b=removal_rate_other_5b_cold_conv,Concentration_init=Init_con,Concentration_inlet=Inlet_con,Ventilation_1=Ventilation_cold_conv_1,Ventilation_2=Ventilation_cold_conv_2,Risk_cold_conv=Risk_cold_conv_overall[duration/time_step+1,])
Spearman_cold_rain<-data.frame(DE1=DE1,DE2=DE2,DE3=DE3,DE4=DE4,DE5=DE5,Inhalation_rate=Inhalation,Dose_response=r_inf,Generation_rate_1=Generation_cold_rain_1,Generation_rate_2=Generation_cold_rain_2,Generation_rate_3=Generation_cold_rain_3,Generation_rate_4=Generation_cold_rain_4,Generation_rate_5=Generation_cold_rain_5,Removal_rate_1a=removal_rate_other_1a_cold_rain,Removal_rate_2a=removal_rate_other_2a_cold_rain,Removal_rate_3a=removal_rate_other_3a_cold_rain,Removal_rate_4a=removal_rate_other_4a_cold_rain,Removal_rate_5a=removal_rate_other_5a_cold_rain,Removal_rate_1b=removal_rate_other_1b_cold_rain,Removal_rate_2b=removal_rate_other_2b_cold_rain,Removal_rate_3b=removal_rate_other_3b_cold_rain,Removal_rate_4b=removal_rate_other_4b_cold_rain,Removal_rate_5b=removal_rate_other_5b_cold_rain,Concentration_init=Init_con,Concentration_inlet=Inlet_con,Ventilation_1=Ventilation_cold_rain_1,Ventilation_2=Ventilation_cold_rain_2,Risk_cold_rain=Risk_cold_rain_overall[duration/time_step+1,])
Spearman_hot_conv<-data.frame(DE1=DE1,DE2=DE2,DE3=DE3,DE4=DE4,DE5=DE5,DE6=DE6,DE7=DE7,DE8=DE8,DE9=DE9,Inhalation_rate=Inhalation,Dose_response=r_inf,Generation_rate_1=Generation_hot_conv_1b,Generation_rate_2=Generation_hot_conv_2b,Generation_rate_3=Generation_hot_conv_3b,Generation_rate_4=Generation_hot_conv_4b,Generation_rate_5=Generation_hot_conv_5b,Generation_rate_6=Generation_hot_conv_6b,Generation_rate_7=Generation_hot_conv_7b,Generation_rate_8=Generation_hot_conv_8b,Generation_rate_9=Generation_hot_conv_9b,Removal_rate_1a=removal_rate_other_1a_hot_conv,Removal_rate_2a=removal_rate_other_2a_hot_conv,Removal_rate_3a=removal_rate_other_3a_hot_conv,Removal_rate_4a=removal_rate_other_4a_hot_conv,Removal_rate_5a=removal_rate_other_5a_hot_conv,Removal_rate_6a=removal_rate_other_6a_hot_conv,Removal_rate_7a=removal_rate_other_7a_hot_conv,Removal_rate_8a=removal_rate_other_8a_hot_conv,Removal_rate_9a=removal_rate_other_9a_hot_conv,Removal_rate_1b=removal_rate_other_1b_hot_conv,Removal_rate_2b=removal_rate_other_2b_hot_conv,Removal_rate_3b=removal_rate_other_3b_hot_conv,Removal_rate_4b=removal_rate_other_4b_hot_conv,Removal_rate_5b=removal_rate_other_5b_hot_conv,Removal_rate_6b=removal_rate_other_6b_hot_conv,Removal_rate_7b=removal_rate_other_7b_hot_conv,Removal_rate_8b=removal_rate_other_8b_hot_conv,Removal_rate_9b=removal_rate_other_9b_hot_conv,Concentration_init=Init_con,Concentration_inlet=Inlet_con,Ventilation_1=Ventilation_hot_conv_1,Ventilation_2=Ventilation_hot_conv_2,Risk_hot_conv=Risk_hot_conv_overall[duration/time_step+1,])
Spearman_hot_rain<-data.frame(DE1=DE1,DE2=DE2,DE3=DE3,DE4=DE4,DE5=DE5,DE6=DE6,DE7=DE7,DE8=DE8,DE9=DE9,Inhalation_rate=Inhalation,Dose_response=r_inf,Generation_rate_1=Generation_hot_rain_1b,Generation_rate_2=Generation_hot_rain_2b,Generation_rate_3=Generation_hot_rain_3b,Generation_rate_4=Generation_hot_rain_4b,Generation_rate_5=Generation_hot_rain_5b,Generation_rate_6=Generation_hot_rain_6b,Generation_rate_7=Generation_hot_rain_7b,Generation_rate_8=Generation_hot_rain_8b,Generation_rate_9=Generation_hot_rain_9b,Removal_rate_1a=removal_rate_other_1a_hot_rain,Removal_rate_2a=removal_rate_other_2a_hot_rain,Removal_rate_3a=removal_rate_other_3a_hot_rain,Removal_rate_4a=removal_rate_other_4a_hot_rain,Removal_rate_5a=removal_rate_other_5a_hot_rain,Removal_rate_6a=removal_rate_other_6a_hot_rain,Removal_rate_7a=removal_rate_other_7a_hot_rain,Removal_rate_8a=removal_rate_other_8a_hot_rain,Removal_rate_9a=removal_rate_other_9a_hot_rain,Removal_rate_1b=removal_rate_other_1b_hot_rain,Removal_rate_2b=removal_rate_other_2b_hot_rain,Removal_rate_3b=removal_rate_other_3b_hot_rain,Removal_rate_4b=removal_rate_other_4b_hot_rain,Removal_rate_5b=removal_rate_other_5b_hot_rain,Removal_rate_6b=removal_rate_other_6b_hot_rain,Removal_rate_7b=removal_rate_other_7b_hot_rain,Removal_rate_8b=removal_rate_other_8b_hot_rain,Removal_rate_9b=removal_rate_other_9b_hot_rain,Concentration_init=Init_con,Concentration_inlet=Inlet_con,Ventilation_1=Ventilation_hot_rain_1,Ventilation_2=Ventilation_hot_rain_2,Risk_hot_conv=Risk_hot_rain_overall[duration/time_step+1,])

Data_spearman_cold_conv<-numeric()
for (i in 1:26){
  corr<-abs(cor(x=Spearman_cold_conv[,i],y=Spearman_cold_conv[,27],method="spearman"))
  Data_spearman_cold_conv<-rbind(Data_spearman_cold_conv,corr)
}
colnames(Data_spearman_cold_conv)<-c("Spearman")
Data_spearman_cold_rain<-numeric()
for (i in 1:26){
  corr<-abs(cor(x=Spearman_cold_rain[,i],y=Spearman_cold_rain[,27],method="spearman"))
  Data_spearman_cold_rain<-rbind(Data_spearman_cold_rain,corr)
}
colnames(Data_spearman_cold_rain)<-c("Spearman")

Data_spearman_hot_conv<-numeric()
for (i in 1:42){
  corr<-abs(cor(x=Spearman_hot_conv[,i],y=Spearman_hot_conv[,43],method="spearman"))
  Data_spearman_hot_conv<-rbind(Data_spearman_hot_conv,corr)
}
colnames(Data_spearman_hot_conv)<-c("Spearman")
Data_spearman_hot_rain<-numeric()
for (i in 1:42){
  corr<-abs(cor(x=Spearman_hot_rain[,i],y=Spearman_hot_rain[,43],method="spearman"))
  Data_spearman_hot_rain<-rbind(Data_spearman_hot_rain,corr)
}
colnames(Data_spearman_hot_rain)<-c("Spearman")

Data_spearman_cold_conv<-as.data.frame(Data_spearman_cold_conv)
Data_spearman_cold_rain<-as.data.frame(Data_spearman_cold_rain)
Data_spearman_hot_conv<-as.data.frame(Data_spearman_hot_conv)
Data_spearman_hot_rain<-as.data.frame(Data_spearman_hot_rain)

Data_spearman_cold_conv$Parameter<-c("DE1","DE2","DE3","DE4","DE5","Inhalation rate","Dose-response parameter","Aerosol generation rate (1-2 ??m)","Aerosol generation rate (2-3 ??m)","Aerosol generation rate (3-4 ??m)","Aerosol generation rate (4-5 ??m)","Aerosol generation rate (5-6 ??m)","Aerosol decay rate (1-2 ??m,shower on)","Aerosol decay rate (2-3 ??m,shower on)","Aerosol decay rate (3-4 ??m,shower on)","Aerosol decay rate (4-5 ??m,shower on)","Aerosol decay rate (5-6 ??m,shower on)","Aerosol decay rate (1-2 ??m,shower off)","Aerosol decay rate (2-3 ??m,shower off)","Aerosol decay rate (3-4 ??m,shower off)","Aerosol decay rate (4-5 ??m,shower off)","Aerosol decay rate (5-6 ??m,shower off)","Initial concentration","Steady state concentration","Ventilation rate (shower on)","Ventilation rate (shower off)")
Data_spearman_cold_conv$Condition<-rep("Cold water_conventional showerhead",26)
Data_spearman_cold_rain$Parameter<-c("DE1","DE2","DE3","DE4","DE5","Inhalation rate","Dose-response parameter","Aerosol generation rate (1-2 ??m)","Aerosol generation rate (2-3 ??m)","Aerosol generation rate (3-4 ??m)","Aerosol generation rate (4-5 ??m)","Aerosol generation rate (5-6 ??m)","Aerosol decay rate (1-2 ??m,shower on)","Aerosol decay rate (2-3 ??m,shower on)","Aerosol decay rate (3-4 ??m,shower on)","Aerosol decay rate (4-5 ??m,shower on)","Aerosol decay rate (5-6 ??m,shower on)","Aerosol decay rate (1-2 ??m,shower off)","Aerosol decay rate (2-3 ??m,shower off)","Aerosol decay rate (3-4 ??m,shower off)","Aerosol decay rate (4-5 ??m,shower off)","Aerosol decay rate (5-6 ??m,shower off)","Initial concentration","Steady state concentration","Ventilation rate (shower on)","Ventilation rate (shower off)")
Data_spearman_cold_rain$Condition<-rep("Cold water_rain showerhead",26)
Data_spearman_hot_conv$Parameter<-c("DE1","DE2","DE3","DE4","DE5","DE6","DE7","DE8","DE9","Inhalation rate","Dose-response parameter","Aerosol generation rate (1-2 ??m)","Aerosol generation rate (2-3 ??m)","Aerosol generation rate (3-4 ??m)","Aerosol generation rate (4-5 ??m)","Aerosol generation rate (5-6 ??m)","Aerosol generation rate (6-7 ??m)","Aerosol generation rate (7-8 ??m)","Aerosol generation rate (8-9 ??m)","Aerosol generation rate (9-10 ??m)","Aerosol decay rate (1-2 ??m,shower on)","Aerosol decay rate (2-3 ??m,shower on)","Aerosol decay rate (3-4 ??m,shower on)","Aerosol decay rate (4-5 ??m,shower on)","Aerosol decay rate (5-6 ??m,shower on)","Aerosol decay rate (6-7 ??m,shower on)","Aerosol decay rate (7-8 ??m,shower on)","Aerosol decay rate (8-9 ??m,shower on)","Aerosol decay rate (9-10 ??m,shower on)","Aerosol decay rate (1-2 ??m,shower off)","Aerosol decay rate (2-3 ??m,shower off)","Aerosol decay rate (3-4 ??m,shower off)","Aerosol decay rate (4-5 ??m,shower off)","Aerosol decay rate (5-6 ??m,shower off)","Aerosol decay rate (6-7 ??m,shower off)","Aerosol decay rate (7-8 ??m,shower off)","Aerosol decay rate (8-9 ??m,shower off)","Aerosol decay rate (9-10 ??m,shower off)","Initial concentration","Steady state concentration","Ventilation rate (shower on)","Ventilation rate (shower off)")
Data_spearman_hot_conv$Condition<-rep("Hot water_conventional showerhead",42)
Data_spearman_hot_rain$Parameter<-c("DE1","DE2","DE3","DE4","DE5","DE6","DE7","DE8","DE9","Inhalation rate","Dose-response parameter","Aerosol generation rate (1-2 ??m)","Aerosol generation rate (2-3 ??m)","Aerosol generation rate (3-4 ??m)","Aerosol generation rate (4-5 ??m)","Aerosol generation rate (5-6 ??m)","Aerosol generation rate (6-7 ??m)","Aerosol generation rate (7-8 ??m)","Aerosol generation rate (8-9 ??m)","Aerosol generation rate (9-10 ??m)","Aerosol decay rate (1-2 ??m,shower on)","Aerosol decay rate (2-3 ??m,shower on)","Aerosol decay rate (3-4 ??m,shower on)","Aerosol decay rate (4-5 ??m,shower on)","Aerosol decay rate (5-6 ??m,shower on)","Aerosol decay rate (6-7 ??m,shower on)","Aerosol decay rate (7-8 ??m,shower on)","Aerosol decay rate (8-9 ??m,shower on)","Aerosol decay rate (9-10 ??m,shower on)","Aerosol decay rate (1-2 ??m,shower off)","Aerosol decay rate (2-3 ??m,shower off)","Aerosol decay rate (3-4 ??m,shower off)","Aerosol decay rate (4-5 ??m,shower off)","Aerosol decay rate (5-6 ??m,shower off)","Aerosol decay rate (6-7 ??m,shower off)","Aerosol decay rate (7-8 ??m,shower off)","Aerosol decay rate (8-9 ??m,shower off)","Aerosol decay rate (9-10 ??m,shower off)","Initial concentration","Steady state concentration","Ventilation rate (shower on)","Ventilation rate (shower off)")
Data_spearman_hot_rain$Condition<-rep("Hot water_rain showerhead",42)
Data_spearman_cold<-rbind(Data_spearman_cold_conv,Data_spearman_cold_rain)
Data_spearman_hot<-rbind(Data_spearman_hot_conv,Data_spearman_hot_rain)


ggplot(Data_spearman_cold,aes(x=Parameter,y=Spearman,fill=Condition))+geom_bar(stat="identity",position=position_dodge())+ylab("Spearman rank correlation coefficient")+coord_flip()+theme(panel.grid.major = element_blank(), 
                                                                                                                                                                                           panel.grid.minor = element_blank(),panel.background = element_blank(),
                                                                                                                                                                                           panel.border = element_rect(colour = "black", fill=NA, size=2),text=element_text(size=15))
ggplot(Data_spearman_hot,aes(x=Parameter,y=Spearman,fill=Condition))+geom_bar(stat="identity",position=position_dodge())+ylab("Spearman rank correlation coefficient")+coord_flip()+theme(panel.grid.major = element_blank(), 
                                                                                                                                                                                          panel.grid.minor = element_blank(),panel.background = element_blank(),
                                                                                                                                                                                          panel.border = element_rect(colour = "black", fill=NA, size=2),text=element_text(size=15))

