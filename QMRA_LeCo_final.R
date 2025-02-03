library(deSolve)
library(ggplot2)
library(ReacTran)
library(reshape2)
library(mc2d)
library(MASS)
scientific_10 <- function(x) {
  parse(text=gsub("e", " %*% 10^", scales::scientific_format()(x)))
}

Tem<-303 #Temperature
Length<-1.2 #Length of shower stall(m)
Width<-0.9 #Width of shower stall(m)
Height<-2.3 #Height of shower stall(m)
Aera<-Length*Width #Cross sectional aera of shower stall(m)
Volume<-Aera*Height #Volume of shower stall(m3)
Flow<-0.015 # flow rate of APS (m3/min)
Air_flow<-0 #Vertical natural air flow rate (m3/min)
path<-0.066e-6 #mean free path of air (m)
vis<-1.57e-5 #Air viscosity
v_air<-(Air_flow+Flow)/Aera #Mean air velocity (m/min)
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
release_constant<-5
######### Import the dataset named Con_leg here ##########
glm.leg<-glm.nb(Flushed~First_draw,data=Con_leg,link='identity')
betas<-coef(glm.leg)
vcov<-vcov(glm.leg)
betas_simulated<-MASS::mvrnorm(1000, betas, vcov)
coe_leg<-unname(betas_simulated[,2])
Intercept_leg<-unname(betas_simulated[,1])
Init_con_cold<-1000*rgamma(1000,shape=0.7,rate=5.912e-5)
Init_con_hot<-1000*rgamma(1000,shape=0.361,rate=1.09e-4)
steady_con_cold<-Init_con_cold*coe_leg+Intercept_leg
steady_con_hot<-Init_con_hot*coe_leg+Intercept_leg


Inhalation<-runif(1000,min=0.013,max=0.017)
Fraction<-c(0.175,0.1639,0.1556,0.0667,0.0389,0.025,0.0278,0.05,0.0528)
DE1<-runif(1000,min=0.23,max=0.53)
DE2<-runif(1000,min=0.36,max=0.62)
DE3<-runif(1000,min=0.29,max=0.62)
DE4<-runif(1000,min=0.19,max=0.61)
DE5<-runif(1000,min=0.1,max=0.52)
DE6<-runif(1000,min=0.03,max=0.4)
DE7<-runif(1000,min=0.03,max=0.29)
DE8<-runif(1000,min=0.01,max=0.19)
DE9<-runif(1000,min=0.01,max=0.12)
DE<-cbind(DE1=DE1,DE2=DE2,DE3=DE3,DE4=DE4,DE5=DE5,DE6=DE6,DE7=DE7,DE8=DE8,DE9=DE9)
r_inf<-rlnorm(1000,meanlog=-2.93,sdlog=0.49)

duration<-15
time_step<-0.01
time <- seq(from=0, to=duration, by = time_step)

#risk characterization parameter
Morbidity_child<-0.000322
Morbidity_adult<-0.012147
Morbidity_elderly<-0.037531
Disease_severity_Pontiac<-0.1
Disease_severity_LD<-0.3
Disease_duration_Pontiac<-3.5
Disease_duration_LD<-21
Pontiac_percentage<-0.95
LD_percentage<-0.05
Mortality_LD<-0.004
Mean_age_paitent<-62.2
Life_expenctancy<-84
YLD<-(Disease_severity_Pontiac*Disease_duration_Pontiac/365*Pontiac_percentage)+(Disease_severity_LD*Disease_duration_LD/365*LD_percentage)
YLL<-(Life_expenctancy-Mean_age_paitent)*Mortality_LD
DALY<-YLL+YLD


#cold water conventional showerhead
Generation_cold_conv_1<-runif(1000,min=2e-11,max=2.3e-11)
Generation_cold_conv_2<-runif(1000,min=2.69e-11,max=4.41e-11)
Generation_cold_conv_3<-runif(1000,min=3.13e-11,max=8.68e-11)
Generation_cold_conv_4<-runif(1000,min=2.47e-11,max=6.83e-11)
Generation_cold_conv_5<-runif(1000,min=1.98e-11,max=1.21e-10)

removal_rate_other_1a_cold_conv<-runif(1000,min=-0.03392,max=0.15142)
removal_rate_other_2a_cold_conv<-runif(1000,min=-0.09078,max=0.00687)
removal_rate_other_3a_cold_conv<-runif(1000,min=-0.02253,max=0.02335)
removal_rate_other_4a_cold_conv<-runif(1000,min=-0.26586,max=0.03944)
removal_rate_other_5a_cold_conv<-runif(1000,min=-0.21247,max=0.06478)

removal_rate_other_1b_cold_conv<-runif(1000,min=0.03294,max=0.45037)
removal_rate_other_2b_cold_conv<-runif(1000,min=0.11459,max=0.72107)
removal_rate_other_3b_cold_conv<-runif(1000,min=0.20259,max=0.91004)
removal_rate_other_4b_cold_conv<-runif(1000,min=0.20903,max=1.10935)
removal_rate_other_5b_cold_conv<-runif(1000,min=0.2875,max=1.15)

Ventilation_cold_conv_1<-runif(1000,min=0.24425,max=0.54668)
Ventilation_cold_conv_2<-runif(1000,min=0.1,max=0.54668)

Dose_cold_conv_overall<-numeric()
for (i in 1:1000) {
  state <- c(C_aer_1 = 0,C_aer_2=0,C_aer_3=0,C_aer_4=0,C_aer_5=0,C_leg=Init_con_cold[i],Dose=0)
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
    dC_leg<-release_constant*(steady_con_cold[i]-C_leg)
    dDose<-(C_aer_1*Fraction[1]*DE1[i]+C_aer_2*Fraction[2]*DE2[i]+C_aer_3*Fraction[3]*DE3[i]+C_aer_4*Fraction[4]*DE4[i]+C_aer_5*Fraction[5]*DE5[i])*C_leg*Inhalation[i]
    return(list(c(dC_aer_1,dC_aer_2,dC_aer_3,dC_aer_4,dC_aer_5,dC_leg,dDose)))
  }
  
  out_cold_conv <- ode.1D(y = state, times = time, func = differetial, parms = NULL,nspec=1)
  Dose_cold_conv<-out_cold_conv[,"Dose"]
  Dose_cold_conv_overall<-cbind(Dose_cold_conv_overall,Dose_cold_conv)
}

Dose_cold_conv_1<-Dose_cold_conv_overall[101,]
Dose_cold_conv_2<-Dose_cold_conv_overall[201,]-Dose_cold_conv_overall[101,]
Dose_cold_conv_3<-Dose_cold_conv_overall[301,]-Dose_cold_conv_overall[201,]
Dose_cold_conv_4<-Dose_cold_conv_overall[401,]-Dose_cold_conv_overall[301,]
Dose_cold_conv_5<-Dose_cold_conv_overall[501,]-Dose_cold_conv_overall[401,]
Dose_cold_conv_6<-Dose_cold_conv_overall[601,]-Dose_cold_conv_overall[501,]
Dose_cold_conv_7<-Dose_cold_conv_overall[701,]-Dose_cold_conv_overall[601,]
Dose_cold_conv_8<-Dose_cold_conv_overall[801,]-Dose_cold_conv_overall[701,]
Dose_cold_conv_9<-Dose_cold_conv_overall[901,]-Dose_cold_conv_overall[801,]
Dose_cold_conv_10<-Dose_cold_conv_overall[1001,]-Dose_cold_conv_overall[901,]
Dose_cold_conv_11<-Dose_cold_conv_overall[1101,]-Dose_cold_conv_overall[1001,]
Dose_cold_conv_12<-Dose_cold_conv_overall[1201,]-Dose_cold_conv_overall[1101,]
Dose_cold_conv_13<-Dose_cold_conv_overall[1301,]-Dose_cold_conv_overall[1201,]
Dose_cold_conv_14<-Dose_cold_conv_overall[1401,]-Dose_cold_conv_overall[1301,]
Dose_cold_conv_15<-Dose_cold_conv_overall[1501,]-Dose_cold_conv_overall[1401,]

Risk_cold_conv_t<-1-exp(-r_inf*t(Dose_cold_conv_overall))
Risk_cold_conv<-t(Risk_cold_conv_t)
Risk_cold_conv_child<-Morbidity_child*Risk_cold_conv[duration/time_step+1,]
Risk_cold_conv_adult<-Morbidity_adult*Risk_cold_conv[duration/time_step+1,]
Risk_cold_conv_elderly<-Morbidity_elderly*Risk_cold_conv[duration/time_step+1,]
Risk_cold_conv_annual<-1-(1-Risk_cold_conv)^365
Risk_cold_conv_child_annual<-1-(1-Risk_cold_conv_child)^365
Risk_cold_conv_adult_annual<-1-(1-Risk_cold_conv_adult)^365
Risk_cold_conv_elderly_annual<-1-(1-Risk_cold_conv_elderly)^365
Risk_cold_conv_DALY<-Risk_cold_conv_annual[duration/time_step+1,]*DALY
Dose_cold_conv_discrete<-rbind(Dose_cold_conv_1,Dose_cold_conv_2,Dose_cold_conv_3,Dose_cold_conv_4,Dose_cold_conv_5,Dose_cold_conv_6,Dose_cold_conv_7,Dose_cold_conv_8,Dose_cold_conv_9,Dose_cold_conv_10,Dose_cold_conv_11,Dose_cold_conv_12,Dose_cold_conv_13,Dose_cold_conv_14,Dose_cold_conv_15)
Risk_cold_conv_discrete_t<-1-exp(-r_inf*t(Dose_cold_conv_discrete))
Risk_cold_conv_discrete<-t(Risk_cold_conv_discrete_t)
Risk_cold_conv_discrete_annual<-1-(1-Risk_cold_conv_discrete)^365

median_risk_cold_conv<-apply(Risk_cold_conv_annual,1,median)
percentile_risk_low_cold_conv<-apply(Risk_cold_conv_annual,1,quantile,probs=c(.25))
percentile_risk_high_cold_conv<-apply(Risk_cold_conv_annual,1,quantile,probs=c(.75))
risk_overall_cold_conv<-data.frame(Time=seq(from=0,to=duration,by=time_step),risk_median=median_risk_cold_conv,risk_lower=percentile_risk_low_cold_conv,risk_higher=percentile_risk_high_cold_conv,Condition=rep("Cold water",1501),Type=rep("Conventional showerhead",1501))
Final_risk_cold_conv_child<-data.frame(Risk=as.vector(Risk_cold_conv_child_annual),Group=rep("Child",1000),Condition=rep("Cold water_conventional showerhead",1000))
Final_risk_cold_conv_adult<-data.frame(Risk=as.vector(Risk_cold_conv_adult_annual),Group=rep("Adult",1000),Condition=rep("Cold water_conventional showerhead",1000))
Final_risk_cold_conv_elderly<-data.frame(Risk=as.vector(Risk_cold_conv_elderly_annual),Group=rep("Elderly",1000),Condition=rep("Cold water_conventional showerhead",1000))
Final_risk_cold_conv_DALY<-data.frame(Risk=as.vector(Risk_cold_conv_DALY),Group=rep("DALY",1000),Condition=rep("Cold water_conventional showerhead",1000))
Risk_cold_conv_discrete_1<-data.frame(Risk=Risk_cold_conv_discrete_annual[1,],Time=rep("A",1000),Condition=rep("Cold water",1000),Type=rep("Conventional showerhead",1000))
Risk_cold_conv_discrete_2<-data.frame(Risk=Risk_cold_conv_discrete_annual[2,],Time=rep("B",1000),Condition=rep("Cold water",1000),Type=rep("Conventional showerhead",1000))
Risk_cold_conv_discrete_3<-data.frame(Risk=Risk_cold_conv_discrete_annual[3,],Time=rep("C",1000),Condition=rep("Cold water",1000),Type=rep("Conventional showerhead",1000))
Risk_cold_conv_discrete_4<-data.frame(Risk=Risk_cold_conv_discrete_annual[4,],Time=rep("D",1000),Condition=rep("Cold water",1000),Type=rep("Conventional showerhead",1000))
Risk_cold_conv_discrete_5<-data.frame(Risk=Risk_cold_conv_discrete_annual[5,],Time=rep("E",1000),Condition=rep("Cold water",1000),Type=rep("Conventional showerhead",1000))
Risk_cold_conv_discrete_6<-data.frame(Risk=Risk_cold_conv_discrete_annual[6,],Time=rep("F",1000),Condition=rep("Cold water",1000),Type=rep("Conventional showerhead",1000))
Risk_cold_conv_discrete_7<-data.frame(Risk=Risk_cold_conv_discrete_annual[7,],Time=rep("G",1000),Condition=rep("Cold water",1000),Type=rep("Conventional showerhead",1000))
Risk_cold_conv_discrete_8<-data.frame(Risk=Risk_cold_conv_discrete_annual[8,],Time=rep("H",1000),Condition=rep("Cold water",1000),Type=rep("Conventional showerhead",1000))
Risk_cold_conv_discrete_9<-data.frame(Risk=Risk_cold_conv_discrete_annual[9,],Time=rep("I",1000),Condition=rep("Cold water",1000),Type=rep("Conventional showerhead",1000))
Risk_cold_conv_discrete_10<-data.frame(Risk=Risk_cold_conv_discrete_annual[10,],Time=rep("J",1000),Condition=rep("Cold water",1000),Type=rep("Conventional showerhead",1000))
Risk_cold_conv_discrete_11<-data.frame(Risk=Risk_cold_conv_discrete_annual[11,],Time=rep("K",1000),Condition=rep("Cold water",1000),Type=rep("Conventional showerhead",1000))
Risk_cold_conv_discrete_12<-data.frame(Risk=Risk_cold_conv_discrete_annual[12,],Time=rep("L",1000),Condition=rep("Cold water",1000),Type=rep("Conventional showerhead",1000))
Risk_cold_conv_discrete_13<-data.frame(Risk=Risk_cold_conv_discrete_annual[13,],Time=rep("M",1000),Condition=rep("Cold water",1000),Type=rep("Conventional showerhead",1000))
Risk_cold_conv_discrete_14<-data.frame(Risk=Risk_cold_conv_discrete_annual[14,],Time=rep("N",1000),Condition=rep("Cold water",1000),Type=rep("Conventional showerhead",1000))
Risk_cold_conv_discrete_15<-data.frame(Risk=Risk_cold_conv_discrete_annual[15,],Time=rep("O",1000),Condition=rep("Cold water",1000),Type=rep("Conventional showerhead",1000))
Risk_cold_conv_discrete_overall<-rbind(Risk_cold_conv_discrete_1,Risk_cold_conv_discrete_2,Risk_cold_conv_discrete_3,Risk_cold_conv_discrete_4,Risk_cold_conv_discrete_5,Risk_cold_conv_discrete_6,Risk_cold_conv_discrete_7,Risk_cold_conv_discrete_8,Risk_cold_conv_discrete_9,Risk_cold_conv_discrete_10,Risk_cold_conv_discrete_11,Risk_cold_conv_discrete_12,Risk_cold_conv_discrete_13,Risk_cold_conv_discrete_14,Risk_cold_conv_discrete_15)


#cold water rain showerhead
Generation_cold_rain_1<-runif(1000,min=1.695E-11,max=2.233E-11)
Generation_cold_rain_2<-runif(1000,min=2.2407E-11,max=3.021E-11)
Generation_cold_rain_3<-runif(1000,min=2.301E-11,max=3.5048E-11)
Generation_cold_rain_4<-runif(1000,min=1.5375E-11,max=2.7E-11)
Generation_cold_rain_5<-runif(1000,min=7.3938E-12,max=3.1116E-11)

removal_rate_other_1a_cold_rain<-runif(1000,min=-0.10211,max=0.16524)
removal_rate_other_2a_cold_rain<-runif(1000,min=-0.17806,max=0.1288)
removal_rate_other_3a_cold_rain<-runif(1000,min=-0.20758,max=0.15475)
removal_rate_other_4a_cold_rain<-runif(1000,min=-0.19983,max=0.27813)
removal_rate_other_5a_cold_rain<-runif(1000,min=-0.18641,max=1.07692)

removal_rate_other_1b_cold_rain<-runif(1000,min=-0.05278,max=0.28631)
removal_rate_other_2b_cold_rain<-runif(1000,min=0.03164,max=0.33589)
removal_rate_other_3b_cold_rain<-runif(1000,min=0.02437,max=0.31603)
removal_rate_other_4b_cold_rain<-runif(1000,min=-0.03146,max=0.39477)
removal_rate_other_5b_cold_rain<-runif(1000,min=-0.02947,max=0.30382)

Ventilation_cold_rain_1<-runif(1000,min=0.67487,max=1.50561)
Ventilation_cold_rain_2<-runif(1000,min=0,max=0.65927)

Dose_cold_rain_overall<-numeric()
for (i in 1:1000) {
  state <- c(C_aer_1 = 0,C_aer_2=0,C_aer_3=0,C_aer_4=0,C_aer_5=0,C_leg=Init_con_cold[i],Dose=0)
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
    dC_leg<-release_constant*(steady_con_cold[i]-C_leg)
    dDose<-(C_aer_1*Fraction[1]*DE1[i]+C_aer_2*Fraction[2]*DE2[i]+C_aer_3*Fraction[3]*DE3[i]+C_aer_4*Fraction[4]*DE4[i]+C_aer_5*Fraction[5]*DE5[i])*C_leg*Inhalation[i]
    return(list(c(dC_aer_1,dC_aer_2,dC_aer_3,dC_aer_4,dC_aer_5,dC_leg,dDose)))
  }
  
  out_cold_rain <- ode.1D(y = state, times = time, func = differetial, parms = NULL,nspec=1)
  Dose_cold_rain<-out_cold_rain[,"Dose"]
  Dose_cold_rain_overall<-cbind(Dose_cold_rain_overall,Dose_cold_rain)
}
Dose_cold_rain_1<-Dose_cold_rain_overall[101,]
Dose_cold_rain_2<-Dose_cold_rain_overall[201,]-Dose_cold_rain_overall[101,]
Dose_cold_rain_3<-Dose_cold_rain_overall[301,]-Dose_cold_rain_overall[201,]
Dose_cold_rain_4<-Dose_cold_rain_overall[401,]-Dose_cold_rain_overall[301,]
Dose_cold_rain_5<-Dose_cold_rain_overall[501,]-Dose_cold_rain_overall[401,]
Dose_cold_rain_6<-Dose_cold_rain_overall[601,]-Dose_cold_rain_overall[501,]
Dose_cold_rain_7<-Dose_cold_rain_overall[701,]-Dose_cold_rain_overall[601,]
Dose_cold_rain_8<-Dose_cold_rain_overall[801,]-Dose_cold_rain_overall[701,]
Dose_cold_rain_9<-Dose_cold_rain_overall[901,]-Dose_cold_rain_overall[801,]
Dose_cold_rain_10<-Dose_cold_rain_overall[1001,]-Dose_cold_rain_overall[901,]
Dose_cold_rain_11<-Dose_cold_rain_overall[1101,]-Dose_cold_rain_overall[1001,]
Dose_cold_rain_12<-Dose_cold_rain_overall[1201,]-Dose_cold_rain_overall[1101,]
Dose_cold_rain_13<-Dose_cold_rain_overall[1301,]-Dose_cold_rain_overall[1201,]
Dose_cold_rain_14<-Dose_cold_rain_overall[1401,]-Dose_cold_rain_overall[1301,]
Dose_cold_rain_15<-Dose_cold_rain_overall[1501,]-Dose_cold_rain_overall[1401,]

Risk_cold_rain_t<-1-exp(-r_inf*t(Dose_cold_rain_overall))
Risk_cold_rain<-t(Risk_cold_rain_t)
Risk_cold_rain_child<-Morbidity_child*Risk_cold_rain[duration/time_step+1,]
Risk_cold_rain_adult<-Morbidity_adult*Risk_cold_rain[duration/time_step+1,]
Risk_cold_rain_elderly<-Morbidity_elderly*Risk_cold_rain[duration/time_step+1,]
Risk_cold_rain_annual<-1-(1-Risk_cold_rain)^365
Risk_cold_rain_child_annual<-1-(1-Risk_cold_rain_child)^365
Risk_cold_rain_adult_annual<-1-(1-Risk_cold_rain_adult)^365
Risk_cold_rain_elderly_annual<-1-(1-Risk_cold_rain_elderly)^365
Risk_cold_rain_DALY<-Risk_cold_rain_annual[duration/time_step+1,]*DALY
Dose_cold_rain_discrete<-rbind(Dose_cold_rain_1,Dose_cold_rain_2,Dose_cold_rain_3,Dose_cold_rain_4,Dose_cold_rain_5,Dose_cold_rain_6,Dose_cold_rain_7,Dose_cold_rain_8,Dose_cold_rain_9,Dose_cold_rain_10,Dose_cold_rain_11,Dose_cold_rain_12,Dose_cold_rain_13,Dose_cold_rain_14,Dose_cold_rain_15)
Risk_cold_rain_discrete_t<-1-exp(-r_inf*t(Dose_cold_rain_discrete))
Risk_cold_rain_discrete<-t(Risk_cold_rain_discrete_t)
Risk_cold_rain_discrete_annual<-1-(1-Risk_cold_rain_discrete)^365

median_risk_cold_rain<-apply(Risk_cold_rain_annual,1,median)
percentile_risk_low_cold_rain<-apply(Risk_cold_rain_annual,1,quantile,probs=c(.25))
percentile_risk_high_cold_rain<-apply(Risk_cold_rain_annual,1,quantile,probs=c(.75))
risk_overall_cold_rain<-data.frame(Time=seq(from=0,to=duration,by=time_step),risk_median=median_risk_cold_rain,risk_lower=percentile_risk_low_cold_rain,risk_higher=percentile_risk_high_cold_rain,Condition=rep("Cold water",1501),Type=rep("Rain showerhead",1501))
Final_risk_cold_rain_child<-data.frame(Risk=as.vector(Risk_cold_rain_child_annual),Group=rep("Child",1000),Condition=rep("Cold water_rain showerhead",1000))
Final_risk_cold_rain_adult<-data.frame(Risk=as.vector(Risk_cold_rain_adult_annual),Group=rep("Adult",1000),Condition=rep("Cold water_rain showerhead",1000))
Final_risk_cold_rain_elderly<-data.frame(Risk=as.vector(Risk_cold_rain_elderly_annual),Group=rep("Elderly",1000),Condition=rep("Cold water_rain showerhead",1000))
Final_risk_cold_rain_DALY<-data.frame(Risk=as.vector(Risk_cold_rain_DALY),Group=rep("DALY",1000),Condition=rep("Cold water_rain showerhead",1000))
Risk_cold_rain_discrete_1<-data.frame(Risk=Risk_cold_rain_discrete_annual[1,],Time=rep("A",1000),Condition=rep("Cold water",1000),Type=rep("Rain showerhead",1000))
Risk_cold_rain_discrete_2<-data.frame(Risk=Risk_cold_rain_discrete_annual[2,],Time=rep("B",1000),Condition=rep("Cold water",1000),Type=rep("Rain showerhead",1000))
Risk_cold_rain_discrete_3<-data.frame(Risk=Risk_cold_rain_discrete_annual[3,],Time=rep("C",1000),Condition=rep("Cold water",1000),Type=rep("Rain showerhead",1000))
Risk_cold_rain_discrete_4<-data.frame(Risk=Risk_cold_rain_discrete_annual[4,],Time=rep("D",1000),Condition=rep("Cold water",1000),Type=rep("Rain showerhead",1000))
Risk_cold_rain_discrete_5<-data.frame(Risk=Risk_cold_rain_discrete_annual[5,],Time=rep("E",1000),Condition=rep("Cold water",1000),Type=rep("Rain showerhead",1000))
Risk_cold_rain_discrete_6<-data.frame(Risk=Risk_cold_rain_discrete_annual[6,],Time=rep("F",1000),Condition=rep("Cold water",1000),Type=rep("Rain showerhead",1000))
Risk_cold_rain_discrete_7<-data.frame(Risk=Risk_cold_rain_discrete_annual[7,],Time=rep("G",1000),Condition=rep("Cold water",1000),Type=rep("Rain showerhead",1000))
Risk_cold_rain_discrete_8<-data.frame(Risk=Risk_cold_rain_discrete_annual[8,],Time=rep("H",1000),Condition=rep("Cold water",1000),Type=rep("Rain showerhead",1000))
Risk_cold_rain_discrete_9<-data.frame(Risk=Risk_cold_rain_discrete_annual[9,],Time=rep("I",1000),Condition=rep("Cold water",1000),Type=rep("Rain showerhead",1000))
Risk_cold_rain_discrete_10<-data.frame(Risk=Risk_cold_rain_discrete_annual[10,],Time=rep("J",1000),Condition=rep("Cold water",1000),Type=rep("Rain showerhead",1000))
Risk_cold_rain_discrete_11<-data.frame(Risk=Risk_cold_rain_discrete_annual[11,],Time=rep("K",1000),Condition=rep("Cold water",1000),Type=rep("Rain showerhead",1000))
Risk_cold_rain_discrete_12<-data.frame(Risk=Risk_cold_rain_discrete_annual[12,],Time=rep("L",1000),Condition=rep("Cold water",1000),Type=rep("Rain showerhead",1000))
Risk_cold_rain_discrete_13<-data.frame(Risk=Risk_cold_rain_discrete_annual[13,],Time=rep("M",1000),Condition=rep("Cold water",1000),Type=rep("Rain showerhead",1000))
Risk_cold_rain_discrete_14<-data.frame(Risk=Risk_cold_rain_discrete_annual[14,],Time=rep("N",1000),Condition=rep("Cold water",1000),Type=rep("Rain showerhead",1000))
Risk_cold_rain_discrete_15<-data.frame(Risk=Risk_cold_rain_discrete_annual[15,],Time=rep("O",1000),Condition=rep("Cold water",1000),Type=rep("Rain showerhead",1000))
Risk_cold_rain_discrete_overall<-rbind(Risk_cold_rain_discrete_1,Risk_cold_rain_discrete_2,Risk_cold_rain_discrete_3,Risk_cold_rain_discrete_4,Risk_cold_rain_discrete_5,Risk_cold_rain_discrete_6,Risk_cold_rain_discrete_7,Risk_cold_rain_discrete_8,Risk_cold_rain_discrete_9,Risk_cold_rain_discrete_10,Risk_cold_rain_discrete_11,Risk_cold_rain_discrete_12,Risk_cold_rain_discrete_13,Risk_cold_rain_discrete_14,Risk_cold_rain_discrete_15)

#hot water conventional showerhead
Generation_hot_conv_1a<-runif(1000,min=2.2104E-10,max=1.02E-09)
Generation_hot_conv_2a<-runif(1000,min=2.0104E-10,max=6.58E-09)
Generation_hot_conv_3a<-runif(1000,min=4.1136E-10,max=3.07E-08)
Generation_hot_conv_4a<-runif(1000,min=6.2326E-10,max=7.64E-08)
Generation_hot_conv_5a<-runif(1000,min=7.8596E-10,max=2.1178E-07)
Generation_hot_conv_6a<-runif(1000,min=1.4247E-10,max=2.1223E-07)
Generation_hot_conv_7a<-runif(1000,min=1.1186E-10,max=1.9265E-07)
Generation_hot_conv_8a<-runif(1000,min=1.2061E-10,max=1.5419E-07)
Generation_hot_conv_9a<-runif(1000,min=7.1876E-11,max=1.6936E-07)

Generation_hot_conv_1b<-runif(1000,min=3.4275E-11,max=5.604E-10)
Generation_hot_conv_2b<-runif(1000,min=5.282E-11,max=2.25E-09)
Generation_hot_conv_3b<-runif(1000,min=8.0189E-11,max=9.44E-09)
Generation_hot_conv_4b<-runif(1000,min=1.0306E-10,max=1.7616E-08)
Generation_hot_conv_5b<-runif(1000,min=9.5484E-11,max=7.5146E-08)
Generation_hot_conv_6b<-runif(1000,min=2.4337E-11,max=3.035E-08)
Generation_hot_conv_7b<-runif(1000,min=1.109E-11,max=5.287E-08)
Generation_hot_conv_8b<-runif(1000,min=7.1793E-12,max=7.1795E-08)
Generation_hot_conv_9b<-runif(1000,min=4.5555E-12,max=4.9834E-08)

removal_rate_other_1a_hot_conv<-runif(1000,min=0.50755,max=1.21256)
removal_rate_other_2a_hot_conv<-runif(1000,min=0.46896,max=0.88404)
removal_rate_other_3a_hot_conv<-runif(1000,min=0.40365,max=0.63496)
removal_rate_other_4a_hot_conv<-runif(1000,min=0.4675,max=0.50088)
removal_rate_other_5a_hot_conv<-runif(1000,min=0.16738,max=0.93573)
removal_rate_other_6a_hot_conv<-runif(1000,min=0.33759,max=0.80233)
removal_rate_other_7a_hot_conv<-runif(1000,min=0.07428,max=0.6891)
removal_rate_other_8a_hot_conv<-runif(1000,min=0.02779,max=0.60142)
removal_rate_other_9a_hot_conv<-runif(1000,min=-0.00834,max=0.41798)


removal_rate_other_1b_hot_conv<-runif(1000,min=1.15664,max=1.9894)
removal_rate_other_2b_hot_conv<-runif(1000,min=1.22065,max=2.82191)
removal_rate_other_3b_hot_conv<-runif(1000,min=1.67278,max=3.06179)
removal_rate_other_4b_hot_conv<-runif(1000,min=1.90191,max=3.29054)
removal_rate_other_5b_hot_conv<-runif(1000,min=1.26614,max=2.99778)
removal_rate_other_6b_hot_conv<-runif(1000,min=0.89358,max=3.05415)
removal_rate_other_7b_hot_conv<-runif(1000,min=0.86361,max=2.89174)
removal_rate_other_8b_hot_conv<-runif(1000,min=0.79677,max=2.79785)
removal_rate_other_9b_hot_conv<-runif(1000,min=0.67654,max=2.47257)


Ventilation_hot_conv_1<-runif(1000,min=0.25,max=0.63)
Ventilation_hot_conv_2<-runif(1000,min=0.6167,max=1.16698)

Dose_hot_conv_overall<-numeric()
for (i in 1:1000) {
  state <- c(C_aer_1 = 0,C_aer_2=0,C_aer_3=0,C_aer_4=0,C_aer_5=0,C_aer_6=0,C_aer_7=0,C_aer_8=0,C_aer_9=0,C_leg=Init_con_hot[i],Dose=0)
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
    
    dC_leg<-release_constant*(steady_con_hot[i]-C_leg)
    dDose<-(C_aer_1*Fraction[1]*DE1[i]+C_aer_2*Fraction[2]*DE2[i]+C_aer_3*Fraction[3]*DE3[i]+C_aer_4*Fraction[4]*DE4[i]+C_aer_5*Fraction[5]*DE5[i]+C_aer_6*Fraction[6]*DE6[i]+C_aer_7*Fraction[7]*DE7[i]+C_aer_8*Fraction[8]*DE8[i]+C_aer_9*Fraction[9]*DE9[i])*C_leg*Inhalation[i]
    return(list(c(dC_aer_1,dC_aer_2,dC_aer_3,dC_aer_4,dC_aer_5,dC_aer_6,dC_aer_7,dC_aer_8,dC_aer_9,dC_leg,dDose)))
  }
  
  out_hot_conv <- ode.1D(y = state, times = time, func = differetial, parms = NULL,nspec=1)
  Dose_hot_conv<-out_hot_conv[,"Dose"]
  Dose_hot_conv_overall<-cbind(Dose_hot_conv_overall,Dose_hot_conv)
}
Dose_hot_conv_1<-Dose_hot_conv_overall[101,]
Dose_hot_conv_2<-Dose_hot_conv_overall[201,]-Dose_hot_conv_overall[101,]
Dose_hot_conv_3<-Dose_hot_conv_overall[301,]-Dose_hot_conv_overall[201,]
Dose_hot_conv_4<-Dose_hot_conv_overall[401,]-Dose_hot_conv_overall[301,]
Dose_hot_conv_5<-Dose_hot_conv_overall[501,]-Dose_hot_conv_overall[401,]
Dose_hot_conv_6<-Dose_hot_conv_overall[601,]-Dose_hot_conv_overall[501,]
Dose_hot_conv_7<-Dose_hot_conv_overall[701,]-Dose_hot_conv_overall[601,]
Dose_hot_conv_8<-Dose_hot_conv_overall[801,]-Dose_hot_conv_overall[701,]
Dose_hot_conv_9<-Dose_hot_conv_overall[901,]-Dose_hot_conv_overall[801,]
Dose_hot_conv_10<-Dose_hot_conv_overall[1001,]-Dose_hot_conv_overall[901,]
Dose_hot_conv_11<-Dose_hot_conv_overall[1101,]-Dose_hot_conv_overall[1001,]
Dose_hot_conv_12<-Dose_hot_conv_overall[1201,]-Dose_hot_conv_overall[1101,]
Dose_hot_conv_13<-Dose_hot_conv_overall[1301,]-Dose_hot_conv_overall[1201,]
Dose_hot_conv_14<-Dose_hot_conv_overall[1401,]-Dose_hot_conv_overall[1301,]
Dose_hot_conv_15<-Dose_hot_conv_overall[1501,]-Dose_hot_conv_overall[1401,]

Risk_hot_conv_t<-1-exp(-r_inf*t(Dose_hot_conv_overall))
Risk_hot_conv<-t(Risk_hot_conv_t)
Risk_hot_conv_child<-Morbidity_child*Risk_hot_conv[duration/time_step+1,]
Risk_hot_conv_adult<-Morbidity_adult*Risk_hot_conv[duration/time_step+1,]
Risk_hot_conv_elderly<-Morbidity_elderly*Risk_hot_conv[duration/time_step+1,]
Risk_hot_conv_annual<-1-(1-Risk_hot_conv)^365
Risk_hot_conv_child_annual<-1-(1-Risk_hot_conv_child)^365
Risk_hot_conv_adult_annual<-1-(1-Risk_hot_conv_adult)^365
Risk_hot_conv_elderly_annual<-1-(1-Risk_hot_conv_elderly)^365
Risk_hot_conv_DALY<-Risk_hot_conv_annual[duration/time_step+1,]*DALY
Dose_hot_conv_discrete<-rbind(Dose_hot_conv_1,Dose_hot_conv_2,Dose_hot_conv_3,Dose_hot_conv_4,Dose_hot_conv_5,Dose_hot_conv_6,Dose_hot_conv_7,Dose_hot_conv_8,Dose_hot_conv_9,Dose_hot_conv_10,Dose_hot_conv_11,Dose_hot_conv_12,Dose_hot_conv_13,Dose_hot_conv_14,Dose_hot_conv_15)
Risk_hot_conv_discrete_t<-1-exp(-r_inf*t(Dose_hot_conv_discrete))
Risk_hot_conv_discrete<-t(Risk_hot_conv_discrete_t)
Risk_hot_conv_discrete_annual<-1-(1-Risk_hot_conv_discrete)^365

median_risk_hot_conv<-apply(Risk_hot_conv_annual,1,median)
percentile_risk_low_hot_conv<-apply(Risk_hot_conv_annual,1,quantile,probs=c(.25))
percentile_risk_high_hot_conv<-apply(Risk_hot_conv_annual,1,quantile,probs=c(.75))
risk_overall_hot_conv<-data.frame(Time=seq(from=0,to=duration,by=time_step),risk_median=median_risk_hot_conv,risk_lower=percentile_risk_low_hot_conv,risk_higher=percentile_risk_high_hot_conv,Condition=rep("Hot water",1501),Type=rep("Conventional showerhead",1501))
Final_risk_hot_conv_child<-data.frame(Risk=as.vector(Risk_hot_conv_child_annual),Group=rep("Child",1000),Condition=rep("Hot water_conventional showerhead",1000))
Final_risk_hot_conv_adult<-data.frame(Risk=as.vector(Risk_hot_conv_adult_annual),Group=rep("Adult",1000),Condition=rep("Hot water_conventional showerhead",1000))
Final_risk_hot_conv_elderly<-data.frame(Risk=as.vector(Risk_hot_conv_elderly_annual),Group=rep("Elderly",1000),Condition=rep("Hot water_conventional showerhead",1000))
Final_risk_hot_conv_DALY<-data.frame(Risk=as.vector(Risk_hot_conv_DALY),Group=rep("DALY",1000),Condition=rep("Hot water_conventional showerhead",1000))
Risk_hot_conv_discrete_1<-data.frame(Risk=Risk_hot_conv_discrete_annual[1,],Time=rep("A",1000),Condition=rep("Hot water",1000),Type=rep("Conventional showerhead",1000))
Risk_hot_conv_discrete_2<-data.frame(Risk=Risk_hot_conv_discrete_annual[2,],Time=rep("B",1000),Condition=rep("Hot water",1000),Type=rep("Conventional showerhead",1000))
Risk_hot_conv_discrete_3<-data.frame(Risk=Risk_hot_conv_discrete_annual[3,],Time=rep("C",1000),Condition=rep("Hot water",1000),Type=rep("Conventional showerhead",1000))
Risk_hot_conv_discrete_4<-data.frame(Risk=Risk_hot_conv_discrete_annual[4,],Time=rep("D",1000),Condition=rep("Hot water",1000),Type=rep("Conventional showerhead",1000))
Risk_hot_conv_discrete_5<-data.frame(Risk=Risk_hot_conv_discrete_annual[5,],Time=rep("E",1000),Condition=rep("Hot water",1000),Type=rep("Conventional showerhead",1000))
Risk_hot_conv_discrete_6<-data.frame(Risk=Risk_hot_conv_discrete_annual[6,],Time=rep("F",1000),Condition=rep("Hot water",1000),Type=rep("Conventional showerhead",1000))
Risk_hot_conv_discrete_7<-data.frame(Risk=Risk_hot_conv_discrete_annual[7,],Time=rep("G",1000),Condition=rep("Hot water",1000),Type=rep("Conventional showerhead",1000))
Risk_hot_conv_discrete_8<-data.frame(Risk=Risk_hot_conv_discrete_annual[8,],Time=rep("H",1000),Condition=rep("Hot water",1000),Type=rep("Conventional showerhead",1000))
Risk_hot_conv_discrete_9<-data.frame(Risk=Risk_hot_conv_discrete_annual[9,],Time=rep("I",1000),Condition=rep("Hot water",1000),Type=rep("Conventional showerhead",1000))
Risk_hot_conv_discrete_10<-data.frame(Risk=Risk_hot_conv_discrete_annual[10,],Time=rep("J",1000),Condition=rep("Hot water",1000),Type=rep("Conventional showerhead",1000))
Risk_hot_conv_discrete_11<-data.frame(Risk=Risk_hot_conv_discrete_annual[11,],Time=rep("K",1000),Condition=rep("Hot water",1000),Type=rep("Conventional showerhead",1000))
Risk_hot_conv_discrete_12<-data.frame(Risk=Risk_hot_conv_discrete_annual[12,],Time=rep("L",1000),Condition=rep("Hot water",1000),Type=rep("Conventional showerhead",1000))
Risk_hot_conv_discrete_13<-data.frame(Risk=Risk_hot_conv_discrete_annual[13,],Time=rep("M",1000),Condition=rep("Hot water",1000),Type=rep("Conventional showerhead",1000))
Risk_hot_conv_discrete_14<-data.frame(Risk=Risk_hot_conv_discrete_annual[14,],Time=rep("N",1000),Condition=rep("Hot water",1000),Type=rep("Conventional showerhead",1000))
Risk_hot_conv_discrete_15<-data.frame(Risk=Risk_hot_conv_discrete_annual[15,],Time=rep("O",1000),Condition=rep("Hot water",1000),Type=rep("Conventional showerhead",1000))
Risk_hot_conv_discrete_overall<-rbind(Risk_hot_conv_discrete_1,Risk_hot_conv_discrete_2,Risk_hot_conv_discrete_3,Risk_hot_conv_discrete_4,Risk_hot_conv_discrete_5,Risk_hot_conv_discrete_6,Risk_hot_conv_discrete_7,Risk_hot_conv_discrete_8,Risk_hot_conv_discrete_9,Risk_hot_conv_discrete_10,Risk_hot_conv_discrete_11,Risk_hot_conv_discrete_12,Risk_hot_conv_discrete_13,Risk_hot_conv_discrete_14,Risk_hot_conv_discrete_15)

#hot water rain showerhead
Generation_hot_rain_1a<-runif(1000,min=1.7321E-10,max=1.3288E-09)
Generation_hot_rain_2a<-runif(1000,min=9.5476E-10,max=9.0608E-09)
Generation_hot_rain_3a<-runif(1000,min=3.2477E-09,max=5.0365E-08)
Generation_hot_rain_4a<-runif(1000,min=8.2171E-09,max=1.3585E-07)
Generation_hot_rain_5a<-runif(1000,min=1.4606E-08,max=3.1964E-07)
Generation_hot_rain_6a<-runif(1000,min=1.0996E-08,max=1.7492E-07)
Generation_hot_rain_7a<-runif(1000,min=1.7714E-08,max=1.5338E-07)
Generation_hot_rain_8a<-runif(1000,min=2.9283E-08,max=1.7254E-07)
Generation_hot_rain_9a<-runif(1000,min=3.1438E-08,max=1.6774E-07)

Generation_hot_rain_1b<-runif(1000,min=2.5915E-11,max=5.7689E-10)
Generation_hot_rain_2b<-runif(1000,min=1.5191E-10,max=4.6506E-09)
Generation_hot_rain_3b<-runif(1000,min=8.1961E-10,max=2.2612E-08)
Generation_hot_rain_4b<-runif(1000,min=3.0104E-09,max=3.643E-08)
Generation_hot_rain_5b<-runif(1000,min=9.7552E-09,max=1.9205E-07)
Generation_hot_rain_6b<-runif(1000,min=8.8716E-09,max=1.8682E-08)
Generation_hot_rain_7b<-runif(1000,min=1.6535E-08,max=5.0233E-08)
Generation_hot_rain_8b<-runif(1000,min=2.1309E-08,max=8.3122E-08)
Generation_hot_rain_9b<-runif(1000,min=2.3829E-08,max=5.4671E-08)

removal_rate_other_1a_hot_rain<-runif(1000,min=1.18591,max=5.26714)
removal_rate_other_2a_hot_rain<-runif(1000,min=1.09574,max=8.33159)
removal_rate_other_3a_hot_rain<-runif(1000,min=1.07225,max=6.55227)
removal_rate_other_4a_hot_rain<-runif(1000,min=1.00537,max=5.30007)
removal_rate_other_5a_hot_rain<-runif(1000,min=0.88638,max=3.23186)
removal_rate_other_6a_hot_rain<-runif(1000,min=0.15746,max=2.27229)
removal_rate_other_7a_hot_rain<-runif(1000,min=0.434,max=2.17553)
removal_rate_other_8a_hot_rain<-runif(1000,min=0.55297,max=2.55994)
removal_rate_other_9a_hot_rain<-runif(1000,min=0.39294,max=2.16557)


removal_rate_other_1b_hot_rain<-runif(1000,min=0.21866,max=1.13216)
removal_rate_other_2b_hot_rain<-runif(1000,min=1.2167,max=1.57934)
removal_rate_other_3b_hot_rain<-runif(1000,min=1.59455,max=3.93018)
removal_rate_other_4b_hot_rain<-runif(1000,min=1.6216,max=6.98146)
removal_rate_other_5b_hot_rain<-runif(1000,min=1.51361,max=12.57266)
removal_rate_other_6b_hot_rain<-runif(1000,min=1.64361,max=13.13833)
removal_rate_other_7b_hot_rain<-runif(1000,min=1.64411,max=12.83658)
removal_rate_other_8b_hot_rain<-runif(1000,min=1.59212,max=12.52764)
removal_rate_other_9b_hot_rain<-runif(1000,min=1.56964,max=13.39847)


Ventilation_hot_rain_1<-runif(1000,min=0.3715,max=1.13751)
Ventilation_hot_rain_2<-runif(1000,min=0.03451,max=1.13751)

Dose_hot_rain_overall<-numeric()
for (i in 1:1000) {
  state <- c(C_aer_1 = 0,C_aer_2=0,C_aer_3=0,C_aer_4=0,C_aer_5=0,C_aer_6=0,C_aer_7=0,C_aer_8=0,C_aer_9=0,C_leg=Init_con_hot[i],Dose=0)
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
    G_1<-if (t<1) {Generation_hot_rain_1a[i]} else if (t>1 && t<10) {Generation_hot_rain_1b[i]} else {0}
    numda_1<-if (t<1) {0} else if (t>1 && t<10) {removal_rate_other_1a_hot_rain[i]} else {removal_rate_other_1b_hot_rain[i]}
    dC_aer_1<-G_1/Volume-Flow/Volume*C_aer_1-Q_air/Volume*C_aer_1-v_floor[1]*Aera/Volume*C_aer_1-v_vertical[1]*(Length*Height*2+Width*Height*2)/Volume*C_aer_1-numda_1*C_aer_1
    G_2<-if (t<1) {Generation_hot_rain_2a[i]} else if (t>1 && t<10) {Generation_hot_rain_2b[i]} else {0}
    numda_2<-if (t<1) {0} else if (t>1 && t<10) {removal_rate_other_2a_hot_rain[i]} else {removal_rate_other_2b_hot_rain[i]}
    dC_aer_2<-G_2/Volume-Flow/Volume*C_aer_2-Q_air/Volume*C_aer_2-v_floor[2]*Aera/Volume*C_aer_2-v_vertical[2]*(Length*Height*2+Width*Height*2)/Volume*C_aer_2-numda_2*C_aer_2
    G_3<-if (t<1) {Generation_hot_rain_3a[i]} else if (t>1 && t<10) {Generation_hot_rain_3b[i]} else {0}
    numda_3<-if (t<1) {0} else if (t>1 && t<10) {removal_rate_other_3a_hot_rain[i]} else {removal_rate_other_3b_hot_rain[i]}
    dC_aer_3<-G_3/Volume-Flow/Volume*C_aer_3-Q_air/Volume*C_aer_3-v_floor[3]*Aera/Volume*C_aer_3-v_vertical[3]*(Length*Height*2+Width*Height*2)/Volume*C_aer_3-numda_3*C_aer_3
    G_4<-if (t<1) {Generation_hot_rain_4a[i]} else if (t>1 && t<10) {Generation_hot_rain_4b[i]} else {0}
    numda_4<-if (t<1) {0} else if (t>1 && t<10) {removal_rate_other_4a_hot_rain[i]} else {removal_rate_other_4b_hot_rain[i]}
    dC_aer_4<-G_4/Volume-Flow/Volume*C_aer_4-Q_air/Volume*C_aer_4-v_floor[4]*Aera/Volume*C_aer_3-v_vertical[4]*(Length*Height*2+Width*Height*2)/Volume*C_aer_4-numda_4*C_aer_4
    G_5<-if (t<1) {Generation_hot_rain_5a[i]} else if (t>1 && t<10) {Generation_hot_rain_5b[i]} else {0}
    numda_5<-if (t<1) {0} else if (t>1 && t<10) {removal_rate_other_5a_hot_rain[i]} else {removal_rate_other_5b_hot_rain[i]}
    dC_aer_5<-G_5/Volume-Flow/Volume*C_aer_5-Q_air/Volume*C_aer_5-v_floor[5]*Aera/Volume*C_aer_5-v_vertical[5]*(Length*Height*2+Width*Height*2)/Volume*C_aer_5-numda_5*C_aer_5
    G_6<-if (t<1) {Generation_hot_rain_6a[i]} else if (t>1 && t<10) {Generation_hot_rain_6b[i]} else {0}
    numda_6<-if (t<1) {0} else if (t>1 && t<10) {removal_rate_other_6a_hot_rain[i]} else {removal_rate_other_6b_hot_rain[i]}
    dC_aer_6<-G_6/Volume-Flow/Volume*C_aer_6-Q_air/Volume*C_aer_6-v_floor[6]*Aera/Volume*C_aer_6-v_vertical[6]*(Length*Height*2+Width*Height*2)/Volume*C_aer_6-numda_6*C_aer_6
    G_7<-if (t<1) {Generation_hot_rain_7a[i]} else if (t>1 && t<10) {Generation_hot_rain_7b[i]} else {0}
    numda_7<-if (t<1) {0} else if (t>1 && t<10) {removal_rate_other_7a_hot_rain[i]} else {removal_rate_other_7b_hot_rain[i]}
    dC_aer_7<-G_7/Volume-Flow/Volume*C_aer_7-Q_air/Volume*C_aer_7-v_floor[7]*Aera/Volume*C_aer_7-v_vertical[7]*(Length*Height*2+Width*Height*2)/Volume*C_aer_7-numda_7*C_aer_7
    G_8<-if (t<1) {Generation_hot_rain_8a[i]} else if (t>1 && t<10) {Generation_hot_rain_8b[i]} else {0}
    numda_8<-if (t<1) {0} else if (t>1 && t<10) {removal_rate_other_8a_hot_rain[i]} else {removal_rate_other_8b_hot_rain[i]}
    dC_aer_8<-G_8/Volume-Flow/Volume*C_aer_8-Q_air/Volume*C_aer_8-v_floor[8]*Aera/Volume*C_aer_8-v_vertical[8]*(Length*Height*2+Width*Height*2)/Volume*C_aer_8-numda_8*C_aer_8
    G_9<-if (t<1) {Generation_hot_rain_9a[i]} else if (t>1 && t<10) {Generation_hot_rain_9b[i]} else {0}
    numda_9<-if (t<1) {0} else if (t>1 && t<10) {removal_rate_other_9a_hot_rain[i]} else {removal_rate_other_9b_hot_rain[i]}
    dC_aer_9<-G_9/Volume-Flow/Volume*C_aer_9-Q_air/Volume*C_aer_9-v_floor[9]*Aera/Volume*C_aer_9-v_vertical[9]*(Length*Height*2+Width*Height*2)/Volume*C_aer_9-numda_9*C_aer_9
    
    dC_leg<-release_constant*(steady_con_hot[i]-C_leg)
    dDose<-(C_aer_1*Fraction[1]*DE1[i]+C_aer_2*Fraction[2]*DE2[i]+C_aer_3*Fraction[3]*DE3[i]+C_aer_4*Fraction[4]*DE4[i]+C_aer_5*Fraction[5]*DE5[i]+C_aer_6*Fraction[6]*DE6[i]+C_aer_7*Fraction[7]*DE7[i]+C_aer_8*Fraction[8]*DE8[i]+C_aer_9*Fraction[9]*DE9[i])*C_leg*Inhalation[i]
    return(list(c(dC_aer_1,dC_aer_2,dC_aer_3,dC_aer_4,dC_aer_5,dC_aer_6,dC_aer_7,dC_aer_8,dC_aer_9,dC_leg,dDose)))
  }
  
  out_hot_rain <- ode.1D(y = state, times = time, func = differetial, parms = NULL,nspec=1)
  Dose_hot_rain<-out_hot_rain[,"Dose"]
  Dose_hot_rain_overall<-cbind(Dose_hot_rain_overall,Dose_hot_rain)
  
}
Dose_hot_rain_1<-Dose_hot_rain_overall[101,]
Dose_hot_rain_2<-Dose_hot_rain_overall[201,]-Dose_hot_rain_overall[101,]
Dose_hot_rain_3<-Dose_hot_rain_overall[301,]-Dose_hot_rain_overall[201,]
Dose_hot_rain_4<-Dose_hot_rain_overall[401,]-Dose_hot_rain_overall[301,]
Dose_hot_rain_5<-Dose_hot_rain_overall[501,]-Dose_hot_rain_overall[401,]
Dose_hot_rain_6<-Dose_hot_rain_overall[601,]-Dose_hot_rain_overall[501,]
Dose_hot_rain_7<-Dose_hot_rain_overall[701,]-Dose_hot_rain_overall[601,]
Dose_hot_rain_8<-Dose_hot_rain_overall[801,]-Dose_hot_rain_overall[701,]
Dose_hot_rain_9<-Dose_hot_rain_overall[901,]-Dose_hot_rain_overall[801,]
Dose_hot_rain_10<-Dose_hot_rain_overall[1001,]-Dose_hot_rain_overall[901,]
Dose_hot_rain_11<-Dose_hot_rain_overall[1101,]-Dose_hot_rain_overall[1001,]
Dose_hot_rain_12<-Dose_hot_rain_overall[1201,]-Dose_hot_rain_overall[1101,]
Dose_hot_rain_13<-Dose_hot_rain_overall[1301,]-Dose_hot_rain_overall[1201,]
Dose_hot_rain_14<-Dose_hot_rain_overall[1401,]-Dose_hot_rain_overall[1301,]
Dose_hot_rain_15<-Dose_hot_rain_overall[1501,]-Dose_hot_rain_overall[1401,]

Risk_hot_rain_t<-1-exp(-r_inf*t(Dose_hot_rain_overall))
Risk_hot_rain<-t(Risk_hot_rain_t)
Risk_hot_rain_child<-Morbidity_child*Risk_hot_rain[duration/time_step+1,]
Risk_hot_rain_adult<-Morbidity_adult*Risk_hot_rain[duration/time_step+1,]
Risk_hot_rain_elderly<-Morbidity_elderly*Risk_hot_rain[duration/time_step+1,]
Risk_hot_rain_annual<-1-(1-Risk_hot_rain)^365
Risk_hot_rain_child_annual<-1-(1-Risk_hot_rain_child)^365
Risk_hot_rain_adult_annual<-1-(1-Risk_hot_rain_adult)^365
Risk_hot_rain_elderly_annual<-1-(1-Risk_hot_rain_elderly)^365
Risk_hot_rain_DALY<-Risk_hot_rain_annual[duration/time_step+1,]*DALY
Dose_hot_rain_discrete<-rbind(Dose_hot_rain_1,Dose_hot_rain_2,Dose_hot_rain_3,Dose_hot_rain_4,Dose_hot_rain_5,Dose_hot_rain_6,Dose_hot_rain_7,Dose_hot_rain_8,Dose_hot_rain_9,Dose_hot_rain_10,Dose_hot_rain_11,Dose_hot_rain_12,Dose_hot_rain_13,Dose_hot_rain_14,Dose_hot_rain_15)
Risk_hot_rain_discrete_t<-1-exp(-r_inf*t(Dose_hot_rain_discrete))
Risk_hot_rain_discrete<-t(Risk_hot_rain_discrete_t)
Risk_hot_rain_discrete_annual<-1-(1-Risk_hot_rain_discrete)^365

median_risk_hot_rain<-apply(Risk_hot_rain_annual,1,median)
percentile_risk_low_hot_rain<-apply(Risk_hot_rain_annual,1,quantile,probs=c(.25))
percentile_risk_high_hot_rain<-apply(Risk_hot_rain_annual,1,quantile,probs=c(.75))
risk_overall_hot_rain<-data.frame(Time=seq(from=0,to=duration,by=time_step),risk_median=median_risk_hot_rain,risk_lower=percentile_risk_low_hot_rain,risk_higher=percentile_risk_high_hot_rain,Condition=rep("Hot water",1501),Type=rep("Rain showerhead",1501))
Final_risk_hot_rain_child<-data.frame(Risk=as.vector(Risk_hot_rain_child_annual),Group=rep("Child",1000),Condition=rep("Hot water_rain showerhead",1000))
Final_risk_hot_rain_adult<-data.frame(Risk=as.vector(Risk_hot_rain_adult_annual),Group=rep("Adult",1000),Condition=rep("Hot water_rain showerhead",1000))
Final_risk_hot_rain_elderly<-data.frame(Risk=as.vector(Risk_hot_rain_elderly_annual),Group=rep("Elderly",1000),Condition=rep("Hot water_rain showerhead",1000))
Final_risk_hot_rain_DALY<-data.frame(Risk=as.vector(Risk_hot_rain_DALY),Group=rep("DALY",1000),Condition=rep("Hot water_rain showerhead",1000))
Risk_hot_rain_discrete_1<-data.frame(Risk=Risk_hot_rain_discrete_annual[1,],Time=rep("A",1000),Condition=rep("Hot water",1000),Type=rep("Rain showerhead",1000))
Risk_hot_rain_discrete_2<-data.frame(Risk=Risk_hot_rain_discrete_annual[2,],Time=rep("B",1000),Condition=rep("Hot water",1000),Type=rep("Rain showerhead",1000))
Risk_hot_rain_discrete_3<-data.frame(Risk=Risk_hot_rain_discrete_annual[3,],Time=rep("C",1000),Condition=rep("Hot water",1000),Type=rep("Rain showerhead",1000))
Risk_hot_rain_discrete_4<-data.frame(Risk=Risk_hot_rain_discrete_annual[4,],Time=rep("D",1000),Condition=rep("Hot water",1000),Type=rep("Rain showerhead",1000))
Risk_hot_rain_discrete_5<-data.frame(Risk=Risk_hot_rain_discrete_annual[5,],Time=rep("E",1000),Condition=rep("Hot water",1000),Type=rep("Rain showerhead",1000))
Risk_hot_rain_discrete_6<-data.frame(Risk=Risk_hot_rain_discrete_annual[6,],Time=rep("F",1000),Condition=rep("Hot water",1000),Type=rep("Rain showerhead",1000))
Risk_hot_rain_discrete_7<-data.frame(Risk=Risk_hot_rain_discrete_annual[7,],Time=rep("G",1000),Condition=rep("Hot water",1000),Type=rep("Rain showerhead",1000))
Risk_hot_rain_discrete_8<-data.frame(Risk=Risk_hot_rain_discrete_annual[8,],Time=rep("H",1000),Condition=rep("Hot water",1000),Type=rep("Rain showerhead",1000))
Risk_hot_rain_discrete_9<-data.frame(Risk=Risk_hot_rain_discrete_annual[9,],Time=rep("I",1000),Condition=rep("Hot water",1000),Type=rep("Rain showerhead",1000))
Risk_hot_rain_discrete_10<-data.frame(Risk=Risk_hot_rain_discrete_annual[10,],Time=rep("J",1000),Condition=rep("Hot water",1000),Type=rep("Rain showerhead",1000))
Risk_hot_rain_discrete_11<-data.frame(Risk=Risk_hot_rain_discrete_annual[11,],Time=rep("K",1000),Condition=rep("Hot water",1000),Type=rep("Rain showerhead",1000))
Risk_hot_rain_discrete_12<-data.frame(Risk=Risk_hot_rain_discrete_annual[12,],Time=rep("L",1000),Condition=rep("Hot water",1000),Type=rep("Rain showerhead",1000))
Risk_hot_rain_discrete_13<-data.frame(Risk=Risk_hot_rain_discrete_annual[13,],Time=rep("M",1000),Condition=rep("Hot water",1000),Type=rep("Rain showerhead",1000))
Risk_hot_rain_discrete_14<-data.frame(Risk=Risk_hot_rain_discrete_annual[14,],Time=rep("N",1000),Condition=rep("Hot water",1000),Type=rep("Rain showerhead",1000))
Risk_hot_rain_discrete_15<-data.frame(Risk=Risk_hot_rain_discrete_annual[15,],Time=rep("O",1000),Condition=rep("Hot water",1000),Type=rep("Rain showerhead",1000))
Risk_hot_rain_discrete_overall<-rbind(Risk_hot_rain_discrete_1,Risk_hot_rain_discrete_2,Risk_hot_rain_discrete_3,Risk_hot_rain_discrete_4,Risk_hot_rain_discrete_5,Risk_hot_rain_discrete_6,Risk_hot_rain_discrete_7,Risk_hot_rain_discrete_8,Risk_hot_rain_discrete_9,Risk_hot_rain_discrete_10,Risk_hot_rain_discrete_11,Risk_hot_rain_discrete_12,Risk_hot_rain_discrete_13,Risk_hot_rain_discrete_14,Risk_hot_rain_discrete_15)



Risk_over_time<-rbind(risk_overall_cold_conv,risk_overall_hot_conv,risk_overall_cold_rain,risk_overall_hot_rain)

Final_risk_overall<-rbind(Final_risk_cold_rain_child,Final_risk_cold_rain_adult,Final_risk_cold_rain_elderly,Final_risk_hot_rain_child,Final_risk_hot_rain_adult,Final_risk_hot_rain_elderly,Final_risk_cold_conv_child,Final_risk_cold_conv_adult,Final_risk_cold_conv_elderly,Final_risk_hot_conv_child,Final_risk_hot_conv_adult,Final_risk_hot_conv_elderly)
Final_risk_overall$Group<-factor(Final_risk_overall$Group,levels = c("Child","Adult","Elderly"))
Final_risk_DALY<-rbind(Final_risk_cold_conv_DALY,Final_risk_cold_rain_DALY,Final_risk_hot_conv_DALY,Final_risk_hot_rain_DALY)
Risk_discrete<-rbind(Risk_cold_conv_discrete_overall,Risk_cold_rain_discrete_overall,Risk_hot_conv_discrete_overall,Risk_hot_rain_discrete_overall)
Risk_threshold<-data.frame(Time=seq(from=0,to=duration,by=time_step),Risk=rep(1e-4,duration/time_step+1))
DALY_threshould<-data.frame(Time=seq(from=0,to=duration,by=time_step),Risk=rep(1e-6,duration/time_step+1))
ggplot(Risk_over_time,aes(x=Time,y=risk_median))+geom_line()+scale_y_continuous(trans='log10',breaks=c(1e-10,1e-9,1e-8,1e-7,1e-6,1e-5,1e-4,1e-3,1e-2,1e-1,1),label=scientific_10)+scale_x_continuous(breaks = seq(0,duration,b=1))+xlab("Time (min)")+ylab("Cumulative risk of infection")+geom_ribbon(aes(x=Time,ymax=risk_higher,ymin=risk_lower,fill=Condition),linetype="dashed",alpha=0.5,col="black",show.legend = FALSE)+scale_fill_brewer(palette="Pastel1",direction=-1)+facet_grid(Condition~Type)+theme(panel.grid.major = element_blank(), 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                          panel.grid.minor = element_blank(),panel.background = element_blank(),
                                                                                                                                                                                                                                                                                                                                                                                                                                                                          panel.border = element_rect(colour = "black", fill=NA, size=2),text=element_text(size = 18))+geom_line(Risk_threshold,mapping=aes(Time,Risk),colour="red",linetype="dashed")

ggplot(Final_risk_overall,mapping=aes(Group,Risk,fill=Condition))+stat_boxplot(geom = "errorbar", width = 0.4,position=position_dodge(0.8))+geom_boxplot(outlier.shape = NA,width=0.6,position =position_dodge(0.8))+ylab("Risk of illness")+scale_y_continuous(trans='log10',breaks = c(1e-10,1e-9,1e-8,1e-7,1e-6,1e-5,1e-4,1e-3,1e-2,1e-1,1),label=scientific_10)+theme(panel.grid.major = element_blank(), 
                                                                                                                                                                                                                                                                                                                                               panel.grid.minor = element_blank(),panel.background = element_blank(),
                                                                                                                                                                                                                                                                                                                                               panel.border = element_rect(colour = "black", fill=NA, size=2),text = element_text(size = 18),axis.title.x = element_blank())+scale_fill_brewer(palette="RdBu",direction=-1)
ggplot(Final_risk_DALY,aes(Condition,Risk))+stat_boxplot(geom = "errorbar", width = 0.2,position=position_dodge(0.6))+geom_boxplot(outlier.shape = NA,width=0.4,position =position_dodge(0.6),fill="#FFCC99")+ylab("DALY")+scale_y_continuous(trans='log10',breaks = c(1e-10,1e-9,1e-8,1e-7,1e-6,1e-5,1e-4,1e-3,1e-2,1e-1,1),label=scientific_10)+theme(panel.grid.major = element_blank(), 
                                                                                                                                                                                                                                                                                                                                                                  panel.grid.minor = element_blank(),panel.background = element_blank(),
                                                                                                                                                                                                                                                                                                                                                                  panel.border = element_rect(colour = "black", fill=NA, size=2),text = element_text(size=18),axis.title.x = element_blank())
ggplot(Risk_discrete,aes(Time,Risk,fill=Condition))+stat_boxplot(geom = "errorbar", width = 0.4,position=position_dodge(2))+geom_boxplot(outlier.shape = NA,width=0.6,position =position_dodge(2))+xlab("Time (min)")+ylab("Risk of infection")+scale_x_discrete(labels=c("0-1","1-2","2-3","3-4","4-5","5-6","6-7","7-8","8-9","9-10","10-11","11-12","12-13","13-14","14-15"))+scale_y_continuous(trans='log10',breaks = c(1e-12,1e-10,1e-8,1e-6,1e-4,1e-2,1),label=scientific_10)+facet_grid(Condition~Type)+theme(panel.grid.major = element_blank(), 
                                                                                                                                                                                                                                                                                                                                               panel.grid.minor = element_blank(),panel.background = element_blank(),
                                                                                                                                                                                                                                                                                                                                               panel.border = element_rect(colour = "black", fill=NA, size=2),text = element_text(size=18))+scale_fill_brewer(palette="RdBu",direction=-1)
# sensitivity analysis
Spearman_cold_conv<-data.frame(DE1=DE1,DE2=DE2,DE3=DE3,DE4=DE4,DE5=DE5,Inhalation_rate=Inhalation,Dose_response=r_inf,Generation_rate_1=Generation_cold_conv_1,Generation_rate_2=Generation_cold_conv_2,Generation_rate_3=Generation_cold_conv_3,Generation_rate_4=Generation_cold_conv_4,Generation_rate_5=Generation_cold_conv_5,Removal_rate_1a=removal_rate_other_1a_cold_conv,Removal_rate_2a=removal_rate_other_2a_cold_conv,Removal_rate_3a=removal_rate_other_3a_cold_conv,Removal_rate_4a=removal_rate_other_4a_cold_conv,Removal_rate_5a=removal_rate_other_5a_cold_conv,Removal_rate_1b=removal_rate_other_1b_cold_conv,Removal_rate_2b=removal_rate_other_2b_cold_conv,Removal_rate_3b=removal_rate_other_3b_cold_conv,Removal_rate_4b=removal_rate_other_4b_cold_conv,Removal_rate_5b=removal_rate_other_5b_cold_conv,Concentration_init=Init_con_cold,Ventilation_1=Ventilation_cold_conv_1,Ventilation_2=Ventilation_cold_conv_2,Risk_cold_conv=Risk_cold_conv_annual[duration/time_step+1,])
Spearman_cold_rain<-data.frame(DE1=DE1,DE2=DE2,DE3=DE3,DE4=DE4,DE5=DE5,Inhalation_rate=Inhalation,Dose_response=r_inf,Generation_rate_1=Generation_cold_rain_1,Generation_rate_2=Generation_cold_rain_2,Generation_rate_3=Generation_cold_rain_3,Generation_rate_4=Generation_cold_rain_4,Generation_rate_5=Generation_cold_rain_5,Removal_rate_1a=removal_rate_other_1a_cold_rain,Removal_rate_2a=removal_rate_other_2a_cold_rain,Removal_rate_3a=removal_rate_other_3a_cold_rain,Removal_rate_4a=removal_rate_other_4a_cold_rain,Removal_rate_5a=removal_rate_other_5a_cold_rain,Removal_rate_1b=removal_rate_other_1b_cold_rain,Removal_rate_2b=removal_rate_other_2b_cold_rain,Removal_rate_3b=removal_rate_other_3b_cold_rain,Removal_rate_4b=removal_rate_other_4b_cold_rain,Removal_rate_5b=removal_rate_other_5b_cold_rain,Concentration_init=Init_con_cold,Ventilation_1=Ventilation_cold_rain_1,Ventilation_2=Ventilation_cold_rain_2,Risk_cold_rain=Risk_cold_rain_annual[duration/time_step+1,])
Spearman_hot_conv<-data.frame(DE1=DE1,DE2=DE2,DE3=DE3,DE4=DE4,DE5=DE5,DE6=DE6,DE7=DE7,DE8=DE8,DE9=DE9,Inhalation_rate=Inhalation,Dose_response=r_inf,Generation_rate_1=Generation_hot_conv_1b,Generation_rate_2=Generation_hot_conv_2b,Generation_rate_3=Generation_hot_conv_3b,Generation_rate_4=Generation_hot_conv_4b,Generation_rate_5=Generation_hot_conv_5b,Generation_rate_6=Generation_hot_conv_6b,Generation_rate_7=Generation_hot_conv_7b,Generation_rate_8=Generation_hot_conv_8b,Generation_rate_9=Generation_hot_conv_9b,Removal_rate_1a=removal_rate_other_1a_hot_conv,Removal_rate_2a=removal_rate_other_2a_hot_conv,Removal_rate_3a=removal_rate_other_3a_hot_conv,Removal_rate_4a=removal_rate_other_4a_hot_conv,Removal_rate_5a=removal_rate_other_5a_hot_conv,Removal_rate_6a=removal_rate_other_6a_hot_conv,Removal_rate_7a=removal_rate_other_7a_hot_conv,Removal_rate_8a=removal_rate_other_8a_hot_conv,Removal_rate_9a=removal_rate_other_9a_hot_conv,Removal_rate_1b=removal_rate_other_1b_hot_conv,Removal_rate_2b=removal_rate_other_2b_hot_conv,Removal_rate_3b=removal_rate_other_3b_hot_conv,Removal_rate_4b=removal_rate_other_4b_hot_conv,Removal_rate_5b=removal_rate_other_5b_hot_conv,Removal_rate_6b=removal_rate_other_6b_hot_conv,Removal_rate_7b=removal_rate_other_7b_hot_conv,Removal_rate_8b=removal_rate_other_8b_hot_conv,Removal_rate_9b=removal_rate_other_9b_hot_conv,Concentration_init=Init_con_hot,Ventilation_1=Ventilation_hot_conv_1,Ventilation_2=Ventilation_hot_conv_2,Risk_hot_conv=Risk_hot_conv_annual[duration/time_step+1,])
Spearman_hot_rain<-data.frame(DE1=DE1,DE2=DE2,DE3=DE3,DE4=DE4,DE5=DE5,DE6=DE6,DE7=DE7,DE8=DE8,DE9=DE9,Inhalation_rate=Inhalation,Dose_response=r_inf,Generation_rate_1=Generation_hot_rain_1b,Generation_rate_2=Generation_hot_rain_2b,Generation_rate_3=Generation_hot_rain_3b,Generation_rate_4=Generation_hot_rain_4b,Generation_rate_5=Generation_hot_rain_5b,Generation_rate_6=Generation_hot_rain_6b,Generation_rate_7=Generation_hot_rain_7b,Generation_rate_8=Generation_hot_rain_8b,Generation_rate_9=Generation_hot_rain_9b,Removal_rate_1a=removal_rate_other_1a_hot_rain,Removal_rate_2a=removal_rate_other_2a_hot_rain,Removal_rate_3a=removal_rate_other_3a_hot_rain,Removal_rate_4a=removal_rate_other_4a_hot_rain,Removal_rate_5a=removal_rate_other_5a_hot_rain,Removal_rate_6a=removal_rate_other_6a_hot_rain,Removal_rate_7a=removal_rate_other_7a_hot_rain,Removal_rate_8a=removal_rate_other_8a_hot_rain,Removal_rate_9a=removal_rate_other_9a_hot_rain,Removal_rate_1b=removal_rate_other_1b_hot_rain,Removal_rate_2b=removal_rate_other_2b_hot_rain,Removal_rate_3b=removal_rate_other_3b_hot_rain,Removal_rate_4b=removal_rate_other_4b_hot_rain,Removal_rate_5b=removal_rate_other_5b_hot_rain,Removal_rate_6b=removal_rate_other_6b_hot_rain,Removal_rate_7b=removal_rate_other_7b_hot_rain,Removal_rate_8b=removal_rate_other_8b_hot_rain,Removal_rate_9b=removal_rate_other_9b_hot_rain,Concentration_init=Init_con_hot,Ventilation_1=Ventilation_hot_rain_1,Ventilation_2=Ventilation_hot_rain_2,Risk_hot_conv=Risk_hot_rain_annual[duration/time_step+1,])

Data_spearman_cold_conv<-numeric()
for (i in 1:25){
  corr<-abs(cor(x=Spearman_cold_conv[,i],y=Spearman_cold_conv[,26],method="spearman"))
  Data_spearman_cold_conv<-rbind(Data_spearman_cold_conv,corr)
}
colnames(Data_spearman_cold_conv)<-c("Spearman")
Data_spearman_cold_rain<-numeric()
for (i in 1:25){
  corr<-abs(cor(x=Spearman_cold_rain[,i],y=Spearman_cold_rain[,26],method="spearman"))
  Data_spearman_cold_rain<-rbind(Data_spearman_cold_rain,corr)
}
colnames(Data_spearman_cold_rain)<-c("Spearman")

Data_spearman_hot_conv<-numeric()
for (i in 1:41){
  corr<-abs(cor(x=Spearman_hot_conv[,i],y=Spearman_hot_conv[,42],method="spearman"))
  Data_spearman_hot_conv<-rbind(Data_spearman_hot_conv,corr)
}
colnames(Data_spearman_hot_conv)<-c("Spearman")
Data_spearman_hot_rain<-numeric()
for (i in 1:41){
  corr<-abs(cor(x=Spearman_hot_rain[,i],y=Spearman_hot_rain[,42],method="spearman"))
  Data_spearman_hot_rain<-rbind(Data_spearman_hot_rain,corr)
}
colnames(Data_spearman_hot_rain)<-c("Spearman")

Data_spearman_cold_conv<-as.data.frame(Data_spearman_cold_conv)
Data_spearman_cold_rain<-as.data.frame(Data_spearman_cold_rain)
Data_spearman_hot_conv<-as.data.frame(Data_spearman_hot_conv)
Data_spearman_hot_rain<-as.data.frame(Data_spearman_hot_rain)

Data_spearman_cold_conv$Parameter<-c("DE1","DE2","DE3","DE4","DE5","Inhalation rate","Dose-response parameter","Aerosol generation rate (1-2 μm)","Aerosol generation rate (2-3 μm)","Aerosol generation rate (3-4 μm)","Aerosol generation rate (4-5 μm)","Aerosol generation rate (5-6 μm)","Aerosol decay rate (1-2 μm,shower on)","Aerosol decay rate (2-3 μm,shower on)","Aerosol decay rate (3-4 μm,shower on)","Aerosol decay rate (4-5 μm,shower on)","Aerosol decay rate (5-6 μm,shower on)","Aerosol decay rate (1-2 μm,shower off)","Aerosol decay rate (2-3 μm,shower off)","Aerosol decay rate (3-4 μm,shower off)","Aerosol decay rate (4-5 μm,shower off)","Aerosol decay rate (5-6 μm,shower off)","Concentration of pathogen in water","Ventilation rate (shower on)","Ventilation rate (shower off)")
Data_spearman_cold_conv$Condition<-rep("Cold water_conventional showerhead",25)
Data_spearman_cold_rain$Parameter<-c("DE1","DE2","DE3","DE4","DE5","Inhalation rate","Dose-response parameter","Aerosol generation rate (1-2 μm)","Aerosol generation rate (2-3 μm)","Aerosol generation rate (3-4 μm)","Aerosol generation rate (4-5 μm)","Aerosol generation rate (5-6 μm)","Aerosol decay rate (1-2 μm,shower on)","Aerosol decay rate (2-3 μm,shower on)","Aerosol decay rate (3-4 μm,shower on)","Aerosol decay rate (4-5 μm,shower on)","Aerosol decay rate (5-6 μm,shower on)","Aerosol decay rate (1-2 μm,shower off)","Aerosol decay rate (2-3 μm,shower off)","Aerosol decay rate (3-4 μm,shower off)","Aerosol decay rate (4-5 μm,shower off)","Aerosol decay rate (5-6 μm,shower off)","Concentration of pathogen in water","Ventilation rate (shower on)","Ventilation rate (shower off)")
Data_spearman_cold_rain$Condition<-rep("Cold water_rain showerhead",25)
Data_spearman_hot_conv$Parameter<-c("DE1","DE2","DE3","DE4","DE5","DE6","DE7","DE8","DE9","Inhalation rate","Dose-response parameter","Aerosol generation rate (1-2 μm)","Aerosol generation rate (2-3 μm)","Aerosol generation rate (3-4 μm)","Aerosol generation rate (4-5 μm)","Aerosol generation rate (5-6 μm)","Aerosol generation rate (6-7 μm)","Aerosol generation rate (7-8 μm)","Aerosol generation rate (8-9 μm)","Aerosol generation rate (9-10 μm)","Aerosol decay rate (1-2 μm,shower on)","Aerosol decay rate (2-3 μm,shower on)","Aerosol decay rate (3-4 μm,shower on)","Aerosol decay rate (4-5 μm,shower on)","Aerosol decay rate (5-6 μm,shower on)","Aerosol decay rate (6-7 μm,shower on)","Aerosol decay rate (7-8 μm,shower on)","Aerosol decay rate (8-9 μm,shower on)","Aerosol decay rate (9-10 μm,shower on)","Aerosol decay rate (1-2 μm,shower off)","Aerosol decay rate (2-3 μm,shower off)","Aerosol decay rate (3-4 μm,shower off)","Aerosol decay rate (4-5 μm,shower off)","Aerosol decay rate (5-6 μm,shower off)","Aerosol decay rate (6-7 μm,shower off)","Aerosol decay rate (7-8 μm,shower off)","Aerosol decay rate (8-9 μm,shower off)","Aerosol decay rate (9-10 μm,shower off)","Concentration of pathogen in water","Ventilation rate (shower on)","Ventilation rate (shower off)")
Data_spearman_hot_conv$Condition<-rep("Hot water_conventional showerhead",41)
Data_spearman_hot_rain$Parameter<-c("DE1","DE2","DE3","DE4","DE5","DE6","DE7","DE8","DE9","Inhalation rate","Dose-response parameter","Aerosol generation rate (1-2 μm)","Aerosol generation rate (2-3 μm)","Aerosol generation rate (3-4 μm)","Aerosol generation rate (4-5 μm)","Aerosol generation rate (5-6 μm)","Aerosol generation rate (6-7 μm)","Aerosol generation rate (7-8 μm)","Aerosol generation rate (8-9 μm)","Aerosol generation rate (9-10 μm)","Aerosol decay rate (1-2 μm,shower on)","Aerosol decay rate (2-3 μm,shower on)","Aerosol decay rate (3-4 μm,shower on)","Aerosol decay rate (4-5 μm,shower on)","Aerosol decay rate (5-6 μm,shower on)","Aerosol decay rate (6-7 μm,shower on)","Aerosol decay rate (7-8 μm,shower on)","Aerosol decay rate (8-9 μm,shower on)","Aerosol decay rate (9-10 μm,shower on)","Aerosol decay rate (1-2 μm,shower off)","Aerosol decay rate (2-3 μm,shower off)","Aerosol decay rate (3-4 μm,shower off)","Aerosol decay rate (4-5 μm,shower off)","Aerosol decay rate (5-6 μm,shower off)","Aerosol decay rate (6-7 μm,shower off)","Aerosol decay rate (7-8 μm,shower off)","Aerosol decay rate (8-9 μm,shower off)","Aerosol decay rate (9-10 μm,shower off)","Concentration of pathogen in water","Ventilation rate (shower on)","Ventilation rate (shower off)")
Data_spearman_hot_rain$Condition<-rep("Hot water_rain showerhead",41)
Data_spearman_cold<-rbind(Data_spearman_cold_conv,Data_spearman_cold_rain)
Data_spearman_hot<-rbind(Data_spearman_hot_conv,Data_spearman_hot_rain)


ggplot(Data_spearman_cold,aes(x=Parameter,y=Spearman,fill=Condition))+geom_bar(stat="identity",position=position_dodge())+ylab("Spearman rank correlation coefficient")+coord_flip()+theme(panel.grid.major = element_blank(), 
                                                                                                                                                                                           panel.grid.minor = element_blank(),panel.background = element_blank(),
                                                                                                                                                                                           panel.border = element_rect(colour = "black", fill=NA, size=2),text=element_text(size=15))
ggplot(Data_spearman_hot,aes(x=Parameter,y=Spearman,fill=Condition))+geom_bar(stat="identity",position=position_dodge())+ylab("Spearman rank correlation coefficient")+coord_flip()+theme(panel.grid.major = element_blank(), 
                                                                                                                                                                                          panel.grid.minor = element_blank(),panel.background = element_blank(),
                                                                                                                                                                                          panel.border = element_rect(colour = "black", fill=NA, size=2),text=element_text(size=15))

