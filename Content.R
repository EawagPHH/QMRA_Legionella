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

set.seed(100)
# Parameter bacteria
release_constant<-5

#demographic parameter
Morbidity_child<-0.0077
Morbidity_adult<-0.097
Morbidity_elderly<-0.75
Disease_severity_Pontiac<-0.1
Disease_severity_LD<-0.3
Disease_duration_Pontiac<-3.5
Disease_duration_LD<-21
Pontiac_percentage<-0.95
LD_percentage<-0.05
Mortality_LD<-0.005
Mean_age_paitent<-62.2
Life_expenctancy<-82.5
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

Ventilation_cold_conv_1<-runif(1000,min=0.24425,max=0.54668)
Ventilation_cold_conv_2<-runif(1000,min=0.1,max=0.54668)

Ventilation_hot_conv_1<-runif(1000,min=0.25,max=0.63)
Ventilation_hot_conv_2<-runif(1000,min=0.6167,max=1.16698)

Generation_rate_cold_conv<-c(mean(Generation_cold_conv_1),mean(Generation_cold_conv_2),mean(Generation_cold_conv_3),mean(Generation_cold_conv_4),mean(Generation_cold_conv_5))
Generation_rate_cold_rain<-c(mean(Generation_cold_rain_1),mean(Generation_cold_rain_2),mean(Generation_cold_rain_3),mean(Generation_cold_rain_4),mean(Generation_cold_rain_5))
Generation_rate_hot_conv_a<-c(mean(Generation_hot_conv_1a),mean(Generation_hot_conv_2a),mean(Generation_hot_conv_3a),mean(Generation_hot_conv_4a),mean(Generation_hot_conv_5a),mean(Generation_hot_conv_6a),mean(Generation_hot_conv_7a),mean(Generation_hot_conv_8a),mean(Generation_hot_conv_9a))
Generation_rate_hot_conv_b<-c(mean(Generation_hot_conv_1b),mean(Generation_hot_conv_2b),mean(Generation_hot_conv_3b),mean(Generation_hot_conv_4b),mean(Generation_hot_conv_5b),mean(Generation_hot_conv_6b),mean(Generation_hot_conv_7b),mean(Generation_hot_conv_8b),mean(Generation_hot_conv_9b))
Generation_rate_hot_rain_a<-c(mean(Generation_hot_rain_1a),mean(Generation_hot_rain_2a),mean(Generation_hot_rain_3a),mean(Generation_hot_rain_4a),mean(Generation_hot_rain_5a),mean(Generation_hot_rain_6a),mean(Generation_hot_rain_7a),mean(Generation_hot_rain_8a),mean(Generation_hot_rain_9a))
Generation_rate_hot_rain_b<-c(mean(Generation_hot_rain_1b),mean(Generation_hot_rain_2b),mean(Generation_hot_rain_3b),mean(Generation_hot_rain_4b),mean(Generation_hot_rain_5b),mean(Generation_hot_rain_6b),mean(Generation_hot_rain_7b),mean(Generation_hot_rain_8b),mean(Generation_hot_rain_9b))
Decay_cold_conv_a<-c(mean(removal_rate_other_1a_cold_conv),mean(removal_rate_other_2a_cold_conv),mean(removal_rate_other_3a_cold_conv),mean(removal_rate_other_4a_cold_conv),mean(removal_rate_other_5a_cold_conv))
Decay_cold_conv_b<-c(mean(removal_rate_other_1b_cold_conv),mean(removal_rate_other_2b_cold_conv),mean(removal_rate_other_3b_cold_conv),mean(removal_rate_other_4b_cold_conv),mean(removal_rate_other_5b_cold_conv))
Decay_cold_rain_a<-c(mean(removal_rate_other_1a_cold_rain),mean(removal_rate_other_2a_cold_rain),mean(removal_rate_other_3a_cold_rain),mean(removal_rate_other_4a_cold_rain),mean(removal_rate_other_5a_cold_rain))
Decay_cold_rain_b<-c(mean(removal_rate_other_1b_cold_rain),mean(removal_rate_other_2b_cold_rain),mean(removal_rate_other_3b_cold_rain),mean(removal_rate_other_4b_cold_rain),mean(removal_rate_other_5b_cold_rain))
Decay_hot_conv_a<-c(mean(removal_rate_other_1a_hot_conv),mean(removal_rate_other_2a_hot_conv),mean(removal_rate_other_3a_hot_conv),mean(removal_rate_other_4a_hot_conv),mean(removal_rate_other_5a_hot_conv),mean(removal_rate_other_6a_hot_conv),mean(removal_rate_other_7a_hot_conv),mean(removal_rate_other_8a_hot_conv),mean(removal_rate_other_9a_hot_conv))
Decay_hot_conv_b<-c(mean(removal_rate_other_1b_hot_conv),mean(removal_rate_other_2b_hot_conv),mean(removal_rate_other_3b_hot_conv),mean(removal_rate_other_4b_hot_conv),mean(removal_rate_other_5b_hot_conv),mean(removal_rate_other_6b_hot_conv),mean(removal_rate_other_7b_hot_conv),mean(removal_rate_other_8b_hot_conv),mean(removal_rate_other_9b_hot_conv))
Decay_hot_rain_a<-c(mean(removal_rate_other_1a_hot_rain),mean(removal_rate_other_2a_hot_rain),mean(removal_rate_other_3a_hot_rain),mean(removal_rate_other_4a_hot_rain),mean(removal_rate_other_5a_hot_rain),mean(removal_rate_other_6a_hot_rain),mean(removal_rate_other_7a_hot_rain),mean(removal_rate_other_8a_hot_rain),mean(removal_rate_other_9a_hot_rain))
Decay_hot_rain_b<-c(mean(removal_rate_other_1b_hot_rain),mean(removal_rate_other_2b_hot_rain),mean(removal_rate_other_3b_hot_rain),mean(removal_rate_other_4b_hot_rain),mean(removal_rate_other_5b_hot_rain),mean(removal_rate_other_6b_hot_rain),mean(removal_rate_other_7b_hot_rain),mean(removal_rate_other_8b_hot_rain),mean(removal_rate_other_9b_hot_rain))

Generation_rate_cold_conv_sd<-c(sd(Generation_cold_conv_1),sd(Generation_cold_conv_2),sd(Generation_cold_conv_3),sd(Generation_cold_conv_4),sd(Generation_cold_conv_5))
Generation_rate_cold_rain_sd<-c(sd(Generation_cold_rain_1),sd(Generation_cold_rain_2),sd(Generation_cold_rain_3),sd(Generation_cold_rain_4),sd(Generation_cold_rain_5))
Generation_rate_hot_conv_a_sd<-c(sd(Generation_hot_conv_1a),sd(Generation_hot_conv_2a),sd(Generation_hot_conv_3a),sd(Generation_hot_conv_4a),sd(Generation_hot_conv_5a),sd(Generation_hot_conv_6a),sd(Generation_hot_conv_7a),sd(Generation_hot_conv_8a),sd(Generation_hot_conv_9a))
Generation_rate_hot_conv_b_sd<-c(sd(Generation_hot_conv_1b),sd(Generation_hot_conv_2b),sd(Generation_hot_conv_3b),sd(Generation_hot_conv_4b),sd(Generation_hot_conv_5b),sd(Generation_hot_conv_6b),sd(Generation_hot_conv_7b),sd(Generation_hot_conv_8b),sd(Generation_hot_conv_9b))
Generation_rate_hot_rain_a_sd<-c(sd(Generation_hot_rain_1a),sd(Generation_hot_rain_2a),sd(Generation_hot_rain_3a),sd(Generation_hot_rain_4a),sd(Generation_hot_rain_5a),sd(Generation_hot_rain_6a),sd(Generation_hot_rain_7a),sd(Generation_hot_rain_8a),sd(Generation_hot_rain_9a))
Generation_rate_hot_rain_b_sd<-c(mean(Generation_hot_rain_1b),mean(Generation_hot_rain_2b),mean(Generation_hot_rain_3b),mean(Generation_hot_rain_4b),mean(Generation_hot_rain_5b),mean(Generation_hot_rain_6b),mean(Generation_hot_rain_7b),mean(Generation_hot_rain_8b),mean(Generation_hot_rain_9b))
Decay_cold_conv_a_sd<-c(sd(removal_rate_other_1a_cold_conv),sd(removal_rate_other_2a_cold_conv),sd(removal_rate_other_3a_cold_conv),sd(removal_rate_other_4a_cold_conv),sd(removal_rate_other_5a_cold_conv))
Decay_cold_conv_b_sd<-c(sd(removal_rate_other_1b_cold_conv),sd(removal_rate_other_2b_cold_conv),sd(removal_rate_other_3b_cold_conv),sd(removal_rate_other_4b_cold_conv),sd(removal_rate_other_5b_cold_conv))
Decay_cold_rain_a_sd<-c(sd(removal_rate_other_1a_cold_rain),sd(removal_rate_other_2a_cold_rain),sd(removal_rate_other_3a_cold_rain),sd(removal_rate_other_4a_cold_rain),sd(removal_rate_other_5a_cold_rain))
Decay_cold_rain_b_sd<-c(sd(removal_rate_other_1b_cold_rain),sd(removal_rate_other_2b_cold_rain),sd(removal_rate_other_3b_cold_rain),sd(removal_rate_other_4b_cold_rain),sd(removal_rate_other_5b_cold_rain))
Decay_hot_conv_a_sd<-c(sd(removal_rate_other_1a_hot_conv),sd(removal_rate_other_2a_hot_conv),sd(removal_rate_other_3a_hot_conv),sd(removal_rate_other_4a_hot_conv),sd(removal_rate_other_5a_hot_conv),sd(removal_rate_other_6a_hot_conv),sd(removal_rate_other_7a_hot_conv),sd(removal_rate_other_8a_hot_conv),sd(removal_rate_other_9a_hot_conv))
Decay_hot_conv_b_sd<-c(sd(removal_rate_other_1b_hot_conv),sd(removal_rate_other_2b_hot_conv),sd(removal_rate_other_3b_hot_conv),sd(removal_rate_other_4b_hot_conv),sd(removal_rate_other_5b_hot_conv),sd(removal_rate_other_6b_hot_conv),sd(removal_rate_other_7b_hot_conv),sd(removal_rate_other_8b_hot_conv),sd(removal_rate_other_9b_hot_conv))
Decay_hot_rain_a_sd<-c(sd(removal_rate_other_1a_hot_rain),sd(removal_rate_other_2a_hot_rain),sd(removal_rate_other_3a_hot_rain),sd(removal_rate_other_4a_hot_rain),sd(removal_rate_other_5a_hot_rain),sd(removal_rate_other_6a_hot_rain),sd(removal_rate_other_7a_hot_rain),sd(removal_rate_other_8a_hot_rain),sd(removal_rate_other_9a_hot_rain))
Decay_hot_rain_b_sd<-c(sd(removal_rate_other_1b_hot_rain),sd(removal_rate_other_2b_hot_rain),sd(removal_rate_other_3b_hot_rain),sd(removal_rate_other_4b_hot_rain),sd(removal_rate_other_5b_hot_rain),sd(removal_rate_other_6b_hot_rain),sd(removal_rate_other_7b_hot_rain),sd(removal_rate_other_8b_hot_rain),sd(removal_rate_other_9b_hot_rain))
#dose estimation for cold water scenario
aerosol_cold<-function(shower_time,L,W,H,aerosol_size_class, Generation_rate,Decay_rate_during_shower,Decay_rate_after_shower,Ventilation,Legionella_initial,Legionella_final,Fraction){
  time_step<-0.01
  time<-seq(from=0, to=shower_time, by = time_step)
  state <- c(C_aer=0,C_leg=Legionella_initial,Dose=0)
  differetial <- function(t, state, parms,G,Q_air,numda){
    C_aer<-state[1]
    C_leg<-state[2]
    Dose<-state[3]
    Q_air<-Ventilation
    G<-if (t<shower_time-5) {Generation_rate} else {0}
    numda<-if (t<shower_time-5) {Decay_rate_during_shower} else {Decay_rate_during_shower}
    dC_aer<-G/(L*W*H)-Q_air/(L*W*H)*C_aer-v_floor[aerosol_size_class]*1/H*C_aer-v_vertical[aerosol_size_class]*(L*H*2+W*H*2)/(L*W*H)*C_aer-numda*C_aer
    dC_leg<-release_constant*(Legionella_final-C_leg)
    dDose<-C_aer*C_leg*Fraction
    return(list(c(dC_aer,dC_leg,dDose)))
  }
  output<- ode.1D(y = state, times = time, func = differetial, parms = NULL,nspec=1)
  return(output[,"Dose"])
}

#dose estimation for hot water scenario
aerosol_hot<-function(shower_time,L,W,H,aerosol_size_class, Generation_rate_initial,Generation_rate_steady_state,Decay_rate_during_shower,Decay_rate_after_shower,Ventilation,Legionella_initial,Legionella_final,Fraction){
  time_step<-0.01
  time<-seq(from=0, to=shower_time, by = time_step)
  state <- c(C_aer=0,C_leg=Legionella_initial,Dose=0)
  differetial <- function(t, state, parms,G,Q_air,numda){
    C_aer<-state[1]
    C_leg<-state[2]
    Dose<-state[3]
    Q_air<-Ventilation
    G<-if (t<1) {Generation_rate_initial} else if (t>1 && t<shower_time-5) {Generation_rate_steady_state} else {0}
    numda<-if (t<1) {0} else if (t>1 && t<shower_time-5) {Decay_rate_during_shower} else {Decay_rate_after_shower}
    dC_aer<-G/(L*W*H)-Q_air/(L*W*H)*C_aer-v_floor[aerosol_size_class]*1/H*C_aer-v_vertical[aerosol_size_class]*(L*H*2+W*H*2)/(L*W*H)*C_aer-numda*C_aer
    dC_leg<-release_constant*(Legionella_final-C_leg)
    dDose<-C_aer*C_leg*Fraction
    return(list(c(dC_aer,dC_leg,dDose)))
  }
  output<- ode.1D(y = state, times = time, func = differetial, parms = NULL,nspec=1)
  return(output[,"Dose"])
}

#estimation of risk for cold water
risk_cold<-function(dataframe,Deposition_efficiency,Inhalation,dose_response,exposure_frequency,endpoint){
  Dose_1<-outer(dataframe[,1],Deposition_efficiency[,1])+outer(dataframe[,2],Deposition_efficiency[,2])+outer(dataframe[,3],Deposition_efficiency[,3])+outer(dataframe[,4],Deposition_efficiency[,4])+outer(dataframe[,5],Deposition_efficiency[,5])
  Dose_2<-t(Dose_1)*Inhalation
  risk_1<-(1-exp(-dose_response*Dose_2))*endpoint
  risk_2<-1-(1-risk_1)^exposure_frequency
  risk_annual<-t(risk_2)
  return(risk_annual)
}
#estimation of risk for hot water
risk_hot<-function(dataframe,Deposition_efficiency,Inhalation,dose_response,exposure_frequency,endpoint){
  Dose_1<-outer(dataframe[,1],Deposition_efficiency[,1])+outer(dataframe[,2],Deposition_efficiency[,2])+outer(dataframe[,3],Deposition_efficiency[,3])+outer(dataframe[,4],Deposition_efficiency[,4])+outer(dataframe[,5],Deposition_efficiency[,5])+outer(dataframe[,6],Deposition_efficiency[,6])+outer(dataframe[,7],Deposition_efficiency[,7])+outer(dataframe[,8],Deposition_efficiency[,8])+outer(dataframe[,9],Deposition_efficiency[,9])
  Dose_2<-t(Dose_1)*Inhalation
  risk_1<-(1-exp(-dose_response*Dose_2))*endpoint
  risk_2<-1-(1-risk_1)^exposure_frequency
  risk_annual<-t(risk_2)
  return(risk_annual)
}