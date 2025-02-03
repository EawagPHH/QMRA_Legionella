source("C:/Users/lizha/Desktop/Content.R")
library(gridExtra)


#alternative scenario modelling (Hot Water)
#Change Initial Legionella concentration
risk_initial_legionella<-numeric()
for(j in 1:6){
Dose_hot_combine<-numeric()
for (i in 1:9){
  Dose_hot_combine<-cbind(Dose_hot_combine,aerosol_hot(15,Length,Width,Height,i,Generation_rate_hot_rain_a[i],Generation_rate_hot_rain_b[i],Decay_hot_rain_a[i],Decay_hot_rain_b[i],0.7,1e6-(j-1)*((1E6-0)/5),1e5,Fraction[i]))
  
}
c_risk<-risk_hot(Dose_hot_combine,DE,Inhalation,r_inf,365,1)
median_risk<-apply(c_risk,1,median)
risk_initial_legionella<-cbind(risk_initial_legionella,median_risk[1501])
}
risk_initial_legionella<-data.frame(Change=c(0,-1,-2,-3,-4,-5)*(1E6/5)/1e6,Risk=as.vector(risk_initial_legionella),Parameter=rep("Initial pathogen concentration in water",6))
#Change steady state Legionella concentration
risk_steady_legionella<-numeric()
for(j in 1:6){
  Dose_hot_combine<-numeric()
  for (i in 1:9){
    Dose_hot_combine<-cbind(Dose_hot_combine,aerosol_hot(15,Length,Width,Height,i,Generation_rate_hot_rain_a[i],Generation_rate_hot_rain_b[i],Decay_hot_rain_a[i],Decay_hot_rain_b[i],0.7,1e6,1e5-(j-1)*((1e5-0)/5),Fraction[i]))
    
  }
  c_risk<-risk_hot(Dose_hot_combine,DE,Inhalation,r_inf,365,1)
  median_risk<-apply(c_risk,1,median)
  risk_steady_legionella<-cbind(risk_steady_legionella,median_risk[1501])
}
risk_steady_legionella<-data.frame(Change=c(0,-1,-2,-3,-4,-5)*(1e5/5)/1e5,Risk=as.vector(risk_steady_legionella),Parameter=rep("Steady state pathogen concentration in water",6))

#Change initial aerosol generation rate
risk_aerosol_generation_init<-numeric()
for(j in 1:6){
  Dose_hot_combine<-numeric()
  for (i in 1:9){
    Dose_hot_combine<-cbind(Dose_hot_combine,aerosol_hot(15,Length,Width,Height,i,Generation_rate_hot_rain_a[i]-(j-1)*((Generation_rate_hot_rain_a[i]-0)/5),Generation_rate_hot_rain_b[i],Decay_hot_rain_a[i],Decay_hot_rain_b[i],0.7,1e6,1e5,Fraction[i]))
    
  }
  c_risk<-risk_hot(Dose_hot_combine,DE,Inhalation,r_inf,365,1)
  median_risk<-apply(c_risk,1,median)
  risk_aerosol_generation_init<-cbind(risk_aerosol_generation_init,median_risk[1501])
}
risk_aerosol_generation_init<-data.frame(Change=c(0,-1,-2,-3,-4,-5)*0.2,Risk=as.vector(risk_aerosol_generation_init),Parameter=rep("Initial aerosol generation rate",6))

#Change steady state aerosol generation rate
risk_aerosol_generation_steady<-numeric()
for(j in 1:6){
  Dose_hot_combine<-numeric()
  for (i in 1:9){
    Dose_hot_combine<-cbind(Dose_hot_combine,aerosol_hot(15,Length,Width,Height,i,Generation_rate_hot_rain_a[i],Generation_rate_hot_rain_b[i]-(j-1)*((Generation_rate_hot_rain_b[i]-0)/5),Decay_hot_rain_a[i],Decay_hot_rain_b[i],0.7,1e6,1e5,Fraction[i]))
    
  }
  c_risk<-risk_hot(Dose_hot_combine,DE,Inhalation,r_inf,365,1)
  median_risk<-apply(c_risk,1,median)
  risk_aerosol_generation_steady<-cbind(risk_aerosol_generation_steady,median_risk[1501])
}
risk_aerosol_generation_steady<-data.frame(Change=c(0,-1,-2,-3,-4,-5)*0.2,Risk=as.vector(risk_aerosol_generation_steady),Parameter=rep("Steady state aerosol generation rate",6))

#Change ventilation
risk_ventilation<-numeric()
for(j in 1:6){
  Dose_hot_combine<-numeric()
  for (i in 1:9){
    Dose_hot_combine<-cbind(Dose_hot_combine,aerosol_hot(15,Length,Width,Height,i,Generation_rate_hot_rain_a[i],Generation_rate_hot_rain_b[i],Decay_hot_rain_a[i],Decay_hot_rain_b[i],0.7+(j-1)*((1.4-0.7)/5),1e6,1e5,Fraction[i]))
    
  }
  c_risk<-risk_hot(Dose_hot_combine,DE,Inhalation,r_inf,365,1)
  median_risk<-apply(c_risk,1,median)
  risk_ventilation<-cbind(risk_ventilation,median_risk[1501])
}
risk_ventilation<-data.frame(Change=c(0,1,2,3,4,5)*((1.4-0.7)/5)/0.7,Risk=as.vector(risk_ventilation),Parameter=rep("Ventilation rate",6))

#Change shower stall
risk_stall<-numeric()
for(j in 1:6){
  Dose_hot_combine<-numeric()
  for (i in 1:9){
    Dose_hot_combine<-cbind(Dose_hot_combine,aerosol_hot(15,Length+(j-1)*(2.4-Length)/5,Width+(j-1)*(1.8-Width)/5,Height+(j-1)*(1.8-Height)/5,i,Generation_rate_hot_rain_a[i],Generation_rate_hot_rain_b[i],Decay_hot_rain_a[i],Decay_hot_rain_b[i],0.7,1e6,1e5,Fraction[i]))
    
  }
  c_risk<-risk_hot(Dose_hot_combine,DE,Inhalation,r_inf,365,1)
  median_risk<-apply(c_risk,1,median)
  risk_stall<-cbind(risk_stall,median_risk[1501])
}
risk_stall<-data.frame(Change=c(0,1,2,3,4,5)*(2.4-Length)/5/Length,Risk=as.vector(risk_stall),Parameter=rep("Shower stall volume",6))



#Change exposure frequency
risk_exposure_frequency<-numeric()
for(j in 1:6){
  Dose_hot_combine<-numeric()
  for (i in 1:9){
    Dose_hot_combine<-cbind(Dose_hot_combine,aerosol_hot(15,Length,Width,Height,i,Generation_rate_hot_rain_a[i],Generation_rate_hot_rain_b[i],Decay_hot_rain_a[i],Decay_hot_rain_b[i],0.7,1e6,1e5,Fraction[i]))
    
  }
  c_risk<-risk_hot(Dose_hot_combine,DE,Inhalation,r_inf,365-(j-1)*(182)/5,1)
  median_risk<-apply(c_risk,1,median)
  risk_exposure_frequency<-cbind(risk_exposure_frequency,median_risk[1501])
}
risk_exposure_frequency<-data.frame(Change=c(0,-1,-2,-3,-4,-5)*182/5/365,Risk=as.vector(risk_exposure_frequency),Parameter=rep("Shower frequency",6))




data_risk<-rbind(risk_initial_legionella,risk_steady_legionella,risk_aerosol_generation_init,risk_aerosol_generation_steady,risk_exposure_frequency,risk_stall,risk_ventilation)
data_risk$Change<-data_risk$Change*100
ggplot(data=data_risk,aes(Change,Risk,colour = Parameter,shape = Parameter))+xlab("Relative change (%)")+ylab("Annual risk of infection")+scale_y_continuous(trans = "log10",breaks = c(5e-5,1e-4,1.5e-4,2e-4,2.5e-4),limits = c(5e-5,3e-4),labels = scientific_10)+scale_x_continuous(breaks=seq(from=-100,to=100,by=20))+scale_color_manual(values = c("#CC0000","#FFCC33","#339900","#99CCCC","#336699","#9933CC","#FF33FF"))+scale_shape_manual(values=c(15,16,17,18,0,1,2))+geom_point(size=2.5)+geom_line(size=1)+theme(panel.grid.major = element_blank(), 
                                                                                                                                                                                    panel.grid.minor = element_blank(),panel.background = element_blank(),
                                                                                                                                                                                    panel.border = element_rect(colour = "black", fill=NA, size=2),text=element_text(size=20),axis.title.y = element_text( vjust =2.5))
#alternative scenario modelling (Cold Water)
#Change Initial Legionella concentration
risk_initial_legionella_cold<-numeric()
for(j in 1:6){
  Dose_cold_combine<-numeric()
  for (i in 1:5){
    Dose_cold_combine<-cbind(Dose_cold_combine,aerosol_cold(15,Length,Width,Height,i,Generation_rate_cold_rain[i],Decay_cold_rain_a[i],Decay_cold_rain_b[i],0.7,1e8-(j-1)*((1E8-0)/5),1e7,Fraction[i]))
    
  }
  c_risk_cold<-risk_cold(Dose_cold_combine,DE,Inhalation,r_inf,365,1)
  median_risk_cold<-apply(c_risk_cold,1,median)
  risk_initial_legionella_cold<-cbind(risk_initial_legionella_cold,median_risk_cold[1501])
}
risk_initial_legionella_cold<-data.frame(Change=c(0,-1,-2,-3,-4,-5)*(1E7/5)/1e7,Risk=as.vector(risk_initial_legionella_cold),Parameter=rep("Initial pathogen concentration in water",6))
#Change steady state Legionella concentration
risk_steady_legionella_cold<-numeric()
for(j in 1:6){
  Dose_cold_combine<-numeric()
  for (i in 1:5){
    Dose_cold_combine<-cbind(Dose_cold_combine,aerosol_cold(15,Length,Width,Height,i,Generation_rate_cold_rain[i],Decay_cold_rain_a[i],Decay_cold_rain_b[i],0.7,1e8,1e7-(j-1)*((1e7-0)/5),Fraction[i]))
    
  }
  c_risk_cold<-risk_cold(Dose_cold_combine,DE,Inhalation,r_inf,365,1)
  median_risk_cold<-apply(c_risk_cold,1,median)
  risk_steady_legionella_cold<-cbind(risk_steady_legionella_cold,median_risk_cold[1501])
}
risk_steady_legionella_cold<-data.frame(Change=c(0,-1,-2,-3,-4,-5)*(1e7/5)/1e7,Risk=as.vector(risk_steady_legionella_cold),Parameter=rep("Steady state pathogen concentration in water",6))


#Change steady state aerosol generation rate
risk_aerosol_generation_steady_cold<-numeric()
for(j in 1:6){
  Dose_cold_combine<-numeric()
  for (i in 1:5){
    Dose_cold_combine<-cbind(Dose_cold_combine,aerosol_cold(15,Length,Width,Height,i,Generation_rate_cold_rain[i]-(j-1)*((Generation_rate_cold_rain[i]-0)/5),Decay_cold_rain_a[i],Decay_cold_rain_b[i],0.7,1e8,1e7,Fraction[i]))
    
  }
  c_risk_cold<-risk_cold(Dose_cold_combine,DE,Inhalation,r_inf,365,1)
  median_risk_cold<-apply(c_risk_cold,1,median)
  risk_aerosol_generation_steady_cold<-cbind(risk_aerosol_generation_steady_cold,median_risk_cold[1501])
}
risk_aerosol_generation_steady_cold<-data.frame(Change=c(0,-1,-2,-3,-4,-5)*0.2,Risk=as.vector(risk_aerosol_generation_steady_cold),Parameter=rep("Steady state aerosol generation rate",6))

#Change ventilation
risk_ventilation_cold<-numeric()
for(j in 1:6){
  Dose_cold_combine<-numeric()
  for (i in 1:5){
    Dose_cold_combine<-cbind(Dose_cold_combine,aerosol_cold(15,Length,Width,Height,i,Generation_rate_cold_rain[i],Decay_cold_rain_a[i],Decay_cold_rain_b[i],0.7+(j-1)*((1.4-0.7)/5),1e8,1e7,Fraction[i]))
    
  }
  c_risk_cold<-risk_cold(Dose_cold_combine,DE,Inhalation,r_inf,365,1)
  median_risk_cold<-apply(c_risk_cold,1,median)
  risk_ventilation_cold<-cbind(risk_ventilation_cold,median_risk_cold[1501])
}
risk_ventilation_cold<-data.frame(Change=c(0,1,2,3,4,5)*((1.4-0.7)/5)/0.7,Risk=as.vector(risk_ventilation_cold),Parameter=rep("Ventilation rate",6))

#Change shower stall
risk_stall_cold<-numeric()
for(j in 1:6){
  Dose_cold_combine<-numeric()
  for (i in 1:5){
    Dose_cold_combine<-cbind(Dose_cold_combine,aerosol_cold(15,Length+(j-1)*(2.4-Length)/5,Width+(j-1)*(1.8-Width)/5,Height+(j-1)*(1.8-Height)/5,i,Generation_rate_cold_rain[i],Decay_cold_rain_a[i],Decay_cold_rain_b[i],0.7,1e8,1e7,Fraction[i]))
    
  }
  c_risk_cold<-risk_cold(Dose_cold_combine,DE,Inhalation,r_inf,365,1)
  median_risk_cold<-apply(c_risk_cold,1,median)
  risk_stall_cold<-cbind(risk_stall_cold,median_risk_cold[1501])
}
risk_stall_cold<-data.frame(Change=c(0,1,2,3,4,5)*(2.4-Length)/5/Length,Risk=as.vector(risk_stall_cold),Parameter=rep("Shower stall volume",6))



#Change exposure frequency
risk_exposure_frequency_cold<-numeric()
for(j in 1:6){
  Dose_cold_combine<-numeric()
  for (i in 1:5){
    Dose_cold_combine<-cbind(Dose_cold_combine,aerosol_cold(15,Length,Width,Height,i,Generation_rate_cold_rain[i],Decay_cold_rain_a[i],Decay_cold_rain_b[i],0.7,1e8,1e7,Fraction[i]))
    
  }
  c_risk_cold<-risk_cold(Dose_cold_combine,DE,Inhalation,r_inf,365-(j-1)*(182)/5,1)
  median_risk_cold<-apply(c_risk_cold,1,median)
  risk_exposure_frequency_cold<-cbind(risk_exposure_frequency_cold,median_risk_cold[1501])
}
risk_exposure_frequency_cold<-data.frame(Change=c(0,-1,-2,-3,-4,-5)*182/5/365,Risk=as.vector(risk_exposure_frequency_cold),Parameter=rep("Shower frequency",6))




data_risk_cold<-rbind(risk_initial_legionella_cold,risk_steady_legionella_cold,risk_aerosol_generation_steady_cold,risk_exposure_frequency_cold,risk_stall_cold,risk_ventilation_cold)
data_risk_cold$Change<-data_risk_cold$Change*100
ggplot(data=data_risk_cold,aes(Change,Risk,colour = Parameter,shape = Parameter))+xlab("Relative change (%)")+ylab("Annual risk of infection")+scale_y_continuous(trans = "log10",limits = c(5e-5,3e-4),breaks = c(5e-5,1e-4,1.5e-4,2e-4,2.5e-4,3e-4),labels = scientific_10)+scale_x_continuous(breaks=seq(from=-100,to=100,by=20))+scale_shape_manual(values=c(16,17,18,0,1,2))+geom_point(size=2.5)+geom_line(size=1)+theme(panel.grid.major = element_blank(), 
                                                                                                                                                                                                                                                                                                                                                    panel.grid.minor = element_blank(),panel.background = element_blank(),
                                                                                                                                                                                                                                                                                                                                                    panel.border = element_rect(colour = "black", fill=NA, size=2),text=element_text(size=20),axis.title.y = element_text( vjust =2.5))+scale_color_manual(values = c("#FFCC33","#339900","#99CCCC","#336699","#9933CC","#FF33FF"))
