library(shiny)
library(shiny.fluent)
library(bslib)
library(gridlayout)


#build shiny spp
ui <- fluidPage(
  tags$style(HTML("
    .card-header {
      background-color: #ADD8E6; 
      color: black; 
      padding: 10px; 
      border-radius: 5px; 
      font-weight: bold; 
    }
  ")),
  tags$style(HTML("
    .card-heat {
      background-color: #FFDAB9; 
      color: black; 
      padding: 10px; 
      border-radius: 5px; 
      font-weight: bold; 
    }
  ")),
  tabsetPanel(
    id = "tabs",
    tabPanel("About the app",
             grid_card("description", "Description",
                       HTML("<p>This app is a developed QMRA framework evaluating the infection risks of L. pneumophila for showers.&nbsp;To calculate the risk, the user needs to input values for the model parameters.&nbsp;The detailed explanation of each model parameter is provided below")),
             grid_card("parameter", "Parameter", 
                       HTML("<p>Parameters related to environmental condition </p>
                              <li>1.&nbsp;Ventilation:&nbsp;ventilation rate for the shower stall or bathroom.</li>
                               <li>2.&nbsp;Shower stall dimension:&nbsp;Length, width and height of the shower stall.</li>
                               <li>3.&nbsp;Water temperature:&nbsp;either cold water (25℃) or hot water (40℃)  condition.</li>
                               <li>4.&nbsp;Showerhead:&nbsp;either conventional showerhead or rain showerhead.</li>
                             <p>(Based on the choice of showerhead types and water temperatures, aerosol generation rates obtained from our empirical study will be provided for the risk estimation. Mean and standard deviations of the aerosol generation rates are available. However, the data we provided are based on public shower stalls in Switzerland. Users can also input the size-resolved aerosol generation rates specfic to their scenarios.)</p>
                               <p>Water quality </p>
                               <li>1.&nbsp;First draw sample concentration:&nbsp;the concentration of L. pneumophila in the first liter of shower water.</li>
                               <li>2.&nbsp;Flushed sample concentration:&nbsp;the average concentration of L. pneumophila in composite shower water samples. Steady state concentration of L. pneumophila during flushing can provide more accurate risk estimate if such data are available. </li>
                               <p>Parameters related to human and human activities  </p>
                               <li>1.&nbsp;Shower time:&nbsp;the time spent for a shower (from turning on showers to turning off showers).</li>
                               <li>2.&nbsp;Shower frequency:&nbsp;how many times of showers do the person take per day?</li>
                               <li>3.&nbsp;Morbidity:&nbsp;users can chose different demographic groups and input the specific morbidity ratios based on epidemiology data available</li>
                            "))
             
    ),
    tabPanel("QMRA",
             grid_page(
               layout = c(
                 "      300px   1fr    1fr  1fr ",
                 "100px header  header header header",
                 "1fr   sidebar_a plot_a plot_b plot_c",
                 "1fr   sidebar_a plot_d plot_d sidebar_b"
               ),
               grid_card_text("header", "Quantitative microbial risk assessment for LeCo project"),
               grid_card("sidebar_a", "Exposure parameters",
                         class="card-header",
                         numericInput("Ventilation", "Ventilation (m3/min)", value = 0.7),
                         numericInput("Con_leg_init", "First draw sample concentration (CFU/L)", value = 1000),
                         numericInput("Con_leg_steady", "Flushed sample concentration (CFU/L)", value = 100),
                         numericInput("Shower_length", "Length of shower stall (m)", value = 1.2),
                         numericInput("Shower_width", "Width of shower stall (m)", value = 0.9),
                         numericInput("Shower_height", "Height of shower stall (m)", value = 2.3),
                         numericInput("Shower_frequency", "Shower frequency (#/day)", value = 1),
                         numericInput("Shower_duration", "Shower time (min)", value = 15),
                         checkboxGroupInput("Water_tem", "Water temperature", choices = c("Hot water" = 1, "Cold water" = 2), selected = 1),
                         checkboxGroupInput("showerhead", "Showerhead type", choices = c("Conventional showerhead" = 1, "Rain showerhead" = 2), selected = 1),
                         selectInput("Generation", "Aerosol generation rate", choices = c("Mean" = 1, "Mean + standard deviation" = 2, "Mean - standard deviation" = 3), selected = 1),selectInput("aerosol", "Do you want to use your own aerosol generation rates?", choices = c("No" = 1, "Yes" = 2), selected = 1),uiOutput("dynamic")),
               grid_card("sidebar_b", "Parameters for demographic groups",
                         class="card-heat",
                         checkboxGroupInput("Target_group", "Targeted group of people", choices = c("Adult" = 1, "Elderly" = 2, "Child" = 3), selected = 1),
                         uiOutput("demographic")),
               grid_card("plot_a", "Conventional showerhead",plotOutput("cumulative_risk_conventional")),
               grid_card("plot_b", "Rain showerhead", plotOutput("cumulative_risk_rain")),
               grid_card("plot_c", "User defined scenario", plotOutput("cumulative_risk_user")),
               grid_card("plot_d", "Risk for different demographic groups", plotOutput("risk_group"))
             )
    )
  )
)

server <- function(input, output) {
source("C:/Users/lizha/Desktop/Content.R")
  Generation_rate_hot_conv_init<-reactive({if (input$Generation==1) {Generation_rate_hot_conv_a} else if (input$Generation==2) {Generation_rate_hot_conv_a+Generation_rate_hot_conv_a_sd} else if (input$Generation==3) {Generation_rate_hot_conv_a-Generation_rate_hot_conv_a_sd}}) 
  Generation_rate_hot_conv_steady<-reactive({if (input$Generation==1) {Generation_rate_hot_conv_b} else if (input$Generation==2) {Generation_rate_hot_conv_b+Generation_rate_hot_conv_b_sd} else if (input$Generation==3) {Generation_rate_hot_conv_b-Generation_rate_hot_conv_b_sd}})   
  Generation_rate_hot_rain_init<-reactive({if (input$Generation==1) {Generation_rate_hot_rain_a} else if (input$Generation==2) {Generation_rate_hot_rain_a+Generation_rate_hot_rain_a_sd} else if (input$Generation==3) {Generation_rate_hot_rain_a-Generation_rate_hot_rain_a_sd}}) 
  Generation_rate_hot_rain_steady<-reactive({if (input$Generation==1) {Generation_rate_hot_rain_b} else if (input$Generation==2) {Generation_rate_hot_rain_b+Generation_rate_hot_rain_b_sd} else if (input$Generation==3) {Generation_rate_hot_rain_b-Generation_rate_hot_rain_b_sd}})   
  Generation_rate_cold_conv_steady<-reactive({if (input$Generation==1) {Generation_rate_cold_conv} else if (input$Generation==2) {Generation_rate_cold_conv+Generation_rate_cold_conv_sd} else if (input$Generation==3) {Generation_rate_cold_conv-Generation_rate_cold_conv_sd}}) 
  Generation_rate_cold_rain_steady<-reactive({if (input$Generation==1) {Generation_rate_cold_rain} else if (input$Generation==2) {Generation_rate_cold_rain+Generation_rate_cold_rain_sd} else if (input$Generation==3) {Generation_rate_cold_rain-Generation_rate_cold_rain_sd}}) 
  Morbidity_adult_define<-reactive({if (1 %in% input$Target_group) {input$Morbidity_adult} else {Morbidity_adult}})
  Morbidity_elderly_define<-reactive({if (2 %in% input$Target_group) {input$Morbidity_elderly} else {Morbidity_elderly}})
  Morbidity_child_define<-reactive({if (3 %in% input$Target_group) {input$Morbidity_child} else {Morbidity_child}})
  
  #cold water conventional showerhead
  Dose_cold_conv_1<-reactive({aerosol_cold(input$Shower_duration,input$Shower_length,input$Shower_width,input$Shower_height,1,Generation_rate_cold_conv_steady()[1],Decay_cold_conv_a[1],Decay_cold_conv_b[1],input$Ventilation,input$Con_leg_init*1000,input$Con_leg_steady*1000,Fraction[1])})
  Dose_cold_conv_2<-reactive({aerosol_cold(input$Shower_duration,input$Shower_length,input$Shower_width,input$Shower_height,2,Generation_rate_cold_conv_steady()[2],Decay_cold_conv_a[2],Decay_cold_conv_b[2],input$Ventilation,input$Con_leg_init*1000,input$Con_leg_steady*1000,Fraction[2])})
  Dose_cold_conv_3<-reactive({aerosol_cold(input$Shower_duration,input$Shower_length,input$Shower_width,input$Shower_height,3,Generation_rate_cold_conv_steady()[3],Decay_cold_conv_a[3],Decay_cold_conv_b[3],input$Ventilation,input$Con_leg_init*1000,input$Con_leg_steady*1000,Fraction[3])})
  Dose_cold_conv_4<-reactive({aerosol_cold(input$Shower_duration,input$Shower_length,input$Shower_width,input$Shower_height,4,Generation_rate_cold_conv_steady()[4],Decay_cold_conv_a[4],Decay_cold_conv_b[4],input$Ventilation,input$Con_leg_init*1000,input$Con_leg_steady*1000,Fraction[4])})
  Dose_cold_conv_5<-reactive({aerosol_cold(input$Shower_duration,input$Shower_length,input$Shower_width,input$Shower_height,5,Generation_rate_cold_conv_steady()[5],Decay_cold_conv_a[5],Decay_cold_conv_b[5],input$Ventilation,input$Con_leg_init*1000,input$Con_leg_steady*1000,Fraction[5])})
  Dose_cold_conv_combine <- reactive({data.frame(Dose_cold_conv_1(),Dose_cold_conv_2(),Dose_cold_conv_3(),Dose_cold_conv_4(),Dose_cold_conv_5())})
  risk_cold_conv_general<-reactive({risk_cold(Dose_cold_conv_combine(),DE,Inhalation,r_inf,input$Shower_frequency*365,1)})
  risk_cold_conv_adult<-reactive({risk_cold(Dose_cold_conv_combine(),DE,Inhalation,r_inf,input$Shower_frequency*365,Morbidity_adult_define())[input$Shower_duration/0.01+1,]})
  risk_cold_conv_elderly<-reactive({risk_cold(Dose_cold_conv_combine(),DE,Inhalation,r_inf,input$Shower_frequency*365,Morbidity_elderly_define())[input$Shower_duration/0.01+1,]})
  risk_cold_conv_child<-reactive({risk_cold(Dose_cold_conv_combine(),DE,Inhalation,r_inf,input$Shower_frequency*365,Morbidity_child_define())[input$Shower_duration/0.01+1,]})
  risk_cold_child_conv_data<-reactive({data.frame(Risk=as.vector(risk_cold_conv_child()),Group=rep("Child",1000),Condition=rep("Cold water_conventional showerhead",1000),Temperature=rep("Cold water",1000),Type=rep("Conventional showerhead",1000))})
  risk_cold_adult_conv_data<-reactive({data.frame(Risk=as.vector(risk_cold_conv_adult()),Group=rep("Adult",1000),Condition=rep("Cold water_conventional showerhead",1000),Temperature=rep("Cold water",1000),Type=rep("Conventional showerhead",1000))})
  risk_cold_elderly_conv_data<-reactive({data.frame(Risk=as.vector(risk_cold_conv_elderly()),Group=rep("Elderly",1000),Condition=rep("Cold water_conventional showerhead",1000),Temperature=rep("Cold water",1000),Type=rep("Conventional showerhead",1000))})
  median_risk_cold_conv<-reactive({apply(risk_cold_conv_general(),1,median)})
  lower_risk_cold_conv<-reactive({apply(risk_cold_conv_general(),1,quantile,probs=c(.25))})
  higher_risk_cold_conv<-reactive({apply(risk_cold_conv_general(),1,quantile,probs=c(.75))})
  risk_cold_conv_time<-reactive({data.frame(Time=seq(from=0,to=input$Shower_duration,by=0.01),risk_median=median_risk_cold_conv(),risk_lower=lower_risk_cold_conv(),risk_higher=higher_risk_cold_conv(),Condition=rep("Cold water",input$Shower_duration/0.01+1))})
  
  #cold water rain showerhead
  Dose_cold_rain_1<-reactive({aerosol_cold(input$Shower_duration,input$Shower_length,input$Shower_width,input$Shower_height,1,Generation_rate_cold_rain_steady()[1],Decay_cold_rain_a[1],Decay_cold_rain_b[1],input$Ventilation,input$Con_leg_init*1000,input$Con_leg_steady*1000,Fraction[1])})
  Dose_cold_rain_2<-reactive({aerosol_cold(input$Shower_duration,input$Shower_length,input$Shower_width,input$Shower_height,2,Generation_rate_cold_rain_steady()[2],Decay_cold_rain_a[2],Decay_cold_rain_b[2],input$Ventilation,input$Con_leg_init*1000,input$Con_leg_steady*1000,Fraction[2])})
  Dose_cold_rain_3<-reactive({aerosol_cold(input$Shower_duration,input$Shower_length,input$Shower_width,input$Shower_height,3,Generation_rate_cold_rain_steady()[3],Decay_cold_rain_a[3],Decay_cold_rain_b[3],input$Ventilation,input$Con_leg_init*1000,input$Con_leg_steady*1000,Fraction[3])})
  Dose_cold_rain_4<-reactive({aerosol_cold(input$Shower_duration,input$Shower_length,input$Shower_width,input$Shower_height,4,Generation_rate_cold_rain_steady()[4],Decay_cold_rain_a[4],Decay_cold_rain_b[4],input$Ventilation,input$Con_leg_init*1000,input$Con_leg_steady*1000,Fraction[4])})
  Dose_cold_rain_5<-reactive({aerosol_cold(input$Shower_duration,input$Shower_length,input$Shower_width,input$Shower_height,5,Generation_rate_cold_rain_steady()[5],Decay_cold_rain_a[5],Decay_cold_rain_b[5],input$Ventilation,input$Con_leg_init*1000,input$Con_leg_steady*1000,Fraction[5])})
  Dose_cold_rain_combine <- reactive({data.frame(Dose_cold_rain_1(),Dose_cold_rain_2(),Dose_cold_rain_3(),Dose_cold_rain_4(),Dose_cold_rain_5())})
  risk_cold_rain_general<-reactive({risk_cold(Dose_cold_rain_combine(),DE,Inhalation,r_inf,input$Shower_frequency*365,1)})
  risk_cold_rain_adult<-reactive({risk_cold(Dose_cold_rain_combine(),DE,Inhalation,r_inf,input$Shower_frequency*365,Morbidity_adult_define())[input$Shower_duration/0.01+1,]})
  risk_cold_rain_elderly<-reactive({risk_cold(Dose_cold_rain_combine(),DE,Inhalation,r_inf,input$Shower_frequency*365,Morbidity_elderly_define())[input$Shower_duration/0.01+1,]})
  risk_cold_rain_child<-reactive({risk_cold(Dose_cold_rain_combine(),DE,Inhalation,r_inf,input$Shower_frequency*365,Morbidity_child_define())[input$Shower_duration/0.01+1,]})
  risk_cold_child_rain_data<-reactive({data.frame(Risk=as.vector(risk_cold_rain_child()),Group=rep("Child",1000),Condition=rep("Cold water_rain showerhead",1000),Temperature=rep("Cold water",1000),Type=rep("Rain showerhead",1000))})
  risk_cold_adult_rain_data<-reactive({data.frame(Risk=as.vector(risk_cold_rain_adult()),Group=rep("Adult",1000),Condition=rep("Cold water_rain showerhead",1000),Temperature=rep("Cold water",1000),Type=rep("Rain showerhead",1000))})
  risk_cold_elderly_rain_data<-reactive({data.frame(Risk=as.vector(risk_cold_rain_elderly()),Group=rep("Elderly",1000),Condition=rep("Cold water_rain showerhead",1000),Temperature=rep("Cold water",1000),Type=rep("Rain showerhead",1000))})
  median_risk_cold_rain<-reactive({apply(risk_cold_rain_general(),1,median)})
  lower_risk_cold_rain<-reactive({apply(risk_cold_rain_general(),1,quantile,probs=c(.25))})
  higher_risk_cold_rain<-reactive({apply(risk_cold_rain_general(),1,quantile,probs=c(.75))})
  risk_cold_rain_time<-reactive({data.frame(Time=seq(from=0,to=input$Shower_duration,by=0.01),risk_median=median_risk_cold_rain(),risk_lower=lower_risk_cold_rain(),risk_higher=higher_risk_cold_rain(),Condition=rep("Cold water",input$Shower_duration/0.01+1))})
  
  #hot water conventional showerhead
  Dose_hot_conv_1<-reactive({aerosol_hot(input$Shower_duration,input$Shower_length,input$Shower_width,input$Shower_height,1,Generation_rate_hot_conv_init()[1],Generation_rate_hot_conv_steady()[1],Decay_hot_conv_a[1],Decay_hot_conv_b[1],input$Ventilation,input$Con_leg_init*1000,input$Con_leg_steady*1000,Fraction[1])})
  Dose_hot_conv_2<-reactive({aerosol_hot(input$Shower_duration,input$Shower_length,input$Shower_width,input$Shower_height,2,Generation_rate_hot_conv_init()[2],Generation_rate_hot_conv_steady()[2],Decay_hot_conv_a[2],Decay_hot_conv_b[2],input$Ventilation,input$Con_leg_init*1000,input$Con_leg_steady*1000,Fraction[2])})
  Dose_hot_conv_3<-reactive({aerosol_hot(input$Shower_duration,input$Shower_length,input$Shower_width,input$Shower_height,3,Generation_rate_hot_conv_init()[3],Generation_rate_hot_conv_steady()[2],Decay_hot_conv_a[3],Decay_hot_conv_b[3],input$Ventilation,input$Con_leg_init*1000,input$Con_leg_steady*1000,Fraction[3])})
  Dose_hot_conv_4<-reactive({aerosol_hot(input$Shower_duration,input$Shower_length,input$Shower_width,input$Shower_height,4,Generation_rate_hot_conv_init()[4],Generation_rate_hot_conv_steady()[4],Decay_hot_conv_a[4],Decay_hot_conv_b[4],input$Ventilation,input$Con_leg_init*1000,input$Con_leg_steady*1000,Fraction[4])})
  Dose_hot_conv_5<-reactive({aerosol_hot(input$Shower_duration,input$Shower_length,input$Shower_width,input$Shower_height,5,Generation_rate_hot_conv_init()[5],Generation_rate_hot_conv_steady()[5],Decay_hot_conv_a[5],Decay_hot_conv_b[5],input$Ventilation,input$Con_leg_init*1000,input$Con_leg_steady*1000,Fraction[5])})
  Dose_hot_conv_6<-reactive({aerosol_hot(input$Shower_duration,input$Shower_length,input$Shower_width,input$Shower_height,6,Generation_rate_hot_conv_init()[6],Generation_rate_hot_conv_steady()[6],Decay_hot_conv_a[6],Decay_hot_conv_b[6],input$Ventilation,input$Con_leg_init*1000,input$Con_leg_steady*1000,Fraction[6])})
  Dose_hot_conv_7<-reactive({aerosol_hot(input$Shower_duration,input$Shower_length,input$Shower_width,input$Shower_height,7,Generation_rate_hot_conv_init()[7],Generation_rate_hot_conv_steady()[7],Decay_hot_conv_a[7],Decay_hot_conv_b[7],input$Ventilation,input$Con_leg_init*1000,input$Con_leg_steady*1000,Fraction[7])})
  Dose_hot_conv_8<-reactive({aerosol_hot(input$Shower_duration,input$Shower_length,input$Shower_width,input$Shower_height,8,Generation_rate_hot_conv_init()[8],Generation_rate_hot_conv_steady()[8],Decay_hot_conv_a[8],Decay_hot_conv_b[8],input$Ventilation,input$Con_leg_init*1000,input$Con_leg_steady*1000,Fraction[8])})
  Dose_hot_conv_9<-reactive({aerosol_hot(input$Shower_duration,input$Shower_length,input$Shower_width,input$Shower_height,9,Generation_rate_hot_conv_init()[9],Generation_rate_hot_conv_steady()[9],Decay_hot_conv_a[9],Decay_hot_conv_b[9],input$Ventilation,input$Con_leg_init*1000,input$Con_leg_steady*1000,Fraction[9])})
  Dose_hot_conv_combine <- reactive({data.frame(Dose_hot_conv_1(),Dose_hot_conv_2(),Dose_hot_conv_3(),Dose_hot_conv_4(),Dose_hot_conv_5(),Dose_hot_conv_6(),Dose_hot_conv_7(),Dose_hot_conv_8(),Dose_hot_conv_9())})
  risk_hot_conv_general<-reactive({risk_hot(Dose_hot_conv_combine(),DE,Inhalation,r_inf,input$Shower_frequency*365,1)})
  risk_hot_conv_adult<-reactive({risk_hot(Dose_hot_conv_combine(),DE,Inhalation,r_inf,input$Shower_frequency*365,Morbidity_adult_define())[input$Shower_duration/0.01+1,]})
  risk_hot_conv_elderly<-reactive({risk_hot(Dose_hot_conv_combine(),DE,Inhalation,r_inf,input$Shower_frequency*365,Morbidity_elderly_define())[input$Shower_duration/0.01+1,]})
  risk_hot_conv_child<-reactive({risk_hot(Dose_hot_conv_combine(),DE,Inhalation,r_inf,input$Shower_frequency*365,Morbidity_child_define())[input$Shower_duration/0.01+1,]})
  risk_hot_conv_child_data<-reactive({data.frame(Risk=as.vector(risk_hot_conv_child()),Group=rep("Child",1000),Condition=rep("Hot water_conventional showerhead",1000),Temperature=rep("Hot water",1000),Type=rep("Conventional showerhead",1000))})
  risk_hot_conv_adult_data<-reactive({data.frame(Risk=as.vector(risk_hot_conv_adult()),Group=rep("Adult",1000),Condition=rep("Hot water_conventional showerhead",1000),Temperature=rep("Hot water",1000),Type=rep("Conventional showerhead",1000))})
  risk_hot_conv_elderly_data<-reactive({data.frame(Risk=as.vector(risk_hot_conv_elderly()),Group=rep("Elderly",1000),Condition=rep("Hot water_conventional showerhead",1000),Temperature=rep("Hot water",1000),Type=rep("Conventional showerhead",1000))})
  median_risk_hot_conv<-reactive({apply(risk_hot_conv_general(),1,median)})
  lower_risk_hot_conv<-reactive({apply(risk_hot_conv_general(),1,quantile,probs=c(.25))})
  higher_risk_hot_conv<-reactive({apply(risk_hot_conv_general(),1,quantile,probs=c(.75))})
  risk_hot_conv_time<-reactive({data.frame(Time=seq(from=0,to=input$Shower_duration,by=0.01),risk_median=median_risk_hot_conv(),risk_lower=lower_risk_hot_conv(),risk_higher=higher_risk_hot_conv(),Condition=rep("Hot water",input$Shower_duration/0.01+1))})
  
  #hot water rain showerhead
  Dose_hot_rain_1<-reactive({aerosol_hot(input$Shower_duration,input$Shower_length,input$Shower_width,input$Shower_height,1,Generation_rate_hot_rain_init()[1],Generation_rate_hot_rain_steady()[1],Decay_hot_rain_a[1],Decay_hot_rain_b[1],input$Ventilation,input$Con_leg_init*1000,input$Con_leg_steady*1000,Fraction[1])})
  Dose_hot_rain_2<-reactive({aerosol_hot(input$Shower_duration,input$Shower_length,input$Shower_width,input$Shower_height,2,Generation_rate_hot_rain_init()[2],Generation_rate_hot_rain_steady()[2],Decay_hot_rain_a[2],Decay_hot_rain_b[2],input$Ventilation,input$Con_leg_init*1000,input$Con_leg_steady*1000,Fraction[2])})
  Dose_hot_rain_3<-reactive({aerosol_hot(input$Shower_duration,input$Shower_length,input$Shower_width,input$Shower_height,3,Generation_rate_hot_rain_init()[3],Generation_rate_hot_rain_steady()[2],Decay_hot_rain_a[3],Decay_hot_rain_b[3],input$Ventilation,input$Con_leg_init*1000,input$Con_leg_steady*1000,Fraction[3])})
  Dose_hot_rain_4<-reactive({aerosol_hot(input$Shower_duration,input$Shower_length,input$Shower_width,input$Shower_height,4,Generation_rate_hot_rain_init()[4],Generation_rate_hot_rain_steady()[4],Decay_hot_rain_a[4],Decay_hot_rain_b[4],input$Ventilation,input$Con_leg_init*1000,input$Con_leg_steady*1000,Fraction[4])})
  Dose_hot_rain_5<-reactive({aerosol_hot(input$Shower_duration,input$Shower_length,input$Shower_width,input$Shower_height,5,Generation_rate_hot_rain_init()[5],Generation_rate_hot_rain_steady()[5],Decay_hot_rain_a[5],Decay_hot_rain_b[5],input$Ventilation,input$Con_leg_init*1000,input$Con_leg_steady*1000,Fraction[5])})
  Dose_hot_rain_6<-reactive({aerosol_hot(input$Shower_duration,input$Shower_length,input$Shower_width,input$Shower_height,6,Generation_rate_hot_rain_init()[6],Generation_rate_hot_rain_steady()[6],Decay_hot_rain_a[6],Decay_hot_rain_b[6],input$Ventilation,input$Con_leg_init*1000,input$Con_leg_steady*1000,Fraction[6])})
  Dose_hot_rain_7<-reactive({aerosol_hot(input$Shower_duration,input$Shower_length,input$Shower_width,input$Shower_height,7,Generation_rate_hot_rain_init()[7],Generation_rate_hot_rain_steady()[7],Decay_hot_rain_a[7],Decay_hot_rain_b[7],input$Ventilation,input$Con_leg_init*1000,input$Con_leg_steady*1000,Fraction[7])})
  Dose_hot_rain_8<-reactive({aerosol_hot(input$Shower_duration,input$Shower_length,input$Shower_width,input$Shower_height,8,Generation_rate_hot_rain_init()[8],Generation_rate_hot_rain_steady()[8],Decay_hot_rain_a[8],Decay_hot_rain_b[8],input$Ventilation,input$Con_leg_init*1000,input$Con_leg_steady*1000,Fraction[8])})
  Dose_hot_rain_9<-reactive({aerosol_hot(input$Shower_duration,input$Shower_length,input$Shower_width,input$Shower_height,9,Generation_rate_hot_rain_init()[9],Generation_rate_hot_rain_steady()[9],Decay_hot_rain_a[9],Decay_hot_rain_b[9],input$Ventilation,input$Con_leg_init*1000,input$Con_leg_steady*1000,Fraction[9])})
  Dose_hot_rain_combine <- reactive({data.frame(Dose_hot_rain_1(),Dose_hot_rain_2(),Dose_hot_rain_3(),Dose_hot_rain_4(),Dose_hot_rain_5(),Dose_hot_rain_6(),Dose_hot_rain_7(),Dose_hot_rain_8(),Dose_hot_rain_9())})
  risk_hot_rain_general<-reactive({risk_hot(Dose_hot_rain_combine(),DE,Inhalation,r_inf,input$Shower_frequency*365,1)})
  risk_hot_rain_adult<-reactive({risk_hot(Dose_hot_rain_combine(),DE,Inhalation,r_inf,input$Shower_frequency*365,Morbidity_adult_define())[input$Shower_duration/0.01+1,]})
  risk_hot_rain_elderly<-reactive({risk_hot(Dose_hot_rain_combine(),DE,Inhalation,r_inf,input$Shower_frequency*365,Morbidity_elderly_define())[input$Shower_duration/0.01+1,]})
  risk_hot_rain_child<-reactive({risk_hot(Dose_hot_rain_combine(),DE,Inhalation,r_inf,input$Shower_frequency*365,Morbidity_child_define())[input$Shower_duration/0.01+1,]})
  risk_hot_rain_child_data<-reactive({data.frame(Risk=as.vector(risk_hot_rain_child()),Group=rep("Child",1000),Condition=rep("Hot water_rain showerhead",1000),Temperature=rep("Hot water",1000),Type=rep("Rain showerhead",1000))})
  risk_hot_rain_adult_data<-reactive({data.frame(Risk=as.vector(risk_hot_rain_adult()),Group=rep("Adult",1000),Condition=rep("Hot water_rain showerhead",1000),Temperature=rep("Hot water",1000),Type=rep("Rain showerhead",1000))})
  risk_hot_rain_elderly_data<-reactive({data.frame(Risk=as.vector(risk_hot_rain_elderly()),Group=rep("Elderly",1000),Condition=rep("Hot water_rain showerhead",1000),Temperature=rep("Hot water",1000),Type=rep("Rain showerhead",1000))})
  median_risk_hot_rain<-reactive({apply(risk_hot_rain_general(),1,median)})
  lower_risk_hot_rain<-reactive({apply(risk_hot_rain_general(),1,quantile,probs=c(.25))})
  higher_risk_hot_rain<-reactive({apply(risk_hot_rain_general(),1,quantile,probs=c(.75))})
  risk_hot_rain_time<-reactive({data.frame(Time=seq(from=0,to=input$Shower_duration,by=0.01),risk_median=median_risk_hot_rain(),risk_lower=lower_risk_hot_rain(),risk_higher=higher_risk_hot_rain(),Condition=rep("Hot water",input$Shower_duration/0.01+1))})
  
  #User defined scenario
  Dose_define_1<-reactive({aerosol_cold(input$Shower_duration,input$Shower_length,input$Shower_width,input$Shower_height,1,1e-6*input$Generation_1,0,0,input$Ventilation,input$Con_leg_init*1000,input$Con_leg_steady*1000,Fraction[1])})
  Dose_define_2<-reactive({aerosol_cold(input$Shower_duration,input$Shower_length,input$Shower_width,input$Shower_height,2,1e-6*input$Generation_2,0,0,input$Ventilation,input$Con_leg_init*1000,input$Con_leg_steady*1000,Fraction[2])})
  Dose_define_3<-reactive({aerosol_cold(input$Shower_duration,input$Shower_length,input$Shower_width,input$Shower_height,3,1e-6*input$Generation_3,0,0,input$Ventilation,input$Con_leg_init*1000,input$Con_leg_steady*1000,Fraction[3])})
  Dose_define_4<-reactive({aerosol_cold(input$Shower_duration,input$Shower_length,input$Shower_width,input$Shower_height,4,1e-6*input$Generation_4,0,0,input$Ventilation,input$Con_leg_init*1000,input$Con_leg_steady*1000,Fraction[4])})
  Dose_define_5<-reactive({aerosol_cold(input$Shower_duration,input$Shower_length,input$Shower_width,input$Shower_height,5,1e-6*input$Generation_5,0,0,input$Ventilation,input$Con_leg_init*1000,input$Con_leg_steady*1000,Fraction[5])})
  Dose_define_6<-reactive({aerosol_cold(input$Shower_duration,input$Shower_length,input$Shower_width,input$Shower_height,6,1e-6*input$Generation_6,0,0,input$Ventilation,input$Con_leg_init*1000,input$Con_leg_steady*1000,Fraction[6])})
  Dose_define_7<-reactive({aerosol_cold(input$Shower_duration,input$Shower_length,input$Shower_width,input$Shower_height,7,1e-6*input$Generation_7,0,0,input$Ventilation,input$Con_leg_init*1000,input$Con_leg_steady*1000,Fraction[7])})
  Dose_define_8<-reactive({aerosol_cold(input$Shower_duration,input$Shower_length,input$Shower_width,input$Shower_height,8,1e-6*input$Generation_8,0,0,input$Ventilation,input$Con_leg_init*1000,input$Con_leg_steady*1000,Fraction[8])})
  Dose_define_9<-reactive({aerosol_cold(input$Shower_duration,input$Shower_length,input$Shower_width,input$Shower_height,9,1e-6*input$Generation_9,0,0,input$Ventilation,input$Con_leg_init*1000,input$Con_leg_steady*1000,Fraction[9])})
  Dose_define_combine <- reactive({data.frame(Dose_define_1(),Dose_define_2(),Dose_define_3(),Dose_define_4(),Dose_define_5(),Dose_define_6(),Dose_define_7(),Dose_define_8(),Dose_define_9())})
  risk_define_general<-reactive({risk_hot(Dose_define_combine(),DE,Inhalation,r_inf,input$Shower_frequency*365,1)})
  risk_define_adult<-reactive({risk_hot(Dose_define_combine(),DE,Inhalation,r_inf,input$Shower_frequency*365,Morbidity_adult_define())[input$Shower_duration/0.01+1,]})
  risk_define_elderly<-reactive({risk_hot(Dose_define_combine(),DE,Inhalation,r_inf,input$Shower_frequency*365,Morbidity_elderly_define())[input$Shower_duration/0.01+1,]})
  risk_define_child<-reactive({risk_hot(Dose_define_combine(),DE,Inhalation,r_inf,input$Shower_frequency*365,Morbidity_child_define())[input$Shower_duration/0.01+1,]})
  risk_define_child_data<-reactive({data.frame(Risk=as.vector(risk_define_child()),Group=rep("Child",1000),Condition=rep("User defined scenario",1000),Temperature=rep("User defined",1000),Type=rep("User defined",1000))})
  risk_define_adult_data<-reactive({data.frame(Risk=as.vector(risk_define_adult()),Group=rep("Adult",1000),Condition=rep("User defined scenario",1000),Temperature=rep("User defined",1000),Type=rep("User defined",1000))})
  risk_define_elderly_data<-reactive({data.frame(Risk=as.vector(risk_define_elderly()),Group=rep("Elderly",1000),Condition=rep("User defined scenario",1000),Temperature=rep("User defined",1000),Type=rep("User defined",1000))})
  median_risk_define<-reactive({apply(risk_define_general(),1,median)})
  lower_risk_define<-reactive({apply(risk_define_general(),1,quantile,probs=c(.25))})
  higher_risk_define<-reactive({apply(risk_define_general(),1,quantile,probs=c(.75))})
  risk_define_time<-reactive({data.frame(Time=seq(from=0,to=input$Shower_duration,by=0.01),risk_median=median_risk_define(),risk_lower=lower_risk_define(),risk_higher=higher_risk_define(),Condition=rep("User defined scenario",input$Shower_duration/0.01+1))})
  
  risk_group_combine<-reactive({rbind(risk_cold_elderly_conv_data(),risk_cold_adult_conv_data(),risk_cold_child_conv_data(),risk_cold_adult_rain_data(),risk_cold_elderly_rain_data(),risk_cold_child_rain_data(),risk_hot_conv_adult_data(),risk_hot_conv_elderly_data(),risk_hot_conv_child_data(),risk_hot_rain_adult_data(),risk_hot_rain_child_data(),risk_hot_rain_elderly_data())})
  data_blank<-data.frame(Risk=c(0),Group=c("Demographic group"),Condition=c("Water temperature"))
  data_blank_2<-reactive({data.frame(Time=seq(from=0,to=input$Shower_duration,by=0.01),risk_median=rep(0,input$Shower_duration/0.01+1),risk_lower=rep(0,input$Shower_duration/0.01+1),risk_higher=rep(0,input$Shower_duration/0.01+1),Condition=rep("Water temperature",input$Shower_duration/0.01+1))})
  #output
  output$dynamic<-renderUI({
    if (input$aerosol==2) {
      tagList(numericInput("Generation_1","Generation rate for aerosol of 1-2 μm (g/min)",value = 1e-5),
      numericInput("Generation_2","Generation rate for aerosol of 2-3 μm (g/min)",value = 1e-5),
      numericInput("Generation_3","Generation rate for aerosol of 3-4 μm (g/min)",value = 1e-5),
      numericInput("Generation_4","Generation rate for aerosol of 4-5 μm (g/min)",value = 1e-5),
      numericInput("Generation_5","Generation rate for aerosol of 5-6 μm (g/min)",value = 1e-5),
      numericInput("Generation_6","Generation rate for aerosol of 6-7 μm (g/min)",value = 1e-5),
      numericInput("Generation_7","Generation rate for aerosol of 7-8 μm (g/min)",value = 1e-5),
      numericInput("Generation_8","Generation rate for aerosol of 8-9 μm (g/min)",value = 1e-5),
      numericInput("Generation_9","Generation rate for aerosol of 9-10 μm (g/min)",value = 1e-5))
    }
    })
  output$demographic<-renderUI({
    if(identical(input$Target_group,c("1","2","3"))){
      tagList(numericInput("Morbidity_adult","Morbidity ratio for adults",value = Morbidity_adult),
              numericInput("Morbidity_elderly","Morbidity ratio for the elderly",value = Morbidity_elderly),
              numericInput("Morbidity_child","Morbidity ratio for children",value = Morbidity_child)
              )
    } else if (identical(input$Target_group,c("1","2"))) {
      tagList(numericInput("Morbidity_adult","Morbidity ratio for adults",value = Morbidity_adult),
              numericInput("Morbidity_elderly","Morbidity ratio for the elderly",value = Morbidity_elderly))
    } else if (identical(input$Target_group,c("1","3"))) {
      tagList(numericInput("Morbidity_adult","Morbidity ratio for adults",value = Morbidity_adult),
              numericInput("Morbidity_child","Morbidity ratio for children",value = Morbidity_child))
    } else if (identical(input$Target_group,c("2","3"))) {
      tagList(numericInput("Morbidity_elderly","Morbidity ratio for the elderly",value = Morbidity_elderly),
              numericInput("Morbidity_child","Morbidity ratio for children",value = Morbidity_child))
    } else if (identical(input$Target_group,c("1"))) {
      tagList(numericInput("Morbidity_adult","Morbidity ratio for adults",value = Morbidity_adult))
    } else if (identical(input$Target_group,c("2"))) {
      tagList(numericInput("Morbidity_elderly","Morbidity ratio for the elderly",value = Morbidity_elderly))
    } else if (identical(input$Target_group,c("3"))) {
      tagList(numericInput("Morbidity_child","Morbidity ratio for children",value = Morbidity_child))
    }
  })
  
  
  output$cumulative_risk_conventional <- renderPlot({
    risk_time_conv <- if (identical(input$Water_tem,c("1")) && 1 %in% input$showerhead) { risk_hot_conv_time() } else if (identical(input$Water_tem,c("2")) && 1 %in% input$showerhead) { risk_cold_conv_time() } else if (identical(input$Water_tem,c("1","2")) && 1 %in% input$showerhead) {rbind(risk_cold_conv_time(),risk_hot_conv_time())} else {data_blank_2()}
    ggplot(risk_time_conv,aes(x=Time,y=risk_median,group=Condition))+geom_line()+scale_y_continuous(trans='log10',label=scientific_10)+scale_x_continuous(breaks = seq(0,input$Shower_duration,b=1))+xlab("Time (min)")+ylab("Cumulative risk of infection")+geom_ribbon(aes(x=Time,ymax=risk_higher,ymin=risk_lower,fill=Condition),show.legend = TRUE,linetype="dashed",alpha=0.5,col="black")+scale_fill_brewer(palette="Pastel1",direction=-1)+theme(panel.grid.major = element_blank(), 
                                                                                                                                                                                                                                                                                                                                                                                                                                  panel.grid.minor = element_blank(),panel.background = element_blank(),
                                                                                                                                                                                                                                                                                                                                                                                                                                  panel.border = element_rect(colour = "black", fill=NA, size=2),text=element_text(size=20))
  })
  output$cumulative_risk_rain <- renderPlot({
    risk_time_rain <- if (identical(input$Water_tem,c("1")) && 2 %in% input$showerhead) { risk_hot_rain_time() } else if (identical(input$Water_tem,c("2")) && 2 %in% input$showerhead) { risk_cold_rain_time() } else if (identical(input$Water_tem,c("1","2")) && 2 %in% input$showerhead) {rbind(risk_cold_rain_time(),risk_hot_rain_time())} else {data_blank_2()}
    ggplot(risk_time_rain,aes(x=Time,y=risk_median,group=Condition))+geom_line()+scale_y_continuous(trans='log10',label=scientific_10)+scale_x_continuous(breaks = seq(0,input$Shower_duration,b=1))+xlab("Time (min)")+ylab("Cumulative risk of infection")+geom_ribbon(aes(x=Time,ymax=risk_higher,ymin=risk_lower,fill=Condition),show.legend = TRUE,linetype="dashed",alpha=0.5,col="black")+scale_fill_brewer(palette="Pastel1",direction=-1)+theme(panel.grid.major = element_blank(), 
                                                                                                                                                                                                                                                                                                                                                                                                                                                     panel.grid.minor = element_blank(),panel.background = element_blank(),
                                                                                                                                                                                                                                                                                                                                                                                                                                                     panel.border = element_rect(colour = "black", fill=NA, size=2),text=element_text(size=20))
  })
  output$cumulative_risk_user <- renderPlot({
    risk_time_user <- if (input$aerosol==2) { risk_define_time() } else {data_blank_2()} 
    ggplot(risk_time_user,aes(x=Time,y=risk_median))+geom_line()+scale_y_continuous(trans='log10',label=scientific_10)+scale_x_continuous(breaks = seq(0,input$Shower_duration,b=1))+xlab("Time (min)")+ylab("Cumulative risk of infection")+geom_ribbon(aes(x=Time,ymax=risk_higher,ymin=risk_lower,fill=Condition),show.legend = FALSE,linetype="dashed",alpha=0.5,col="black")+scale_fill_brewer(palette="Greys",direction=1)+theme(panel.grid.major = element_blank(), 
                                                                                                                                                                                                                                                                                                                                                                                                                                                          panel.grid.minor = element_blank(),panel.background = element_blank(),
                                                                                                                                                                                                                                                                                                                                                                                                                                                          panel.border = element_rect(colour = "black", fill=NA, size=2),text=element_text(size=20))
  })
  output$risk_group <- renderPlot({
    risk_people<-if (input$aerosol==1) {risk_group_combine()} else {rbind(risk_group_combine(),risk_define_adult_data(),risk_define_child_data(),risk_define_elderly_data())}
    risk_people_1 <-if (identical(input$Water_tem,c("1"))) {dplyr::filter(risk_people, Temperature %in% c("User defined","Hot water"))} else if (identical(input$Water_tem,c("2"))) {dplyr::filter(risk_people, Temperature %in% c("User defined","Cold water"))} else if (identical(input$Water_tem,c("1","2"))) {risk_people} else if (input$aerosol==2 && identical(input$Water_tem,c())){dplyr::filter(risk_people, Temperature %in% c("User defined"))} else {data_blank}
    risk_people_2<-if (identical(input$showerhead,c("1"))) {dplyr::filter(risk_people_1, Type %in% c("User defined","Conventional showerhead"))} else if (identical(input$showerhead,c("2"))) {dplyr::filter(risk_people_1, Type %in% c("User defined","Rain showerhead"))} else if (identical(input$showerhead,c("1","2"))) {risk_people_1} else if (input$aerosol==2 && identical(input$showerhead,c())){dplyr::filter(risk_people_1, Type %in% c("User defined"))} else {data_blank}
    risk_people_3<-if (identical(input$Target_group,c("1"))) {dplyr::filter(risk_people_2, Group=="Adult")} else if (identical(input$Target_group,c("2"))) {dplyr::filter(risk_people_2, Group=="Elderly")} else if (identical(input$Target_group,c("3"))) {dplyr::filter(risk_people_2, Group=="Child")} else if (identical(input$Target_group,c("1","2"))) {dplyr::filter(risk_people_2, Group %in% c("Adult","Elderly"))} else if (identical(input$Target_group,c("1","3"))) {dplyr::filter(risk_people_2, Group %in% c("Adult","Child"))} else if (identical(input$Target_group,c("2","3"))) {dplyr::filter(risk_people_2, Group %in% c("Child","Elderly"))} else if (identical(input$Target_group,c("1","2","3"))) {risk_people_2} else {data_blank} 
    
    ggplot(risk_people_3,mapping=aes(Group,Risk,fill=Condition))+stat_boxplot(geom = "errorbar", width = 0.4,position=position_dodge(0.8))+geom_boxplot(outlier.shape = NA,width=0.6,position =position_dodge(0.8))+ylab("Risk of illness")+scale_y_continuous(trans='log10',breaks = c(1e-10,1e-9,1e-8,1e-7,1e-6,1e-5,1e-4,1e-3,1e-2,1e-1,1),label=scientific_10)+theme(panel.grid.major = element_blank(), 
                                                                                                                                                                                                                                                                                                                                                                              panel.grid.minor = element_blank(),panel.background = element_blank(),
                                                                                                                                                                                                                                                                                                                                                                              panel.border = element_rect(colour = "black", fill=NA, size=2),text = element_text(size = 18))+scale_fill_brewer(palette="RdBu",direction=-1)
    
})
}

shinyApp(ui = ui, server = server)