library(MASS)
library(pscl)

#Import the dataset named Con_leg

glm.leg<-glm.nb(Flushed ~ log(First_draw+1), data = Con_leg,
                control = glm.control(maxit = 1000, epsilon = 1e-8))
predictions <- predict(glm.leg,type = "response",se.fit = TRUE)

Con_leg$predicted <- predictions$fit
Con_leg$lower <- predictions$fit - 1.96 * predictions$se.fit
Con_leg$upper <- predictions$fit + 1.96 * predictions$se.fit



ggplot(Con_leg, aes(x = First_draw, y = Flushed)) +
  
  geom_point(alpha = 0.6, size = 2, color = "steelblue") +
  
  
  geom_ribbon(aes(ymin = lower, ymax = upper),
              fill = "skyblue", 
              color = "black",        
              size = 0.5,            
              linetype = "dashed",     
              alpha = 0.3) +          
  
  geom_line(aes(y = predicted), 
            color = "red", 
            size = 1.2,
            lineend = "round") +
  
  labs(
    x = "First draw sample count (CFU)",
    y = "Flushed sample count (CFU)"
  ) +
  theme_bw() +  
  theme(
    panel.border = element_rect(color = "black", size = 1),  
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  scale_x_continuous(trans="log10") +
  scale_y_continuous(trans="log10")



