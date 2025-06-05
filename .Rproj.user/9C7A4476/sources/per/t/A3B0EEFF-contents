#-------------------------------------
# Work space preparation
#-------------------------------------

library(ggplot2)
library(MASS)
library(fitdistrplus)
library(metRology)
library(DAAG)
library(ggdist)
library(gridExtra)
library(lme4)
library(mgcv)


#-------------------------------------
# Data import and carpentry
#-------------------------------------

data <- read.csv("nmc_speeds2.csv")
data$WLH_ID <- as.factor(data$WLH_ID)

#Define when they are moving/resting
data$Active <- 0
data[which(data$est > 0),"Active"] <- 1



#-------------------------------------
# Descriptive statistics
#-------------------------------------

sum(data$Active == 1)
sum(data$Active == 0)

# Fit a scaled t distribution to all of the temperature data
all_temps <- fitdist(data$temp_c,
                     "t.scaled",
                     start=list(df=3,
                                mean = 5,
                                sd = 4))

coef(all_temps)

#Median and 95% quantiles
qt.scaled(c(0.025,0.5,0.975),
          mean = coef(all_temps)[2],
          sd = coef(all_temps)[3],
          df = coef(all_temps)[1])


temps <- seq(-40, 30, length.out = 1000)
temp_ests <- all_temps$estimate
y_temp <- do.call(dt.scaled, c(list(x = temps), as.list(ests)))
temp_dist <- data.frame(temp_c = temps,
                        y = y_temp)

#Subset estimates where they are moving
move <- data[which(data$est > 0),]

# Fit a scaled t distribution to the temperature data when animals are moving
move_temps <- fitdist(move$temp_c,
                      "t.scaled",
                      start=list(df=3,
                                 mean = 5,
                                 sd = 4))

coef(move_temps)

#Median and 95% quantiles
qt.scaled(c(0.025,0.5,0.975),
          mean = coef(move_temps)[2],
          sd = coef(move_temps)[3],
          df = coef(move_temps)[1])


temps <- seq(-40, 30, length.out = 1000)
move_ests <- move_temps$estimate
y_move <- do.call(dt.scaled, c(list(x = temps), as.list(move_ests)))
move_dist <- data.frame(temp_c = temps,
                        y = y_move)

FIG <-
  ggplot() +
  geom_histogram(data = data, aes(x = temp_c, y = ..density..), binwidth = 1, fill = "#CCDBDC", alpha = 0.3) +
  geom_path(data = temp_dist, aes(y=y, x=temp_c), linewidth = 0.4, col = "#CCDBDC") + 
  geom_vline(aes(xintercept = coef(all_temps)[2]), linewidth = 0.4, col = "#CCDBDC") +
  
  geom_histogram(data = move, aes(x = temp_c, y = ..density..), binwidth = 1, fill = "#19535F", alpha = 0.3) +
  geom_path(data = move_dist, aes(y=y, x=temp_c), linewidth = 0.4, col = "#19535F") +
  geom_vline(aes(xintercept = coef(move_temps)[2]), linewidth = 0.4, col = "#19535F") +
  
  ylab(expression(bold(Density)))+
  xlab(expression(bold(Temperature~(degree*C))))+
  
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(size=8, family = "sans", face = "bold"),
        axis.title.x = element_text(size=8, family = "sans", face = "bold"),
        axis.text.y = element_text(size=6, family = "sans"),
        axis.text.x  = element_text(size=6, family = "sans"),
        plot.title = element_text(hjust = -0.05, size = 8, family = "sans", face = "bold"),
        legend.position = "none",
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", color = NA)) +
  scale_x_continuous(limits = c(-40, 30), expand = c(0,0.01)) +
  scale_y_continuous(limits = c(0, 0.07), expand = c(0,0))

#Save the figures
ggsave(FIG,
       width = 5.75, height = 2.5, units = "in",
       dpi = 600,
       bg = "transparent",
       file="NMC_Temps.png")

#-------------------------------------
# Can temperature predict caribou activity
#-------------------------------------

#Can temperature predict if caribou are active or not
null_mod <- glmer(Active ~ 1 + (1 | WLH_ID), family = binomial, data = data)
FIT <- glmer(Active ~ temp_c + (1 | WLH_ID), family = binomial, data = data)
FIT2 <- glmer(Active ~ temp_c+ I(temp_c^2) + (1 | WLH_ID), family = binomial, data = data)
FIT3 <- gam(Active ~ s(temp_c, bs = "tp", k = 5) + s(WLH_ID, bs = 're'),
            family = binomial(link = "logit"),
            data = data,
            method = "REML")
MuMIn::AICc(FIT3,FIT2,FIT,null_mod)


DATA <- data.frame(temp_c = seq(-40,30, 0.01),
                   WLH_ID = as.factor("new"))

Pred <- predict(FIT3, newdata = DATA, se = TRUE, type = "link")

DATA$Active <- exp(Pred$fit)/(1+exp(Pred$fit))
DATA$Active.max <- exp(Pred$fit+(1.96*Pred$se.fit))/(1+exp(Pred$fit+(1.96*Pred$se.fit)))
DATA$Active.min <- exp(Pred$fit-(1.96*Pred$se.fit))/(1+exp(Pred$fit-(1.96*Pred$se.fit)))




#Figure of the speeds
speed_fig <- 
  ggplot(data=data, aes(x=temp_c, y=est)) +
  ggtitle("A") +
  geom_point(size = 0.8, alpha = 0.3,stroke = 0,shape=16, colour = "#3471bc") +
  
  ylab(expression(bold(Movement~speed~(m/s))))+
  xlab(expression(bold(Temperature~(degree*C))))+
  
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(size=8, family = "sans", face = "bold"),
        axis.title.x = element_text(size=8, family = "sans", face = "bold"),
        axis.text.y = element_text(size=6, family = "sans"),
        axis.text.x  = element_text(size=6, family = "sans"),
        plot.title = element_text(hjust = -0.05, size = 8, family = "sans", face = "bold"),
        legend.position = "none",
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", color = NA)) +
  scale_x_continuous(limits = c(-40, 30), expand = c(0,0.01)) +
  scale_y_continuous(limits = c(0, 1), expand = c(0,0.01))



FIT_fig <-
  ggplot(data = data, aes(x = temp_c)) +
  ggtitle("B") +
  
  geom_ribbon(data = DATA, aes(x=temp_c, ymin = Active.min, ymax = Active.max), fill = "grey70", alpha = 0.3) + 
  geom_path(data = DATA, aes(y=Active, x=temp_c), linewidth = 0.4, col = "black") +
  #geom_path(data = DATA, aes(y=Active.max, x=temp_c), linewidth = 0.4, col = "grey80", linetype = "dashed") +
  #geom_path(data = DATA, aes(y=Active.min, x=temp_c), linewidth = 0.4, col = "grey80", linetype = "dashed") +
  
  stat_slab(data = data,
            aes(x = temp_c, y = Active, 
                side = ifelse(Active == 0, "top", "bottom"),
                fill = Active),
            density = "histogram",
            scale = 0.4, breaks = 100, size = 1/2, alpha = 0.6, col = NA) +
  scale_fill_gradient(high = "#3471bc", low = "#e36b20") +
  
  ylab(expression(bold(Probability~of~moving)))+
  xlab(expression(bold(Temperature~(degree*C))))+
  
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(size=8, family = "sans", face = "bold"),
        axis.title.x = element_text(size=8, family = "sans", face = "bold"),
        axis.text.y = element_text(size=6, family = "sans"),
        axis.text.x  = element_text(size=6, family = "sans"),
        plot.title = element_text(hjust = -0.05, size = 8, family = "sans", face = "bold"),
        legend.position = "none",
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", color = NA)) +
  scale_x_continuous(limits = c(-40, 30), expand = c(0,0.01)) +
  scale_y_continuous(limits = c(0, 1), expand = c(0,0.01))



FIG <-
  grid.arrange(speed_fig,
               FIT_fig,
               ncol=2,
               nrow=1)

#Save the figures
ggsave(FIG,
       width = 5.75, height = 2.5, units = "in",
       dpi = 600,
       bg = "transparent",
       file="NMC_Speeds.png")




