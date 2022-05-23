
# This code is supplement to: J. Schwarzer et al. A theoretical model for CO2 fluxes of accurate incubation measurements, in prep.
# further information: https://ghg-in-permafrost.awi.de/



#setwd()

library(ggplot2) #library for plotting
library(dplyr)  # data munging library
library(tidyr)  # data munging library
library(deSolve) # library for solving differential equations
library(minpack.lm) # library for least squares fit using levenberg-marquart algorithm


# variables names ---------------------------------------------------------------

# SOC = degradable Soil organic Carbon in mol/l
#       for the model the SOC is simplified to organic Carbon in the volume of the aqueous phase
# Caq = CO2 in Water phase of the soil sample in mol/l
# Cgas = CO2 in gas headspace in mol/l
# k1 = rate constant for microbial degradation
# kh1 = rate constant representing CO2 going from water to gas-phase
# KH = Henry Constant
# Vw = volume of water phase in l
# Vg = volume of gas phase in l


# constants -------------------------------------------------------------------
# general gas constant [(kg * m²)/(s² * mol * K)]
R = 8.31446261815324

# Temperatur from C to K
CinK = 273.15
# experiment temperature [K]
Temp = 20 + CinK 

# convert Henry constant  from Tnorm (298.15K) to experiment T (293.15K)
# published value for KHcp at normT (Sander, 2015): 3.4e-2 in mol/l*atm
# Sander, R. Compilation of Henry’s Law Constants (Version 4.0) for Water as Solvent. Atmos. Chem. Phys. 2015, 15 (8), 4399–4981. https://doi.org/10.5194/acp-15-4399-2015.
KHcpT = 3.4e-2
# published conversion factor see Sander (2015): 2400
Cfactor= 2400
KHcp=KHcpT * exp(Cfactor*((1/Temp)-(1/298.15)))

# for the model the inverse KH with units in mol/l is needed:
# Henry Constant in p_g/c_aq [l*atm/mol]
KH = 1/KHcp
# Henry Constant in c_g/c_aq [no unit]
KH =round( (KH * 0.001 * 101325) / (R*Temp), 3)

# Gas volume of measurement set up
Vg = 0.22 # [l]
# Water volume of measurement set up
Vw <- 0.018 # [l]



# get the data --------------------------------------------------------------
df <- read.csv('./day01_01')

# subset df for time vector and CO2 concentration in mol/l
df <- df %>% select(CO2_moll, DiffTime)
colnames(df) <- c("CO2_moll", "time")

#gather concentration into key-value pairs
dfg = df %>% gather(species, conc, -time)


# have a look at the concentration development measured:
ggplot(data=dfg, aes(x=time, y=conc, color=species)) +
  geom_point(size=1) +
  labs(x = "time [s]", y= bquote(~CO[2]~"[mol/l]"))



# kinetic model --------------------------------------------------------------------

# rxnrate- function with a set of three ODEs that describe the system:

# c - concentration of species SOC, Caq and Cgas
# t - time
# parms- list of rate constants k1 and kh1 

rxnrate=function(t,c,parms){
  
  k1=parms$k1
  kh1=parms$kh1
  
  # derivatives dc/dt are computed below
  r=rep(0,length(c))
  r[1]=-k1*c["SOC"] #dcA/dt
  r[2]= k1*c["SOC"]-kh1*c["Caq"]+kh1/KH*c["Cgas"] #dcB/dt
  r[3]= kh1*c["Caq"]*(Vw/Vg)-kh1/KH*c["Cgas"]*(Vw/Vg)#dcC/dt
  
  return(list(r))
  
}



# run the model with some made up values for initial concentrations and rate constants:
# start values for initial concentrations
cinit <- c(SOC=1,Caq=0.75e-06,Cgas=5e-05)
parms <- list(k1=1.5e-09,kh1=0.04)
t <- df$time

# perform ode
out <- ode(y=cinit, times=t, func=rxnrate, parms=parms)
head(out)

# plot of model and measured CO2 concentrations
outdf=data.frame(out[,c(1,3,4)])

#plot for modeled vs. measured
ggplot()+
  geom_line(data=df, aes(x=time, y=CO2_moll, color="measured"),size=1)+
  geom_line(data=outdf, aes(x=time, y = Cgas, color="modelled"))+
  scale_color_manual(name = "CO2 in gasphase", values = c("measured" = "#F8766D", "modelled" = "#00BFC4"))




# parameter fitting----------------------------------------------------------

# ssq - function to calculate residuals of model vs. measured data

# fitparams - parameters to be estimated: k1, kh1 and the 3 initial concentrations of SOC, Caq, and Cgas
# df - dataframe containing experimental data

ssq=function(fitparms, df){
  
  # inital concentration
  cinit = fitparms[3:5]
  names(cinit) = c("SOC", "Caq", "Cgas")
  
  # time points for which conc. is reported
  t=c(seq(0,length(df$time),1),df$time)
  t=sort(unique(t))
  
  # parameters from the parameter estimation routine
  k1=fitparms[1]
  kh1=fitparms[2]
  
  # solve ODE for a given set of parameters
  out=ode(y=cinit, times=t, func=rxnrate, parms=list(k1=k1,kh1=kh1))
  
  # Filter data that contains time points where data is available
  outdf=data.frame(out[,c(1,4)])
  outdf=outdf[outdf$time %in% df$time,]
  
  # Evaluate predicted vs experimental residual
  preddf = outdf %>% gather(species, conc, -time)
  expdf = df %>% gather(species, conc, -time)
  
  ssqres=preddf$conc-expdf$conc
  
  # return predicted vs experimental residual
  return(ssqres)
  
}




# initial guess for parameters
fitparms=c(k1=1.5e-09,kh1=0.04,
           SOC=1,Caq=0.75e-06,Cgas=mean(df$CO2_moll[1:3]))


# fitting using levenberg marquart algorithm
# lower bound of 0 imposed for fitting parameters since rate constants and concentrations are non-negative
fitval=nls.lm(par=fitparms, fn=ssq, lower = rep(0, 5),
              control = nls.lm.control(nprint = 1), df = df)

# fittet parameters:
fitval$par

# extract fitet parameters to perform ODE
parms <- list(k1=as.numeric(fitval$par[1]),
              kh1=as.numeric(fitval$par[2]))
cinit <- c(SOC=as.numeric(fitval$par[3]),
           Caq=as.numeric(fitval$par[4]),
           Cgas=as.numeric(fitval$par[5]))

# perform ode with fittet params
out <- ode(y=cinit, times=t, func=rxnrate, parms=parms)
outdf=data.frame(out[,c(1,3,4)])

# plot data and nls.lm fittet model
ggplot()+
  geom_point(data=df, aes(x=time, y=CO2_moll, color="measured"),size=1)+
  geom_line(data=outdf, aes(x=time, y = Cgas, color="modelled"))+
  scale_color_manual(name = bquote(~CO[2]~"in gasphase"), values = c("measured" = "#F8766D", "modelled" = "#00BFC4"))+
  labs(x = "time [s]", y= bquote(~CO[2]~"[mol/l]"))

# linear model -------------------

# measurement length 
l <- length(df$time)
# a third of ther measurement length
thrd <- l/3
# linear model for last third of measurement
lin <- lm(CO2_moll~time, data = df[(l-thrd):l,])
# save coefficients into list of fits
lin$coefficients[1] # intercept
lin$coefficients[2] # slope

# add correction factor from analytic solution
corr <- 1+ (Vw/(Vg*KH))
lin$coefficients[2] * corr


                           
# plot linear fit of the last third         
ggplot()+
  geom_point(data=df, aes(x=time, y=CO2_moll),color="#CC7987", size=1.5)+
  geom_abline(intercept = lin$coefficients[1], slope = lin$coefficients[2],
              color="#A8DEC8", lwd=1.2)+
  geom_abline(intercept = lin$coefficients[1], slope = lin$coefficients[2]*corr, 
              color="#425BCE", lwd=1.2, linetype="dashed")+
  labs(x = "time [s]", y= bquote(~CO[2]~"[mol/l]"))

ggplot()+
  geom_point(data=df, aes(x=time, y=CO2_moll,color="measured"), size=1.5)+
  geom_abline(aes(color="linear regression", 
                  intercept = lin$coefficients[1], slope = lin$coefficients[2]),
              size=1.5)+
  geom_abline(aes(color="corrected linear regression",
                  intercept = lin$coefficients[1], slope = lin$coefficients[2]*corr),
              size=1.5, linetype="dashed")+
  labs(x = "time [s]", y= bquote(~CO[2]~"[mol/l]"))+
  scale_color_manual(name='',
                     breaks=c('measured', 'linear regression', 'corrected linear regression'),
                     values=c('measured'='#CC7987', 'linear regression'='#A8DEC8', 'corrected linear regression'='#425BCE'))



