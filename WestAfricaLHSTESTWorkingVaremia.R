
#World Bank Project code for Non-Commercial Sex
rm(list=ls())
library(FME)
#library(ggplot2)
pars=list(T_f=500000, T_m=500000, beta_FM=0.001, beta_MF=0.002, theta=25, circum=0.42, 
          f_FGP=0.02, f_FBB=0.64, f_FNB=0.56,f_FTS=0.2,f_mf=0.31, f_mc=0.16, gamma_V=2, gamma=0.125,
          mu = 0.03,p_F=0.03, p_M=0.03, p_F1=0.03, p_M1=0.03,P_CBB=0.005,P_CNB=0.005,
          P_CTS=0.01,P_FBB=0.005,P_FNB=0.005,P_FTS=0.03, pmf=0.23,pmc=0.38,df_BB=2,df_NB=2, df_TS=2, 
          dm_BB=2,dm_NB=2,dm_TS=2,Z=0.6,CFBB_CBB1=48,CFNB_CNB1=48,CFTS_CTS1=2,CCBB_FBB=295,CCNB_FNB=73, CMGP_FGP=0.83,CMGP_FGP1=0.83,
          CCTS_FTS=6,CCBB_FTS=6,CCNB_FTS=6,eta_T=0.16, t_gamma=0.031)
FSW.ode<- function(pars,times=seq(0,60,by=20))
{derivs <- function(t,state,pars)
  
{ # returns rate of change
  with(as.list(c(state, pars)), 
{
  
 #Treatment rates
   
  if(t<20){
    t_cov =0
  } else if(t>=20&t<=27){
    t_cov=(t-20)*0.02
   }else {
     t_cov=0.14
   }
  
    
  #Duration 
  #t_gamma=1/t_du #treatment rate
  alpha_BB=1/df_BB
  alpha_NB=1/df_NB
  alpha_TS=1/df_TS
  delta_BB=1/dm_BB
  delta_NB=1/dm_NB
  delta_TS=1/dm_TS
  
  #Population sizes
  #initial sizes
  N0_FGP=T_f-(P_FBB*T_f+P_FNB*T_f+P_FTS*T_f) 
  N0_MGP=abs(T_m-(CCBB_FBB*P_FBB*T_f/CFBB_CBB1+CCNB_FNB*P_FNB*T_f/CFNB_CNB1+(1-Z)*CCTS_FTS*P_FTS*T_f/CFTS_CTS1))
  N0_MGPint=T_m-(CCBB_FBB*P_FBB*T_f/CFBB_CBB1+CCNB_FNB*P_FNB*T_f/CFNB_CNB1+(1-Z)*CCTS_FTS*P_FTS*T_f/CFTS_CTS1)
  N0_MGP=(T_f-(P_FBB*T_f+P_FNB*T_f+P_FTS*T_f))
  N0_CBB=CCBB_FBB*P_FBB*T_f/CFBB_CBB1
  N0_CNB=CCNB_FNB*P_FNB*T_f/CFNB_CNB1
  N0_CTS=(1-Z)*CCTS_FTS*P_FTS*T_f/CFTS_CTS1
  N0_FBB=P_FBB*T_f
  N0_FNB=P_FNB*T_f
  N0_FTS=P_FTS*T_f
  
  #variable population sizes
  N_FGP=S_FGP+I_FGV+I_FGP
  N_MGP=S_MGP+I_MGV+I_MGP
  
  N_FBB=S_FBB+I_FBV+I_FBB
  N_FNB=S_FNB+I_FNV+I_FNB
  N_FTS=S_FTS+I_FTV+I_FTS
  
  N_CBB=S_CBB+I_CBV+I_CBB
  N_CNB=S_CNB+I_CNV+I_CNB
  N_CTS=S_CTS+I_CTV+I_CTS
  
  N_FSW=N_FBB+N_FNB+N_FTS
  N_C=N_CBB+N_CNB+N_CTS
  N_F=N_FGP+N_FSW
  N_M=N_MGP+N_C
  
  #PI
  
    
  PI_CNB=N_CNB/(N_CBB+N_CNB)
  PI_CBB=N_CBB/(N_CBB+N_CNB)
  #PI_CTS=N_CTS/(N_CTS+N_CBB+N_CNB)
  

   
  P_MGPFGP=N_FGP/(pmf*N_FBB+pmf*N_FNB+N_FGP)
  P_MGPFBB=pmf*N_FBB/(pmf*N_FBB+pmf*N_FNB+N_FGP)
  P_MGPFNB=pmf*N_FNB/(pmf*N_FBB+pmf*N_FNB+N_FGP)
  P_FGPMGP=N_MGP/(pmc*N_CBB+pmc*N_CNB+N_MGP)
  P_FGPCBB=pmc*N_CBB/(pmc*N_CBB+pmc*N_CNB+N_MGP)
  P_FGPCNB=pmc*N_CNB/(pmc*N_CBB+pmc*N_CNB+N_MGP)
  
  #contact constraints
  CFBB_CBB=CCBB_FBB*N_FBB/N_CBB
  CFNB_CNB=CCNB_FNB*N_FNB/N_CNB
  CFTS_CTS=(1-Z)*(CCTS_FTS*N_FTS)/N_CTS
  CFTS_CBB=PI_CBB*Z*CCTS_FTS*N_FTS/N_CBB
  CFTS_CNB=PI_CNB*Z*CCTS_FTS*N_FTS/N_CNB
  CFGP_MGP=P_FGPMGP*(CMGP_FGP*N_FGP)/(P_MGPFGP*N_MGP)
  CMGP_FBB=P_MGPFBB*(CFGP_MGP*N_MGP)/(pmf*N_FBB)
  CMGP_FNB=P_MGPFNB*(CFGP_MGP*N_MGP)/(pmf*N_FNB)
  CFGP_CBB=P_FGPCBB*(CMGP_FGP*N_FGP)/(pmc*N_CBB)
  CFGP_CNB=P_FGPCNB*(CMGP_FGP*N_FGP)/(pmc*N_CNB)
  
  #General population
  lambda_MGP=beta_FM*CFGP_MGP*52*(1-circum)*((1-f_mf)*P_MGPFBB*(theta*I_FBV+((1-t_cov)+t_cov*eta_T)*I_FBB)/N_FBB+(1-f_mf)*P_MGPFNB*(theta*I_FNV+((1-t_cov)+t_cov*eta_T)*I_FNB)/N_FNB+(1-f_FGP)*P_MGPFGP*(theta*I_FGV+((1-t_cov)+t_cov*eta_T)*I_FGP)/N_FGP)
  lambda_FGP=beta_MF*CMGP_FGP*52*((1-f_mc)*P_FGPCBB*(theta*I_CBV+((1-t_cov)+t_cov*eta_T)*I_CBB)/N_CBB+(1-f_mc)*P_FGPCNB*(theta*I_CNV+((1-t_cov)+t_cov*eta_T)*I_CNB)/N_CNB+(1-f_FGP)*P_FGPMGP*(theta*I_MGV+((1-t_cov)+t_cov*eta_T)*I_MGP)/N_MGP)
 
    
  #clients
  lambda_CBB=beta_FM*(1-circum)*(CFBB_CBB*2*(1-f_FBB)*(theta*I_FBV+((1-t_cov)+t_cov*eta_T)*I_FBB)/N_FBB + CFTS_CBB*20*(1-f_FTS)*(theta*I_FTV+((1-t_cov)+t_cov*eta_T)*I_FTS)/N_FTS+CFGP_CBB*52*(1-f_FGP)*(theta*I_FGV+((1-t_cov)+t_cov*eta_T)*I_FGP)/N_FGP)
  lambda_CNB=beta_FM*(1-circum)*(CFNB_CNB*4*(1-f_FNB)*(theta*I_FNV+((1-t_cov)+t_cov*eta_T)*I_FNB)/N_FNB + CFTS_CNB*20*(1-f_FTS)*(theta*I_FTV+((1-t_cov)+t_cov*eta_T)*I_FTS)/N_FTS+CFGP_CNB*52*(1-f_FGP)*(theta*I_FGV+((1-t_cov)+t_cov*eta_T)*I_FGP)/N_FGP)
  lambda_CTS=beta_FM*(1-circum)*CFTS_CTS*20*(1-f_FTS)*(theta*I_FTV+((1-t_cov)+t_cov*eta_T)*I_FTS)/N_FTS
  
  #sex workers
  lambda_FBB=beta_MF*(CCBB_FBB*2*(1-f_FBB)*(theta*I_CBV+((1-t_cov)+t_cov*eta_T)*I_CBB)/N_CBB+CMGP_FBB*52*(1-f_FGP)*(theta*I_MGV+((1-t_cov)+t_cov*eta_T)*I_MGP)/N_MGP)
  lambda_FNB=beta_MF*(CCNB_FNB*4*(1-f_FNB)*(theta*I_CNV+((1-t_cov)+t_cov*eta_T)*I_CNB)/N_CNB+CMGP_FNB*52*(1-f_FGP)*(theta*I_MGV+((1-t_cov)+t_cov*eta_T)*I_MGP)/N_MGP)
  lambda_FTS=beta_MF*CCTS_FTS*20*(1-f_FTS)*(PI_CBB*Z*(theta*I_CBV+((1-t_cov)+t_cov*eta_T)*I_CBB)/N_CBB + PI_CNB*Z*(theta*I_CNV+((1-t_cov)+t_cov*eta_T)*I_CNB)/N_CNB + (1-Z)*P_FTS*(theta*I_CTV+((1-t_cov)+t_cov*eta_T)*I_CTS)/N_CTS)
    
  #General population
  
  
  #Females
  dS_FGP = p_F*N0_FGP+ gamma*I_FGP+alpha_BB*S_FBB+alpha_NB*S_FNB+alpha_TS*S_FTS-lambda_FGP*S_FGP-mu*S_FGP-(alpha_BB*(S_FBB+I_FBB)+alpha_NB*(S_FNB+I_FNB)+alpha_TS*(S_FTS+I_FTS))*(S_FGP/(S_FGP+I_FGP))
  dI_FGV=lambda_FGP*S_FGP-(gamma_V+mu)*I_FGV
  dI_FGP=gamma_V*I_FGV+alpha_BB*I_FBB+alpha_NB*I_FNB+alpha_TS*I_FTS-(mu+(1-t_cov)*gamma+t_cov*t_gamma)*I_FGP-(alpha_BB*(S_FBB+I_FBB)+alpha_NB*(S_FNB+I_FNB)+alpha_TS*(S_FTS+I_FTS))*(I_FGP/(S_FGP+I_FGP))
  
  #Males
  dS_MGP = p_M*N0_MGP+gamma*I_MGP+delta_BB*S_CBB+delta_NB*S_CNB+delta_TS*S_CTS-lambda_MGP*S_MGP-mu*S_MGP-(delta_BB*(S_CBB+I_CBB)+delta_NB*(S_CNB+I_CNB)+delta_TS*(S_CTS+I_CTS))*(S_MGP/(S_MGP+I_MGP))
  dI_MGV = lambda_MGP*S_MGP-(gamma_V+mu)*I_MGV
  dI_MGP= gamma_V*I_MGV+delta_BB*I_CBB+delta_NB*I_CNB+delta_TS*I_CTS-(mu+(1-t_cov)*gamma+t_cov*t_gamma)*I_MGP-(delta_BB*(S_CBB+I_CBB)+delta_NB*(S_CNB+I_CNB)+delta_TS*(S_CTS+I_CTS))*(I_MGP/(S_MGP+I_MGP))
  
  #Clients
  
  dS_CBB=p_M1*N0_CBB+ gamma*I_CBB+delta_BB*(S_CBB+I_CBB)*(S_MGP/(S_MGP+I_MGP))-lambda_CBB*S_CBB-(delta_BB+mu)*S_CBB
  dI_CBV=lambda_CBB*S_CBB-(gamma_V+mu)*I_CBV
  dI_CBB=gamma_V*I_CBV+delta_BB*(S_CBB+I_CBB)*(I_MGP/(S_MGP+I_MGP))-(mu+(1-t_cov)*gamma+t_cov*t_gamma+delta_BB)*I_CBB
  
  dS_CNB=p_M1*N0_CNB+ gamma*I_CNB+delta_NB*(S_CNB+I_CNB)*(S_MGP/(S_MGP+I_MGP))-lambda_CNB*S_CNB-(delta_NB+mu)*S_CNB
  dI_CNV=lambda_CNB*S_CNB-(gamma_V+mu)*I_CNV
  dI_CNB=gamma_V*I_CNV+delta_NB*(S_CNB+I_CNB)*(I_MGP/(S_MGP+I_MGP))-(mu+(1-t_cov)*gamma+t_cov*t_gamma+delta_NB)*I_CNB
  
  dS_CTS=p_M1*N0_CTS+gamma*I_CTS+delta_TS*(S_CTS+I_CTS)*(S_MGP/(S_MGP+I_MGP))-lambda_CTS*S_CTS-(delta_TS+mu)*S_CTS
  dI_CTV=lambda_CTS*S_CTS-(gamma_V+mu)*I_CTV
  dI_CTS=gamma_V*I_CTV+delta_TS*(S_CTS+I_CTS)*(I_MGP/(S_MGP+I_MGP))-(mu+(1-t_cov)*gamma+t_cov*t_gamma+delta_TS)*I_CTS
  
  #FSWs
  dS_FBB=p_F1*N0_FBB+ gamma*I_FBB+alpha_BB*(S_FBB+I_FBB)*(S_FGP/(S_FGP+I_FGP))-lambda_FBB*S_FBB-(alpha_BB+mu)*S_FBB
  dI_FBV=lambda_FBB*S_FBB-(gamma_V+mu)*I_FBV
  dI_FBB=gamma_V*I_FBV+alpha_BB*(S_FBB+I_FBB)*(I_FGP/(S_FGP+I_FGP))-(mu+(1-t_cov)*gamma+t_cov*t_gamma+alpha_BB)*I_FBB
  
  dS_FNB=p_F1*N0_FNB+ gamma*I_FNB+alpha_NB*(S_FNB+I_FNB)*(S_FGP/(S_FGP+I_FGP))-lambda_FNB*S_FNB-(alpha_NB+mu)*S_FNB
  dI_FNV=lambda_FNB*S_FNB-(gamma_V+mu)*I_FNV
  dI_FNB=gamma_V*I_FNV+alpha_NB*(S_FNB+I_FNB)*(I_FGP/(S_FGP+I_FGP))-(mu+(1-t_cov)*gamma+t_cov*t_gamma+alpha_NB)*I_FNB
  
  dS_FTS=p_F1*N0_FTS+gamma*I_FTS+alpha_TS*(S_FTS+I_FTS)-lambda_FTS*S_FTS-(alpha_TS+mu)*S_FTS
  dI_FTV=lambda_FTS*S_FTS-(gamma_V+mu)*I_FTV
  dI_FTS=gamma_V*I_FTV-(mu+(1-t_cov)*gamma+t_cov*t_gamma+alpha_TS)*I_FTS
  
  return(list(c(dS_FGP,dI_FGV, dI_FGP,dS_MGP,dI_MGV,dI_MGP,dS_CBB,dI_CBV,dI_CBB,dS_CNB,dI_CNV,dI_CNB,dS_CTS,dI_CTV,dI_CTS,dS_FBB,dI_FBV,dI_FBB,dS_FNB,dI_FNV, dI_FNB,dS_FTS,dI_FTV,dI_FTS), 
              N_F=N_FGP+N_FSW,N_M=N_MGP+N_C,
              N0_MGPint=T_m-(CCBB_FBB*P_FBB*T_f/CFBB_CBB1+CCNB_FNB*P_FNB*T_f/CFNB_CNB1+(1-Z)*CCTS_FTS*P_FTS*T_f/CFTS_CTS1),
              N0_CBB=CCBB_FBB*P_FBB*T_f/CFBB_CBB1,
              N0_CNB=CCNB_FNB*P_FNB*T_f/CFNB_CNB1,
              N0_CTS=(1-Z)*CCTS_FTS*P_FTS*T_f/CFTS_CTS1,
              Prev_F=((I_FGP+I_FGV+I_FBB+I_FBV+I_FNB+I_FNV+I_FTS+I_FTV)/N_F)*100,
              Prev_M=((I_MGP+I_MGV+I_CBB+I_CBV+I_CNB+I_CNV+I_CTS+I_CTV)/N_M)*100,
              Prev_CBB=((I_CBB+I_CBV)/N_CBB)*100,Prev_CNB=((I_CNB+I_CNV)/N_CNB)*100,
              Prev_CTS=((I_CTS+I_CTV)/N_CTS)*100,Prev_FBB=((I_FBB+I_FBV)/N_FBB)*100,
              Prev_FNB=((I_FNB+I_FNV)/N_FNB)*100,Prev_FTS=((I_FTS+I_FTV)/N_FTS)*100,
              CFBB_CBB=CCBB_FBB*N_FBB/N_CBB,
              CFNB_CNB=CCNB_FNB*N_FNB/N_CNB,
              CFTS_CTS=(1-Z)*(CCTS_FTS*N_FTS)/N_CTS, 
              CFTS_CBB=PI_CBB*Z*CCTS_FTS*N_FTS/N_CBB,
              CFTS_CNB=PI_CNB*Z*CCTS_FTS*N_FTS/N_CNB,
              CFGP_MGP=P_FGPMGP*(CMGP_FGP*N_FGP)/(N_MGP),
              CMGP_FBB=P_MGPFBB*(CFGP_MGP*N_MGP)/(pmf*N_FBB),
              CMGP_FNB=P_MGPFNB*(CFGP_MGP*N_MGP)/(pmf*N_FNB),
              CFGP_CBB=P_FGPCBB*(CMGP_FGP*N_FGP)/(pmc*N_CBB),
              CFGP_CNB=P_FGPCNB*(CMGP_FGP*N_FGP)/(pmc*N_CNB),
              N_CBB=S_CBB+I_CBV+I_CBB,
              N_CNB=S_CNB+I_CNV+I_CNB,
              PI_CNB=N_CNB/(N_CBB+N_CNB), 
              PI_CBB=N_CBB/(N_CBB+N_CNB),
              P_MGPFGP=N_FGP/(pmf*N_FBB+pmf*N_FNB+N_FGP),
              P_MGPFBB=pmf*N_FBB/(pmf*N_FBB+pmf*N_FNB+N_FGP),
              P_MGPFNB=pmf*N_FNB/(pmf*N_FBB+pmf*N_FNB+N_FGP),
              P_FGPMGP=N_MGP/(pmc*N_CBB+pmc*N_CNB+N_MGP),
              P_FGPCBB=pmc*N_CBB/(pmc*N_CBB+pmc*N_CNB+N_MGP),
              P_FGPCNB=pmc*N_CNB/(pmc*N_CBB+pmc*N_CNB+N_MGP),t_cov,
              N_FBB=S_FBB+I_FBV+I_FBB,
              N_FNB=S_FNB+I_FNV+I_FNB,
              N_FTS=S_FTS+I_FTV+I_FTS,              
              N_CTS=S_CTS+I_CTV+I_CTS))})
}
 #initial conditions
  state<-c(S_FGP=pars$T_f-(pars$P_FBB*pars$T_f+pars$P_FNB*pars$T_f+pars$P_FTS*pars$T_f), I_FGV=0,I_FGP=0,S_MGP=abs(pars$T_m-(pars$CCBB_FBB*pars$P_FBB*pars$T_f/pars$CFBB_CBB1+pars$CCNB_FNB*pars$P_FNB*pars$T_f/pars$CFNB_CNB1+(1-pars$Z)*pars$CCTS_FTS*pars$P_FTS*pars$T_f/pars$CFTS_CTS1)),
          I_MGV=0,I_MGP=0,S_CBB=pars$CCBB_FBB*pars$P_FBB*pars$T_f/pars$CFBB_CBB1,I_CBV=0,I_CBB=0,S_CNB=pars$CCNB_FNB*pars$P_FNB*pars$T_f/pars$CFNB_CNB1,
          I_CNV=0,I_CNB=0,S_CTS=(1-pars$Z)*pars$CCTS_FTS*pars$P_FTS*pars$T_f/pars$CFTS_CTS1,I_CTV=0,I_CTS=0,S_FBB=pars$P_FBB*pars$T_f,I_FBV=0,I_FBB=1,S_FNB=pars$P_FNB*pars$T_f,I_FNV=0,I_FNB=0,S_FTS=pars$P_FTS*pars$T_f,I_FTV=0,I_FTS=0) 
# ode solves the model by integration
 return(as.data.frame(ode(y = state,times = times, func = derivs, parms = pars)))
}

out=FSW.ode(pars)

#Time series plots
R0<-cbind(out$N_F,out$N_M)
R1<-cbind(out$S_FGP,out$I_FGP,out$S_MGP,out$I_MGP)
R2<-cbind(out$S_CBB,out$I_CBB,out$S_CNB,out$I_CNB,out$S_CTS,out$I_CTS)
R3<-cbind(out$S_FBB,out$I_FBB,out$S_FNB,out$I_FNB,out$S_FTS,out$I_FTS)


#Prevalences
R4<-cbind(out$Prev_F,out$Prev_M)
R5<-cbind(out$Prev_CBB,out$Prev_CNB,out$Prev_CTS)
R6<-cbind(out$Prev_FBB,out$Prev_FNB,out$Prev_FTS)

#High Viraemia
R7<-cbind(out$S_FGP,out$I_FGV,out$S_MGP,out$I_MGV)
R8<-cbind(out$S_CBB,out$I_CBV,out$S_CNB,out$I_CNV,out$S_CTS,out$I_CTV)
R9<-cbind(out$S_FBB,out$I_FBV,out$S_FNB,out$I_FNV,out$S_FTS,out$I_FTV)


t<-out[,1]

par(mfrow=c(1,1))
matplot(t,R0,type="l",lwd = c(2,2), col =c("blue","red"), xlab = "Time (years)", ylab = "Population size")
legend("topright", c("N_F","N_M"),lty = 1:2,bty="n",col =c("blue","red"), lwd = c(2,2),cex = 0.6)


par(mfrow=c(1,1))
matplot(t,R1,type="l",lwd = c(2,2,2,2), col =c("blue","red","black","brown"), xlab = "Time (years)", ylab = "Population size")
legend("center", c("S_FGP","I_FGP", "S_MGP","I_MGP"),lty = 1:4,bty="n",col =c("blue","red","black","yellow"), lwd = c(2,2,2,2),cex = 0.6)

par(mfrow=c(1,2))
matplot(t,R2,type="l",lwd = c(2,2,2,2,2,2), col =c("forestgreen","green","blue","red","black","brown"), xlab = "Time (years)", ylab = "Populations size")
legend("topright", c("S_CBB","I_CBB","S_CNB", "I_CNB","S_CTS", "I_CTS"),lty = 1:6,bty="n",col=c("forestgreen","green","blue","red","black","brown"), lwd =c(2,2,2,2,2,2),cex = 0.6)

matplot(t,R3,type="l",lwd = c(2,2,2,2,2,2), col =c("forestgreen","green","blue","red","black", "brown"), xlab = "Time (years)", ylab = "Population size")
legend("topright", c("S_FBB","I_FBB","S_FNB", "I_FNB","S_FTS", "I_FTS"),lty = 1:6,bty="n",col=c("forestgreen","green","blue","red","black","brown"), lwd = c(2,2,2,2,2,2),cex = 0.6)
#Prevtes
par(mfrow=c(1,1))
matplot(t,R4,type="l",lwd = c(2,2), col =c("blue","red"), xlab = "Time (years)", ylab = "Prevalence")
legend("topright", c("Prev_F","Prev_M"),lty = 1:2,bty="n",col =c("blue","red"), lwd = c(2,2),cex = 0.6)

par(mfrow=c(1,2))
matplot(t,R5, type="l",lwd = c(2,2,2), col =c("forestgreen","green","blue"), xlab = "Time (years)", ylab = "Prevalence")
legend("topright", c("Prev_CBB","Prev_CNB","Prev_CTS"),lty = 1:3,bty="n",col=c("forestgreen","green","blue"), lwd =c(2,2,2),cex = 0.6)

matplot(t,R6,type="l",lwd = c(2,2,2), col =c("forestgreen","green","blue"), xlab = "Time (years)", ylab = "Prevalence")
legend("topright", c("Prev_FBB","Prev_FNB","Prev_FTS"),lty = 1:3,bty="n", col=c("forestgreen","green","blue"), lwd = c(2,2,2),cex = 0.6)
#write.table(out,"out.txt",sep="\t")

# parameter ranges
#v     "beta_FM","beta_MF","P_FBB", "P_FNB",  "P_FTS",  "pmf"    "pmc"  "f_mf"  "f_mc"   "df_BB", "df_NB",  "df_TS",  "dm_BB",    "dm_NB", "dm_TS"  "Z" , CCBB_FBB,CCNB_FNB, "CFBB_CBB1" "CFNB_CNB1" "CCTS_FTS" "f_FBB" "f_FNB" "f_FTS", t_gamma,  eta_T)
# min<-c(0.0006,   0.0006,   0.001,    0.0015,     0.001,  0.12,    0.28,  0.26,   0.8  ,    1,       1,       5,        5,         5,       1 ,     0,      252   ,   42,       12,         12 ,         3,      0.41,    0.41,   0.10,   0.023,      0)
# max<-c(0.001,    0.002,    0.011,    0.03,       0.20,   0.49,    0.58,  0.47,   0.32 ,    6,       6,       15,       10,        10,      10 ,    1,      1092  ,  336,       48,         48 ,         24,     0.83,    0.68,   0.39,    0.04,    0.16)
#v     "beta_FM","beta_MF","P_FBB", "P_FNB",  "P_FTS",  "pmf"    "pmc"  "f_mf"  "f_mc"   "df_BB", "df_NB",  "df_TS",  "dm_BB",    "dm_NB", "dm_TS"  "Z",CCBB_FBB,CCNB_FNB, "CFBB_CBB1" "CFNB_CNB1" "CCTS_FTS" "f_FBB" "f_FNB" "f_FTS", t_gamma, eta_T)
min<-c(0.0009,   0.0015,   0.006,    0.01575,  0.001,   0.305,    0.43,  0.365,   0.56  ,    4,      4,        15,       7.5,       7.5,     5.5 ,   0,    798,    189,        30,       30,         13.5,        0.62,   0.545,  0.245,   0.0315, 0.08)
max<-c(0.0009,   0.0015,   0.006,    0.01575,  0.20,    0.305,    0.43,  0.365,   0.56 ,     4,      4,        15,       7.5,       7.5,     5.5 ,   1,    798,    189,        30,       30,         13.5,       0.62,   0.545,  0.245,   0.0315, 0.08)


parRanges<-cbind(min,max)
rownames(parRanges)= c("beta_FM","beta_MF","P_FBB","P_FNB","P_FTS","pmf","pmc","f_mf","f_mc","df_BB","df_NB","df_TS","dm_BB","dm_NB","dm_TS","Z","CCBB_FBB","CCNB_FNB","CFBB_CBB1","CFNB_CNB1","CCTS_FTS","f_FBB","f_FNB","f_FTS","t_gamma","eta_T")
parRanges

# Sensitivity analysis of all parameters (LHS)
#SA0=sensRange(func =FSW.ode, parms = pars, dist = "latin",sensvar = c("N0_MGPint","N0_CBB","N0_CNB","N0_CTS","CFBB_CBB", "CFGP_CBB","CFTS_CBB", "CFNB_CNB","CFGP_CNB","CFTS_CNB","PI_CNB", "PI_CBB","P_MGPFGP", "P_MGPFBB","P_MGPFNB", "P_FGPMGP", "P_FGPCBB", "P_FGPCNB","Prev_F","Prev_M","Prev_CBB","Prev_CNB","Prev_CTS","Prev_FBB","Prev_FNB","Prev_FTS"),parRange = parRanges, num=1000)
SA0=sensRange(func =FSW.ode, parms = pars, dist = "latin",sensvar = c("N0_MGPint","N0_CBB","N0_CNB","N0_CTS","CFTS_CBB","CFTS_CNB","P_FGPCBB", "P_FGPCNB","P_MGPFBB","P_MGPFNB","Prev_F","Prev_M","Prev_CBB","Prev_CNB","Prev_CTS","Prev_FBB","Prev_FNB","Prev_FTS","N_FBB","N_FNB","N_FTS","N_CBB","N_CNB","N_CTS"),parRange = parRanges, num=2000)
write.table(SA0,"SA04.txt",sep="\t")

#selecting output parameters

d = read.table("SA04.txt", sep="\t",header = T)
d1=subset(d, select=c(beta_FM,beta_MF,P_FBB,P_FNB,P_FTS,df_BB,
                      df_NB,df_TS,dm_BB,dm_NB,dm_TS, Z,N0_MGPint60,N0_CBB60,N0_CNB60,N0_CTS60,N_FBB60,N_FNB60,N_FTS60,N_CBB60,N_CNB60,N_CTS60,CCBB_FBB,CCNB_FNB,CFTS_CBB60,CFTS_CNB60, CFBB_CBB1,CFNB_CNB1,CCTS_FTS,P_FGPCBB60,P_FGPCNB60,P_MGPFBB60,P_MGPFNB60,f_FBB,f_FNB,f_FTS,Prev_F60,Prev_M60,Prev_CBB60,Prev_CNB60, Prev_CTS60,
                      Prev_FBB60,Prev_FNB60,Prev_FTS60))
#All scenarios output

# d2=subset(d1,(Prev_F60>=0&Prev_F60<=7) &(Prev_M60>=0 &Prev_M60<=3.5)&(Prev_CBB60>=1 & Prev_CBB60<=7)&(Prev_CNB60>=1&Prev_CNB60<=7)& (Prev_CTS60>=0&Prev_CTS60<=4)&(Prev_FBB60>=15&Prev_FBB60<=35)&
#             (Prev_FNB60>=10&Prev_FNB60<=25)&(Z>=0&Z<=1)&(N0_MGPint60>0))
# write.table(d2,"ScenarioALL.txt",sep="\t")
d2z=subset(d1,(Prev_F60>=0&Prev_F60<=7) &(Prev_M60>=0 &Prev_M60<=4)&(Prev_FBB60>=15&Prev_FBB60<=48)&(Prev_FNB60>=10&Prev_FNB60<=25)&(N0_MGPint60>0))
write.table(d2z,"ScenarioALLrecent.txt",sep="\t")

#Scenario1 output

d3=subset(d1,(Prev_F60>=0&Prev_F60<=2) &(Prev_M60>=0 &Prev_M60<=1.5)&(Prev_FBB60>=15&Prev_FBB60<=48)&(Prev_FNB60>=10&Prev_FNB60<=25)&(N0_MGPint60>0))
write.table(d3,"Scenario1.txt",sep="\t")

#Scenario2 output

d4=subset(d1,(Prev_F60>=4&Prev_F60<=7)&(Prev_M60>=2.5 &Prev_M60<=4)&(Prev_FBB60>=15&Prev_FBB60<=48)&(Prev_FNB60>=10&Prev_FNB60<=25)&(N0_MGPint60>0))
write.table(d4,"Scenario2.txt",sep="\t")

d5=subset(d1,(Prev_F60>=0&Prev_F60<=25)&(Prev_M60>=0 &Prev_M60<=25)&(N0_MGPint60>0))
write.table(d5,"Scenario25prev.txt",sep="\t")
#Scenario3 output
d6=subset(d1,(Prev_F60>=0&Prev_F60<=7)&(Prev_M60>=0 &Prev_M60<=4)&(N0_MGPint60>0))
write.table(d6,"Scenarioprev07.txt",sep="\t")

d7=subset(d1,(Prev_F60>=0&Prev_F60<=2) &(Prev_M60>=0 &Prev_M60<=1.5)&(N0_MGPint60>0)&(Prev_FBB60>=15&Prev_FBB60<=48))
write.table(d7,"Scenario3.txt",sep="\t")

#Scenario4 output

d8=subset(d1,(Prev_F60>=4&Prev_F60<=7)&(Prev_M60>=2.5 &Prev_M60<=4)&(N0_MGPint60>0))
write.table(d8,"Scenario4.txt",sep="\t")


#Scenario5 output

#d7=subset(d1,(Prev_F30>=0&Prev_F30<=3) &(Prev_M30>=0 &Prev_M30<=1)&(Prev_CBB30>=1 & Prev_CBB30<=7)&(Prev_CNB30>=1&Prev_CNB30<=7)&(Prev_FBB30>=15&Prev_FBB30<=35)&(Prev_FNB30>=10&Prev_FNB30<=25)&(Z>=0&Z<=1))
#write.table(d7,"Scenario5.txt",sep="\t")

#Scenario6 output

#d8=subset(d1,(Prev_F30>=3&Prev_F30<=7)&(Prev_M30>=1 &Prev_M30<=3.5)&(Prev_CNB30>=1&Prev_CNB30<=7)&(Prev_FNB30>=10&Prev_FNB30<=25))
#write.table(d8,"Scenario6.txt",sep="\t")


SA1 <- summary(SA0)
par(mfrow=c(2,4))
plot(SA1, xlab = "time (years)",ylab = "HIV prevalence", mfrow = NULL,quant = F, col = c("lightblue","darkblue"), legpos = "topright")
mtext(outer = TRUE, line = -1.5, side = 1, "",cex = 1.25)








