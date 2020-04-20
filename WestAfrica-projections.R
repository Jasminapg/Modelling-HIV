
#World Bank Project code for Non-Commercial Sex
rm(list=ls())
library(FME)

derivs <- function(t,state,pars)
  
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
  lambda_FTS=beta_MF*CCTS_FTS*20*(1-f_FTS)*(PI_CBB*Z*(theta*I_CBV+((1-t_cov)+t_cov*eta_T)*I_CBB)/N_CBB + PI_CNB*Z*(theta*I_CNV+((1-t_cov)+t_cov*eta_T)*I_CNB)/N_CNB + (1-Z)*(theta*I_CTV+((1-t_cov)+t_cov*eta_T)*I_CTS)/N_CTS)
  
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
              N_CTS=S_CTS+I_CTV+I_CTS))})}

d = read.table("parameters.txt", sep="\t",header = T)
pr1=subset(d,(P_FTS>=0&P_FTS<0.07))
#pr1=subset(d,(P_FTS>=0.07&P_FTS<=0.20))

  for(i in 1:5){
    beta_FM=pr1[i,1]
    beta_MF=pr1[i,2]
    P_FBB =pr1[i,3] 
    P_FNB=pr1[i,4]  
    P_FTS=pr1[i,5]  
    pmf=pr1[i,6]    
    pmc=pr1[i,7]    
    f_mf=pr1[i,8]  
    f_mc=pr1[i,9]   
    df_BB=pr1[i,10]  
    df_NB=pr1[i,11]  
    df_TS=pr1[i,12]  
    dm_BB=pr1[i,13]  
    dm_NB=pr1[i,14]  
    dm_TS=pr1[i,15]  
    Z=pr1[i,16]      
    CCBB_FBB=pr1[i,17]
    CCNB_FNB=pr1[i,18]
    CFBB_CBB1=pr1[i,19]
    CFNB_CNB1=pr1[i,20]
    CCTS_FTS=pr1[i,21] 
    f_FBB =pr1[i,22]   
    f_FNB=pr1[i,23]    
    f_FTS=pr1[i,24]    
    t_gamma=pr1[i,25]  
    eta_T=pr1[i,26]    
    
    
    t=seq(1975,2080,by=1)
    
    pars=list(T_f=500000, T_m=500000, beta_FM, beta_MF, theta=25, circum=0.42,f_FGP=0.02, f_FBB, f_FNB,f_FTS,f_mf,
              f_mc, gamma_V=2, gamma=0.125,mu=0.03,p_F=0.03, p_M=0.03, p_F1=0.03, p_M1=0.03,P_FBB=P_FBB,P_FNB=P_FNB,
              P_FTS=P_FTS, pmf,pmc, df_BB,df_NB, df_TS, dm_BB,dm_NB,dm_TS,Z=Z,CFBB_CBB1=CFBB_CBB1,CFNB_CNB1=CFNB_CNB1,CFTS_CTS1=2,
              CCBB_FBB=CCBB_FBB,CCNB_FNB=CCNB_FNB, CMGP_FGP=0.83,CMGP_FGP1=0.83,CCTS_FTS=CCTS_FTS,eta_T, t_gamma)
    
    pars2=list(T_f=500000, T_m=500000, beta_FM, beta_MF, theta=25, circum=0.42,f_FGP=0.02, f_FBB, f_FNB,f_FTS,f_mf,
              f_mc, gamma_V=2, gamma=0.125,mu=0.03,p_F=0.03, p_M=0.03, p_F1=0.03, p_M1=0.03,P_FBB=P_FBB,P_FNB=P_FNB,
              P_FTS=P_FTS, pmf,pmc, df_BB,df_NB, df_TS, dm_BB,dm_NB,dm_TS,Z=0,CFBB_CBB1=CFBB_CBB1,CFNB_CNB1=CFNB_CNB1,CFTS_CTS1=2,
              CCBB_FBB=CCBB_FBB,CCNB_FNB=CCNB_FNB, CMGP_FGP=0.83,CMGP_FGP1=0.83,CCTS_FTS=CCTS_FTS,eta_T, t_gamma)
    
    pars3=list(T_f=500000, T_m=500000, beta_FM, beta_MF, theta=25, circum=0.0,f_FGP=0.02, f_FBB, f_FNB,f_FTS,f_mf,
               f_mc, gamma_V=2, gamma=0.125,mu=0.03,p_F=0.03, p_M=0.03, p_F1=0.03, p_M1=0.03,P_FBB=P_FBB,P_FNB=P_FNB,
               P_FTS=P_FTS, pmf,pmc, df_BB,df_NB, df_TS, dm_BB,dm_NB,dm_TS,Z=Z,CFBB_CBB1=CFBB_CBB1,CFNB_CNB1=CFNB_CNB1,CFTS_CTS1=2,
               CCBB_FBB=CCBB_FBB,CCNB_FNB=CCNB_FNB, CMGP_FGP=0.83,CMGP_FGP1=0.83,CCTS_FTS=CCTS_FTS,eta_T, t_gamma)
    
    
    #initial conditions
    state<-c(S_FGP=pars$T_f-(pars$P_FBB*pars$T_f+pars$P_FNB*pars$T_f+pars$P_FTS*pars$T_f),I_FGV=0,I_FGP=0,
             S_MGP=abs(pars$T_m-(pars$CCBB_FBB*pars$P_FBB*pars$T_f/pars$CFBB_CBB1+pars$CCNB_FNB*pars$P_FNB*pars$T_f/pars$CFNB_CNB1+(1-pars$Z)*pars$CCTS_FTS*pars$P_FTS*pars$T_f/pars$CFTS_CTS1)),
             I_MGV=0,I_MGP=0,S_CBB=pars$CCBB_FBB*pars$P_FBB*pars$T_f/pars$CFBB_CBB1,I_CBV=0,I_CBB=0,S_CNB=pars$CCNB_FNB*pars$P_FNB*pars$T_f/pars$CFNB_CNB1,I_CNV=0,I_CNB=0,
             S_CTS=(1-pars$Z)*pars$CCTS_FTS*pars$P_FTS*pars$T_f/pars$CFTS_CTS1,I_CTV=0,I_CTS=0,S_FBB=pars$P_FBB*pars$T_f,I_FBV=0,I_FBB=20,S_FNB=pars$P_FNB*pars$T_f,I_FNV=0,
             I_FNB=20,S_FTS=pars$P_FTS*pars$T_f,I_FTV=0,I_FTS=5) 
    
    out1=lsoda(state,t,derivs,pars)
    out2=lsoda(state,t,derivs,pars2)
    out3=lsoda(state,t,derivs,pars3)
    
                     #c("t","Prev_F","Prev_M", "Prev_CBB", "Prev_CNB", "Prev_CTS", "Prev_FBB", "Prev_FNB", "Prev_FTS")  
groupprevalence=cbind(t,out1[,32], out1[,33], out1[,34],out1[,35],out1[,36], out1[,37],out1[,38],out1[,39])       
write.table(groupprevalence,"groupprevalence.csv",sep=",",append = T, col.names =F)

#Baseline
prevf=rbind(out1[,32])
prevm=rbind(out1[,33])
prevcbb=rbind(out1[,34])
prevcnb=rbind(out1[,35])
prevcts=rbind(out1[,36])
prevfbb=rbind(out1[,37])
prevfnb=rbind(out1[,38])
prevfts=rbind(out1[,39])

write.table(prevf,"prevf.csv",sep=",",append = T, col.names =F)  
write.table(prevm,"prevm.csv",sep=",",append = T, col.names =F) 
write.table(prevcbb,"prevcbb.csv",sep=",",append = T, col.names =F) 
write.table(prevcnb,"prevcnb.csv",sep=",",append = T, col.names =F) 
write.table(prevcts,"prevcts.csv",sep=",",append = T, col.names =F) 
write.table(prevfbb,"prevfbb.csv",sep=",",append = T, col.names =F) 
write.table(prevfnb,"prevfnb.csv",sep=",",append = T, col.names =F) 
write.table(prevfts,"prevfts.csv",sep=",",append = T, col.names =F) 


#Zeta=0
prevf1=rbind(out2[,32])
prevm1=rbind(out2[,33])
prevcbb1=rbind(out2[,34])
prevcnb1=rbind(out2[,35])
prevcts1=rbind(out2[,36])
prevfbb1=rbind(out2[,37])
prevfnb1=rbind(out2[,38])
prevfts1=rbind(out2[,39])

write.table(prevf1,"prevf1.csv",sep=",",append = T, col.names =F)  
write.table(prevm1,"prevm1.csv",sep=",",append = T, col.names =F) 
write.table(prevcbb1,"prevcbb1.csv",sep=",",append = T, col.names =F) 
write.table(prevcnb1,"prevcnb1.csv",sep=",",append = T, col.names =F) 
write.table(prevcts1,"prevcts1.csv",sep=",",append = T, col.names =F) 
write.table(prevfbb1,"prevfbb1.csv",sep=",",append = T, col.names =F) 
write.table(prevfnb1,"prevfnb1.csv",sep=",",append = T, col.names =F) 
write.table(prevfts1,"prevfts1.csv",sep=",",append = T, col.names =F) 

#Circumcision=0

prevf2=rbind(out3[,32])
prevm2=rbind(out3[,33])
prevcbb2=rbind(out3[,34])
prevcnb2=rbind(out3[,35])
prevcts2=rbind(out3[,36])
prevfbb2=rbind(out3[,37])
prevfnb2=rbind(out3[,38])
prevfts2=rbind(out3[,39])

write.table(prevf2,"prevf2.csv",sep=",",append = T, col.names =F)  
write.table(prevm2,"prevm2.csv",sep=",",append = T, col.names =F) 
write.table(prevcbb2,"prevcbb2.csv",sep=",",append = T, col.names =F) 
write.table(prevcnb2,"prevcnb2.csv",sep=",",append = T, col.names =F) 
write.table(prevcts2,"prevcts2.csv",sep=",",append = T, col.names =F) 
write.table(prevfbb2,"prevfbb2.csv",sep=",",append = T, col.names =F) 
write.table(prevfnb2,"prevfnb2.csv",sep=",",append = T, col.names =F) 
write.table(prevfts2,"prevfts2.csv",sep=",",append = T, col.names =F) 


  }

d1=read.table("prevf.csv",sep=",",header = F)[,-1]
d2=read.table("prevm.csv",sep=",",header = F)[,-1]
d3=read.table("prevcbb.csv",sep=",",header = F)[,-1]
d4=read.table("prevcnb.csv",sep=",",header = F)[,-1]
d5=read.table("prevcts.csv",sep=",",header = F)[,-1]
d6=read.table("prevfbb.csv",sep=",",header = F)[,-1]
d7=read.table("prevfnb.csv",sep=",",header = F)[,-1]
d8=read.table("prevfts.csv",sep=",",header = F)[,-1]

d11=read.table("prevf1.csv",sep=",",header = F)[,-1]
d21=read.table("prevm1.csv",sep=",",header = F)[,-1]
d31=read.table("prevcbb1.csv",sep=",",header = F)[,-1]
d41=read.table("prevcnb1.csv",sep=",",header = F)[,-1]
d51=read.table("prevcts1.csv",sep=",",header = F)[,-1]
d61=read.table("prevfbb1.csv",sep=",",header = F)[,-1]
d71=read.table("prevfnb1.csv",sep=",",header = F)[,-1]
d81=read.table("prevfts1.csv",sep=",",header = F)[,-1]

d12=read.table("prevf2.csv",sep=",",header = F)[,-1]
d22=read.table("prevm2.csv",sep=",",header = F)[,-1]
d32=read.table("prevcbb2.csv",sep=",",header = F)[,-1]
d42=read.table("prevcnb2.csv",sep=",",header = F)[,-1]
d52=read.table("prevcts2.csv",sep=",",header = F)[,-1]
d62=read.table("prevfbb2.csv",sep=",",header = F)[,-1]
d72=read.table("prevfnb2.csv",sep=",",header = F)[,-1]
d82=read.table("prevfts2.csv",sep=",",header = F)[,-1]

par(mfrow=c(2,4))
plot(t,d1[1,], main = "females", type = "l", ylim = c(0, 10),  xlab = "Time (years)", ylab = "HIV prevalence for females",col = 1) 
for (i in 1:100){
  lines(t,d1[i,], col = topo.colors(100)[i]) 
} 
plot(t,d2[1,], main = "males", type = "l", ylim = c(0, 10),  xlab = "Time (years)", ylab = "HIV prevalence for males",col = 1) 
for (i in 1:100){
  lines(t,d2[i,], col = topo.colors(100)[i]) 
} 
plot(t,d3[1,],main = "CBB", type = "l", ylim = c(0, 50),  xlab = "Time (years)", ylab = "HIV prevalence for CBB",col = 1) 
for (i in 1:100){
  lines(t,d3[i,], col = topo.colors(100)[i]) 
} 
plot(t,d4[1,],main = "CNB", type = "l", ylim = c(0, 50),  xlab = "Time (years)", ylab = "HIV prevalence for CNB",col = 1) 
for (i in 1:100){
  lines(t,d4[i,], col = topo.colors(100)[i]) 
} 
plot(t,d5[1,],main = "CTS", type = "l", ylim = c(0, 50),  xlab = "Time (years)", ylab = "HIV prevalence for CTS",col = 1) 
for (i in 1:100){
  lines(t,d5[i,], col = topo.colors(100)[i]) 
} 
plot(t,d6[1,],main = "FBB", type = "l", ylim = c(0, 50),  xlab = "Time (years)", ylab = "HIV prevalence for FBB",col = 1) 
for (i in 1:100){
  lines(t,d6[i,], col = topo.colors(100)[i]) 
} 
plot(t,d7[1,],main = "FNB", type = "l", ylim = c(0, 50),  xlab = "Time (years)", ylab = "HIV prevalence for FNB",col = 1) 
for (i in 1:100){
  lines(t,d7[i,], col = topo.colors(100)[i]) 
} 
plot(t,d8[1,],main = "FTS",type = "l", ylim = c(0, 50),  xlab = "Time (years)", ylab = "HIV prevalence for FTS",col = 1) 
for (i in 1:100){
  lines(t,d8[i,], col = topo.colors(100)[i]) 
} 


par(mfrow=c(2,4))


plot(t,sapply(d1[-1, ], quantile, 0.50, names=FALSE), main = "females",type = "l", ylim = c(0, 50),  xlab = "Time (years)", ylab = "HIV prevalence for females",col = "blue") 
lines(t,sapply(d1[-1, ], quantile, 0.10, names=FALSE), col = "green") 
lines(t,sapply(d1[-1, ], quantile, 0.25, names=FALSE), col = "brown") 
lines(t,sapply(d1[-1, ], quantile, 0.75, names=FALSE), col = "yellow") 
lines(t,sapply(d1[-1, ], quantile, 0.90, names=FALSE), col = "red") 
legend("topleft", c("10 percentile","25 percentile","50 percentile", "75 percentile", "90 percentile"),lty =1,bty="n", col=c("green","brown","blue","yellow","red"), lwd = c(1,1,1,1,1),cex = 0.6)

plot(t,sapply(d2[-1, ], quantile, 0.50, names=FALSE), main = "males", type = "l", ylim = c(0, 50),  xlab = "Time (years)", ylab = "HIV prevalence for males",col = "blue") 
lines(t,sapply(d2[-1, ], quantile, 0.10, names=FALSE), col = "green") 
lines(t,sapply(d2[-1, ], quantile, 0.25, names=FALSE), col = "brown") 
lines(t,sapply(d2[-1, ], quantile, 0.75, names=FALSE), col = "yellow") 
lines(t,sapply(d2[-1, ], quantile, 0.90, names=FALSE), col = "red") 
legend("topleft", c("10 percentile","25 percentile","50 percentile", "75 percentile", "90 percentile"),lty = 1,bty="n", col=c("green","brown","blue","yellow","red"), lwd = c(1,1,1,1,1),cex = 0.6)

plot(t,sapply(d3[-1, ], quantile, 0.50, names=FALSE), main = "CBB", type = "l", ylim = c(0, 50),  xlab = "Time (years)", ylab = "HIV prevalence of CBB",col = "blue") 
lines(t,sapply(d3[-1, ], quantile, 0.10, names=FALSE), col = "green") 
lines(t,sapply(d3[-1, ], quantile, 0.25, names=FALSE), col = "brown") 
lines(t,sapply(d3[-1, ], quantile, 0.75, names=FALSE), col = "yellow") 
lines(t,sapply(d3[-1, ], quantile, 0.90, names=FALSE), col = "red") 
legend("topleft", c("10 percentile","25 percentile","50 percentile", "75 percentile", "90 percentile"),lty = 1,bty="n", col=c("green","brown","blue","yellow","red"), lwd = c(1,1,1,1,1),cex = 0.6)

plot(t,sapply(d4[-1, ], quantile, 0.50, names=FALSE), main = "CNB", type = "l", ylim = c(0, 50),  xlab = "Time (years)", ylab = "HIV prevalence for CNB",col = "blue") 
lines(t,sapply(d4[-1, ], quantile, 0.10, names=FALSE), col = "green") 
lines(t,sapply(d4[-1, ], quantile, 0.25, names=FALSE), col = "brown") 
lines(t,sapply(d4[-1, ], quantile, 0.75, names=FALSE), col = "yellow") 
lines(t,sapply(d4[-1, ], quantile, 0.90, names=FALSE), col = "red") 
legend("topleft", c("10 percentile","25 percentile","50 percentile", "75 percentile", "90 percentile"),lty = 1,bty="n", col=c("green","brown","blue","yellow","red"), lwd = c(1,1,1,1,1),cex = 0.6)

plot(t,sapply(d5[-1, ], quantile, 0.50, names=FALSE), main = "CTS", type = "l", ylim = c(0, 50),  xlab = "Time (years)", ylab = "HIV prevalence CTS",col = "blue") 
lines(t,sapply(d5[-1, ], quantile, 0.10, names=FALSE), col = "green") 
lines(t,sapply(d5[-1, ], quantile, 0.25, names=FALSE), col = "brown") 
lines(t,sapply(d5[-1, ], quantile, 0.75, names=FALSE), col = "yellow") 
lines(t,sapply(d5[-1, ], quantile, 0.90, names=FALSE), col = "red") 
legend("topleft", c("10 percentile","25 percentile","50 percentile", "75 percentile", "90 percentile"),lty = 1,bty="n", col=c("green","brown","blue","yellow","red"), lwd = c(1,1,1,1,1),cex = 0.6)

plot(t,sapply(d6[-1, ], quantile, 0.50, names=FALSE), main = "FBB", type = "l", ylim = c(0, 50),  xlab = "Time (years)", ylab = "HIV prevalence FBB",col = "blue") 
lines(t,sapply(d6[-1, ], quantile, 0.10, names=FALSE), col = "green") 
lines(t,sapply(d6[-1, ], quantile, 0.25, names=FALSE), col = "brown") 
lines(t,sapply(d6[-1, ], quantile, 0.75, names=FALSE), col = "yellow") 
lines(t,sapply(d6[-1, ], quantile, 0.90, names=FALSE), col = "red") 
legend("topleft", c("10 percentile","25 percentile","50 percentile", "75 percentile", "90 percentile"),lty = 1,bty="n", col=c("green","brown","blue","yellow","red"), lwd = c(1,1,1,1,1),cex = 0.6)

plot(t,sapply(d7[-1, ], quantile, 0.50, names=FALSE),main = "FNB", type = "l", ylim = c(0, 50),  xlab = "Time (years)", ylab = "HIV prevalence FNB",col = "blue") 
lines(t,sapply(d7[-1, ], quantile, 0.10, names=FALSE), col = "green") 
lines(t,sapply(d7[-1, ], quantile, 0.25, names=FALSE), col = "brown") 
lines(t,sapply(d7[-1, ], quantile, 0.75, names=FALSE), col = "yellow") 
lines(t,sapply(d7[-1, ], quantile, 0.90, names=FALSE), col = "red") 
legend("topleft", c("10 percentile","25 percentile","50 percentile", "75 percentile", "90 percentile"),lty = 1,bty="n", col=c("green","brown","blue","yellow","red"), lwd = c(1,1,1,1,1),cex = 0.6)

plot(t,sapply(d8[-1, ], quantile, 0.50, names=FALSE), main = "FTS",type = "l", ylim = c(0, 50),  xlab = "Time (years)", ylab = "HIV prevalence FTS",col = "blue") 
lines(t,sapply(d8[-1, ], quantile, 0.10, names=FALSE), col = "green") 
lines(t,sapply(d8[-1, ], quantile, 0.25, names=FALSE), col = "brown") 
lines(t,sapply(d8[-1, ], quantile, 0.75, names=FALSE), col = "yellow") 
lines(t,sapply(d8[-1, ], quantile, 0.90, names=FALSE), col = "red") 
legend("topleft", c("10 percentile","25 percentile","50 percentile", "75 percentile", "90 percentile"),lty = 1,bty="n", col=c("green","brown","blue","yellow","red"), lwd = c(1,1,1,1,1),cex = 0.6)

#Impact Zeta & circumcision
par(mfrow=c(2,2))

plot(t,sapply(d2[-1, ], quantile, 0.50, names=FALSE), main = "males",lwd=2,type = "l", ylim = c(0, 100),  xlab = "Time (years)", ylab = "HIV prevalence for males",col = "black") 
lines(t,sapply(d21[-1, ], quantile, 0.50, names=FALSE),lwd=2, col = "blue") 
lines(t,sapply(d22[-1, ], quantile, 0.50, names=FALSE),lwd=2, col = "red") 
legend("topleft", c("Baseline","No mixing with FTS","No circumcision"),lty =1,bty="n", col=c("black","blue","red"), lwd = c(2,2,2),cex = 0.8)

plot(t,sapply(d3[-1, ], quantile, 0.50, names=FALSE), main = "CBB",lwd=2,type = "l", ylim = c(0, 100),  xlab = "Time (years)", ylab = "HIV prevalence for CBB",col = "black") 
lines(t,sapply(d31[-1, ], quantile, 0.50, names=FALSE), lwd=2,col = "blue") 
lines(t,sapply(d32[-1, ], quantile, 0.50, names=FALSE), lwd=2,col = "red") 
legend("topleft", c("Baseline","No mixing with FTS","No circumcision"),lty =1,bty="n", col=c("black","blue","red"), lwd = c(2,2,2),cex = 0.8)

plot(t,sapply(d4[-1, ], quantile, 0.50, names=FALSE), main = "CNB",lwd=2,type = "l", ylim = c(0, 100),  xlab = "Time (years)", ylab = "HIV prevalence for CNB",col = "black") 
lines(t,sapply(d41[-1, ], quantile, 0.50, names=FALSE),lwd=2, col = "blue") 
lines(t,sapply(d42[-1, ], quantile, 0.50, names=FALSE),lwd=2, col = "red") 
legend("topleft", c("Baseline","No mixing with FTS","No circumcision"),lty =1,bty="n", col=c("black","blue","red"), lwd = c(2,2,2),cex =0.8)

plot(t,sapply(d5[-1, ], quantile, 0.50, names=FALSE), main = "CTS",lwd=2,type = "l", ylim = c(0, 100),  xlab = "Time (years)", ylab = "HIV prevalence for CTS",col = "black") 
lines(t,sapply(d51[-1, ], quantile, 0.50, names=FALSE),lwd=2, col = "blue") 
lines(t,sapply(d52[-1, ], quantile, 0.50, names=FALSE),lwd=2, col = "red") 
legend("topleft", c("Baseline","No mixing with FTS","No circumcision"),lty =1,bty="n", col=c("black","blue","red"), lwd = c(2,2,2),cex = 0.8)

plot(t,sapply(d1[-1, ], quantile, 0.50, names=FALSE), main = "females",type = "l",lwd=2, ylim = c(0, 100),  xlab = "Time (years)", ylab = "HIV prevalence for females",col = "black") 
lines(t,sapply(d11[-1, ], quantile, 0.50, names=FALSE),lwd=2, col = "blue") 
lines(t,sapply(d12[-1, ], quantile, 0.50, names=FALSE),lwd=2, col = "red") 
legend("topleft", c("Baseline","No mixing with FTS","No circumcision"),lty =1,bty="n", col=c("black","blue","red"), lwd = c(2,2,2),cex = 0.8)

plot(t,sapply(d6[-1, ], quantile, 0.50, names=FALSE), main = "FBB",lwd=2,type = "l", ylim = c(0, 100),  xlab = "Time (years)", ylab = "HIV prevalence for FBB",col = "black") 
lines(t,sapply(d61[-1, ], quantile, 0.50, names=FALSE), lwd=2,col = "blue") 
lines(t,sapply(d62[-1, ], quantile, 0.50, names=FALSE), lwd=2,col = "red") 
legend("topleft", c("Baseline","No mixing with FTS","No circumcision"),lty =1,bty="n", col=c("black","blue","red"), lwd = c(2,2,2),cex = 0.8)

plot(t,sapply(d7[-1, ], quantile, 0.50, names=FALSE), main = "FNB",lwd=2,type = "l", ylim = c(0, 100),  xlab = "Time (years)", ylab = "HIV prevalence for FNB",col = "black") 
lines(t,sapply(d71[-1, ], quantile, 0.50, names=FALSE),lwd=2, col = "blue") 
lines(t,sapply(d72[-1, ], quantile, 0.50, names=FALSE),lwd=2, col = "red") 
legend("topleft", c("Baseline","No mixing with FTS","No circumcision"),lty =1,bty="n", col=c("black","blue","red"), lwd = c(2,2,2),cex = 0.8)

plot(t,sapply(d8[-1, ], quantile, 0.50, names=FALSE), main = "FTS",lwd=2,type = "l", ylim = c(0, 100),  xlab = "Time (years)", ylab = "HIV prevalence for FTS",col = "black") 
lines(t,sapply(d81[-1, ], quantile, 0.50, names=FALSE), lwd=2,col = "blue") 
lines(t,sapply(d82[-1, ], quantile, 0.50, names=FALSE), lwd=2,col = "red") 
legend("topleft", c("Baseline","No mixing with FTS","No circumcision"),lty =1,bty="n", col=c("black","blue","red"), lwd = c(2,2,2),cex = 0.8)



par(mfrow=c(2,2))


boxplot(d2,range=0,col="black",main = "males",type = "l", xlab = "Time (years)", ylab = "HIV prevalence for males",ylim = c(0, 100), xaxt="n",border=c("black"))
boxplot(d21, add=TRUE,range=0, col="blue",ylim = c(0, 100), xaxt="n",border=c("blue"))
boxplot(d22, add=TRUE,range=0, col="red", xaxt="n",border=c("red"))
axis(side = 1, at=1:106, labels=1975:2080, tck =-.02, las=1, cex.axis=1)
legend("topright", c("Baseline","No mixing with FTS","No circumcision"), fill = c("black", "blue","red"),bty="n", cex = 0.5)

boxplot(d3,range=0,col="black",main = "CBB",type = "l", xlab = "Time (years)", ylab = "HIV prevalence for CBB",ylim = c(0, 100), xaxt="n",border=c("black"))
boxplot(d31, add=TRUE,range=0, col="blue",ylim = c(0, 100), xaxt="n",border=c("blue"))
boxplot(d32, add=TRUE,range=0, col="red", xaxt="n",border=c("red"))
axis(side = 1, at=1:106, labels=1975:2080, tck =-.02, las=1, cex.axis=1)
legend("topright", c("Baseline","No mixing with FTS","No circumcision"), fill = c("black", "blue","red"),bty="n", cex = 0.5)

boxplot(d4,range=0,col="black",main = "CNB",type = "l", xlab = "Time (years)", ylab = "HIV prevalence for CNB",ylim = c(0, 100), xaxt="n",border=c("black"))
boxplot(d41, add=TRUE,range=0, col="blue",ylim = c(0, 100), xaxt="n",border=c("blue"))
boxplot(d42, add=TRUE,range=0, col="red", xaxt="n",border=c("red"))
axis(side = 1, at=1:106, labels=1975:2080, tck =-.02, las=1, cex.axis=1)
legend("topright", c("Baseline","No mixing with FTS","No circumcision"), fill = c("black", "blue","red"),bty="n", cex = 0.5)

boxplot(d5,range=0,col="black",main = "CTS",type = "l", xlab = "Time (years)", ylab = "HIV prevalence for CTS",ylim = c(0, 100), xaxt="n",border=c("black"))
boxplot(d51, add=TRUE,range=0, col="blue",ylim = c(0, 100), xaxt="n",border=c("blue"))
boxplot(d52, add=TRUE,range=0, col="red", xaxt="n",border=c("red"))
axis(side = 1, at=1:106, labels=1975:2080, tck =-.02, las=1, cex.axis=1)
legend("topright", c("Baseline","No mixing with FTS","No circumcision"), fill = c("black", "blue","red"),bty="n", cex = 0.5)

boxplot(d1,range=0,col="black",main = "females",type = "l", xlab = "Time (years)", ylab = "HIV prevalence for females",ylim = c(0, 100), xaxt="n",border=c("black"))
boxplot(d11, add=TRUE,range=0, col="blue",ylim = c(0, 100), xaxt="n",border=c("blue"))
boxplot(d12, add=TRUE,range=0, col="red", xaxt="n",border=c("red"))
axis(side = 1, at=1:106, labels=1975:2080, tck =-.02, las=1, cex.axis=1)
legend("topright", c("Baseline","No mixing with FTS","No circumcision"), fill = c("black", "blue","red"),bty="n", cex = 0.5)


boxplot(d6,range=0,col="black",main = "FBB",type = "l", xlab = "Time (years)", ylab = "HIV prevalence for FBB",ylim = c(0, 100), xaxt="n",border=c("black"))
boxplot(d61, add=TRUE,range=0, col="blue",ylim = c(0, 100), xaxt="n",border=c("blue"))
boxplot(d62, add=TRUE,range=0, col="red", xaxt="n",border=c("red"))
axis(side = 1, at=1:106, labels=1975:2080, tck =-.02, las=1, cex.axis=1)
legend("topright", c("Baseline","No mixing with FTS","No circumcision"), fill = c("black", "blue","red"),bty="n", cex = 0.5)

boxplot(d7,range=0,col="black",main = "FNB",type = "l", xlab = "Time (years)", ylab = "HIV prevalence for FNB",ylim = c(0, 100), xaxt="n",border=c("black"))
boxplot(d71, add=TRUE,range=0, col="blue",ylim = c(0, 100), xaxt="n",border=c("blue"))
boxplot(d72, add=TRUE,range=0, col="red", xaxt="n",border=c("red"))
axis(side = 1, at=1:106, labels=1975:2080, tck =-.02, las=1, cex.axis=1)
legend("topright", c("Baseline","No mixing with FTS","No circumcision"), fill = c("black", "blue","red"),bty="n", cex = 0.5)

boxplot(d8,range=0,col="black",main = "FTS",type = "l", xlab = "Time (years)", ylab = "HIV prevalence for FTS",ylim = c(0, 100), xaxt="n",border=c("black"))
boxplot(d81, add=TRUE,range=0, col="blue",ylim = c(0, 100), xaxt="n",border=c("blue"))
boxplot(d82, add=TRUE,range=0, col="red", xaxt="n",border=c("red"))
axis(side = 1, at=1:106, labels=1975:2080, tck =-.02, las=1, cex.axis=1)
legend("topright", c("Baseline","No mixing with FTS","No circumcision"), fill = c("black", "blue","red"),bty="n", cex = 0.5)





par(mfrow=c(2,4))
plot(t,d1[1,], main = "females", pch = 1, ylim = c(0, 100),  xlab = "Time (years)", ylab = "HIV prevalence for females",col = "black") 
for (i in 1:100){
  points(t,d1[i,], col = "black")
  points(t,d11[i,], col="blue")
  points(t,d12[i,], col="red")
} 
legend("topleft", c("Baseline","No mixing with FTS","No circumcision"),pch=c(1,1,1),bty="n", col=c("black","blue","red"),cex = 0.8)


plot(t,d2[1,], main = "males", pch = 1, ylim = c(0, 100),  xlab = "Time (years)", ylab = "HIV prevalence for males",col = "black") 
for (i in 1:100){
  points(t,d2[i,], col = "black")
  points(t,d21[i,], col="blue")
  points(t,d22[i,], col="red")
} 
legend("topleft", c("Baseline","No mixing with FTS","No circumcision"),pch=c(1,1,1),bty="n", col=c("black","blue","red"),cex = 0.8)

plot(t,d3[1,], main = "CBB", pch = 1, ylim = c(0, 100),  xlab = "Time (years)", ylab = "HIV prevalence for CBB",col = "black") 
for (i in 1:100){
  points(t,d3[i,], col = "black")
  points(t,d31[i,], col="blue")
  points(t,d32[i,], col="red")
} 
legend("topleft", c("Baseline","No mixing with FTS","No circumcision"),pch=c(1,1,1),bty="n", col=c("black","blue","red"),cex = 0.8)

plot(t,d4[1,], main = "CNB", pch = 1, ylim = c(0, 100),  xlab = "Time (years)", ylab = "HIV prevalence for CTS",col = "black") 
for (i in 1:100){
  points(t,d4[i,], col = "black")
  points(t,d41[i,], col="blue")
  points(t,d42[i,], col="red")
} 
legend("topleft", c("Baseline","No mixing with FTS","No circumcision"),pch=c(1,1,1),bty="n", col=c("black","blue","red"),cex = 0.8)

plot(t,d5[1,], main = "CTS", pch = 1, ylim = c(0, 100),  xlab = "Time (years)", ylab = "HIV prevalence for CTS",col = "black") 
for (i in 1:100){
  points(t,d5[i,], col = "black")
  points(t,d51[i,], col="blue")
  points(t,d52[i,], col="red")
} 
legend("topleft", c("Baseline","No mixing with FTS","No circumcision"),pch=c(1,1,1),bty="n", col=c("black","blue","red"),cex = 0.8)

plot(t,d6[1,], main = "FBB", pch = 1, ylim = c(0, 100),  xlab = "Time (years)", ylab = "HIV prevalence for FBB",col = "black") 
for (i in 1:100){
  points(t,d6[i,], col = "black")
  points(t,d61[i,], col="blue")
  points(t,d62[i,], col="red")
} 
legend("topleft", c("Baseline","No mixing with FTS","No circumcision"),pch=c(1,1,1),bty="n", col=c("black","blue","red"),cex = 0.8)

plot(t,d7[1,], main = "FNB", pch = 1, ylim = c(0, 100),  xlab = "Time (years)", ylab = "HIV prevalence for FNB",col = "black") 
for (i in 1:100){
  points(t,d7[i,], col = "black")
  points(t,d71[i,], col="blue")
  points(t,d72[i,], col="red")
} 
legend("topleft", c("Baseline","No mixing with FTS","No circumcision"),pch=c(1,1,1),bty="n", col=c("black","blue","red"),cex = 0.8)

plot(t,d8[1,], main = "FTS", pch = 1, ylim = c(0, 100),  xlab = "Time (years)", ylab = "HIV prevalence for FTS",col = "black") 
for (i in 1:100){
  points(t,d8[i,], col = "black")
  points(t,d81[i,], col="blue")
  points(t,d82[i,], col="red")
} 
legend("topleft", c("Baseline","No mixing with FTS","No circumcision"),pch=c(1,1,1),bty="n", col=c("black","blue","red"),cex = 0.8)







