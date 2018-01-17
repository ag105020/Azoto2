'''
Created on 2013/01/05


#----------------------------------------------------------------------------------------------
*About 317s
In this version, we make 605_06_01_01 like 314_41 for updating parameters for the second
paper
#----------------------------------------------------------------------------------------------
#note

In this file, we will use the advanced data file that Eric helped
We are referring to 314_59 for this purpose.
#----------------------------------------------------------------------------------------------
@author: ag105020
'''
from pylab import *
from Buhler87_data_07 import *
import time
rc('text', usetex=True)
rc('font', family='serif')
rc('font', serif='Times New Roman')


##########################################################################################################################################
def az2(Nin):
    #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    #1.Parameters----------------------------
    #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    O2sat=0.225         #(mol m-3) Saturated oxygen concentration (Post et al 1982)
    O2satratio=0.05         #(diemnsionless) The saturation ratio of oxygen
    global O2percent
    O2percent=O2satratio*100
    O2=O2sat*O2satratio     #(mol m-3) oxygen concentration of the medium
    
    #plotcolor=buhler_data.plotcolor
    #OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    #defining basic parameters
    #OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    mu25=8.903*10**(-4)           #(kg/m/s) Dynamic viscosity (dounenseikeisuu) of water at T=25C degrees Kestin 1978 (15 Computation of diffusion coefficient for different temperature.xlsx)
    mu30=7.975*10**(-4)          #(kg/m/s) Dynamic viscosity (dounenseikeisuu) of water at T=30C degrees Kestin 1978 (15 Computation of diffusion coefficient for different temperature.xlsx)
    Kelvin=273.15
    TeffectonD25to30=(30+Kelvin)/(25+Kelvin)/mu30*mu25  #(dimensionless) converting diffusion coefficient from 25C to 30C
#    O2=kuhla_data.O2concentration
#    plotcolor=kuhla_data.plotcolor
    #manage later
    
    Dh=0.15 #(h-1) dilution rate
    #Dh=0.3
    D=Dh/3600 #(s-1) dilution rate
    
    CNstep=0.1
    CNwidth=15
    CN=arange(CNstep,CNwidth+CNstep,CNstep) #(molSucrose molN-1) C to N ratio in the incoming medium
    CNnumber=size(CN)
    U=arange(0,CNnumber,1)
    global CHin
    CHin=30*CN      #(molC m-3): CH concentration in the poring media
    #Nin=2.5         #(molN m-3): N concentration in the poring media
    #print(CHin)
#    Umax=Dhnumber
#    U=arange(0,Umax,1)    #for loop
    D=Dh/3600          #(s-1) growth rate  
    Rv=(18.547*log(Dh)+123.3)/100    #(dimensionless) vital ratio (the ratio of vital cells) (181-11) Equation from
                                        # from "Postgate 1973.xlsx"
#     for i in U:
#         if Rv<0:
#             Rv=nan
#         else:
#             break
#     Rv[Rv<0]=nan
    
    Dhgrow=Dh/Rv                        #(h-1) Growth of actually growing cells (181-11)
    Dgrow=D/Rv                          #(s-1) Growth of actually growing cells (181-11) 
    O2sat=0.225             #(mol m-3) 100% oxygen from post 1982
    O2min=0.225/100         #(mol m-3) 1% oxygen from post 1982
    O2percent=O2/O2sat*100
    Llinear=(0.0124*O2percent+2.1009)*10**(-6)    #(m) cell length (from "cell size.xlsx")
    r2linear=(0.0061*O2percent+1.5224)*10**(-6)   #(m) cell width (from "cell size.xlsx")
    RLg=(Llinear/2+r2linear)/3
    F=0
    hh=0.00079
    h=hh*O2/O2          #diffusivity of glycolipid layer compared to water (132-23)
    b7=hh*20      #(159-16) to make y intercept of l =1
    c7=1.2787679502151037            #(159-16) same
    d7=0.0018083821866941114
    a7=b7**(1/c7)*d7          #(159-16) same
    l=b7*(d7/(O2*b7**(1/c7)+a7))**c7
    ep3=h/l      #diffusivity of glycolipid layer compared to alginate layer (132-23)
    m=4.695556573362741e-18                #(mol cell-1 s-1) maintenance carbohydrate consumption (separated from respiration) (147-13)
    m=0
    mr=1e-18        #(molO2 cell-1 s-1) maintenance respiration
    Rg=1/29            #ratio of glycolipid layer to the cell radius
    R=RLg/(1+Rg)        #Radius of the cell without glycolipid layer(m)    
    Lg=R*Rg             #(m) thickness of glycolipid layer
    Ra=4/19         #(dimensionless) the ratio of alginate layer to the cell radius
    Ra=0
    La=R*Ra         #(m) thickness of alginate layer (132-23)
    x0=1/h*(1/R-1/(R+Lg))+ep3/h*1/(R+Lg)+(1/(R+Lg+La))*(1-ep3/h)  #effect of glycolipid layer (Kei 133-11, 134-17)
    r5=1/(x0*R)         #diffusivity efficiency of cell membrane
    V=4/3*pi*RLg**3   #Volume of the cell (m3/cell)
    C=6             #number of C in one carbohydrate (ex. C(glucose)=6)
    Kch=1             #half saturation of protein CH for biomass production
    Do2_25=2.12*10**(-9)  #Diffusion coefficient of O2 in the water at 25C (m2/s)
    Do2_30=Do2_25*TeffectonD25to30  #(m2/2) diffusion coefficient of O2 in the water at 30C
    Do2=Do2_30*r5         #Diffusion coefficient of O2 in the water (m2/s)
    Dch_25=6.728*10**(-10)  #Diffusion coefficient of glucose in the water (m2/s)
    Dch_30=Dch_25*TeffectonD25to30  #Diffusion coefficient of glucoe in water at 30C (m2/s)
    Bch=16                #(dimensionless) reduced diffusivity for carbohydrate : $$$$$$$$$$$$$$$$$$$
    Dch=Dch_30*r5/Bch          #Diffusion coefficient of glucose in the water (m2/s)
    Dnh4_25=1.98*10**(-9)  #Diffusion coefficient of ammonium in the water at 25C (m2/s)
    Dnh4_30=Dnh4_25*TeffectonD25to30     #Diffusion coefficient of ammonium in water at 30C (m2/s)
    Dnh4=Dnh4_30*r5        #Diffusion coefficient of ammonium in the water (m2/s)
    VV=arange(0,5,1)    #for loop
    O2cri=0.005195425946179848          #(mol/m3) the critical level of O2 concentration inside the cell
    O2cri=0.00000000000001
    lmax00=1000         #(molCs-1cell-1m-3) maximum biomass production rate per volume
    lmax=lmax00*V*ones(size(Dh))      #(molCs-1cell-1) maximum biomass production rate
  #  lmaxgrow=lmax/Rv
    rho=0.22*10**6/12                #(molCbiomass/m3/cell) biomass density in the cell (105-10)
    Q=rho*V             #(molC/cell) biomass per cell
#    CHin=100    #(mol/m3): CH concentration in the poring water
    n=6             #number of C in one bacterial biomass (BB)
    a=10.8             #number of H in one BB
    b=2.9             #number of O in one BB
    c=1.5             #number of N in one BB
    d=4*n+a-2*b-3*c #inverse number in the coefficient of BB in the half reaction for one e-
    BB=12*n+1*a+16*b+14*c #(g/mol): mass of BB (bacterial biomass)
    y=c/d           #coefficient of NH4+ in the half reaction of BB production
    z=1/4           #coefficient of NH4+ in the half reaction of nitrogen-fixation
    pr=1/1.32      #(dimensionless): protein ration in biomass
    r4=c/n          #(dimensionless): the ratio of N to C in the bacterial biomass    
    KN=Kch*r4       #(mol/m3): half saturation of protein N for biomass production
    #OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    #Computation of fe0, fpr and fn considering material, redox and energy balance
    #OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    fs00=(1/y)/(1/y+1/z)        #The ratio of electron used for protein synthesis
    fn00=(1/z)/(1/y+1/z)        #The ration of electron used for nitrogen fixation
    ep=0.22                   #energy efficiency for the production of energy and the consumption of energy
    dgc0=41.35      #The free energy necessary for the half reaction of glucose production (kJ/e-mol)
    dgATP=50      #(kJ/ATP): energy produced by the reaction of ATP -> ADP (147-19)
    dgn=2*dgATP-dgc0*ep   #The free energy necessary (dg) for the half reaction of nitrogen fixation (kJ/e-mol)
    dgp=35.09-dgc0  #dg for production of pyruvate from glucose (kJ/e-mol))
    dgpc=3.33*1/d*(12*n+1*a+16*b+14*c)  #dg for the production of BB (bacterial biomass) from pyruvate) (147-17)
    dgr=-120.07     #-dg for the energy production pathway (kJ/e-mol)

    if dgp<0:       
        ep1=1/ep    #change ep1 depending on the sign of dgp    
    else:
        ep1=ep
    A=(fn00*dgn+fs00*(dgp/ep1+dgpc/ep))/(-ep*dgr) #A is related to fs0 and fe0
    fe0=A/(1+A)     #the ratio of electron used for energy production
    fs0=1-fe0       #the ratio of electron used for biomass synthesis+nitrogen fixation
    fpr=fs0*fs00    #the ratio of electron used for biomass synthesis
    fn=fs0*fn00     #the ratio of electron used for nitrogen fixation  
    
    epNH4=0.45       #(dimensionless) energy transfer efficiency in ammonium assimilating cells
    #epNH4=0.30
    epNH4=0.54
    
    if dgp<0:       
        epNH4_1=1/epNH4    #change ep1 depending on the sign of dgp    
    else:
        epNH4_1=epNH4
    
    A1=(dgp/epNH4_1+dgpc/epNH4)/(-epNH4*dgr) #A is related to fs0 and fe0
    fe01=A1/(1+A1)     #the ratio of electron used for energy production
    fs01=1-fe01       #the ratio of electron used for biomass synthesis+nitrogen fixation
    fpr1=fs01         #the ratio of electron used for biomass synthesis
                     #the ratio of electron used for nitrogen fixation
      
    #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    #2.--Stoichiometry (to get E1(E for the case O2cri>O2in))------------------
    #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    S1=array([["","CH","H2O","CO2","O2","HCO3-","NH4+","N2","H2","BB","H+","e-"],
              ["'-Rd",1/24,0.25,-0.25,0.,0.,0.,0.,0.,0.,-1.,-1.],
              ["Ra",0.,-0.5,0.,0.25,0.,0.,0.,0.,0.,1.,1.],
              ["Rpr",0.,-(2*n-b+c)/d,(n-c)/d,0.,c/d,c/d,0.,0.,-1/d,1.,1.],
              ["Rn",0.,0.,0.,0.,0.,-0.25,0.125,-0.125,0.,1.25,1.]])
    
    #OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    #for creating S2 (f*R for electron acceptance)
    #OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    Mu=array([[1],[fe0],[fpr],[fn]])  #column of f
    S2=copy(S1[1:5,1:12])           #use copy so S2 does not respond to the change in S1
    
    #convert S2 from string data type to float64 data type
    S2=S2.astype(float64)
    #add numbers for columns and raws for counting 
    S21=arange(1,12,1)
    S22=arange(0,5,1)
    S22=S22.reshape(5,1)
    S2=vstack((S21,S2))
    S2=hstack((S22,S2))
    #calculate f*R
    S2[1:5,1:12]=Mu*S2[1:5,1:12]
    
    #OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    #for creating S3 (f*R for electron donation)
    #OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    S3=vstack((S2[1,1:12],S2[1,1:12],S2[1,1:12]))
    S3=Mu[1:]*S3
    
    #OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    #for creating S4 (f*R for "electron donation + electron acceptance")
    # and RR, which is the entire reaction
    #OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    S4=S2[2:,1:]+S3     #S4 is f*R for "electron donation + electron acceptance"
    p=-S4[1,8]*n                        #(molC/e-mol) CH consumption for protein production
    h0=(S4[1,0]+S4[2,0])*C               #(molC/e-mol) CH consumption for other than energy production
    alp=p/h0                             #(CO2 from (Ra-Rd))/(CH for protein production) (see 72-5)
    beta=-(S4[1,2]+S4[1,4]+S4[2,2])/h0   #(CO2 from (Rpr-Rd) + CO2 from (RN-Rd))/(CH for protein production) (see 72-5)
    RR=S4[0]+S4[1]+S4[2]
    RR1=copy(RR)
    RR1=vstack((S1[0,1:],RR1))
    
    #OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    #calculating r values
    #OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO   
    r1=C*(S4[1,0]+S4[2,0])/(-S4[1,8]*n) #(CH consumption except for energy production)/(biomass synthesis)
    r2=A
    r3=S4[0,3]/(S4[0,0]*C)              #(O2 consumption)/(CH consumption for energy production))
    #OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    #output of each array into CSV files
    #OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
#     savetxt("S1.csv", S1, delimiter=",",fmt='%s')
#     savetxt("S2.csv", S2, delimiter=",",fmt='%2.3f')
#     savetxt("S3.csv", S3, delimiter=",",fmt='%2.8f')
#     savetxt("S4.csv", S4, delimiter=",",fmt='%2.8f')
#     savetxt("RR.csv", RR, delimiter=",",fmt='%2.8f')
#     savetxt("RR1.csv", RR1, delimiter=",",fmt='%s')
 
    #OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    #Getting yield (Y) and the ratio of CO2 production rate to CH consumption (E)
    #OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    Y=-(RR[8]*n)/(C*RR[0])              #Yield
    E1=1/Y-1                            #the ratio of CO2 production rate to CH consumption
    
    #=====================================================================================
    #NH4+ absorption case
    #=====================================================================================
    #OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    #for creating S2 (f*R for electron acceptance)
    #OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
            
    Mua=array([[1],[fe01],[fpr1],[0.0]])  #column of f
    S2a=copy(S1[1:5,1:12])               #use copy so S2 does not respond to the change in S1
    S2a=S2a.astype(float64)
    #add numbers for columns and raws for counting 
    S21a=arange(1,12,1)
    S22a=arange(0,5,1)
    S22a=S22a.reshape(5,1)
    S2a=vstack((S21a,S2a))
    S2a=hstack((S22a,S2a))
    #calculate f*R
    S2a[1:5,1:12]=Mua*S2a[1:5,1:12]
    
    
    #OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    #for creating S3 (f*R for electron donation)
    #OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    S3a=vstack((S2a[1,1:12],S2a[1,1:12],S2a[1,1:12]))
    S3a=Mua[1:]*S3a
    
    #OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    #for creating S4 (f*R for "electron donation + electron acceptance")
    # and RR, which is the entire reaction
    #OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    S4a=S2a[2:,1:]+S3a     #S4 is f*R for "electron donation + electron acceptance"
    RRa=S4a[0]+S4a[1]+S4a[2]
    RR1a=copy(RRa)
    RR1a=vstack((S1[0,1:],RR1a))    
    pa=-S4a[1,8]*n                        #(molC/e-mol) CH consumption for protein production
    ha=(S4a[1,0]+S4a[2,0])*C               #(molC/e-mol) CH consumption for other than energy production
    alpa=pa/ha                             #(CO2 from (Ra-Rd))/(CH for protein production) (see 72-5)
    betaa=-(S4a[1,2]+S4a[1,4]+S4a[2,2])/ha   #(CO2 from (Rpr-Rd) + CO2 from (RN-Rd))/(CH for protein production) (see 72-5)
    #OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    #output of each array into CSV files
    #OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
#     savetxt("S11a.csv", S1, delimiter=",",fmt='%s')
#     savetxt("S21a.csv", S2a, delimiter=",",fmt='%2.3f')
#     savetxt("S31a.csv", S3a, delimiter=",",fmt='%2.8f')
#     savetxt("S41a.csv", S4a, delimiter=",",fmt='%2.8f')
#     savetxt("RR01a.csv", RRa, delimiter=",",fmt='%2.8f')
#     savetxt("RR11a.csv", RR1a, delimiter=",",fmt='%s')
 
    #OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    #Getting yield (Y) and the ratio of CO2 production rate to CH consumption (E)
    #OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO

    Y1=-(RRa[8]*n)/(C*RRa[0])              #Yield
    E3=1/Y1-1 
    
    #OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    #calculating r values
    #OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    r1a=C*(S4a[1,0]+S4a[2,0])/(-S4a[1,8]*n) #(CH consumption except for energy production)/(biomass synthesis)
    r3a=S4a[0,3]/(S4a[0,0]*C) 
    r2a=A1
    
    #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    #Chemostat box model
    #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    #OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    #N (NH4+) limited case (Kei 100)
    #OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    Nc=(KN*Q*Dgrow)/(lmax-Q*Dgrow)                  #(mol/m3): NH4+ concentration in the cell (100-15) (r4 is cancenlled here (158-33)
    lsngrow=lmax*Nc/(Nc+KN)                     #(molC/s/cell): biomass synthesis rate
    N=(lsngrow*r4)/(4*pi*R*Dnh4)+Nc             #(mol/me): NH4+ concentration in the vessel (100-15)
    Vngrow=4*pi*R*Dnh4*(N-Nc)                   #(molN s-1 m-2) nitrogen uptake rate for the growing cells
    Vn=Vngrow*Rv                                #(molN s-1 m-2) nitrogen uptake rate in average
    xn=D*(Nin-N)/Vn       #(cell/m3): number density of the cells (100-16)
    rhog=rho*BB/n                      #(gbiomass/m3/cell): biomass density in the cell
    PRc=rhog*pr                             #(gPR/m3): protein density in the cell
    PRm=PRc*V                               #(gPR/cell): protein mass in one cell
    PRd0n=PRm*xn                            #(gPR/m3): Protein density in the vessel
    PRdn=PRd0n*10**(-3)                     #(mgPR/ml): Protein density in the vessel
    PPP=ones((150,1))                       #used to make PRdn array
    PRdn2=PPP*PRdn                          #array version of PRdn
    
    #OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    #Calculation of potential CH and O2c when nitrogne limiting (new kei 156)
    #OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    
    CHc=(Dgrow*Kch*Q)/(lmax-Dgrow*Q)            #(mol/m3) CH concentration inside the cell
    lsgrow=lmax*CHc/(CHc+Kch)                   #(mol/s/cell): biomass synthesis rate
    ls=lsgrow*Rv                                #(molC s-1 cell-1): Biomass synthesis rate (average)
    CHpot=(D*CHin+xn*4*pi*R*Dch*CHc*Rv)/(D+xn*4*pi*R*Dch*Rv)          #(mol m-3) Highest CH potential during NH4+ uptake only for N source (156-37)(189-7)
    Vchpotgrow=4*pi*R*Dch*(CHpot-CHc)           #(mol s-1 cell-1) Potential carbohydrate uptake rate of the growing cells  (156-56)
    Vchpot=Vchpotgrow*Rv           #(mol s-1 cell-1) Average potential carbohydrate uptake rate (156-56)
    r2pot=(alpa*Vchpotgrow-alpa*lsgrow-lsgrow*betaa-alpa*lsgrow*F-alpa*m)/lsgrow      #(dimensionless) potentially highest r2 (156-36)(189-7)
    O2cpot=(-1)*(lsgrow*r1a*r2pot*r3a)/(4*pi*R*Do2)+O2    #(mol m-3) potentially lowest intracellular O2
    res1grow=4*pi*R*Do2*(O2-O2cpot)             #(mol O2 s-1 cell-)
    res1=res1grow*Rv
    
    #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    #For stack plot of CH consumption (178-1~2)(181-13) Based on 605_06_04_07
    #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    CCHbiomass_a=ls/(-S4a[1,8])*S4a[1,0]                #(molC cell-1 s-1) consumption of CH for biomsas production

    CCHbalancedresp_a=ls/(-S4a[1,8])*S4a[0,0]           #(molC cell-1 s-1) consumption of CH for balanced respiration
                                                    #Consumption of CH for respiration for N2 fixation
    CCHextraresp_a=Vchpot-CCHbiomass_a-CCHbalancedresp_a    #(molC cell-1 s-1) Carbohydrate consumption for extra respiration
    
    #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    
    #------------------------
    #Obtaining CH concentration
    #------------------------
    CH0=(lsgrow*(1+E3+F)+m)/(4*pi*R*Dch)+CHc        #(molC m-3) sucrose concentration in the vessel in ammonium no extra respiration case 

    #------------------------
    #Obtaining Respiration
    #------------------------
    res0grow=lsgrow*r1a*r2a*r3a                 #(mol O2 s-1 cell-1) respiration of growing cells when CH limited
    res0=res0grow*Rv                #(mol O2 s-1 cell-1) the average respiration of all the cells
    
    
    #=====================================================================================
    #NH4+ and N2 mixed case 
    #=====================================================================================    
    df=0.01                            #plese make small enoough value and value that can divide 1 evenly  (e.g. 0.001 etc.) 
    floop=arange(0,1/df,1)          #f for loop
    far=arange(df,1+df,df)          #f arange
    CH=copy(far)                    #This is the only one that I use
    EEi=copy(far)/copy(far)
    xnh4=copy(EEi)
    a7=CNwidth/CNstep
    res2grow=copy(EEi)
    res2=copy(EEi)
    Ncb=copy(EEi)
    N2=copy(EEi)
    CCHbiomass_b=copy(EEi)
    CCHN2fix_b=copy(EEi)
    CCHbalancedresp_b=copy(EEi)
    CCHrespN2fix_b=copy(EEi)
    CCHrespbiomass_b=copy(EEi)
    CCHextraresp_b=copy(EEi)
    
    #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    #$$$$$$$$$start loop to get values of CH$$$$$$$$$$"
    #If you need more values from inside the loop, creat array with the same width as far as I did in CH#
    #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    Ep=copy(far)/copy(far)
    
    #==================
    #For binary search
    #==================
    imax=size(floop)-1
    imin=floop[0]
    center=(imax+imin)//2
    
    for i in floop:
        f=far[i]
        f0a=f*(-1)*RR[0]/RR[8]
        f0b=(1-f)*(-1)*RRa[0]/RRa[8]
        f0=f0a+f0b
        f1=f0a/f0                       #(157-18) conversion of f : (f-1) in nitrogen to carbohydrate
        S4b=f1*S4+(1-f1)*S4a
        RRb=f1*RR+(1-f1)*RRa
        RR1b=copy(RRb)
        RR1b=vstack((S1[0,1:],RR1b))
        Y2=-(RRb[8]*n)/(C*RRb[0])              #Yield
        E5=1/Y2-1 
        r1b=C*(S4b[1,0]+S4b[2,0])/(-S4b[1,8]*n) #(CH consumption except for energy production)/(biomass synthesis)
        r3b=S4b[0,3]/(S4b[0,0]*C) 
        r2b=S4b[0,0]/(S4b[1,0]+S4b[2,0]) 
        pb=-S4b[1,8]*n                        #(molC/e-mol) CH consumption for protein production
        hb=(S4b[1,0]+S4b[2,0])*C               #(molC/e-mol) CH consumption for other than energy production
        alpb=pb/hb                             #(CO2 from (Ra-Rd))/(CH for protein production) (see 72-5)
        betab=-(S4b[1,2]+S4b[1,4]+S4b[2,2])/hb   #(CO2 from (Rpr-Rd) + CO2 from (RN-Rd))/(CH for protein production) (see 72-5)
        #XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
    
        #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
        #3.Calculation to obtain E2 (E when [O2]in>[O2]cri
        #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
        dd=(4*pi*R*Do2*(O2-O2cri))/(alpb*lmax*r1b*r3b)     
        E2=dd*((CHc+Kch)/(CHc))+betab/alpb    #E when [O2]in>[O2]cri
   
        #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
        #4.Chemostat computation
        #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
        
        #OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
        #NH4+ and N2 mixed case
        #OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
        
        O2in=-(lsgrow*r1b*r2b*r3b)/(4*pi*Do2*R)+O2     #(mol/m3): O2 concentration in the cell given there is not Ocri limit (72-6)
        
        
        if O2in<O2cri:
            E=E5
            Ep[i]=E5
        else:
            E=E2
            Ep[i]=E2
        CH[i]=(lsgrow*(1+E+F)+m)/(4*pi*R*Dch)+CHc            #(mol/m3): CH concentration in the environment (146-21)
        xnh4[i]=(D*(Nin-N))/(Rv*lsgrow*r4*(1-f))              #(cell m-3): number density of cell NH4 limited in NH4 and N2 mix (157-37)
        res2grow[i]=lsgrow*r1b*(alpb*E-betab)*r3b       #(mol O2 s-1 cell-1) respiration when dinitrogen and ammonium mixed uptake
        res2[i]=res2grow[i]*Rv
        Ncb[i]=(KN*Q*Dgrow*(1-f))/(lmax-Q*Dgrow*(1-f))                  #(mol/m3): NH4+ concentration in the cell (158-33) (r4 is cancenlled here (158-33)
        N2[i]=(Dgrow*Q*r4*(1-f))/(4*pi*R*Dnh4)+Ncb[i]             #(mol/me): NH4+ concentration in the vessel (100-15)
        
        #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
        #For stack plot of CH consumption (178-1~2)(181-13) Based on 601_00_03
        #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
        CCHbiomass_b[i]=ls/(-S4b[1,8])*S4b[1,0]                #(molC cell-1 s-1) consumption of CH for biomsas production
        CCHN2fix_b[i]=ls/(-S4b[1,8])*S4b[2,0]                  #(molC cell-1 s-1) consumption of CH for N2 fixation
        CCHbalancedresp_b[i]=ls/(-S4b[1,8])*S4b[0,0]           #(molC cell-1 s-1) consumption of CH for balanced respiration
        
        N2fixrespratio=S4[0,0]*f1/(S4[0,0]*f1+S4a[0,0]*(1-f1))      #(dimensionless) respiration for N2 fixing part (including biomass+N2 fixation) for total respiration
        CCHrespN2fixtot=CCHbalancedresp_b[i]*N2fixrespratio     #(molC cell-1 s-1) consumption of CH for N2 fix part of respiraiton (including biomass + N2 fixation)
        CCHrespN2fix_b[i]=CCHrespN2fixtot/(fn00*dgn+fs00*(dgp/ep1+dgpc/ep))*(fn00*dgn)   #(molC cell-1 s-1)
                                                        #Consumption of CH for respiration for N2 fixation
        CCHrespbiomass_b[i]=CCHbalancedresp_b[i]-CCHrespN2fix_b[i]    #(molC cell-1 s-1) Respiration biomass
                                                        #Consumption of CH for respiration for biomass production                                              
        CCHextraresp_b[i]=res2[i]/r3b-CCHbalancedresp_b[i]    #(molC cell-1 s-1) Carbohydrate consumption for extra respiration
    
        #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        #loop end $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
        #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    
    #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    #Optimization
    #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@



    EE=copy(CN)/copy(CN)
    x=EE*xn
    xch1=copy(EE)
    res=copy(EE)*res1           #(mol m-3) respiration rate
    CHp=copy(EE)*CHpot          #(mol m-3) CH concentration in the vessel for plot
    NH4=copy(EE)                #(mol m-3) NH4 concentration in the vessel for plot
    N1=copy(EE)
    fp=copy(EE-EE)
    CCHbiomass=copy(EE-EE)
    CCHN2fix=copy(EE-EE)
    CCHrespbiomass=copy(EE-EE)
    CCHrespN2fix=copy(EE-EE)
    CCHextraresp=copy(EE-EE)
    
    whichone=copy(EE)
    Epp=copy(EE)
    
    center1=copy(EE-EE)
    
    for j in U:
        xch1[j]=(D*(CHin[j]-CH0))/((lsgrow*(1+E3+F)+m)*Rv)  
        N1[j]=(D*Nin-xch1[j]*Rv*lsgrow*r4)/D         #(mol m-3) ammonium concentration in the vessel assuming (158-33)
        
        if O2cpot[j]>O2cri:             #When there is not enough carbon to reduc oxygten concentration when only all ammonium is used
            if xch1[j]>xn:              
                x[j]=xn
                res[j]=res1[j]
                CHp[j]=CHpot[j]
                NH4[j]=N
                whichone[j]=1
                CCHbiomass[j]=CCHbiomass_a
                CCHN2fix[j]=0
                CCHrespbiomass[j]=CCHbalancedresp_a
                CCHrespN2fix[j]=0
                CCHextraresp[j]=CCHextraresp_a[j]
                
            else:
                x[j]=xch1[j]
                res[j]=res0
                CHp[j]=CH0
                NH4[j]=N1[j]
                whichone[j]=2
                CCHbiomass[j]=CCHbiomass_a
                CCHN2fix[j]=0
                CCHrespbiomass[j]=CCHbalancedresp_a
                CCHrespN2fix[j]=0
                CCHextraresp[j]=0
                
        else:                           #this place is to choose the right value when break without going to else part. when there is no extra respiration, res0 and CH0 should be chosen
            if xch1[j]>xn:              #this part is to make sure that the right value is chosen when O2c<O2critical when no extra respiration
                x[j]=xn                 # this case is when extra respiration is done
                res[j]=res1[j] 
                CHp[j]=CHpot[j]
                NH4[j]=N 
                whichone[j]=3
                CCHbiomass[j]=CCHbiomass_a
                CCHN2fix[j]=0
                CCHrespbiomass[j]=CCHbalancedresp_a
                CCHrespN2fix[j]=0
                CCHextraresp[j]=CCHextraresp_a[j]
                
            else:
                x[j]=xch1[j]            #this case is there is no extra respiration is done
                res[j]=res0             #when there is extra carbon, nitrogen fixation kicks in
                CHp[j]=CH0
                NH4[j]=N1[j]
                whichone[j]=4
                CCHbiomass[j]=CCHbiomass_a
                CCHN2fix[j]=0
                CCHrespbiomass[j]=CCHbalancedresp_a
                CCHrespN2fix[j]=0
                CCHextraresp[j]=0
                                
            xch0=D*(CHin[j]-CH[0])/(4*pi*R*Dch*(CH[0]-CHc))/Rv     #(cell/m3): number density of the cells
            if xnh4[0]<xch0:
                while (imax-imin>1):
                    center=(imax+imin)//2
                    #xch=D*(CHin[j]-CH[i])/(4*pi*R*Dch*(CH[i]-CHc))/Rv     #(cell/m3): number density of the cells
                    xch=D*(CHin[j]-CH[center])/(4*pi*R*Dch*(CH[center]-CHc))/Rv     #(cell/m3): number density of the cells
                    xnc=xch-xnh4[center]
                
                    x[j]=xch
                    res[j]=res2[center]
                    CHp[j]=CH[center]
                    NH4[j]=N2[center]
                    fp[j]=far[center]
                    whichone[j]=5
                    Epp[j]=Ep[center]
                    CCHbiomass[j]=CCHbiomass_b[center]
                    CCHN2fix[j]=CCHN2fix_b[center]
                    CCHrespbiomass[j]=CCHrespbiomass_b[center]
                    CCHrespN2fix[j]=CCHrespN2fix_b[center]
                    CCHextraresp[j]=CCHextraresp_b[center]
                   
                    #print(center)
                    if xnc>0:
                        imin=center
                    else:
                        imax=center
            
            center1[j]=center
            #print(center)
            imax=size(floop)-1
            imin=floop[0]
            
                
#                 
#                 
#                 #if xnh4[i]>xch:
#                 if xnc<0:    
#                     break
#                 else:

    #+++++++++++++++++++++++++++++++++++++++++++++++++
    #print(Epp)
    #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    #Calculation of parameters for plotting
    #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    #OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    #Protein and CH
    #OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    
    Vchgrow=4*pi*R*Dch*(CH-CHc)   #(mol/s/cell): CH consumption speed by one growing cell
    Vch=Vchgrow*Rv                 #(mol/s/cell): CH consumption speed by one cell in average
    M=1000*60*60/(0.22*(BB/n)/12*10**6)                  #modifier to get Cs from Vch (105-12)
    Cs0=Vch*M/V                             #(mmolC/(h*gbiomass)): CH consumption rate per biomass
    Cs00=Cs0/pr                             #(mmolC/(h*g protein): CH consumption rate per protein
    Cs=Cs00/12                              #(mmolSucrose/(h*g protein): Sucrose consumption rate per protein 
    PRd0=PRm*x                              #(gPR/m3): Protein density in the vessel
    PRd=PRd0*10**(-3)                       #(mgPR/ml): Protein density in the vessel
    CHs=copy(CHp)/12                #(mol m-3) sucrose

    #OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    #preparation, yield, and N2 fixation for plot
    #OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    resp=res*10**6*60/Q*n/BB/pr/1000
    Yp=lsgrow/(4*pi*R*Dch*(CHp-CHc))*12*BB/n      #(g biomass / mol sucrose) Yield coefficient (158-15) *"Rv values are cancelled"
    N2fix=D*Q*r4*60*fp/Q*n/BB/1000/pr*10**9    #(nmol NH4+/(min*mgprotein))Nitrogen fixation rate convertedto ammonium (see figure 2 of buhler 1987b (160-16)
    fp_to_N2fix=D*Q*r4*60/Q*n/BB/1000/pr*10**9     #(nmol NH4+/(min*mgprotein)) fp to N2fix conversion term
    
    
    #--------------------------------------------------
    #Unit conversion part: For stack plot
    #--------------------------------------------------
    M=1000*60*60/(0.22*(BB/n)/12*10**6)                  #modifier to get Cs from Vch (105-12) 
    conversion1=M/V/pr/12                   #((mmolSucrose/(h*g protein))/(molC cell-1 s-1) Unit conversion term
    CCHbiomass_plot=CCHbiomass*conversion1   #(mmolSucrose/(h*g protein)) CCH for biomass production including respiration for energy
    CCHN2fix_plot=CCHN2fix*conversion1   #(mmolSucrose/(h*g protein)) CCH for N2 fixation including respiration for energy
    CCHrespbiomass_plot=CCHrespbiomass*conversion1     #(mmolSucrose/(h*g protein)) CCH for respiration balanced by biomass production
    CCHrespN2fix_plot=CCHrespN2fix*conversion1      #(mmolSucrose/(h*g protein)) CCH for respiration balanced by N2 fixation
   # CCHbalancedresp_plot=CCHbalancedresp*conversion1     #(mmolSucrose/(h*g protein)) CCH for balanced respiration
    CCHextraresp_plot=CCHextraresp*conversion1   #(mmolSucrose/(h*g protein)) Carbohydrate consumption for extra respiration

    
    return N2fix,PRd,resp,Yp
##########################################################################################################################################
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#6.execution
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

Ninstep=0.1
Ninmax=30
Ninarray=arange(Ninstep,Ninmax+Ninstep,Ninstep) #(molSucrose molN-1) C to N ratio in the incoming medium

for a in Ninarray:
    #start = time.time()   #current time in seconds
    N2fix,PR,Resp,Yp=az2(a)
    #end = time.time()       #current time in seconds
    if a==Ninstep:
        N2fixbox=N2fix
        PRbox=PR
        Respbox=Resp
        Ypbox=Yp
    else:
        N2fixbox=vstack((N2fixbox,N2fix))
        PRbox=vstack((PRbox,PR))
        Respbox=vstack((Respbox,Resp))
        Ypbox=vstack((Ypbox,Yp))
#OOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
#Plot settingff
#OOOOOOOOOOOOOOOOOOOOOOOOOOOOOO    
rcParams.update({'font.size': 24,
                 'lines.markersize':10,
                 'lines.markeredgewidth':0.5})
rcParams.update({'xtick.major.pad': 15})
rcParams.update({'xtick.major.pad': 15})
rcParams.update({'font.serif': 'Times New Roman'})
rcParams.update({'figure.autolayout': True})
rcParams['figure.figsize']=8,6.5
rcParams.update({'figure.facecolor':'W'})
rcParams.update({'lines.linewidth':2.5})
rcParams.update({'text.usetex': True})   #to call real latex
rcParams.update({'text.latex.preamble': ['\\usepackage[greek,english]{babel}']})
rcParams['axes.color_cycle'] = ['blue', 'green', 'red', 'cyan']

#OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
#Plot
#OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
O2pint=int(O2percent)       #Preparing for the title
O2p=str(O2pint)             #Preparing for the title

figure('FIG S1')
pcolormesh(CHin,Ninarray,PRbox,vmax=2,vmin=0)
xlabel('Sucrose (mol C m$^{-3}$)')
ylabel('NH$_4^+$ (mol N m$^{-3}$)')
ylim(ymax=30.00001)
title(O2p+'$\%$ O$_2$')
colorbar().set_label('Protein (mg/ml)')
# savetxt("2DN2fix2.csv", PRbox, delimiter=",",fmt='%2.8f')
# savefig(First_part+Plot2+Last_part) 

show()
    
 
