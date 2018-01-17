'''
Created on 2013/01/05


#----------------------------------------------------------------------------------------------
*About 320
Here I plot the result from 319_03_02

#----------------------------------------------------------------------------------------------
#note
To save time
#----------------------------------------------------------------------------------------------
@author: ag105020
'''
from pylab import *
rc('text', usetex=True)
rc('font', family='serif')
rc('font', serif='Times New Roman')

O2p='30'

CHin=genfromtxt('CHin'+O2p+'.csv',delimiter=',')
Ninarray=genfromtxt('Ninarray'+O2p+'.csv',delimiter=',')
N2fix=genfromtxt('2DN2fix'+O2p+'.csv',delimiter=',')
N2fixV=genfromtxt('2DN2fixV'+O2p+'.csv',delimiter=',')
B12=genfromtxt('B12array'+O2p+'.csv',delimiter=',')
B23=genfromtxt('B23array'+O2p+'.csv',delimiter=',')


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
rcParams.update({'lines.linewidth':7})
rcParams.update({'text.usetex': True})   #to call real latex
rcParams.update({'text.latex.preamble': ['\\usepackage[greek,english]{babel}']})
rcParams['axes.color_cycle'] = ['blue', 'green', 'red', 'cyan']

#OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO

figure(1)
pcolormesh(CHin,Ninarray,N2fix,vmax=32,vmin=0)
xlabel('Sucrose (mol C m$^{-3}$)')
ylabel('NH$_4^+$ (mol N m$^{-3}$)')
title(O2p+'$\%$ O$_2$')
colorbar().set_label('N$_2$ fixation (nmol N min$^{-1}$ mg protein$^{-1}$)')
plot(B12,Ninarray,color='#FF99FF')
plot(B23,Ninarray,color='w')
#savefig("C:\\Users\\Keisuke\\Desktop\\fig1_"+O2p+".png")

figure(5)
pcolormesh(CHin,Ninarray,N2fixV,vmax=50,vmin=0)
xlabel('Sucrose (mol C m$^{-3}$)')
ylabel('NH$_4^+$ (mol N m$^{-3}$)')
title(O2p+'$\%$ O$_2$')
colorbar().set_label('N$_2$ fixation (nmol N min$^{-1}$ ml$^{-1}$)')
plot(B12,Ninarray,color='#FF99FF')
plot(B23,Ninarray,color='w')
#savefig("C:\\Users\\Keisuke\\Desktop\\fig5_"+O2p+".png")


show()