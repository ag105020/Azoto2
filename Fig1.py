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
rc('text', usetex=True)
rc('font', family='serif')
rc('font', serif='Times New Roman')


##########################################################################################################################################
def az2(buhler_data):
    #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    #1.Parameters----------------------------
    #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    O2=buhler_data.O2
    plotcolor=buhler_data.plotcolor
    #OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    #defining basic parameters
    #OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO

    CNstep=0.1
    CNwidth=15
    CN=arange(CNstep,CNwidth+CNstep,CNstep) #(molSucrose molN-1) C to N ratio in the incoming medium
    CNnumber=size(CN)
    U=arange(0,CNnumber,1)
    
    #OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    #Making fitting plot (Based on Buhler 1987, and "15 Data of Nitrogen fixation rate.xlsx"
    #OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    Vmax=33.4
    
    if buhler_data.i==1:
        CNfirst=2.5
        CNhalf=7
    elif buhler_data.i==2:
        CNfirst=4
        CNhalf=10
    else:
        CNfirst=1
        CNhalf=2
       
    N2fix=Vmax*(CN-CNfirst)/(CNhalf+CN-CNfirst*2)
    
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    #5.Plot
    #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    
    #OOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    #Importing experimental data 
    #OOOOOOOOOOOOOOOOOOOOOOOOOOOOOO

    #OOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    #Setting font
    #OOOOOOOOOOOOOOOOOOOOOOOOOOOOOO    
#     rcParams.update({'font.size': 22})
#     rcParams.update({'lines.markersize': 10})
#     rcParams.update({'lines.markeredgewidth': 0.5})
# #    rcParams.update({'lines.markersize': 10})
# #    rcParams.update({'lines.linewidth': 10})
# #    rcParams.update({'marker.size': 10})
    rcParams.update({'font.size': 26,
                     'lines.markersize':10,
                     'lines.markeredgewidth':3.0})
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

    #OOOOOOOOO
    #new plot
    #OOOOOOOOO
    #===============
    #Plot variables
    #===============
    Figlocation="C:\\Users\\Keisuke\\Dropbox (MIT)\\0 Paper\\2 Azotobacter 2\\Diagrams\\02 From 317_03_22_02\\"
    Xlimmax=CNwidth+1e-17
    Folderloc="Azoto2"
   #OOOOOOOOO
    #new plot
    #OOOOOOOOO
    xname='C/N (mol sucrose/mol NH$_4^+$)'
    ##########################
    #OOOOOOOOOOOOOOOOOOOOOOOOO

    
    
     
    if buhler_data.i==1 or buhler_data.i==2:
        fig=figure('FIG 1')
        ax1=fig.add_subplot(111)
      #  ax2=ax1.twinx()
        #ax1.plot(CN, N2fix, linewidth=2.5, color=plotcolor,label='O$_2$=' +str(int(O2/0.225*100))+' ($\%$)')
        ax1.plot(buhler_data.cnNf1, buhler_data.Nf1d, '+', color=plotcolor)
        ax1.plot(buhler_data.cnNf2, buhler_data.Nf2d, 'x', color=plotcolor)
        ax1.set_xlabel(xname)
        ax1.set_ylabel('N$_2$ fixation rate (nmol N min$^{-1}$ mg protein$^{-1}$)',fontsize=22)
        ax1.set_ylim(0,35)
      #  ax2.set_ylim(0,35/fp_to_N2fix*100)
      #  ax2.set_ylabel('N$_2$ fixation ratio ($\%$)')
        xlim(xmax=Xlimmax)
        #ylim(ymax=100)

        legend(loc=2,borderaxespad=0.8, fontsize=25)
    
    
    
    ###################
    #OOOOOOOOOOOOOOOOOO
#     Figurenumber=19
#     StackPlotColors1=('red','cyan')
#     figure(Figurenumber+buhler_data.i)
#     stackplot(CN,fp,1-fp,colors=StackPlotColors1)
#     xlabel(xname)
#     xlim(xmax=Xlimmax)
#     ylabel('N$_2$ : NH$_4^+$')
#     ylim(ymax=1)
#     Savefig(Folderloc,Figurenumber+buhler_data.i)

    #OOOOOOOOOOOOOOOOOO
    ################### 
    
   
    return 
    ##########################################################################################################################################
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#6.execution
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

buhler_data=buhler87()
for a in buhler_data:
    az2(a)
    
show()
    
    
 
