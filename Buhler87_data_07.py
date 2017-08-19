'''
Created on May 18, 2014

@author: Keisuke
'''
from pylab import * 

class Buhler87:    
    def __init__(self,O2percentage,plotcolor,cnpr,prd,cnch,chd,cnre,red,cnY,Yd,cnNf1,cnNf2,Nf1d,Nf2d,imax,i):
        O2saturation=0.225
        self.O2percentage=O2percentage
        self.O2=O2saturation*O2percentage/100
        self.plotcolor=plotcolor
        self.cnpr=cnpr
        self.prd=prd
        self.cnch=cnch
        self.chd=chd
        self.cnre=cnre
        self.red=red
        self.cnY=cnY
        self.Yd=Yd
        self.cnNf1=cnNf1
        self.cnNf2=cnNf2
        self.Nf1d=Nf1d
        self.Nf2d=Nf2d
        self.imax=imax
        self.i=i

def buhler87():    

    a=genfromtxt('From_buhler_1987_2_1.csv', delimiter=',')
    nodata=zeros(a.shape[0])+nan
    O2percentage=[5,15,30,60]
    plotcolor=['red','green','blue','cyan']
    imax=3
    i=arange(0,imax+1,1)
    
    
    buhler_05=Buhler87(O2percentage[0],plotcolor[0],a[:,1],a[:,2],a[:,14],a[:,15],a[:,24],a[:,25],a[:,37],a[:,38],a[:,49],a[:,46],a[:,50],a[:,47],imax,i[0])
    buhler_15=Buhler87(O2percentage[1],plotcolor[1],a[:,4],a[:,5],a[:,17],a[:,18],a[:,27],a[:,28],a[:,40],a[:,41],a[:,55],a[:,52],a[:,56],a[:,53],imax,i[1])
    buhler_30=Buhler87(O2percentage[2],plotcolor[2],a[:,7],a[:,8],nodata,nodata,a[:,30],a[:,31],nodata,nodata,a[:,61],a[:,58],a[:,62],a[:,59],imax,i[2])
    buhler_60=Buhler87(O2percentage[3],plotcolor[3],a[:,10],a[:,11],a[:,20],a[:,21],a[:,33],a[:,34],a[:,43],a[:,44],nodata,a[:,64],nodata,a[:,65],imax,i[3])

    return(buhler_05,buhler_15,buhler_30,buhler_60) 
