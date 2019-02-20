# Code Policy
Please state “Cell Flux Model” and “Keisuke Inomura” in the acknowledgement when your
publication includes the results based on the original/revised code. Or you may consider
including Keisuke Inomura as a co-author depending on the contribution. In either case, the
publication must cite the following paper:
Inomura K, Bragg J, Riemann L Follows MJ. (2018). A quantitative model of nitrogen fixation in the presence of ammonium.
PLOS ONE: e0208282
Keisuke Inomura (University of Washington)
ki24@uw.edu

# Azoto2
For plotting, download all the files and put them in the same folder. The codes are in python 3 and they require Latex and pylab. 

1. Plotting Fig 1: 
Run Fig1.py

2. Plotting Fig 3: 
Run Fig3.py

3. Plotting Fig 4:
Run Fig4.py

4. Plotting Fig 6: 
Run Fig6_calculation.py and then run Fig6_plot.py. To change oxygen concentration, change values of O2satratio in Fig6_calculation.py and O2p in Fig6_plot.py. O2p must be 100 x O2satratio. To change resolution, alter values in tune CNstep and Ninstep in Fig6_calculation.py. 

5. Plotting Fig 7:
Run Fig7.py

6. Plotting S2 Fig:
Run S2_Fig.py

7. Plotting S3 Fig: 
Run S3_Fig.py To change oxygen concentration, change a value of O2satratio. Resolution can be changed in the same way as Fig 6.

