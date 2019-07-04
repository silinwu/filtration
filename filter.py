###################################################################################################################################################################################################
#This is for DATA analysis in the filtration simulation test.
#Special for the initialize the filter.
#After testing testfilter.cpp, using filter.py to get the initial permeability and hydraulitic coeffient of the filter.
#It is strange that the permeability and hydraulitic coeffient of the filter will change with the pressure drop, but not too much. 
#So in testing the permeability of filter, we fix the pressure drop with 5 kPa.
#Built by Silin WU from Hohai University, Nanjing, China. Email: wusilinhhu@126.com
#BE HAPPY WITH THIS DATA ANALYSIS.   -----data update:2019-07-04
###################################################################################################################################################################################################

#####Extension libraries used in this .py##########################################################################################################################################################
import numpy as np
import h5py as h5
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
import sys
import xlwt
import xlrd
import math
########output for the excel ######################################################################################################################################################################
book = xlwt.Workbook()                                                                       #build an excel 
sheet1 = book.add_sheet('filter')                                                            #a sheet in excel named by cakelength 
sheet1.write(0,0,'t(s)')                                                                     #write in the cakelength sheet of excel in row 1 and column 1
sheet1.write(1,0,0) 
sheet1.write(0,5,'Discharge(m3/s)')                                                          #write in row 1 and column 6
sheet1.write(1,5,0) 
sheet1.write(0,7,'Hydraulic conductivity of the filter(m/s)')                                #write in row 1 and column 8
sheet1.write(1,7,np.nan)                                                                     #important : the initial Hydraulic conductivity of the filter need to be tested 
sheet1.write(0,8,'Permittivity of the filter(m/s)')                                          #write in row 1 and column 9
sheet1.write(1,8,np.nan)                                                                     #important : the initial Permittivity of the filter need to be tested                                                                      
#########basic parameters#############################################################################################################################################################################
R =10
RR = 25
Pnx = 1
Pny = 3
pdx = 60
pdy = 60
pny = 3
l  = 3
ratiol = 1e-6
ratiot = 1e-7  
ratiom = 1e-15
ny = ((Pny-1)*(3**0.5*RR+pdy) + 20.0+2*RR + pny*2*R+pny*l)+1
dtout = 1e2
filtertop = int(round((Pny-1)*(3**0.5*RR+pdy) + 2*RR + 20 + 1))                               #the position of the top of filter
filterbottom = 20 + 1                          												  #the position of the bottom of filter
filterlength = (Pny-1)*(3**0.5*RR+pdy) + 2*RR                                                 #the length of the filter
###################read hdf5 ###########################################################################################################################################################################
prefix = "/home/foggy/dlvo-master/testfilter_"
numberhdf5 = 70
for data in range(numberhdf5 - 1):
	data +=1
	name = prefix +str(data).zfill(4) +".h5"   
	f = h5.File(name)
	Np = f['Np'][0]                                                                           #Np:total particles number in Nx *Ny
	Nx = f['Nx'][0]
	Ny = f['Ny'][0]
	isf = np.array(f['PIsFree'])                                                              #-1 represent the filter particles, 1 represent the slurry ones
	Velocity = np.array(f['Velocity_0'])
	Density = np.array(f['Density_0'])
	pos = np.array(f['Pposition'])
	pfc = np.array(f['PForce'])
	pfh = np.array(f['PForceh'])
	posx = pos[0:3*Np-2:3]                                                                    #[start:end:step]
	posy = pos[1:3*Np-1:3]
	pfcx = pfc[0:3*Np-2:3]
	pfcy = pfc[1:3*Np-1:3]
	pfhx = pfh[0:3*Np-2:3]
	pfhy = pfh[1:3*Np-1:3]
###########caculate the reality time(s) and output to the excel built before###############################################################################################################################
	Timer = dtout*data*ratiot			                                                      #output the reality time(s) 
	sheet1.write(data+1,0,Timer)
###########caculate the length of the cake, the pressure drop of the cake and output to the excel built before##########################################################################################
	Density = Density.reshape((Ny,Nx))                                                        #the density distribution
	Densityfiltertop = Density[filtertop-1,:]                                                #the density cake at bottom value  [list]
	Densityfiltertop = np.mean(Densityfiltertop)                                             #the mean density cake at bottom value 
###########caculate the discharge(m3/s) during the filtration and output to the excel built before############################################################################################################
	Nfluent = Nx * Ny
	Area = Nx *1                                                                              #the filtration Area in lattice simulation
	Arear = Area * ratiol * ratiol                                                            #the reality filtration Area(m2)
	Vy = Velocity[1:3*Nfluent:3]                                                              #the velocity distribution in y axis
	Vy = Vy.reshape((Ny,Nx))
	Vy2n = Vy[1,:]                                                                            #the bottom + 1 value of vy
	Density2n = Density[1,:]                                                                  #the bottom + 1 value of density [list]
	Discharge2n = abs(Vy2n.dot(Density2n))/(np.sum(Density[1,:]))                             #the bottom +1 mean velocity of the discharge///////the another calculate method is to use np.mean(Vy1,:), but not precise than this one
	Discharge2nr = Discharge2n*Area*ratiol*ratiol*ratiol/ratiot                               #the reality velocity of the discharge(m3/s)
	sheet1.write(data+1,5,Discharge2nr)                                                       #output the reality discharge(m3/s)                                                     #output the reality cake haudrulic conductivity (m/s)
###########caculate the hydraulic conductivity of filter(m/s), the permittivity(s-1) and output to the excel built before############################################################################################################### 	
	Densityfilterbottom = Density[20+1-1,:]                                                   #the density filter at bottom value  [list]
	Densityfilterbottom = np.mean(Densityfilterbottom)                                        #the mean density filter at bottom value
	flowvel = Discharge2n                                                                     #the bottom + 1 value of vy
	flowvelr = flowvel * ratiol/ratiot                                                        #the reality +1 value of vy for filter(m/s)
	filterlengthr = filterlength * ratiol                                                     #the reality filter length (m)
	Pressuredropfilter = (Densityfiltertop - Densityfilterbottom)/3                           #the pressure drop of the filter 
	Pressuredropfilterr = Pressuredropfilter*ratiom/(ratiol*ratiot*ratiot)                    #the reality pressure drop of the filter(pa)
	filterhcr = flowvelr * filterlengthr/(Pressuredropfilterr/10000)                          #the vertical hydraulic conductivity(m/s)
	filterpmr = flowvelr/(Pressuredropfilterr/10000)                                          #the permittivity of the filter(s-1)
	sheet1.write(data+1,7,filterhcr)                                                          #output the reality filter haudrulic conductivity (m/s)
	sheet1.write(data+1,8,filterpmr)                                                          #output the reality filter permittivity (s-1)
book.save('testfilter.xls')                                                                   #save the data in the excel
######### read the excel data that output last step and plot the figures#################################################################################################################################
wb = xlrd.open_workbook('testfilter.xls')                                                     #open the excel
table = wb.sheet_by_name(u'filter')                                                           #open the sheet
cols1 = table.col_values(0)                                                                   #the data of t(s)
cols6 = table.col_values(5)                                                                   #the data of discharge (m3/s)
cols8 = table.col_values(7)                 												  #the data of Hydraulic conductivity of the filter(m/s)
cols9 = table.col_values(8)                 												  #the data of Permittivity of the filter(s-1)
del cols1[0]                                 											      #delete the first element
del cols6[0]                                 												  #delete the first element
del cols8[0]                                 												  #delete the first element
del cols9[0]                                 												  #delete the first element
#plot t(s) vs cake,filter and system hydraulic conductivity(m/s) in 331
ax = plt.subplot(131)
plt.xlabel("reality time t(s)")
plt.ylabel("reality Permittivity of the filter(s-1)")
plt.plot(cols1,cols9,label="filter",alpha = 0.6)
plt.legend(loc='upper left', bbox_to_anchor=(0.2,0.95))
#plot t(s) vs Hydraulic conductivity of the filter(m/s) in 334
plt.subplot(132)
plt.xlabel("reality time t(s)")
plt.ylabel("Hydraulic conductivity of the filter(m/s)")
plt.plot(cols1,cols8,alpha = 0.6)
#plot t(s) vs discharge(m3/s) in 335
plt.subplot(133)
plt.xlabel("reality time t(s)")
plt.ylabel("reality discharge(m3/s)")
plt.plot(cols1,cols6,alpha = 0.6)
plt.show()
############plot the position of the particle at the last hdf5##################################################################################################################################################################################
fig = plt.figure()
bx = fig.add_subplot(111)  # draw the subgraph, "1",1,1 represent the position of the subgraph; add_subplot(111) is same as add_subplot(1,1,1)
for i in range(Np):
	if(isf[i]<0):
		cir = Circle(xy=(posx[i],posy[i]),radius=RR,alpha=0.2)   #read the line 4; Circle(position,radius,transparent)
	else:
		cir = Circle(xy=(posx[i],posy[i]),radius=R,alpha=0.2)
	bx.add_patch(cir)
#	ax.text(posx[i]-5,posy[i]-2.5,str(i))         # text i at position (posx[i],posx[i])
#ax.quiver(posx,posy,pfcx,pfcy,scale=2)         # the quivers
plt.axis('equal')
plt.hlines(filtertop,0,180)
plt.hlines(21,0,180)
#plt.vlines(0,0,350)
#plt.vlines(180,0,350)
plt.show()
##########################################################################################################################################################################################################