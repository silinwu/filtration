###################################################################################################################################################################################################
#This is for DATA analysis in the filtration simulation test.
#Built by Silin WU from Hohai University, Nanjing, China. Email: wusilinhhu@126.com
#BE HAPPY WITH THIS DATA ANALYSIS.
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
sheet1 = book.add_sheet('test')                                                       
sheet1.write(0,0,'t(s)')                                                                    
sheet1.write(1,0,0) 
sheet1.write(0,1,'Discharge(m3/s)')                                                         
sheet1.write(1,1,0) 
sheet1.write(0,2,'Hydraulic conductivity(m/s)')                                
sheet1.write(1,2,np.nan) 
#########basic parameters#############################################################################################################################################################################
nx = 500
ny = 200
dtout = 1000                                                              #time step ***************
ratiol = 1e-6
ratiot = 1e-7  
ratiom = 1e-15
###############################################################################################################
###################read hdf5 ###########################################################################################################################################################################
prefix = "/media/zhuwei/My Book/filtration/liuquan/NO4/NO.4_"
numberhdf5 = 199                                                            #********************************
filternumber = 0								 							#it is used for caculating the filter porosity
for data in range(numberhdf5 - 1):
	data +=1
	name = prefix +str(data).zfill(4) +".h5"   
	f = h5.File(name)
	Np = f['Np'][0]                                                                           #Np:total particles number in Nx *Ny
	Nx = f['Nx'][0]
	Ny = f['Ny'][0]
	Velocity = np.array(f['Velocity_0'])
	Density = np.array(f['Density_0'])
	pos = np.array(f['Pposition'])
	posx = pos[0:3*Np-2:3]                                                                    #[start:end:step]
	posy = pos[1:3*Np-1:3]
##########caculate the reality time(s) and output to the excel built before###############################################################################################################################
	Timer = dtout*data*ratiot			                                                      #output the reality time(s) 
	sheet1.write(data+1,0,Timer)
###########density######################################################################################################################################################
	Density = Density.reshape((Ny,Nx))                                                   
	Densitytop1 = Density[48+100,:]                                          
	Densitytop1 = np.mean(Densitytop1)                                         
	Densitytop2 = Density[49+100,:]                                            
	Densitytop2 = np.mean(Densitytop2)                                        
	Densitytop3 = Density[50+100,:]                                          
	Densitytop3 = np.mean(Densitytop3)                                        
	Densitytop  = (Densitytop1 + Densitytop2 + Densitytop3)/3 
	Densitybottom1 = Density[48,:]                                                     #the density filter at bottom value  [list]
	Densitybottom1 = np.mean(Densitybottom1)                                        #the mean density filter at bottom value
	Densitybottom2 = Density[49,:]
	Densitybottom2 = np.mean(Densitybottom2)                                        #the mean density filter at bottom value
	Densitybottom3 = Density[50,:]
	Densitybottom3 = np.mean(Densitybottom3)                                        #the mean density filter at bottom value
	Densitybottom  = (Densitybottom1 + Densitybottom2 + Densitybottom3)/3
###########caculate the discharge(m3/s) during the test and output to the excel built before############################################################################################################
	Nfluent = Nx * Ny
	Area = Nx *1                                                                              #the filtration Area in lattice simulation
	Arear = Area * ratiol * ratiol                                                            #the reality filtration Area(m2)
	Vy = Velocity[1:3*Nfluent:3]                                                              #the velocity distribution in y axis
	Vy = Vy.reshape((Ny,Nx))
	Vy2n = Vy[1,:]                                                                            #the bottom + 1 value of vy
	Density2n = Density[1,:]                                                                  #the bottom + 1 value of density [list]
	Discharge2n = abs(Vy2n.dot(Density2n))/(np.sum(Density[1,:]))                             #the bottom +1 mean velocity of the discharge///////the another calculate method is to use np.mean(Vy1,:), but not precise than this one
	Discharge2nr = Discharge2n*Area*ratiol*ratiol*ratiol/ratiot                               #the reality velocity of the discharge(m3/s)
	sheet1.write(data+1,1,Discharge2nr)                                                       #output the reality discharge(m3/s)
##########caculate the hydraulic conductivity (m/s) and output to the excel built before############################################################################################################### 	
	flowvel = Discharge2n                                                                     #the bottom + 1 value of vy
	flowvelr = flowvel * ratiol/ratiot                                                        #the reality +1 value of vy for filter(m/s)
	lengthr = 100 * ratiol                                                     #the reality filter length (m)
	Pressuredrop = (Densitytop - Densitybottom)/3                           #the pressure drop of the filter 
	Pressuredropr = Pressuredrop*ratiom/(ratiol*ratiot*ratiot)                    #the reality pressure drop of the filter(pa)
	hydraulicconductivityr = flowvelr * lengthr/(Pressuredropr/10000)                          #the vertical hydraulic conductivity(m/s)
	sheet1.write(data+1,2,hydraulicconductivityr)                                                          #output the reality filter haudrulic conductivity (m/s)
book.save('hydraulicconductivity NO.4.xls')                                                               #save the data in the excel **************
######### read the excel data that output last step and plot the figures#################################################################################################################################
wb = xlrd.open_workbook('hydraulicconductivity NO.4.xls')                                                 #open the excel   ************
table = wb.sheet_by_name(u'test')                                                       #open the sheet
cols1 = table.col_values(0)                                                                   
cols2 = table.col_values(1)                                                                   
cols3 = table.col_values(2)                                                                   
del cols1[0]                                 											      #delete the first element
del cols2[0]                                										          #delete the first element
del cols3[0]                                 												  #delete the first element
#plot t(s) vs hydraulic conductivity(m/s) in 331
ax = plt.subplot(121)
plt.xlabel("reality time t(s)")
plt.ylabel("reality hydraulic conductivity(m/s)")
plt.plot(cols1,cols3,label="hydraulic conductivity",alpha = 0.6)
plt.legend(loc='upper left', bbox_to_anchor=(0.2,0.95))
#plot t(s) vs cakelength(um) in 332
plt.subplot(122)
plt.xlabel("reality time t(s)")
plt.ylabel("reality discharge(m3/s)")
plt.plot(cols1,cols2,alpha = 0.6)
plt.show()