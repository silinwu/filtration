###################################################################################################################################################################################################
#This is for DATA analysis in the filtration simulation test.
#Built by Silin WU from Hohai University, Nanjing, China. Email: wusilinhhu@126.com
#There are some parameters need to be added into every .py before run .py, list are:
#1 line 42, the initial Hydraulic conductivity of the filter
#2 line 48, the initial Hydraulic conductivity of the system
#3 line 44, the initial Permittivity of the filter
#4 line 52-70, the basic parameters in every test
#5 line 72, the Path of your hdf5 document
#6 line 73, the number of your hdf5 documents
#7 line 193,195, the name of you output excel
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
sheet1 = book.add_sheet('filtration')                                                        #a sheet in excel named by cakelength 
sheet1.write(0,0,'t(s)')                                                                     #write in the cakelength sheet of excel in row 1 and column 1
sheet1.write(1,0,0) 
sheet1.write(0,1,'Cakelength(um)')                                                           #write in row 1 and column 2
sheet1.write(1,1,0) 
sheet1.write(0,2,'Cakelength(m)')                                                            #write in row 1 and column 3
sheet1.write(1,2,0) 
sheet1.write(0,3,'Pressuredrop(Pa)')                                                         #write in row 1 and column 4
sheet1.write(1,3,0) 
sheet1.write(0,4,'Pressuredrop(m)')                                                          #write in row 1 and column 5
sheet1.write(1,4,0) 
sheet1.write(0,5,'Discharge(m3/s)')                                                          #write in row 1 and column 6
sheet1.write(1,5,0) 
sheet1.write(0,6,'Hydraulic conductivity of the cake(m/s)')                                  #write in row 1 and column 7
sheet1.write(1,6,np.nan) 
sheet1.write(0,7,'Hydraulic conductivity of the filter(m/s)')                                #write in row 1 and column 8
sheet1.write(1,7,np.nan)                                                                     #important : the initial Hydraulic conductivity of the filter need to be tested 
sheet1.write(0,8,'Permittivity of the filter(m/s)')                                          #write in row 1 and column 9
sheet1.write(1,8,np.nan)                                                                     #important : the initial Permittivity of the filter need to be tested 
sheet1.write(0,9,'Porosity of the cake')
sheet1.write(1,9,np.nan)
sheet1.write(0,10,'Hydraulic conductivity of the system(m/s)')                               #write in row 1 and column 10
sheet1.write(1,10,np.nan)                                                                    #important : the initial Hydraulic conductivity of the system is equal to the filter
sheet1.write(0,11,'Total filtration volumn(m3)')                                             #write in row 1 and column 11
sheet1.write(1,11,0)                                                                         
#########basic parameters#############################################################################################################################################################################
R = 10
ratio = 2
RR = ratio*R
Pnx = 9 
Pny = 3
pnx = 50
pny = 30
l  = 11
pdx = 38.4
pdy = 38.4
fnx = 7
fny = 15
lll = 100
nfloc = 2
RRR = (((nfloc*2)/3)*3**0.5+1)*R                                              # 3
#RRR = (nfloc*2**0.5+1)*R                                                     # 4
#RRR = (nfloc*2+1)*R                                                          # 6

ratiol = 1e-6
ratiot = 1e-7  
ratiom = 1e-15

ny = ((Pny-1)*(3**0.5*RR+pdy) + 40.0+2*RR + fny*lll)+1
sy = (Pny-1)*(3**0.5*RR+pdy) + 2*RR + 0.5*lll + 40
dtout = 500
cakebottom = int(round(sy - 0.5*lll)) + 1                                                 #the boundary bewteen cake and filter
filterlength = (Pny-1)*(3**0.5*RR+pdy) + 2*RR                                                 #the length of the filter
Filtrationvolumnr = 0
###################read hdf5 ###########################################################################################################################################################################
#prefix = "/home/foggy/dlvo-master/testflocfltitration__"
prefix = "/media/zhuwei/GroupZWTeam/simulation/new/F1/dlvo-master-1.0/testflocfltitrationF1__"
numberhdf5 = 150
for data in range(numberhdf5 - 1):
	data +=1
	name = prefix +str(data).zfill(4) +".h5"   
	f = h5.File(name)
	Np = f['Np'][0]                                                                           #Np:total particles number in Nx *Ny
	Nx = f['Nx'][0]
	Ny = f['Ny'][0]
	isf = np.array(f['PIsFree'])                                                              #-1 represent the filter particles, 1 represent the slurry ones
	Velocity = np.array(f['Velocity_0'])
	PVelocity = np.array(f['PVeloc'])
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
	VPy = PVelocity[1:3*Np-1:3]                                                              #the velocity of particles in y axis
##########caculate the reality time(s) and output to the excel built before###############################################################################################################################
	Timer = dtout*data*ratiot			                                                      #output the reality time(s) 
	sheet1.write(data+1,0,Timer)
###########caculate the length of the cake, the pressure drop of the cake and output to the excel built before##########################################################################################
	cakelengthmaxnum = int((ny-(sy-0.5*lll))/(2*R))                                           #the max piece number of the cake during the simulation test
	cakelength = 0                                                                            #initialize the cakelength
	cakelengthr = 0
	Density = Density.reshape((Ny,Nx))                                                        #the density distribution
	Densitycakebottom = Density[cakebottom-1,:]                                               #the density cake at bottom value  [list]
	Densitycakebottom = np.mean(Densitycakebottom)                                            #the mean density cake at bottom value 
	Pressuredropcaker = 0
	for n in range(cakelengthmaxnum):
		cakelengthBi = sy + 2*R*n                                                             #the length of the first cakelengthB and cakelength is the first piece of the cake
		cakelengthTi = sy + 2*R*(n+1)
		N = 0
		for i in range(Np):
			if posy[i] >= cakelengthBi-0.2*R and posy[i] <= cakelengthTi+0.2*R and abs(VPy[i]) <= 0.01:               
				N +=1
		if N <= fnx * ((RRR/R)-1):                                                            #if the number of the y position of the particles between cakelenthB and cakelenthT is less than pnx +3, then define the length of the cake
			cakelength = n*2*R
			cakelengthr = cakelength*ratiol                                                   #the real cake length(m)
			Densitycaketop = Density[cakebottom + cakelength,:]                               #the density cake at top value  [list]
			Densitycaketop = np.mean(Densitycaketop) 
			Pressuredrop = (Densitycaketop - Densitycakebottom)/3
			Pressuredropcaker = Pressuredrop*ratiom/(ratiol*ratiot*ratiot)                    #the reality pressure drop(pa)
			sheet1.write(data+1,1,cakelengthr*1000000)                                        #output the reality cake length(um)
			sheet1.write(data+1,2,cakelengthr)                                                #output the reality cake length(m)
			sheet1.write(data+1,3,Pressuredropcaker)                                          #output the reality pressure drop(pa)
			sheet1.write(data+1,4,Pressuredropcaker/10000)                                    #output the reality pressure drop(m)
###########caculate the Porosity of the cake and output to the excel built before###################################################################################################################################	
			cakelengthB = cakebottom		                                                  #the bottom of the cake
			cakelengthT = cakebottom + cakelength                                             #the top of the cake
			if cakelengthT == cakelengthB:
				sheet1.write(data+1,9,np.nan)
				break
			particleinsidecakearea = 0                                                        #initialize the total area of the particles inside the cake
			for i in range(Np):	
				if (cakelengthT - posy[i]) >= R and (posy[i] - cakelengthB) >= R:             #particles inside the cake
					particleinsidecakearea += (math.pi)*R*R
				elif (cakelengthT - posy[i]) > R and (posy[i] - cakelengthB) < R and (posy[i] - cakelengthB) >= 0:   #particles at the bottom of the cake, and the centre of the circle is inside the cake
					Dlength = (posy[i] - cakelengthB)
					particleinsidecakearea += (math.pi)*R*R - (R*R*math.acos(Dlength/R)-Dlength*((R*R-Dlength*Dlength)**0.5))
				elif (cakelengthT - posy[i]) < R and (posy[i] - cakelengthB) > R and (cakelengthT - posy[i]) >=0:    #particles at the top of the cake, and the centre of the circle is inside the cake
					Dlength = (cakelengthT - posy[i])
					particleinsidecakearea += (math.pi)*R*R - (R*R*math.acos(Dlength/R)-Dlength*((R*R-Dlength*Dlength)**0.5))
				elif (posy[i] - cakelengthT) < R and (posy[i] - cakelengthT) > 0:
					Dlength = (posy[i] - cakelengthT)
					particleinsidecakearea += R*R*math.acos(Dlength/R)-Dlength*((R*R-Dlength*Dlength)**0.5)
				elif (cakelengthB - posy[i]) < R and (cakelengthB - posy[i]) > 0:
					Dlength = (cakelengthB - posy[i])
					particleinsidecakearea += R*R*math.acos(Dlength/R)-Dlength*((R*R-Dlength*Dlength)**0.5)	
				elif (cakelengthT - posy[i]) < R and (posy[i] - cakelengthB) < R:
					print("error: the wrong particle inside the cake")
				else:
					particleinsidecakearea += 0	
			cakearea = cakelength * Nx
			Porosity = (cakearea-particleinsidecakearea)/cakearea
			sheet1.write(data+1,9,Porosity)
			break
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
	sheet1.write(data+1,5,Discharge2nr)                                                       #output the reality discharge(m3/s)
###########caculate the total filtration volumn during the filtration and output to the excel built before############################################################################################################
	Filtrationvolumnr += dtout*ratiot*Discharge2nr
	sheet1.write(data+1,11,Filtrationvolumnr)
###########caculate the hydaulic conductivity of cake(m/s) and output to the excel built before######################################################################################################################
	if cakelengthr ==0:
		cakehcr = np.nan
	else:
		cakehcr = Discharge2nr*cakelengthr/(Arear*Pressuredropcaker/10000)                    #the Darcy law; the Pressure dropr(pa) was converted to the Unit(m),divived by 10000
	sheet1.write(data+1,6,cakehcr)                                                            #output the reality cake haudrulic conductivity (m/s)
##########caculate the hydraulic conductivity of filter(m/s), the permittivity(s-1) and output to the excel built before############################################################################################################### 	
	Densityfiltertop = Densitycakebottom                                                      #the mean density filter at top value, which is equal to the mean density cake at bottom value
	Densityfilterbottom = Density[40+1-1,:]                                                   #the density filter at bottom value  [list]
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
############caculate the hydraulic conductivity of the system(m/s) and output to the excel built before, the system contains the cake and filter##########################################################################################
	if cakelengthr == 0:
		systemhcr = filterhcr
	else:
		Pressuredropsystemr = Pressuredropfilterr + Pressuredropcaker                         #the reality pressure drop of the system(pa)
		systemlengthr = filterlengthr + cakelengthr                                           #the reality system length(m)
		systemhcr = Discharge2nr*systemlengthr/(Arear*Pressuredropsystemr/10000)              #the reality hydraulic conductivity of the system(m/s)
	sheet1.write(data+1,10,systemhcr)
book.save('test2.xls')                                                                        #save the data in the excel
######### read the excel data that output last step and plot the figures#################################################################################################################################
wb = xlrd.open_workbook('test2.xls')                                                          #open the excel
table = wb.sheet_by_name(u'filtration')                                                       #open the sheet
cols1 = table.col_values(0)                                                                   #the data of t(s)
cols2 = table.col_values(1)                                                                   #the data of cakelength(um)
cols3 = table.col_values(2)                                                                   #the data of cakelength(m)
cols4 = table.col_values(3)                                                                   #the data of Pressuredrop(pa)
cols5 = table.col_values(4)                                                                   #the data of Pressuredrop(m)
cols6 = table.col_values(5)                                                                   #the data of discharge (m3/s)
cols7 = table.col_values(6)                 												  #the data of cake hydraulic conductivity(m/s)
cols8 = table.col_values(7)                 												  #the data of Hydraulic conductivity of the filter(m/s)
cols9 = table.col_values(8)                 												  #the data of Permittivity of the filter(m/s)
cols10 = table.col_values(9)                 												  #the data of Porosity of the cake
cols11 = table.col_values(10)                 												  #the data of Hydraulic conductivity of the system(m/s)
cols12 = table.col_values(11)                 												  #the data of Total volumn during the filtration(m3)
del cols1[0]                                 											      #delete the first element
del cols2[0]                                										          #delete the first element
del cols3[0]                                 												  #delete the first element
del cols4[0]                                											      #delete the first element
del cols5[0]                                 												  #delete the first element
del cols6[0]                                 												  #delete the first element
del cols7[0]                                 												  #delete the first element
del cols8[0]                                 												  #delete the first element
del cols9[0]                                 												  #delete the first element
del cols10[0]                                 												  #delete the first element
del cols11[0]                                                                                 #delete the first element
del cols12[0]                                                                                 #delete the first element
#plot t(s) vs cake,filter and system hydraulic conductivity(m/s) in 331
ax = plt.subplot(331)
plt.xlabel("reality time t(s)")
plt.ylabel("reality hydraulic conductivity(m/s)")
plt.plot(cols1,cols7,label="cake",alpha = 0.6)
plt.plot(cols1,cols8,label="filter",alpha = 0.6)
plt.plot(cols1,cols11,label="system",alpha = 0.6)
plt.legend(loc='upper left', bbox_to_anchor=(0.2,0.95))
ax.set_yscale('log')
#plot t(s) vs cakelength(um) in 332
plt.subplot(332)
plt.xlabel("reality time t(s)")
plt.ylabel("reality cake length(um)")
plt.plot(cols1,cols2,alpha = 0.6)
#plot t(s) vs Porosity of the cake in 333
plt.subplot(333)
plt.xlabel("reality time t(s)")
plt.ylabel("Porosity of the cake")
plt.plot(cols1,cols10,alpha = 0.6)
#plot t(s) vs Hydraulic conductivity of the filter(m/s) in 334
plt.subplot(334)
plt.xlabel("reality time t(s)")
plt.ylabel("Hydraulic conductivity of the filter(m/s)")
plt.plot(cols1,cols8,alpha = 0.6)
#plot t(s) vs discharge(m3/s) in 335
plt.subplot(335)
plt.xlabel("reality time t(s)")
plt.ylabel("reality discharge(m3/s)")
plt.plot(cols1,cols6,alpha = 0.6)
#plot t(s) vs Hydraulic conductivity of the cake(m/s) in 336
plt.subplot(336)
plt.xlabel("reality time t(s)")
plt.ylabel("Hydraulic conductivity of the cake(m/s)")
plt.plot(cols1,cols7,alpha = 0.6)
#plot t(s) vs Total volumn during filtration(m3) in 337
plt.subplot(337)
plt.xlabel("reality time t(s)")
plt.ylabel("Total volumn during filtration(m3)")
plt.plot(cols1,cols12,alpha = 0.6)
plt.show()





###########plot the position of the particle at the last hdf5##################################################################################################################################################################################
fig = plt.figure()
bx = fig.add_subplot(111)  # draw the subgraph, "1",1,1 represent the position of the subgraph; add_subplot(111) is same as add_subplot(1,1,1)

for i in range(Np):
	if(isf[i]<0):
		cir = Circle(xy=(posx[i],posy[i]),radius=RR,alpha=0.2)   #read the line 4; Circle(position,radius,transparent)
	else:
		cir = Circle(xy=(posx[i],posy[i]),radius=R,alpha=0.2)
	bx.add_patch(cir)
	bx.text(posx[i]-5,posy[i]-2.5,str(i))         # text i at position (posx[i],posx[i])
#ax.quiver(posx,posy,pfcx,pfcy,scale=2)         # the quivers
plt.axis('equal')
plt.hlines(cakebottom,0,180)
#plt.hlines(21,0,180)
#plt.vlines(0,0,350)
#plt.vlines(180,0,350)
plt.show()
###########################################################################################################################################################################################################