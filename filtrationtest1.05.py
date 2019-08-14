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
sheet1.write(0,12,'Cake water content')                                                      #write in row 1 and column 12
sheet1.write(1,12,np.nan)  
sheet1.write(0,13,'DE')                                                      				 #write in row 1 and column 13
sheet1.write(1,13,0)
sheet1.write(0,14,'cakerouT')                                                      				 
sheet1.write(1,14,np.nan)   
sheet1.write(0,15,'cakerouB')                                                      				 
sheet1.write(1,15,np.nan)
sheet1.write(0,16,'LOSS')												 				     #particles loss 
sheet1.write(1,16,0)
sheet1.write(0,17,'LOSS(particles/m2)')												 		 #particles loss 
sheet1.write(1,17,0)
sheet1.write(0,18,'Porosity of the filter')
sheet1.write(1,18,np.nan)
#########basic parameters#############################################################################################################################################################################
#######filter##################################################################################################
RR = 20                                                                  # *****************
Pnx = 4                                                                  # *****************
Pny = 3                                                                  # *****************
pdx = 60                                                                 # *****************
pdy = 60                                                                 # *****************
filterlength = (Pny-1)*(3**0.5*RR+pdy) + 2*RR                            #the length of the filter
#######particle################################################################################################
R = 10																	 # *****************
pnx = 10                                                        		 #used when particle filtration  *****************
pny = 5																 #used when particle filtration  *****************
l  = 19.831																 #used when particle filtration  *****************
ny = (Pny-1)*(3**0.5*RR+pdy) + 40.0 + 2*RR + pny*2*R + pny*l + 1   	 #used when particle filtration
sy = (Pny-1)*(3**0.5*RR+pdy) + 40.0 + 2*RR + 0.5*l + R  + 1			 #used when particle filtration
cakebottom = int(round(sy -(0.5*l + R )))								 #used when particle filtration
watercontentinitial = (((l+2*R)*(l+2*R)-3.14*R*R)/(3.14*R*R*2.7))       #used when particle filtration
#######floc####################################################################################################
#R = 10																	 # *****************
#shape = 3																 # *****************
#nfloc = 3																 # *****************
#fnx = 5                                                                  # *****************
#fny = 9                                                                 # *****************
#lll = 126.5                                                                # *****************
#ny = ((Pny-1)*(3**0.5*RR+pdy) + 40.0+2*RR + fny*lll)+1    
#sy = (Pny-1)*(3**0.5*RR+pdy) + 2*RR + 0.5*lll + 40 +1                    #sy is used to find the first floc
#cakebottom = int(round(sy - 0.5*lll))	                                 #the boundary bewteen cake and filter, the bottom of the cake, or the top of the filter
#if shape == 3:
#	RRR = (((nfloc*2)/3)*3**0.5+1)*R
#	watercontentinitial = ((lll*lll-3.14*R*R*((nfloc+2)*(nfloc+1)/2))/(3.14*R*R*((nfloc+2)*(nfloc+1)/2)*2.7))
#elif shape == 4:
#	watercontentinitial = ((lll*lll-3.14*R*R*(nfloc+1)*(nfloc+1))/((3.14*R*R*(nfloc+1)*(nfloc+1))*2.7))
#	RRR = (nfloc*2**0.5+1)*R 
#else:
#	watercontentinitial = ((lll*lll-3.14*R*R*(3*(1+nfloc)*nfloc+1))/(3.14*R*R*(3*(1+nfloc)*nfloc+1)*2.7))
#	RRR = (nfloc*2+1)*R 
#######environment#############################################################################################
dtout = 500                                                              #time step ***************
ratiol = 1e-6
ratiot = 1e-7  
ratiom = 1e-15
Filtrationvolumnr = 0                                                    #initialize the volumn in the filtration
cakelengthmaxnum = int((ny-cakebottom)/(2*R))                            #the max piece number of the cake during the simulation test
###############################################################################################################

###################read hdf5 ###########################################################################################################################################################################
#prefix = "/home/foggy/dlvo-master/test-AN6_"
prefix = "/media/zhuwei/My Book/filtration/0808paper/DN7/test-DN7_"                      #********************************
numberhdf5 = 573                                                            #********************************
filternumber = 0								 							#it is used for caculating the filter porosity
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
	VPy = PVelocity[1:3*Np-1:3]                                                               #the velocity of particles in y axis
	
	if filternumber == 0:
		for i in range(Np):
			if(isf[i]<0):
				filternumber +=1	

##########caculate the reality time(s) and output to the excel built before###############################################################################################################################
	Timer = dtout*data*ratiot			                                                      #output the reality time(s) 
	sheet1.write(data+1,0,Timer)
###########caculate the length of the cake, the pressure drop of the cake and output to the excel built before##########################################################################################
	cakelength = 0                                                                            #initialize the cakelength
	cakelengthr = 0
	Density = Density.reshape((Ny,Nx))                                                        #the density distribution
	Densitycakebottom1 = Density[cakebottom - 1-1,:]                                          #the density cake at bottom value  [list]
	Densitycakebottom1 = np.mean(Densitycakebottom1)                                          #the mean density cake at bottom value 
	Densitycakebottom2 = Density[cakebottom - 1,:]                                            #the density cake at bottom value  [list]
	Densitycakebottom2 = np.mean(Densitycakebottom2)                                          #the mean density cake at bottom value
	Densitycakebottom3 = Density[cakebottom - 1+1,:]                                          #the density cake at bottom value  [list]
	Densitycakebottom3 = np.mean(Densitycakebottom3)                                          #the mean density cake at bottom value
	Densitycakebottom  = (Densitycakebottom1 + Densitycakebottom2 + Densitycakebottom3)/3 
	sheet1.write(data+1,15,Densitycakebottom)
	Pressuredropcaker = 0
	cakewatercontent = 0
	for n in range(cakelengthmaxnum):
		cakelengthBi = sy + 2*R*n                                                             #the length of the first cakelengthB and cakelength is the first piece of the cake
		cakelengthTi = sy + 2*R*(n+1)
		N = 0
		for i in range(Np):
			if posy[i] >= cakelengthBi-0.2*R and posy[i] <= cakelengthTi+0.2*R and abs(VPy[i]) <= 0.02:    # we add a condition that (VPy[i]) <= 0.01, which is mean when meet this condition, the veclosity of the particle is very low, this is convient for testting the cake formed by floc, 0.01 is an experience point.        
				N +=1
		#if N <= fnx * int((RRR/R)):                          #*******floc***************      #if the number of the y position of the particles between cakelenthB and cakelenthT is less than pnx +3, then define the length of the cake
		if N <= pnx + 3:                                      #**************************      #if the number of the y position of the particles between cakelenthB and cakelenthT is less than pnx +3, then define the length of the cake
			cakelength = n*2*R
			cakelengthr = cakelength*ratiol                                                   #the real cake length(m)
			Densitycaketop1 = Density[cakebottom - 1 + int(cakelength)-1,:]               	  #the density cake at top value  [list], 0.8 is an experience value, for the data is more stable
			Densitycaketop1 = np.mean(Densitycaketop1) 	  
			Densitycaketop2 = Density[cakebottom - 1 + int(cakelength),:]                 	  #the density cake at top value  [list], 0.8 is an experience value, for the data is more stable
			Densitycaketop2 = np.mean(Densitycaketop2)	  
			Densitycaketop3 = Density[cakebottom - 1 + int(cakelength)+1,:]               	  #the density cake at top value  [list], 0.8 is an experience value, for the data is more stable
			Densitycaketop3 = np.mean(Densitycaketop3)
			Densitycaketop  = (Densitycaketop1 + Densitycaketop2 + Densitycaketop3)/3
			sheet1.write(data+1,14,Densitycaketop)
			Pressuredrop = (Densitycaketop - Densitycakebottom)/3
			Pressuredropcaker = Pressuredrop*ratiom/(ratiol*ratiot*ratiot)                     #the reality pressure drop(pa)
			sheet1.write(data+1,1,cakelengthr*1000000)                                         #output the reality cake length(um)
			sheet1.write(data+1,2,cakelengthr)                                                 #output the reality cake length(m)
			sheet1.write(data+1,3,Pressuredropcaker)                                           #output the reality pressure drop(pa)
			sheet1.write(data+1,4,Pressuredropcaker/10000)                                     #output the reality pressure drop(m)
###########caculate the Porosity of the cake and output to the excel built before###################################################################################################################################	
			cakelengthB = cakebottom		                                                   #the bottom of the cake
			cakelengthT = cakebottom + cakelength                                              #the top of the cake
			if cakelengthT == cakelengthB:
				sheet1.write(data+1,9,np.nan)
				sheet1.write(data+1,12,np.nan)
				break
			particleinsidecakearea = 0                                                         #initialize the total area of the particles inside the cake
			for i in range(Np):	 
				if (cakelengthT - posy[i]) >= R and (posy[i] - cakelengthB) >= R:              #particles inside the cake
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
			cakewatercontent = ((cakearea-particleinsidecakearea)*1)/(particleinsidecakearea*2.7)
			sheet1.write(data+1,9,Porosity)
			sheet1.write(data+1,12,cakewatercontent)
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
	Densityfilterbottom1 = Density[41-1-1,:]                                                     #the density filter at bottom value  [list]
	Densityfilterbottom1 = np.mean(Densityfilterbottom1)                                        #the mean density filter at bottom value
	Densityfilterbottom2 = Density[41-1,:]
	Densityfilterbottom2 = np.mean(Densityfilterbottom2)                                        #the mean density filter at bottom value
	Densityfilterbottom3 = Density[41-1+1,:]
	Densityfilterbottom3 = np.mean(Densityfilterbottom3)                                        #the mean density filter at bottom value
	Densityfilterbottom  = (Densityfilterbottom1 + Densityfilterbottom2 + Densityfilterbottom3)/3
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
############caculate the DE##########################################################################################################################################################################################################
	M = 0                                               								      #M is used to number the particle inside the cake and filter
	DE = 0                                                                                    #DE is dewatering effiency
	for i in range(Np):
		if posy[i] >= 41 and posy[i] <= cakebottom + cakelength:
			M+=1
	watercontentfinal = max(cakewatercontent, ((Np-M)/Np)*watercontentinitial)
	DE = ((1/(1+watercontentfinal))-(1/(1+watercontentinitial)))/(1/(1+watercontentinitial))
	sheet1.write(data+1,13,DE)
###########caculate the loss particles######################################################################################################################################################
	LOSS = 0
	for i in range(Np):
		if posy[i] == -100:
			LOSS +=1
	LOSSperarea = LOSS / Arear 
	sheet1.write(data+1,16,LOSS)
	sheet1.write(data+1,17,LOSSperarea)
###########caculate the porosity of the filter##########################################################################################################################################################
	particleinsidefilterarea = 0                                                          #initialize the total area of the particles inside the filter
	for i in range(filternumber,Np):	 
		if (cakebottom - posy[i]) >= R and (posy[i] - 41) >= R:                           #particles inside the filter
			particleinsidefilterarea += (math.pi)*R*R
		elif (cakebottom - posy[i]) > R and (posy[i] - 41) < R and (posy[i] - 41) >= 0:   #particles at the bottom of the filter, and the centre of the circle is inside the filter
			Dlength = (posy[i] - 41)
			particleinsidefilterarea += (math.pi)*R*R - (R*R*math.acos(Dlength/R)-Dlength*((R*R-Dlength*Dlength)**0.5))
		elif (cakebottom - posy[i]) < R and (posy[i] - 41) > R and (cakebottom - posy[i]) >=0:    #particles at the top of the filter, and the centre of the circle is inside the filter
			Dlength = (cakebottom - posy[i])
			particleinsidefilterarea += (math.pi)*R*R - (R*R*math.acos(Dlength/R)-Dlength*((R*R-Dlength*Dlength)**0.5))
		elif (posy[i] - cakebottom) < R and (posy[i] - cakebottom) > 0:
			Dlength = (posy[i] - cakebottom)
			particleinsidefilterarea += R*R*math.acos(Dlength/R)-Dlength*((R*R-Dlength*Dlength)**0.5)
		elif (41 - posy[i]) < R and (41 - posy[i]) > 0:
			Dlength = (41 - posy[i])
			particleinsidefilterarea += R*R*math.acos(Dlength/R)-Dlength*((R*R-Dlength*Dlength)**0.5)	
		elif (cakebottom - posy[i]) < R and (posy[i] - 41) < R:
			print("error: the wrong particle inside the filter")
		else:
			particleinsidefilterarea += 0
	filterporosity = (Nx*(cakebottom-41) - filternumber*(math.pi)*RR*RR- particleinsidefilterarea)/(Nx*(cakebottom-41))
	sheet1.write(data+1,18,filterporosity)
book.save('testfiltrationDN7.xls')                                                               #save the data in the excel **************
######### read the excel data that output last step and plot the figures#################################################################################################################################
wb = xlrd.open_workbook('testfiltrationDN7.xls')                                                 #open the excel   ************
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
cols13 = table.col_values(12)																  #cake water content
cols14 = table.col_values(13)                                                                 #DE
cols17 = table.col_values(16)                                                                 #LOSS
cols18 = table.col_values(17)                                                                 #LOSS(particles/m2)
cols19 = table.col_values(18)                                                                 #Porosity of the filter

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
del cols13[0]                                                                                 #delete the first element
del cols14[0]                                                                                 #delete the first element
del cols17[0]                                                                                 #delete the first element
del cols18[0]                                                                                 #delete the first element
del cols19[0]                                                                                 #delete the first element

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
#plot t(s) vs Porosity of the filter in 334
plt.subplot(334)
plt.xlabel("reality time t(s)")
plt.ylabel("Porosity of the filter")
plt.plot(cols1,cols19,alpha = 0.6)
#plot t(s) vs discharge(m3/s) in 335
plt.subplot(335)
plt.xlabel("reality time t(s)")
plt.ylabel("reality discharge(m3/s)")
plt.plot(cols1,cols6,alpha = 0.6)
#plot t(s) vs LOSS in 336
plt.subplot(336)
plt.xlabel("reality time t(s)")
plt.ylabel("LOSS")
plt.plot(cols1,cols17,alpha = 0.6)
#plot t(s) vs Total volumn during filtration(m3) in 337
plt.subplot(337)
plt.xlabel("reality time t(s)")
plt.ylabel("Total volumn during filtration(m3)")
plt.plot(cols1,cols12,alpha = 0.6)
#plot t(s) vs cake water content during filtration in 338
plt.subplot(338)
plt.xlabel("reality time t(s)")
plt.ylabel("cake water content during filtration(m3)")
plt.plot(cols1,cols13,alpha = 0.6)
#plot t(s) vs LOSS(particles/m2) during filtration in 338
plt.subplot(339)
plt.xlabel("reality time t(s)")
plt.ylabel("LOSS(particles/m2)")
plt.plot(cols1,cols18,alpha = 0.6)
plt.show()

############plot the position of the particle at the last hdf5##################################################################################################################################################################################
#fig = plt.figure()
#bx = fig.add_subplot(111)  # draw the subgraph, "1",1,1 represent the position of the subgraph; add_subplot(111) is same as add_subplot(1,1,1)
#
#for i in range(Np):
#	if(isf[i]<0):
#		cir = Circle(xy=(posx[i],posy[i]),radius=RR,alpha=0.2)   #read the line 4; Circle(position,radius,transparent)
#	else:
#		cir = Circle(xy=(posx[i],posy[i]),radius=R,alpha=0.2)
#	bx.add_patch(cir)
#	#bx.text(posx[i]-5,posy[i]-2.5,str(i))         # text i at position (posx[i],posx[i])
##ax.quiver(posx,posy,pfcx,pfcy,scale=2)         # the quivers
#plt.axis('equal')
#plt.hlines(cakebottom,0,180)
##plt.hlines(21,0,180)
##plt.vlines(0,0,350)
##plt.vlines(180,0,350)
#plt.show()
###########################################################################################################################################################################################################