#!/usr/bin/env python
#######################################################################
#    This code is used to generate the real gas property tables based #
#    on the REFPROP data base.                                        #
#                                                                     #
#    Author :   Jianhui Qi        j.qi@uq.edu.au                      #
#    Advisor:   Ingo H. J. Jahn   i.jahn@uq.edu.au                    #
#    Date   :   27-05-2017                                            #
#    version:   1.                                                    #
#######################################################################

# Reference: https://www.cfd-online.com/Forums/openfoam/141821-tabulated-thermophysicalproperties-library.html

import numpy as np 
from pyRefpropMania import *
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy.interpolate import griddata
########################################
# Actaviate the fluid parameter, with REFEPROP:

fluidType  =   'MD4M'                     # Fluid name
myFluid    =   refpropFluid(fluidType)   # Actaviate the REFPROP data base
R          =   8.3144598                 # [J/(K mol)] --> SI unit [kg m2 s-2 K-1 mol-1]
mW         =   458.99328e-3#236.531e-3#44.01e-3                  # [kg/mol] molecure weight 
########################################


#propName  =  'D'  # D: density, P: pressure, H: enthalpy, V: viscosity, 
########################################
# Define the limitation for plotting:
Tmin       =  610. # 300.0                  #[K] The lower limit of the calculation temperature
Tmax       =  670. # 1500.0                #[K] The upper limit of the calculation temperature
            
pmin       =   0.3e6                       #[Pa] The lower limit of the pressure need to be compared
pmax       =   1.7e6                      #[Pa] The upper limit of the pressure need to be compared

#rhomin     =  200.#5                      # [kg/m^3] The lower limit of the desity
#rhomax     =  300.#280                     # [kg/m^3] The upper limit of the pressure

#emin       =  2.96e5#4e5                     # [J/kg]  The lower limit of the internal energy
#emax       =  3.05e5#1.5e6                     # [J/kg]  The upper limit of the internal energy

#hmin       =  3.02e5#2.44e5                  # [J/kg] The lower limit of the enthalpy
#hmax       =  3.05e5#2.0e6                   # [J/kg] The upper limit of the enthalpy


nT         =   100.                     # How many data points for temperature
nP         =   100.                     # How many data points for pressure
nrho       =   100                    # How many data points for density
ne         =   100                     # How many data points for internal energy
nh         =   100                      # How many data points for enthalpy
########################################

# Generate the T, P list to compose the table.
T_list = np.arange(Tmin, Tmax + (Tmax - Tmin)/(nT) , (Tmax - Tmin)/(nT-1.))
p_list = np.arange(pmin, pmax + (pmax - pmin)/(nP) , (pmax - pmin)/(nP-1.))


#D_list = np.arange(rhomin,rhomax + (rhomax - rhomin)/(nrho), (rhomax - rhomin)/(nrho - 1.))

#U_list = np.arange(emin, emax + (emax - emin)/(ne), (emax - emin)/(ne - 1.))

#H_list = np.arange(hmin, hmax + (hmax - hmin)/(nh), (hmax - hmin)/(nh - 1.))
#print e_list


# choose real gas or ideal gas model.

gasEOS  =  "Real"# "Ideal"
table_flag = "ep" #"rhoe"
print "-------------------------------------------------"
print "Now doing calculations..."
print "-------------------------------------------------"


# Initialize all the properties lists.
rho_list     = []  # returning rho value based on T  and p value
mu_list      = []  # returning mu value baded on T and p value
kappa_list   = []  # returning kappa value based on T and p value
Cp_list      = []  # returning Cp value based on T and p value
Cv_list      = []  # returning Cv value based on T and p value
h_list       = []  # returning h value based on T and p value
cpMcv_list   = []  # returning cp minus cv value based on T and p value
A_list       = []  # returning acoustic speed(A) based on T and p value
e_list       = []  # returning e value based on T and p value
h_R_list     = []  # returning h value based on reverse order of T and p value
Thp_list     = []  # returning T based on h and p value
p2_list      = []
T2_list      = []

hrhoe_list    = []
prhoe_list    = []
#Trhoe_list    = []


if gasEOS == "Real":

    
    for i in range(len(T_list)):
        for j in range(len(p_list)):
            rho_list.append(    myFluid.getProps('D','T',T_list[i]-273.15,'P',p_list[j]/1000.)       )
            mu_list.append(     myFluid.getProps('V','T',T_list[i]-273.15,'P',p_list[j]/1000.)/1.e6  )  #[Pa*s] Or [kg/(m s)]
            kappa_list.append(  myFluid.getProps('K','T',T_list[i]-273.15,'P',p_list[j]/1000.)       )
            Cp_list.append(     myFluid.getProps('C','T',T_list[i]-273.15,'P',p_list[j]/1000.)*1000. )  #[J/(kg K)]
            Cv_list.append(     myFluid.getProps('O','T',T_list[i]-273.15,'P',p_list[j]/1000.)*1000. )  #[J/(kg K)]
            A_list.append(      myFluid.getProps('A','T',T_list[i]-273.15,'P',p_list[j]/1000.)       )  #[m/s]
            e_list.append(      myFluid.getProps('U','T',T_list[i]-273.15,'P',p_list[j]/1000.)*1000. )
            h_list.append(      myFluid.getProps('H','T',T_list[i]-273.15,'P',p_list[j]/1000.)*1000. )
            p2_list.append(      p_list[j]) # the p list is used to create p value for prhoe_table
            T2_list.append(      T_list[i]) # the T2 list is used to create T value for Trhoe_table   

    for i in range(len(p_list)):
        for j in range(len(T_list)):
            h_R_list.append( myFluid.getProps('H','T',T_list[j]-273.15,'P',p_list[i]/1000.)*1000. )  #[J/kg]

    for a in range(len(Cp_list)):
            cpMcv_list.append( Cp_list[a] - Cv_list[a])



    for i in range(len(mu_list)):
        if mu_list[i] < 0:
            print "WARNING: viscosity has a negative value"
            print "         1e-06 is given as substitute viscosity"
            mu_list[i] = 1e-06

    
#    for k in range(len(rho_list)):
#        for l in range(len(e_list)):
#            hrhoe_list.append(   myFluid.getProps('H','D',D_list[k],'U',U_list[l]/1000.)*1000)  #[J/kg]
#            prhoe_list.append(   myFluid.getProps('P','D',D_list[k],'U',U_list[l]/1000.)*1000)  #[Pa]


#    for i in range(len(H_list)):
#        for j in range(len(p_list)):
#           Thp_list.append(    myFluid.getProps('T','H',H_list[i]/1000.,'P',p_list[j]/1000.)+273.15) #[K]



else:
    # Ideal gas equation is listed here:

    mu_ideal = 1.0e-06#3.571e-5  # [Pa s]
    Cp_ideal = 5600.0742#1215.25   # [J/(kg K)]
    Ru       = R/mW
    Cv_ideal = Cp_ideal - Ru # [J/(kg K)]
    kappa_ideal= Cp_ideal/Cv_ideal
    cpMcv_ideal = Cp_ideal - Cv_ideal


    
    
    for i in range(len(T_list)):
        for j in range(len(p_list)):
            rho_list.append(  p_list[j]/(Ru*T_list[i]) )
            Cp_list.append(  Cp_ideal )
            mu_list.append(  mu_ideal  )
            kappa_list.append( kappa_ideal )
            cpMcv_list.append( cpMcv_ideal )
            A_list.append( np.sqrt(kappa_ideal*Ru*T_list[i]) )
            e_list.append( Cv_ideal*T_list[i])
            h_list.append(  Cp_ideal*T_list[i]  )
            p2_list.append(      p_list[j]) # the p list is used to create p value for prhoe_table
            T2_list.append(      T_list[i]) # the T2 list is used to create T value for Trhoe_table   


    for i in range(len(p_list)):
        for j in range(len(T_list)):
            h_R_list.append(  Cp_ideal*T_list[j]  ) # reversed h list




if  table_flag == "rhoe":

    print "Creating hrhoe prhoe and Trhoe table ..."

    # Find the limit rho and e:

    rhoMin = np.min(rho_list)
    rhoMax = np.max(rho_list)
    eMin   = np.min(e_list)
    eMax   = np.max(e_list)

    D_list = np.arange(rhoMin,rhoMax + (rhoMax - rhoMin)/(nrho), (rhoMax - rhoMin)/(nrho - 1.))

    U_list = np.arange(eMin, eMax + (eMax - eMin)/(ne), (eMax - eMin)/(ne - 1.)) 

    # Create interpolate table grid
    #grid_rho,grid_e = np.mgrid[rhoMin:rhoMax:100j,eMin:eMax:100j]
    #print D_list
    grid_rho,grid_e = np.meshgrid(D_list,U_list)
    #print grid_rho
    #print grid_e


    # find the avarage value for filling the empty value.
    av_h = np.average(h_list)
    av_p = np.average(p2_list)
    av_A = np.average(A_list)
    av_T = np.average(T2_list)

    # Get h value based on rho ane e value.
    points_rho_e = []
    for i in range(len(rho_list)):
        coordinateTemp = [rho_list[i],e_list[i]]
        points_rho_e.append(coordinateTemp)
    points_rho_e =  np.asarray(points_rho_e) #convert the list to array

    hrhoe_list = griddata(points_rho_e,h_list,(grid_rho,grid_e),method='cubic', fill_value=av_h)
    prhoe_list = griddata(points_rho_e,p2_list,(grid_rho,grid_e),method='cubic', fill_value=av_p)
    Arhoe_list = griddata(points_rho_e,A_list,(grid_rho,grid_e),method='cubic', fill_value=av_A)
    Trhoe_list = griddata(points_rho_e,T2_list,(grid_rho,grid_e),method='cubic', fill_value=av_T)


    #Show the figure
    #fig = plt.figure()
    #ax = fig.add_subplot(111,projection='3d')
    #ax.scatter(grid_rho,grid_e,Trhoe_list)
    #plt.grid()
    #plt.show()



    #=================================#
    print "Creating Thp table ..."


    hMin = np.min(h_list)
    hMax = np.max(h_list)

    print hMin, hMax

    H_list = np.arange(hMin, hMax + (hMax - hMin)/(nh), (hMax - hMin)/(nh - 1.))

    grid_h,grid_p = np.meshgrid(H_list,p2_list)


    points_h_p = []
    for i in range(len(h_list)):
        coordinateTemp = [h_list[i],p2_list[i]]
        points_h_p.append(coordinateTemp)
    points_h_p =  np.asarray(points_h_p) #convert the list to array
    Thp_list = griddata(points_h_p,T2_list,(grid_h,grid_p),method='cubic', fill_value=0)


    #Show the figure
    #fig = plt.figure()
    #ax = fig.add_subplot(111,projection='3d')
    #ax.scatter(grid_h,grid_p,Thp_list)
    #plt.grid()
    #plt.show()







    print 'Starting to write tables. \n'
    #--------------------------------
    print 'Writing densityTable ...'


    rhoF = open("densityTable",'w')
    rhoF.write("( \n")
    for a in range(len(T_list)):
        rhoF.write("(" + str(T_list[a]) + "\n(\n")
        for b in range(len(p_list)):
            sList = ["\t(" + str(p_list[b]) + " " + str(rho_list[int(a*nP+b)]) + ")\n"]
            rhoF.write(" ".join(sList))
        rhoF.write(") ) \n") 
    rhoF.write(");")
    rhoF.close()    






    #--------------------------------
    print 'Writing kappaTable ...'
    kappaF = open("kappaTable",'w')
    kappaF.write("( \n")
    for a in range(len(T_list)):
        kappaF.write("(" + str(T_list[a]) + "\n(\n")
        for b in range(len(p_list)):
            sList = ["\t(" + str(p_list[b]) + " " + str(kappa_list[int(a*nP+b)]) + ")\n"]
            kappaF.write(" ".join(sList))
        kappaF.write(") ) \n") 
    kappaF.write(");")
    kappaF.close()    





    #--------------------------------
    print 'Writing muTable ...'
    muF = open("muTable",'w')
    muF.write("( \n")
    for a in range(len(T_list)):
        muF.write("(" + str(T_list[a]) + "\n(\n")
        for b in range(len(p_list)):
            sList = ["\t(" + str(p_list[b]) + " " + str(mu_list[int(a*nP+b)]) + ")\n"]
            muF.write(" ".join(sList))
        muF.write(") ) \n") 
    muF.write(");")
    muF.close()    



    #--------------------------------
    print 'Writing cpTable ...'
    cpF = open("cpTable",'w')
    cpF.write("( \n")
    for a in range(len(T_list)):
        cpF.write("(" + str(T_list[a]) + "\n(\n")
        for b in range(len(p_list)):
            sList = ["\t(" + str(p_list[b]) + " " + str(Cp_list[int(a*nP+b)]) + ")\n"]
            cpF.write(" ".join(sList))
        cpF.write(") ) \n") 
    cpF.write(");")
    cpF.close()





    #--------------------------------
    #To be noted that the hTable used an invert p,T table.
    print 'Writing hTable ...'
    hF = open("hTable",'w')
    hF.write("( \n")
    for a in range(len(p_list)):
        hF.write("(" + str(p_list[a]) + "\n(\n")
        for b in range(len(T_list)):
            sList = ["\t(" + str(T_list[b]) + " " + str(h_R_list[int(a*nT+b)]) + ")\n"]
            hF.write(" ".join(sList))
        hF.write(") ) \n") 
    hF.write(");")
    hF.close()    



    #--------------------------------
    print 'Writing cpMcvTable ...'
    cpMcvF = open("cpMcvTable",'w')
    cpMcvF.write("( \n")
    for a in range(len(T_list)):
        cpMcvF.write("(" + str(T_list[a]) + "\n(\n")
        for b in range(len(p_list)):
            sList = ["\t(" + str(p_list[b]) + " " + str(cpMcv_list[int(a*nP+b)]) + ")\n"]
            cpMcvF.write(" ".join(sList))
        cpMcvF.write(") ) \n") 
    cpMcvF.write(");")
    cpMcvF.close()    




    #--------------------------------
    print 'Writing ATable ...'
    AF = open("ATable",'w')
    AF.write("( \n")
    for a in range(len(T_list)):
        AF.write("(" + str(T_list[a]) + "\n(\n")
        for b in range(len(p_list)):
            sList = ["\t(" + str(p_list[b]) + " " + str(A_list[int(a*nP+b)]) + ")\n"]
            AF.write(" ".join(sList))
        AF.write(") ) \n") 
    AF.write(");")
    AF.close()    




    #--------------------------------
    print 'Writing hrhoeTable ...'
    hrhoe = open("hrhoeTable",'w')
    hrhoe.write("( \n")
    for a in range(len(D_list)):
        hrhoe.write("(" + str(D_list[a]) + "\n(\n")
        for b in range(len(U_list)):
            sList = ["\t(" + str(U_list[b]) + " " + str(hrhoe_list[b][a]) + ")\n"]
            hrhoe.write(" ".join(sList))
        hrhoe.write(") ) \n") 
    hrhoe.write(");")
    hrhoe.close()    



    #--------------------------------
    print 'Writing prhoeTable ...'
    prhoe = open("prhoeTable",'w')
    prhoe.write("( \n")
    for a in range(len(D_list)):
        prhoe.write("(" + str(D_list[a]) + "\n(\n")
        for b in range(len(U_list)):
            sList = ["\t(" + str(U_list[b]) + " " + str(prhoe_list[b][a]) + ")\n"]
            prhoe.write(" ".join(sList))
        prhoe.write(") ) \n") 
    prhoe.write(");")
    prhoe.close()    


    #--------------------------------
    print 'Writing TrhoeTable ...'
    Trhoe = open("TrhoeTable",'w')
    Trhoe.write("( \n")
    for a in range(len(D_list)):
        Trhoe.write("(" + str(D_list[a]) + "\n(\n")
        for b in range(len(U_list)):
            sList = ["\t(" + str(U_list[b]) + " " + str(Trhoe_list[b][a]) + ")\n"]
            Trhoe.write(" ".join(sList))
        Trhoe.write(") ) \n") 
    Trhoe.write(");")
    Trhoe.close()   


    #--------------------------------
    print 'Writing ArhoeTable ...'
    Arhoe = open("ArhoeTable",'w')
    Arhoe.write("( \n")
    for a in range(len(D_list)):
        Arhoe.write("(" + str(D_list[a]) + "\n(\n")
        for b in range(len(U_list)):
            sList = ["\t(" + str(U_list[b]) + " " + str(Arhoe_list[b][a]) + ")\n"]
            Arhoe.write(" ".join(sList))
        Arhoe.write(") ) \n") 
    Arhoe.write(");")
    Arhoe.close()   




    #--------------------------------
    print 'Writing ThpTable ...'
    ThpF = open("ThpTable",'w')
    ThpF.write("( \n")
    for a in range(len(H_list)):
        ThpF.write("(" + str(H_list[a]) + "\n(\n")
        for b in range(len(p_list)):
            sList = ["\t(" + str(p_list[b]) + " " + str(Thp_list[b][a]) + ")\n"]
            ThpF.write(" ".join(sList))
        ThpF.write(") ) \n") 
    ThpF.write(");")
    ThpF.close()  


elif table_flag == "ep":

    # Find the limit rho and e:

    rhoMin = np.min(rho_list)
    rhoMax = np.max(rho_list)
    eMin   = np.min(e_list)
    eMax   = np.max(e_list)


    D_list = np.arange(rhoMin,rhoMax + (rhoMax - rhoMin)/(nrho), (rhoMax - rhoMin)/(nrho - 1.))

    U_list = np.arange(eMin, eMax + (eMax - eMin)/(ne), (eMax - eMin)/(ne - 1.)) 


    # Create interpolate table grid
    #grid_rho,grid_e = np.mgrid[rhoMin:rhoMax:100j,eMin:eMax:100j]
    #print D_list
    grid_e,grid_p = np.meshgrid(U_list,p_list)
    #print grid_rho
    #print "----",grid_e, grid_p





    # find the avarage value for filling the empty value.
    av_h = 0.#np.average(h_list)
    av_rho = 0# np.average(rho_list)
    av_A = 0. #np.average(A_list)
    av_T = 0. #np.average(T2_list)


    #print e_list
    # Get h value based on rho ane e value.
    points_e_p = []
    for i in range(len(e_list)):
        coordinateTemp = [e_list[i],p2_list[i]]
        points_e_p.append(coordinateTemp)
    points_e_p =  np.asarray(points_e_p) #convert the list to array

    #print U_list,'\n',p_list

   
    hep_list   = griddata(points_e_p,h_list,(grid_e,grid_p),method='cubic', fill_value=av_h)
    rhoep_list = griddata(points_e_p,rho_list,(grid_e,grid_p),method='cubic', fill_value=av_rho)
    Aep_list   = griddata(points_e_p,A_list,(grid_e,grid_p),method='cubic', fill_value=av_A)
    Tep_list   = griddata(points_e_p,T2_list,(grid_e,grid_p),method='cubic', fill_value=av_T)


    #Show the figure
    #fig = plt.figure()
    #ax = fig.add_subplot(111,projection='3d')
    #ax.scatter(grid_rho,grid_e,Trhoe_list)
    #plt.grid()
    #plt.show()

    print "Write prhoe table..."
    grid_rho,grid_e = np.meshgrid(D_list,U_list)
    #print grid_rho
    #print grid_e


    # find the avarage value for filling the empty value.
    av_p = np.average(p2_list)

    # Get h value based on rho ane e value.
    points_rho_e = []
    for i in range(len(rho_list)):
        coordinateTemp = [rho_list[i],e_list[i]]
        points_rho_e.append(coordinateTemp)
    points_rho_e =  np.asarray(points_rho_e) #convert the list to array

    prhoe_list = griddata(points_rho_e,p2_list,(grid_rho,grid_e),method='cubic', fill_value=av_p)









    print 'Starting to write tables. \n'
    #--------------------------------
    print 'Writing densityTable ...'


    rhoF = open("densityTable",'w')
    rhoF.write("( \n")
    for a in range(len(T_list)):
        rhoF.write("(" + str(T_list[a]) + "\n(\n")
        for b in range(len(p_list)):
            sList = ["\t(" + str(p_list[b]) + " " + str(rho_list[int(a*nP+b)]) + ")\n"]
            rhoF.write(" ".join(sList))
        rhoF.write(") ) \n") 
    rhoF.write(");")
    rhoF.close()    






    #--------------------------------
    print 'Writing kappaTable ...'
    kappaF = open("kappaTable",'w')
    kappaF.write("( \n")
    for a in range(len(T_list)):
        kappaF.write("(" + str(T_list[a]) + "\n(\n")
        for b in range(len(p_list)):
            sList = ["\t(" + str(p_list[b]) + " " + str(kappa_list[int(a*nP+b)]) + ")\n"]
            kappaF.write(" ".join(sList))
        kappaF.write(") ) \n") 
    kappaF.write(");")
    kappaF.close()    





    #--------------------------------
    print 'Writing muTable ...'
    muF = open("muTable",'w')
    muF.write("( \n")
    for a in range(len(T_list)):
        muF.write("(" + str(T_list[a]) + "\n(\n")
        for b in range(len(p_list)):
            sList = ["\t(" + str(p_list[b]) + " " + str(mu_list[int(a*nP+b)]) + ")\n"]
            muF.write(" ".join(sList))
        muF.write(") ) \n") 
    muF.write(");")
    muF.close()    



    #--------------------------------
    print 'Writing cpTable ...'
    cpF = open("cpTable",'w')
    cpF.write("( \n")
    for a in range(len(T_list)):
        cpF.write("(" + str(T_list[a]) + "\n(\n")
        for b in range(len(p_list)):
            sList = ["\t(" + str(p_list[b]) + " " + str(Cp_list[int(a*nP+b)]) + ")\n"]
            cpF.write(" ".join(sList))
        cpF.write(") ) \n") 
    cpF.write(");")
    cpF.close()





    #--------------------------------
    #To be noted that the hTable used an invert p,T table.
    print 'Writing hTable ...'
    hF = open("hTable",'w')
    hF.write("( \n")
    for a in range(len(p_list)):
        hF.write("(" + str(p_list[a]) + "\n(\n")
        for b in range(len(T_list)):
            sList = ["\t(" + str(T_list[b]) + " " + str(h_R_list[int(a*nT+b)]) + ")\n"]
            hF.write(" ".join(sList))
        hF.write(") ) \n") 
    hF.write(");")
    hF.close()    



    #--------------------------------
    print 'Writing cpMcvTable ...'
    cpMcvF = open("cpMcvTable",'w')
    cpMcvF.write("( \n")
    for a in range(len(T_list)):
        cpMcvF.write("(" + str(T_list[a]) + "\n(\n")
        for b in range(len(p_list)):
            sList = ["\t(" + str(p_list[b]) + " " + str(cpMcv_list[int(a*nP+b)]) + ")\n"]
            cpMcvF.write(" ".join(sList))
        cpMcvF.write(") ) \n") 
    cpMcvF.write(");")
    cpMcvF.close()    




    #--------------------------------
    print 'Writing ATable ...'
    AF = open("ATable",'w')
    AF.write("( \n")
    for a in range(len(T_list)):
        AF.write("(" + str(T_list[a]) + "\n(\n")
        for b in range(len(p_list)):
            sList = ["\t(" + str(p_list[b]) + " " + str(A_list[int(a*nP+b)]) + ")\n"]
            AF.write(" ".join(sList))
        AF.write(") ) \n") 
    AF.write(");")
    AF.close()    




    #--------------------------------
    print 'Writing hepTable ...'
    hep = open("hepTable",'w')
    hep.write("( \n")
    for a in range(len(U_list)):
        hep.write("(" + str(U_list[a]) + "\n(\n")
        for b in range(len(p_list)):
            sList = ["\t(" + str(p_list[b]) + " " + str(hep_list[b][a]) + ")\n"]
            hep.write(" ".join(sList))
        hep.write(") ) \n") 
    hep.write(");")
    hep.close()    



    #--------------------------------
    print 'Writing prhoeTable ...'
    prhoe = open("prhoeTable",'w')
    prhoe.write("( \n")
    for a in range(len(D_list)):
        prhoe.write("(" + str(D_list[a]) + "\n(\n")
        for b in range(len(U_list)):
            sList = ["\t(" + str(U_list[b]) + " " + str(prhoe_list[b][a]) + ")\n"]
            prhoe.write(" ".join(sList))
        prhoe.write(") ) \n") 
    prhoe.write(");")
    prhoe.close()    


    #--------------------------------
    print 'Writing TepTable ...'
    Tep = open("TepTable",'w')
    Tep.write("( \n")
    for a in range(len(U_list)):
        Tep.write("(" + str(U_list[a]) + "\n(\n")
        for b in range(len(p_list)):
            sList = ["\t(" + str(p_list[b]) + " " + str(Tep_list[b][a]) + ")\n"]
            Tep.write(" ".join(sList))
        Tep.write(") ) \n") 
    Tep.write(");")
    Tep.close()   


    #--------------------------------
    print 'Writing AepTable ...'
    Aep = open("AepTable",'w')
    Aep.write("( \n")
    for a in range(len(U_list)):
        Aep.write("(" + str(U_list[a]) + "\n(\n")
        for b in range(len(p_list)):
            sList = ["\t(" + str(p_list[b]) + " " + str(Aep_list[b][a]) + ")\n"]
            Aep.write(" ".join(sList))
        Aep.write(") ) \n") 
    Aep.write(");")
    Aep.close()   


    #--------------------------------
    print 'Writing rhoepTable ...'
    rhoep = open("rhoepTable",'w')
    rhoep.write("( \n")
    for a in range(len(U_list)):
        rhoep.write("(" + str(U_list[a]) + "\n(\n")
        for b in range(len(p_list)):
            sList = ["\t(" + str(p_list[b]) + " " + str(rhoep_list[b][a]) + ")\n"]
            rhoep.write(" ".join(sList))
        rhoep.write(") ) \n") 
    rhoep.write(");")
    rhoep.close()   
