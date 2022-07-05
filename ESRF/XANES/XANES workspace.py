"""

:platform: Unix, Windows 
:last changed: 2022-07-03

.. moduleauthor:: Tjark Leon Raphael Gröne <tgroene@physnet.uni-hamburg.de>


This Code is used to import and normalize h5 files for XANES spectra
As well rebin, filter or averaging can be applied
Furthermore the data can be exported as .dat file or as .mat file for MCR-ALS application


This little enviroment is based on the DAXS enviroment created by M. Retegan
for users of the ID26 at the ESRF facility for h5 file loading and area normalisation.


"""

import backbone.XANES_class as XANES_class
import pathlib
from time import time

#Suppress logging spam
import logging
logger = logging.getLogger()
logger.setLevel(logging.CRITICAL)

#----------------------------------------------------------------------------------------#
#                                                                                        #
#                                      Configurations                                    #
#                                                                                        #
#----------------------------------------------------------------------------------------#
# Define the counters to use for processing xanes:
counter_XAS = {
    "x": "hdh_energy",            # Define the x-axis, here the incident energy
    "signal": ["det_dtc_apd"],    # Define the data counter, here the main apd detector with dead-time correction
    "monitor": "I02",             # Define the counter used for normalization
    "temperature": "ls335_B"      # Define the temperature chanel of ct scans
}


#Sample Name
sampleName = "ZnAc and S in Oleylamine"

# Setup of the directory tree
rootdir = str(pathlib.Path(__file__).parent.resolve())
outdir = rootdir + "/processed/" + sampleName + "/"



# For the normalisation process
Energys = {
    "reject_first": 70,     	#Define number of rejecte first scans that might be duplicats or clustered with artifacts
    "delta_e0": 0,              #Define difference between measurements e0 and literature e0
    "e_max": 9.798,             #Define last point of a complite fscan to reject unfinished once
    "pre_min": 9.6450,          #Define start point of pre edge fit
    "pre_max": 9.655,           #Define end point of pre edge fit
    "post_min": 9.69,           #Define start point of post edge fit
    "post_max": 9.79827,        #Define end point of post edge fit
}

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#Some tools to work on the plots 
figratio1D =[4,3]   #Define the window size ratio for a 1D plot
figratio2D =[7,4]   #Define the window size ratio for a 2D plot
savefigType = "svg" #Possible formates are svg and png
linewidth = 5       #Define a linewidth for the 1D plots
framewidth = 5      #Define the thickness of the frame (relative to that ticks and offsets will be choosen) 
jump = 10           #Define the step lenght of the scans to show in the 1D plot 
minE = 9.658        #Define the minimum energy to start the plot with, if not used pleas set to None
maxE = 9.677        #Define the maximum energy to end the plot, if not used pleas set to None 
startScan = 42      #Define the first fscan to process for the plot
endScan = 200        #Define the last fscan to process for the plot
timeofset = 0       #Just a visulisation tool to shift the time 0 to over values like the reaktion starting time
raw = False          #Circumvent first normalisation step (area normalisation)
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


# To mark in a 2d plot the area a 1d plot will be take out of
rectangle = {
    "minE": 9.655,      #minimal Energy of the 1d plot
    "maxE": 9.70,       #maximal Energy of the 1d plot
    "mintime": 40,      #starting time of the 1d plot
    "maxtime": 3674,    #end time of the 1d plot
}

#~~~~~~~~~~~~~~~~~~~~~~~~~#
t0 = time() #start time   #
#~~~~~~~~~~~~~~~~~~~~~~~~~#

#----------------------------------------------------------------------------------------#
#                                                                                        #
#                                    Process data set                                    #
#                                                                                        #
#----------------------------------------------------------------------------------------#

#/////////////////\\\\\\\\\\\\\\\\\#
#           Insitu example         #
#/////////////////\\\\\\\\\\\\\\\\\#

filename = '/home/esrf/tjark1998a/Documents/cell_22_xanes_ZnS_olga_0001.h5' #input file directory
sampleName = "ZnAc and S in Oleylamine"
meas3 = XANES_class.XANES(filename,sampleName= sampleName, startScan=startScan, endScan=endScan,
                            counterXAS=counter_XAS, Energys=Energys,timeOfset=timeofset, figratio1D=figratio1D, raw = raw,
                            figratio2D=figratio2D, linewidth=linewidth,framewidth=framewidth )
meas3.norm() #normalise with pre and post edge fit
meas3.Filter("lFilter",7) #Filterer are lFilter or savgol, second argument is strength 
#meas3.saveReload(outputPath=rootdir + "/progressed/Reload") #For insitu very practical to shortan loading times
#meas3.rebin(3) #number of points to rebin
#meas3.average(9) #number of points to average over

#The Plots can all be saved with savefig=True and outputPath=: your output directory
meas3.plot1d(insitu=True,savefig=savefigType,outputPath=outdir + "/XANES/Plot_1D",jump=jump,minE=minE,maxE=maxE) #plot in 1D, insitu=True will define color change and hides ledgend
#meas3.plot2d(vmin=0.88,vmax=1.62,savefig=savefigType,outputPath=rootdir + "HC-4929/processed/Plot_2D") #colormap min and max can be set with vmin, vmax 
meas3.plot2dTemp(vmin=0.88,vmax=1.62,savefig=savefigType,outputPath=outdir + "/XANES/Plot_2D_+temp", rect=rectangle)
#meas3.plot3d(vmin=0.88,vmax=1.62,elev=16,azim=-85) #colormap min and max can be set with vmin and vmax | viewpoint can be set with elev and azim in [°]


#for output files to open in other software
# meas3.outputMatFile(rootdir + "/processed/MAT") #Will output .mat for MCR-ALS, takes output directory
meas3.outputDatFile(outdir + "XANES/Dat") #Will output curent processing result of every spectra as a .dat (meaning for insitu all the filtered, rebind or averaged spectra) 




#/////////////////\\\\\\\\\\\\\\\\\#
#        Some exitu examples       #
#/////////////////\\\\\\\\\\\\\\\\\#


filename = '/home/esrf/tjark1998a/Documents/cell_27_ZnAc_Oamine_XANES_ex-situ-1_0001.h5'
sampleName = "ZnAc in Oleylamine"
meas = XANES_class.XANES(filename,sampleName= sampleName, startScan=1,endScan=20,counterXAS=counter_XAS, Energys=Energys, raw = raw)
meas.norm()
meas.Filter("lFilter",5)
#meas.rebin(3)
meas.average() # For exitu samples very conveniant
#meas.saveReload(outputPath=rootdir + "/processed/Reload") #For insitu very practical to shortan loading times
meas.outputDatFile(outdir + "/XANES/Dat")

filename = '/home/esrf/tjark1998a/Documents/Znac-exsitu-1_0001.h5'
meas2 = XANES_class.XANES(filename,sampleName= "ZnAc Powder", startScan=1,endScan=2,counterXAS=counter_XAS, Energys=Energys, raw = raw)
meas2.norm()
meas2.Filter("lFilter",5)
#meas2.rebin(3)
#meas2.average()
#meas2.saveReload(outputPath=rootdir + "/processed/Reload") #For insitu very practical to shortan loading times
meas2.outputDatFile(outdir + "/XANES/Dat")


# Find e0 of a measurment and output it
meas.find_e0() 
meas2.find_e0()
meas3.find_e0()
print("ZnAc+S+OA: " + str(meas3.e0))
print("ZnAc+OA: " + str(meas.e0))
print("ZnAc Powder: " + str(meas2.e0))


#///////////////////\\\\\\\\\\\\\\\\\\#
#   Plotting multiple spectra in 1D   #
#///////////////////\\\\\\\\\\\\\\\\\\#

x = [meas3.x,meas.x,meas2.x] #energy point should be the samen, but 
matrix = [meas3.matrix,meas.matrix,meas2.matrix] #Intensities
sampleNames = [meas3.sampleName,meas.sampleName,meas2.sampleName] #Sample names: for inhouse include a string with sample name

# The Scan is a list of start scan, end scan and the total number of scans scans 
# The inhouse scan might has only one scan after preprocessing so scan will be [1,1,1] for it
scan = [[meas3.startScan,meas3.endScan,meas3.totalIncludedScans],
        [meas.startScan,meas.endScan,meas.totalIncludedScans],
        [meas2.startScan,meas2.endScan,meas2.totalIncludedScans]
        ]
insitu = [False,False,False] #with scans are insitu meas.


#      Plotting multiple spectra      #
# =====================================#

###### Only thing to change here is plotName and maby outputPath
plotName = "ZnAc in differen environments raw range=90eV"
XANES_class.plot1d_multi(x=x,matrix=matrix,sampleNames=sampleNames,scan=scan,insitu=insitu,
                            linewidth=linewidth, framewidth=framewidth, figratio1D=figratio1D, 
                            minE=minE, maxE=maxE, savefig=savefigType, plotName=plotName, 
                            outputPath=outdir + "XANES/Plot_1D_multi", markWhiteline= True)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
XANES_class.processTime(t0) #determine processing time #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#       #         #     #           #   #########       ########    #       #   ##      #    #########
#       #        # #     #         #    #               #           #       #   # #     #     #     #
#       #       #   #     #       #     #               #           #       #   #  #    #      #   #
#########      #######     #     #      #########       ########    #       #   #   #   #       # #
#       #     #       #     #   #       #               #           #       #   #    #  #        #
#       #    #         #     # #        #               #           #       #   #     # #       ### 
#       #   #           #     #         #########       #           #########   #      ##       ###
