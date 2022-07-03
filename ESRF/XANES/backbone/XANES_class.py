"""

:platform: Unix, Windows
:synopsis: class as backbone for XANES processing 
:last changed: 2022-05-28

.. moduleauthor:: Tjark Leon Raphael Gröne <tgroene@physnet.uni-hamburg.de>

"""
# Import the relevant parts from the DAXS enviroment
import silx
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.colors import ListedColormap

from time import time
from datetime import datetime
from scipy import io as sp_io 
from scipy.signal import savgol_filter, lfilter


#Colorsheme
import backbone.Colormap as cmp

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#S/////////////////////////////////////////////////////////////////////////////////////////////////////S#
#S//                                                                                                 //S#
#S//    |----------------------------------------------------------------------------------------|   //S#
#S//    |                                                                                        |   //S#
#S//    |                                    The XANES Class                                     |   //S#
#S//    |                                                                                        |   //S#
#S//    |----------------------------------------------------------------------------------------|   //S#
#S//                                                                                                 //S#
#S/////////////////////////////////////////////////////////////////////////////////////////////////////S#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#



class XANES:
          
    #----------------------------------------------------------------------------------------#
    #                                           Init                                         #
    #----------------------------------------------------------------------------------------#
    def __init__(self, path: str,sampleName: str,startScan: int,endScan: int,
                counterXAS, Energys, timeOfset=0, figratio1D=[16,9], figratio2D=[7, 4],
                linewidth=3, framewidth=3, reload: str = None, raw = False):
        """
        Arguments
        ----------
        path:           Path to file
        sampleName:     Name of your sample (for Ledgend)
        startScan:      First scan to include 
        endScan:        Last scan to include
        counterXAS:     counters for the Scan measurements    
        timeOfset:      Ofset of the timeax for insitu measurements in s (-10 time axis starts with -10 s)       
        figratio1D:     Ratio of the plot window for 1D plots
        figratio2D:     Ratio of the plot window for 2D plots
        linewidth:      Linewith for every 1D plot
        reload:         Path to a saved MAT file, that was saved with the saveReload() function before
        Energys:        Dictionary containing delta_e0, reject_first, e:max, pre_min, pre_max, post_min, post_max  

        Where   delta_e0 is the ofset of e0 to those in literature
                reject_first is the number of scans at the begining that will be leaft out
                e_max is the last energy point in your scan that needs to be included, 
                    	scans with less energy points will be rejected
                pre_min is the starting energy for the pre edge fit
                pre_max is the ending energy for the pre edge fit
                post_min is the starting energy for the post edge fit
                post_max is the ending energy for the post edge fit 
        """
        self.path = path
        self.sampleName = sampleName
        self.startScan = startScan
        self.endScan = endScan
        self.linewidth = linewidth
        self.framewidth = framewidth

        #The labes, frames andticks are tight to a specific window with, so the ratio musst be normed to that
        ratio1D = figratio1D[0]/figratio1D[1] 
        figratio1D = [ratio1D*9,9]
        self.figratio1D = figratio1D
        ratio2D = figratio2D[0]/figratio2D[1] 
        figratio2D = [ratio2D*4,4]
        self.figratio2D = figratio2D

        self.delta_e0 = Energys["delta_e0"]
        self.reject_first = Energys["reject_first"]
        self.e_max = Energys["e_max"]			    # last energy point of fit
        self.pre_min = Energys["pre_min"]		    # pre edge minimum energy for fit
        self.pre_max = Energys["pre_max"]			    # pre edge maximum energy for fit
        self.post_min = Energys["post_min"]				# post edge minimum energy for fit
        self.post_max = Energys["post_max"]          # post edge maximum energy for fit

   
        self.totalIncludedScans = 1
        self.matrix = []
        self.time_array = []
        self.temp_array = []
        self.temptime_array = []
        self.i0 = []
        self.max_all_int = 0
        self.max_x = 0
        self.min_x = 0
        start_time = 0

        if reload != None:
            loadedmat = sp_io.loadmat(reload)
            self.matrix = loadedmat['data_matrix'][self.startScan:self.endScan]
            self.x = loadedmat['energy_axis'][0]
            self.time_array = loadedmat['time'][0][self.startScan:self.endScan]
            self.temp_array = loadedmat['temp'][self.startScan:self.endScan]
            self.temptime_array = loadedmat['temptime'][0][self.startScan:self.endScan]
        else:
            from daxs.measurements import Source, Xas
            for i in range(self.startScan,self.endScan):

                source = Source(self.path, i, None)
            
                #Reject unfinished and non fscans
                try:
                    self.meas = Xas(source, counters=counterXAS) #only works for fscans
                    with silx.io.open(source.filename) as sf:
                        for selection in source.selections:
                            end_time = sf[selection]["end_time"][()]
                        if start_time == 0:
                            start_time = datetime.strptime(str(end_time)[2:12] + " " + str(end_time)[13:-7], "%Y-%m-%d %H:%M:%S.%f")
                            self.time_array.append(0+timeOfset)
                        else:
                            delta_t = datetime.strptime(str(end_time)[2:12] + " " + str(end_time)[13:-7], "%Y-%m-%d %H:%M:%S.%f") - start_time
                            delta_t = delta_t.total_seconds()
                            self.time_array.append(delta_t+timeOfset)
                except:

                    try:
                        with silx.io.open(source.filename) as sf:
                            for selection in source.selections:
                                end_time = sf[selection]["end_time"][()]
                                measurement = sf[selection]["measurement"]
                                temperature = measurement[counterXAS["temperature"]][()] - 273.15
                                self.temp_array.append(temperature)
                            delta_t = datetime.strptime(str(end_time)[2:12] + " " + str(end_time)[13:-7], "%Y-%m-%d %H:%M:%S.%f") - start_time
                            delta_t = delta_t.total_seconds()
                            self.temptime_array.append(delta_t+timeOfset)
                    except:
                        continue #reject if not fscan nor ct scan
                    continue #reject if not fscan
                if max(self.meas.x) < self.e_max: #only works if fscan was complited to full energy range
                    continue #reject if energy range is incomplete

                #
                if raw != True:
                    self.meas.aggregate() #sums up all scans data
                    self.meas.find_outliers(threshold=9)
                    self.meas.remove_outliers()      
                    self.meas.normalize(mode="area") # normalize to intg. intensity


                signal = self.meas.signal[self.reject_first:]
                self.x = self.meas.x[self.reject_first:]

                self.matrix.append(signal)
                self.i0.append(self.meas.monitor)

    def saveReload(self,outputPath):
        """
        Arguments
        ----------
        outputPath: Path to output directory
        """    
        if outputPath != None:
            ensure_dir(outputPath)
        sp_io.savemat(outputPath + "/" + self.sampleName +'.mat', {'data_matrix': self.matrix, 'energy_axis': self.x,'time': self.time_array ,'temptime': self.temptime_array,'temp': self.temp_array})


    #----------------------------------------------------------------------------------------#
    #                      Linear fit of pre and post edge + finding e0                      #
    #----------------------------------------------------------------------------------------#
    def preedge(self,y):
        """
        Arguments
        ----------
        y:  Intensity array (as a slice of the Intesnsity matrix)
        """ 
        #Find e0
        d2x_y = np.gradient(y,2)
        index_e0 = np.where(d2x_y==max(d2x_y))
        self.e0 = self.x[index_e0]

        preStartindex = None
        preEndindex = None
        postStartindex = None
        postEndindex = None

        for j in self.x:
            if j >= self.pre_min and preStartindex==None:
                preStartindex = int(np.where(self.x==j)[0]) 
            if j > self.pre_max and preEndindex==None:
                preEndindex = int(np.where(self.x==j)[0]) - 1 
            if j >= self.post_min and postStartindex==None:
                postStartindex = int(np.where(self.x==j)[0]) 
            if j > self.post_max and postEndindex==None:
                postEndindex = int(np.where(self.x==j)[0]) - 1  
        pre_x = self.x[preStartindex:preEndindex]
        pre_y = y[preStartindex:preEndindex]
        post_x = self.x[postStartindex:postEndindex]
        post_y = y[postStartindex:postEndindex] 

        self.edgejump = y[int(np.where(self.x==self.e0)[0])]

        pre_m,pre_b = np.polyfit(pre_x, pre_y, 1)
        post_m,post_b = np.polyfit(post_x, post_y, 1)
        self.prefit = (self.x*pre_m + pre_b) + self.edgejump 
        self.postfit = self.x*post_m + post_b

    #----------------------------------------------------------------------------------------#
    #                            Normalization of pre and post edge                          #
    #----------------------------------------------------------------------------------------#
    def norm(self): 
        """
        Arguments
        ----------
        none
        """ 
        #Find e0 and pre and post edge slope
        count = 0
        for i in self.matrix:
            self.preedge(i)
            #normalise spectrum to pre and post edge slope
            index_min_matrix = np.where(i==min(i))
            min_matrix = min(i) / self.prefit[index_min_matrix] 
            for j in range(0,len(self.x)):
                if self.x[j]*1e3 < self.e0:
                    self.matrix[count][j] = (i[j] / self.prefit[j]) - min_matrix[0] 
                elif self.x[j]*1e3 >= self.e0:
                    self.matrix[count][j] = i[j] / self.postfit[j]
            count = count + 1
        self.x = self.x + self.delta_e0

    def find_e0(self):
        for i in self.matrix:
            d2x_y = np.gradient(i,2)
            index_e0 = np.where(d2x_y==max(d2x_y))
            self.e0 = self.x[index_e0]

    #----------------------------------------------------------------------------------------#
    #                                  Export .mat for MCR-ALS                                #
    #----------------------------------------------------------------------------------------#
    def outputMatFile(self,outputPath):
        """
        Arguments
        ----------
        outputPath: Path to output directory
        """    
        if outputPath != None:
            ensure_dir(outputPath)
        sp_io.savemat(outputPath + "/" + self.sampleName +'.mat', {'data_matrix': self.matrix, 'energy_axis': self.x})


    #----------------------------------------------------------------------------------------#
    #                                   Export .dat to export                                #
    #----------------------------------------------------------------------------------------#
    def outputDatFile(self,outputPath):
        """
        Arguments
        ----------
        outputPath: Path to output directory
        """
        if outputPath != None:
            ensure_dir(outputPath)
        #Data output file
        for i, j in enumerate(self.matrix):
            data = np.array([self.x, j])
            outfname = outputPath + "/" + self.path.split("/")[-1].split(".")[0]+ "_" + str(i + self.startScan).zfill(5) + '.txt' #output directory and filename
            np.savetxt(outfname, data.T,header = 'Energy(keV), Intensity (arb.u.) \n' + "Time:  " + str(self.time_array[i]))

    #----------------------------------------------------------------------------------------#
    #                                         Average                                        #
    #----------------------------------------------------------------------------------------#
    def average(self,averagedSteps=None):
        """
        Arguments
        ----------
        averagedSteps: Number of Scans to average over
        """
        pre_array = []
        newmatrix = []

          #get all included scans to averaged over in one step of the data (for exitu includes by default all)

        if averagedSteps == None:
            averagedSteps = self.endScan - self.startScan  
        else:
            self.totalIncludedScans = averagedSteps 

        
        #Loop over all scans from start to end scan
        for i in self.matrix:
            pre_array.append(i)
                
            if len(pre_array) == averagedSteps:
                pre = sum(pre_array) / averagedSteps
                newmatrix.append(pre)
                pre_array = []

        self.matrix = newmatrix


    #----------------------------------------------------------------------------------------#
    #                                         Rebin                                          #
    #----------------------------------------------------------------------------------------#
    def rebin(self,rebin):
        """
        Arguments
        ----------
        rebin: Number of points to rebining to new data point
        """
        count = 0
        x_rebin = []
        pre_rebin = []  

        for i in self.matrix:
            for j in range(0,len(self.x)-rebin,rebin):#
                x_rebin.append(np.array(self.x)[j:j+rebin].sum() / rebin)
                pre_rebin.append(np.array(i[j:j+rebin]).sum() / rebin)  
            
            
            newx = x_rebin
            self.matrix[count] = pre_rebin
            count = count + 1
            x_rebin = []
            pre_rebin = []
        self.x = newx
            
    #----------------------------------------------------------------------------------------#
    #                                          Filter                                        #
    #----------------------------------------------------------------------------------------#
    def Filter(self,Filter,strength):
        """
        Arguments
        ----------
        Filter:     Type of Filter to use (lFilter or savgol)
        strength:   Strength with what the filter is applied 
        """
        count = 0
        for i in self.matrix:
            #applying a filter if set
            if Filter == "savgol":
                    self.matrix[count] = savgol_filter(i,strength,2)
            elif Filter == "lFilter":
                    n = strength  # the larger n is, the smoother curve will be
                    self.matrix[count] = lfilter([1.0 / n] * n,1,i)
            count = count + 1
            
            
    #----------------------------------------------------------------------------------------#
    #                                        plot data                                       #
    #----------------------------------------------------------------------------------------#
    def plot1d(self,startTime = 0, offsetTime = 0, insitu=False,color="black", title: str = None,savefig: str = None,outputPath=None, jump=None, maxE=None, minE=None):
        """
        Arguments
        ----------
        insitu:     Data was colcted insitu
        color:      Difine specific color for plot
        title:      Give plot title
        savefig:    When True saves the plot as a png
        outoutPath: Path to output directory
        startTime:  Time offset by which the first scan is displayed for insitu measurments
        offsetTime: Time the temperature hits its set point
        jump:       Step width for insitu measurments
        minE:       Lower border of the plot in Energy
        maxE:       Upper border of the plot in Energy
        """
        if outputPath != None:
            ensure_dir(outputPath)
        if title == None:
            title = self.sampleName

        plt.figure(figsize=(self.figratio1D[0], self.figratio1D[1]), dpi=50) #
        ax = plt.gca()     

        #define maxima and minima for energy and intensity for ploting purpose
        self.max_x = max(self.x)
        self.min_x = min(self.x)
        for i in self.matrix:
            if max(i) > self.max_all_int:
                self.max_all_int = max(i)


        if insitu == True:
            for i,j in enumerate(self.matrix): #insitu plot
                if i % jump == 0:
                    #loading Colormap for insitu spectra 
                    colormapRB = cmp.cmRB_Time(int(self.time_array[-1])-startTime)
                    newcmp = ListedColormap(colormapRB(np.linspace(0,1,int(self.time_array[-1])-startTime)))
                    colormapRB = colormapRB(np.linspace(0,1,int(self.time_array[-1])-startTime))

                    closest = find_nearest(self.temptime_array, self.time_array[i])
                    tempindex = np.where(self.temptime_array==closest)


                    if len(self.matrix) % 2 ==0:
                        ramp = np.linspace(self.linewidth,self.linewidth/2,int(len(self.matrix)/2)) 
                        thickness = np.concatenate((ramp,ramp[::-1])) 
                    else:
                        ramp = np.linspace(self.linewidth,self.linewidth/2,int(len(self.matrix)/2)) 
                        thickness = np.concatenate(np.concatenate(ramp,np.array([self.linewidth/2])),ramp[::-1]) 


                    if self.time_array[i] < offsetTime:
                        prelabel = "pre "
                    else:
                        prelabel = ""


                    ax.plot(self.x, j, linewidth=thickness[i], color=colormapRB[int(self.time_array[i])], label=prelabel + str(int(self.temp_array[tempindex])) + " °C", alpha=1)
                    ax.legend(prop={'size': 18}, ncol=3,loc='lower right', handlelength=1.6).get_frame().set_linewidth(0.0)

            sm = plt.cm.ScalarMappable(cmap=newcmp, norm=plt.Normalize(vmin=0, vmax=int(self.time_array[-1])-startTime))
            cb = plt.colorbar(sm)
            cb.set_label(label='Time [s]', size=26)
            cb.ax.tick_params(axis='both',which="major", width=self.framewidth, length=self.framewidth*3, labelsize=26)  
            cb.outline.set_linewidth(self.framewidth)
        else:
            #exitu plot
            ax.plot(self.x, self.matrix[0], linewidth=self.linewidth, color=color, label=self.sampleName)
            ax.legend(prop={'size': 18}, ncol=4).get_frame().set_linewidth(0.0)

        ax.set_title(title, fontsize = 26)    
        ax.set_xlabel("Energy [keV]", fontsize = 26)
        ax.set_ylabel("Normalized absorption [edge fraction]", fontsize = 26)
        ax.tick_params(axis='both',which="major", width=self.framewidth, length=self.framewidth*3, labelsize=26)
        ax.spines["top"].set_linewidth(self.framewidth)
        ax.spines["bottom"].set_linewidth(self.framewidth)
        ax.spines["bottom"].set_position(("outward", self.framewidth/2))
        ax.spines["right"].set_linewidth(self.framewidth)
        ax.spines["left"].set_linewidth(self.framewidth)

        if maxE != None and minE !=None:
            plt.xlim(minE,maxE)
        elif maxE != None and minE == None:
            plt.xlim(self.min_x,maxE)
        elif minE != None and maxE == None:
            plt.xlim(minE,self.max_x)
        elif maxE == None and minE == None:
            plt.xlim(self.min_x,self.max_x)
        plt.ylim(bottom=0)
        plt.tight_layout()
        if savefig == "png":
            plt.savefig(outputPath + "/1d_"+ self.sampleName +".png", dpi=600, bbox_inches = "tight")
        if savefig == "svg":
            plt.savefig(outputPath + "/1d_"+ self.sampleName +".svg", dpi=600, bbox_inches = "tight")
        plt.show()
        

    #----------------------------------------------------------------------------------------#
    #                                  Plot insitu XANES in 2D                               #
    #----------------------------------------------------------------------------------------#

    def plot2d(self,vmin=None,vmax=None,savefig: str = None,outputPath=None, rect=None):
        """
        Arguments
        ----------
        vmin:       Intensity the colormap should start
        vmax:       Intensity the colormap should end
        savefig:    When True saves the plot as a png
        outoutPath: Path to output directory
        rect:       Dictionary of the rectangle properties that you can use to mark a certain area in the 2D plot
        """
        if outputPath != None:
            ensure_dir(outputPath)
        #set some plotting settings
        newcmp = cmp.cm2d()

        fig = plt.figure(figsize=(self.figratio2D[0], self.figratio2D[1]), dpi=125) 
        ax = plt.gca() 

        x = self.x
        y = self.time_array
        z = self.matrix

        imap = ax.imshow(z, origin = 'lower', aspect = 'auto', 
                extent = [min(x), max(x), min(y), max(y)], vmin = vmin, vmax = vmax, 
                cmap=newcmp,Interpolation='none')

        ax.set_xticks(np.arange(round(min(x),2),round(max(x),2), 0.02))
        ax.set(xlim = (min(x),max(x)), xlabel = 'Energy [keV]', ylabel="Time [s]")

        ax.spines["top"].set_linewidth(self.framewidth*0.4)
        ax.spines["bottom"].set_linewidth(self.framewidth*0.4)
        ax.spines["right"].set_linewidth(self.framewidth*0.4)
        ax.spines["left"].set_linewidth(self.framewidth*0.4)

        ax.tick_params(axis='both',which="major", width=self.framewidth*0.4, length=self.framewidth*3*0.4, labelsize=10.4)
        ax.set_title(self.sampleName)

        if rect != None:
            rect = patches.Rectangle((rect["minE"], rect["mintime"]), (rect["maxE"]-rect["minE"]), (rect["maxtime"]-rect["mintime"]), linewidth=self.linewidth/2, edgecolor='r', facecolor='none',zorder=100)
            ax.add_patch(rect)

        cb = plt.colorbar(imap)
        cb.set_label(label='Normalized absorption [edge fraction]', size=10.4)
        cb.ax.tick_params(axis='both',which="major", width=self.framewidth*0.4, length=self.framewidth*3*0.4, labelsize=10.4)  
        cb.outline.set_linewidth(self.framewidth*0.4)

        plt.tight_layout()
        if savefig == "png":
            plt.savefig(outputPath + "/2d_"+ self.sampleName +".png", dpi=600, bbox_inches = "tight")
        if savefig == "svg":
            plt.savefig(outputPath + "/2d_"+ self.sampleName +".svg", dpi=600, bbox_inches = "tight")
        plt.show()

    #----------------------------------------------------------------------------------------#
    #                  Plot Insitu data in 2D with temp. curvenext to it                     #
    #----------------------------------------------------------------------------------------#

    def plot2dTemp(self,vmin=None,vmax=None,savefig: str = None,outputPath=None, rect=None):
        """
        Arguments
        ----------
        vmin:       Intensity the colormap should start
        vmax:       Intensity the colormap should end
        savefig:    When True saves the plot as a png
        outoutPath: Path to output directory
        rect:       Dictionary of the rectangle properties that you can use to mark a certain area in the 2D plot
        """
        if outputPath != None:
            ensure_dir(outputPath)
        #set some plotting settings
        newcmp = cmp.cm2d()

        fig, (ax1, ax2) = plt.subplots(1, 2, gridspec_kw={'width_ratios': [1, 4]},figsize=(self.figratio2D[0], self.figratio2D[1]), dpi=125)
        plt.subplots_adjust(wspace = 0.075 )

        time = self.temptime_array
        temp = self.temp_array

        ax1.plot(temp,time, linewidth=self.linewidth, color="red")
        ax1.set(xlabel = 'T [°C]', ylabel="Time [s]",ylim = (min(time),max(time)), xlim=(0,200))
        ax1.tick_params(axis='both',which="major", width=self.framewidth*0.4, length=self.framewidth*3*0.4, labelsize=10.4)

        x = self.x
        y = self.time_array
        z = self.matrix

        ax1.set_ylim(min(y),max(y))

        ax1.spines["top"].set_linewidth(self.framewidth*0.4)
        ax1.spines["bottom"].set_linewidth(self.framewidth*0.4)
        ax1.spines["right"].set_linewidth(self.framewidth*0.4)
        ax1.spines["left"].set_linewidth(self.framewidth*0.4)

        imap = ax2.imshow(z, origin = 'lower', aspect = 'auto', 
                extent = [min(x), max(x), min(y), max(y)], vmin = vmin, vmax = vmax, 
                cmap=newcmp,interpolation='none')

        ax2.set_xticks(np.arange(round(min(x),2),round(max(x),2), 0.02))
        ax2.set(xlim = (min(x),max(x)), xlabel = 'Energy [keV]')

        if rect != None:
            rect = patches.Rectangle((rect["minE"], rect["mintime"]), (rect["maxE"]-rect["minE"]), (rect["maxtime"]-rect["mintime"]), linewidth=self.linewidth/2, edgecolor='r', facecolor='none',zorder=100)
            ax2.add_patch(rect)




        ax2.tick_params(axis='both',which="major", width=self.framewidth*0.4, length=self.framewidth*3*0.4, labelsize=10.4)
        ax2.tick_params(axis="y",left=False,labelleft=False)
        ax2.set_title(self.sampleName)
        
        ax2.spines["top"].set_linewidth(self.framewidth*0.4)
        ax2.spines["bottom"].set_linewidth(self.framewidth*0.4)
        ax2.spines["right"].set_linewidth(self.framewidth*0.4)
        ax2.spines["left"].set_linewidth(self.framewidth*0.4)


        cb = plt.colorbar(imap)
        cb.set_label(label='Normalized absorption [edge fraction]', size=10.4)
        cb.ax.tick_params(axis='both',which="major", width=self.framewidth*0.4, length=self.framewidth*3*0.4, labelsize=10.4)  
        cb.outline.set_linewidth(self.framewidth*0.4)

        plt.tight_layout()
        if savefig == "png":
            plt.savefig(outputPath + "/2d_Temp_" + self.sampleName +".png", dpi=600, bbox_inches = "tight")
        if savefig == "svg":
            plt.savefig(outputPath + "/2d_Temp_" + self.sampleName +".svg", dpi=600, bbox_inches = "tight")
        plt.show()        


    def waterfall3D(self,jump=1,setTime=0):

        fig = plt.figure()
        ax = fig.add_subplot(projection='3d')

        colormapRB = cmp.cmRB_Time(int(self.time_array[-1])-setTime)
        newcmp = ListedColormap(colormapRB(np.linspace(0,1,int(self.time_array[-1])-setTime)))
        colormapRB = colormapRB(np.linspace(0,1,int(self.time_array[-1])-setTime))
        for (i, j) in sorted(enumerate(self.matrix), reverse=True):
        #for i,j in enumerate(self.matrix): 
            if i % jump == 0:
                ax.plot(self.x, j, zs=self.time_array[i], zdir='y',color=colormapRB[int(self.time_array[i])])
                col_obj = ax.fill_between(self.x, min(j), j, step='pre', alpha=0.1, color=colormapRB[int(self.time_array[i])]) 
                ax.add_collection3d(col_obj, zs = self.time_array[i], zdir = 'y')


        ax.set_xlabel('X')
        ax.set_xlim3d(min(self.x), max(self.x))
        ax.set_ylabel('Y')
        ax.set_ylim3d(0,int(self.time_array[-1])-setTime)
        ax.set_zlabel('Z')
        ax.set_zlim3d(0, 1.5)
        ax.set_title("3D Waterfall plot")

        plt.show()


    #----------------------------------------------------------------------------------------#
    #                             Plot insitu data as a 3D plot                              #
    #----------------------------------------------------------------------------------------#

    def plot3d(self, vmin= None, vmax=None, elev=None, azim=None,savefig: str = None,outputPath=None):  
        """
        Arguments
        ----------
        vmin:       Intensity the colormap should start
        vmax:       Intensity the colormap should end
        elev:       Elevation angle 0 equals ground or 1d side perspective
        azim:       Azimutal angel -90 equals 1d front perspective
        savefig:    When True saves the plot as a png
        outoutPath: Path to output directory
        """
        if outputPath != None:
            ensure_dir(outputPath)
        from matplotlib.ticker import LinearLocator, FormatStrFormatter  
        fig = plt.figure()
        ax = fig.gca(projection='3d')
        fig.subplots_adjust(top=1.1, bottom=-.1)


        newcmp = cmp.cm2d()
        # Plot the surface.
        x0, x1 = np.meshgrid(np.array(self.x), np.array(self.time_array))
        surf = ax.plot_surface(x0, x1, np.array(self.matrix), cmap=newcmp,
                            linewidth=-3, antialiased=False, rstride=1, cstride=1, alpha=1, vmin=vmin, vmax=vmax)

        # Customize the z axis.
        # ax.set_zlim(-1.01, 1.01)
        ax.zaxis.set_major_locator(LinearLocator(10))
        ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))

        #ax.set_zticks([])
        ax.set_xlabel("Energy [keV]")
        ax.set_ylabel("Time [s]")#"time [s]")
        ax.set_zlabel("Normalized absorption [edge fraction]")
        ax.view_init(elev=elev, azim=azim)

        # Add a color bar which maps values to colors.
        # fig.colorbar(surf, shrink=0.5, aspect=5)

        plt.tight_layout()
        if savefig == "png":
            plt.savefig(outputPath + "/3d_"+ self.sampleName +".png", dpi=600, bbox_inches = "tight")
        if savefig == "svg":
            plt.savefig(outputPath + "/3d_"+ self.sampleName +".svg", dpi=600, bbox_inches = "tight")
        plt.show()

    #----------------------------------------------------------------------------------------#
    #                                Plot the temp. curve alone                              #
    #----------------------------------------------------------------------------------------#

    def plotTemp(self,savefig=False,outputPath=None):
        """
        Arguments
        ----------
        savefig:    When True saves the plot as a png
        outoutPath: Path to output directory
        """
        if outputPath != None:
            ensure_dir(outputPath)
        plt.figure(figsize=(self.figratio1D[0], self.figratio1D[1]), dpi=50) 
        ax = plt.gca()     


        ax.plot(self.temptime_array, self.temp_array, linewidth=self.linewidth, color="b")

        ax.set_title("Temp", fontsize = 26)    
        ax.set_xlabel("Time [s]", fontsize = 20)
        ax.set_ylabel("Temp [°C]", fontsize = 20)
        ax.tick_params(axis='both',which="major", width=self.framewidth, length=self.framewidth*3, labelsize=20)
        ax.spines["top"].set_linewidth(self.framewidth)
        ax.spines["bottom"].set_linewidth(self.framewidth)
        ax.spines["right"].set_linewidth(self.framewidth)
        ax.spines["left"].set_linewidth(self.framewidth)

        plt.xlim(0,max(self.temptime_array))
        plt.tight_layout()
        if savefig == "png":
            plt.savefig(outputPath + "/temp_"+ self.sampleName +".png", dpi=600, bbox_inches = "tight")
        if savefig == "svg":
            plt.savefig(outputPath + "/temp_"+ self.sampleName +".svg", dpi=600, bbox_inches = "tight")
        plt.show()

    #----------------------------------------------------------------------------------------#
    #               Find the spot the set temp is reached by a linear regressi               #
    #----------------------------------------------------------------------------------------#
    
    def tempRegression(self,Time):
        """
        Arguments
        ----------
        Time:   Dictionary of rampStart, rampEnd, constStart,constEnd

        Where   rampStart is the relative starting time of the ramping
                rampEnd is the relative ending time of the ramping
                constStart is the relative starting time of the region of set constant temperature
                constEnd is the relative ending time of the region of set constant temperature
        """
        self.rampStart = Time["rampStart"]
        self.rampEnd = Time["rampEnd"]
        self.constStart = Time["constStart"]
        self.constEnd = Time["constEnd"]

        rampStartindex=None
        rampEndindex=None
        constStartindex=None
        constEndindex=None

        for i in self.temptime_array:
            if i > self.rampStart and rampStartindex==None:
                rampStartindex = int(self.temptime_array.index(i))
            if i > self.rampEnd and rampEndindex==None:
                rampEndindex = int(self.temptime_array.index(i) - 1)
            if i > self.constStart and constStartindex==None:
                constStartindex = int(self.temptime_array.index(i))
            if i > self.constEnd and constEndindex==None:
                constEndindex = int(self.temptime_array.index(i) - 1)  
        rampTime = self.temptime_array[rampStartindex:rampEndindex]
        rampTemp = self.temp_array[rampStartindex:rampEndindex]
        constTime = self.temptime_array[constStartindex:constEndindex]
        constTemp = self.temp_array[constStartindex:constEndindex]
        ramp_m,ramp_b = np.polyfit(rampTime, rampTemp, 1)
        const_m,const_b = np.polyfit(constTime, constTemp, 1)
        x = np.linspace(self.rampStart,self.constEnd,1000)
        ramp = x*ramp_m + ramp_b
        const = x*const_m + const_b
        self.t0 = x[np.argwhere(np.diff(np.sign(ramp - const))).flatten()]
        print("\n Time the set temperature is reached")
        print("------------------------------------------------------------------")
        print("Time: " + str(self.t0))

    #----------------------------------------------------------------------------------------#
    #                                    Export temp as .xy                                  #
    #----------------------------------------------------------------------------------------#

    def outputTemp(self,outputPath):
        """
        Arguments
        ----------
        outputPath: Path to output directory
        """
        if outputPath != None:
            ensure_dir(outputPath)
        #Data output file
        data = np.array([self.temptime_array, self.temp_array])
        outfname = outputPath + "/" + self.path.split("/")[-1].split(".")[0]+ '.xy' #output directory and filename
        np.savetxt(outfname, data.T,header = 'Time(s) Temperature(°C)')
            


#----------------------------------------------------------------------------------------#
#                                   Get processing time                                  #
#----------------------------------------------------------------------------------------#
def processTime(t0):
    """
    Arguments
    ----------
    t0: Start time
    """
    print("\n Time for processing")
    print("------------------------------------------------------------------")
    print("Time: " + str(time() - t0))

#----------------------------------------------------------------------------------------#
#               Plot multiple exitu or isitu measurments in one 1D plot                  #
#----------------------------------------------------------------------------------------#

def plot1d_multi(x,matrix,sampleNames,scan,insitu=[False],color=None, 
                title: str = None, figratio1D=[16,9], linewidth=3, framewidth=3,
                 maxE=None, minE=None, savefig=None, outputPath=None, plotName=None, markWhiteline=False):
    """
    Arguments
    ----------
    x:              Energy ax   
    matrix:         Intensities 
    sampleNames:    Names of the Samples
    insitu:         Data was colected insitu
    color:          Difine specific color for plot
    title:          Give plot title
    figrati1D:      Ratio of the plot window
    linewidth:      linewith of the plotted graphs
    """
    if title == None:
        title = ""

    if outputPath != None:
        ensure_dir(outputPath)
    
    matrix = np.array(matrix)
    max_all_int = 0
    max_x = 0
    min_x = 1e10
    ratio1D = figratio1D[0]/figratio1D[1] 
    figratio1D = [ratio1D*9,9]
       

    plt.figure(figsize=(figratio1D[0], figratio1D[1]), dpi=50) #
    ax = plt.gca()     

    for i in np.arange(len(matrix)):
        if max(x[i]) > max_x:
            max_x = max(x[i])
        if min(x[i]) < max_x:
            min_x = min(x[i])
        for j in matrix[i]:
            if max(j) > max_all_int:
                max_all_int = max(j)

        if insitu[i] == True:
            #loading Colormap for insitu spectra
            myRange = (scan[i][1] - scan[i][0]) / scan[i][2]
            colormap = cmp.cm(int(myRange)+3)

            count = 0
            for j in matrix[i]: #insitu plot
                ax.plot(x[i], j, linewidth=linewidth, color=colormap[count])   # The x and y values to plot
                count = count + 1
        else:
            #exitu plot
            maxindex = np.where(matrix[i][0]==max(matrix[i][0]))
            color = cmp.color_order[i]
            ax.plot(x[i], matrix[i][0], linewidth=linewidth, color=color, label=sampleNames[i])
            if markWhiteline == True:
                ax.vlines(x=x[i][maxindex], ymin=0, ymax=max(matrix[i][0]), colors=color, ls='--', lw=linewidth*3/5)



    ax.set_title(title, fontsize = 26)    
    ax.set_xlabel("Energy [keV]", fontsize = 26)
    ax.set_ylabel("Normalized absorption [edge fraction]", fontsize = 26)
    ax.tick_params(axis='both',which="major", width=framewidth, length=framewidth*3, labelsize=26)
    ax.spines["top"].set_linewidth(framewidth)
    ax.spines["bottom"].set_linewidth(framewidth)
    ax.spines["right"].set_linewidth(framewidth)
    ax.spines["left"].set_linewidth(framewidth)

    ax.legend(prop={'size': 18,},loc='lower right').get_frame().set_linewidth(0.0)

    if maxE != None and minE !=None:
        plt.xlim(minE,maxE)
    elif maxE != None and minE == None:
        plt.xlim(min_x,maxE)
    elif minE != None and maxE == None:
        plt.xlim(minE,max_x)
    elif maxE == None and minE == None:
        plt.xlim(min_x,max_x)
    plt.ylim(bottom=0)
    plt.tight_layout()
    if savefig == "png":
        plt.savefig(outputPath + "/1d_"+ plotName +".png", dpi=600, bbox_inches = "tight")
    if savefig == "svg":
        plt.savefig(outputPath + "/1d_"+ plotName +".svg", dpi=600, bbox_inches = "tight")
    plt.show()


#----------------------------------------------------------------------------------------#
#                         Finds the preedge for inhouse data                             #
#----------------------------------------------------------------------------------------#

def preedge(x,matrix,pre_min,pre_max,post_min,post_max):
    """
    Arguments
    ----------
    y:  Intensity array (as a slice of the Intesnsity matrix)
    """ 
    #Find e0
    d2x_y = np.gradient(matrix,2)
    index_e0 = np.where(d2x_y==max(d2x_y))
    e0 = x[index_e0]

    preStartindex = None
    preEndindex = None
    postStartindex = None
    postEndindex = None

    for j in x:
        if j >= pre_min and preStartindex==None:
            preStartindex = int(np.where(x==j)[0]) 
        if j > pre_max and preEndindex==None:
            preEndindex = int(np.where(x==j)[0]) - 1 
        if j >= post_min and postStartindex==None:
            postStartindex = int(np.where(x==j)[0]) 
        if j > post_max and postEndindex==None:
            postEndindex = int(np.where(x==j)[0]) - 1  
    pre_x = x[preStartindex:preEndindex]
    pre_y = matrix[preStartindex:preEndindex]
    post_x = x[postStartindex:postEndindex]
    post_y = matrix[postStartindex:postEndindex] 

    edgejump = matrix[int(np.where(x==e0)[0])]

    pre_m,pre_b = np.polyfit(pre_x, pre_y, 1)
    post_m,post_b = np.polyfit(post_x, post_y, 1)
    prefit = (x*pre_m + pre_b) + edgejump 
    postfit = x*post_m + post_b
    return e0, prefit, postfit

#----------------------------------------------------------------------------------------#
#                                 Normalises inhouse data                                #
#----------------------------------------------------------------------------------------#

def norm(x,matrix,Energys): 
    """
    Arguments
    ----------
    x:          Energy   
    matrix:     Intensity
    Energys:    Dictionary containing delta_e0, pre_min, pre_max, post_min, post_max  
    """ 
    #Find e0 and pre and post edge slope
    delta_e0 = Energys["delta_e0"]
    pre_min = Energys["pre_min"]		    # pre edge minimum energy for fit
    pre_max = Energys["pre_max"]			    # pre edge maximum energy for fit
    post_min = Energys["post_min"]				# post edge minimum energy for fit
    post_max = Energys["post_max"]          # post edge maximum energy for fit

    e0, prefit, postfit = preedge(x=x,matrix=matrix,pre_max=pre_max,pre_min=pre_min,post_max=post_max,post_min=post_min)

    #normalise spectrum to pre and post edge slope
    index_min_matrix = np.where(matrix==min(matrix))
    min_matrix = min(matrix) / prefit[index_min_matrix]
    print(min_matrix)
    for j in range(0,len(x)):
        if x[j] < e0:
            matrix[j] = (matrix[j] / prefit[j]) - min_matrix
        elif x[j] >= e0:
            matrix[j] = matrix[j] / postfit[j]
    x = x + delta_e0
    return x, matrix

#----------------------------------------------------------------------------------------#
#                           Set temp regression values in console                        #
#----------------------------------------------------------------------------------------#

def timeInput():
    """
    Arguments
    ----------
    None 
    """ 
    print("\n Set Times for linear regression of reaching set temperature")
    print("------------------------------------------------------------------")
    rampStart = float(input(" rampStart:  "))
    rampEnd = float(input(" rampEnd:    "))
    constStart = float(input(" constStart: "))
    constEnd = float(input(" constEnd:   "))

    Time = {"rampStart": rampStart,
            "rampEnd": rampEnd,
            "constStart": constStart,
            "constEnd": constEnd,}

    return Time

#----------------------------------------------------------------------------------------#
#                  Ensures a output directory will be created on demand                  #
#----------------------------------------------------------------------------------------#

def ensure_dir(dir_path):
    """
    Arguments
    ----------
    dir_path:   Path to directory   
    """ 
    import os
    if not os.path.exists(dir_path):
        os.makedirs(dir_path)


def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx]