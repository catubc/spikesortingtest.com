import struct
import csv
import numpy as np
import os
from math import *
import pandas

from django.http import HttpResponse
import PIL, PIL.Image, StringIO

colors = ['blue','red', 'black']

class Tsf_file(object):

    def __init__(self, fin, sim_dir):
        
        self.fin = fin
        self.sim_dir = sim_dir
        
        #self.read_tsf(fin,sim_dir)

    def read_tsf(self,fin,sim_dir):

        self.header = fin.read(16)
        self.iformat = struct.unpack('i',fin.read(4))[0]
        self.SampleFrequency = struct.unpack('i',fin.read(4))[0]
        self.n_electrodes = struct.unpack('i',fin.read(4))[0]
        self.n_vd_samples = struct.unpack('i',fin.read(4))[0]
        self.vscale_HP = struct.unpack('f',fin.read(4))[0]

        if self.iformat==1001:
            self.Siteloc = np.zeros((2*56), dtype=np.int16)
            self.Siteloc = struct.unpack(str(2*56)+'h', fin.read(2*56*2))
        if self.iformat==1002:
            self.Siteloc = np.zeros((2*self.n_electrodes), dtype=np.int16)
            self.Readloc = np.zeros((self.n_electrodes), dtype=np.int32)
            for i in range(self.n_electrodes):
                self.Siteloc[i*2] = struct.unpack('h', fin.read(2))[0]
                self.Siteloc[i*2+1] = struct.unpack('h', fin.read(2))[0]
                self.Readloc[i] = struct.unpack('i', fin.read(4))[0]

        #print self.n_electrodes, self.n_vd_samples
        
        #self.ec_traces = []  #Save dummy list in case accidentally read elsewhere
        
        #return  #EXIT EARLY - DO NOT LOAD .tsf file unless needed for computing ptp values

        self.ec_traces =  np.fromfile(fin, dtype=np.int16, count=self.n_electrodes*self.n_vd_samples)
        self.ec_traces.shape = self.n_electrodes, self.n_vd_samples

        self.n_cell_spikes = struct.unpack('i',fin.read(4))[0]
        #print "No. ground truth cell spikes: ", self.n_cell_spikes
        if (self.n_cell_spikes>0):
            if (self.iformat==1001):
                self.vertical_site_spacing = struct.unpack('i',fin.read(4))[0]
                self.n_cell_spikes = struct.unpack('i',fin.read(4))[0]

            self.fake_spike_times =  np.fromfile(fin, dtype=np.int32, count=self.n_cell_spikes)
            self.fake_spike_assignment =  np.fromfile(fin, dtype=np.int32, count=self.n_cell_spikes)
            self.fake_spike_channels =  np.fromfile(fin, dtype=np.int32, count=self.n_cell_spikes)

    def compute_ptp(self, units, n_electrodes, ec_traces, vscale_HP):

        import ptps
        #Hacked way of making a numpy array out of uneven list; using trailing zeros
        max_len = 0
        for i in range(len(units)):
            if max_len< len(units[i]): max_len=len(units[i])
        spikes = np.zeros((len(units),max_len), dtype=np.int32)
        for i in range(len(units)):
            spikes[i][0:len(units[i])]=units[i]

        max_ptp = np.zeros(len(units))
        max_ch = np.int32((len(units),0))
        
        #print spikes[27]
        #print len(spikes[27])
        #quit()
        ec_traces = np.int8(ec_traces)  #<--- convert data to int*2 for smaller storage...MAKE SURE OK!!!
        out_array =  ptps.ptps(spikes, n_electrodes, ec_traces, vscale_HP)

        self.ptp = []
        self.maxchan=[]
        for i in range(len(units)):
            self.ptp.append(out_array[0][i])
            self.maxchan.append(out_array[1][i][0])


#class Loadcsv_vispy(object):
    #"""Ground truth file - Espen, Torbjorn, Gaute Vispy data"""
    #def __init__(self, sorted_file, data_dir, tsf):
    
        #self.load_groundtruth(sorted_file, data_dir, tsf)
    
    #def load_groundtruth(self,sorted_file, data_dir, tsf):
        
        ##NB: Ensure that 3 decimal places spiketimes (i.e. ms precision) is sufficient
        #fname = data_dir+sorted_file
        #f = open(fname+'_ground_truth.csv', "r")
        #data_file = np.genfromtxt(f, delimiter=',', dtype=np.float32) #, usecols=(0,))
        #f.close()
        
        ##randomized_filename = fname.replace('/home/catubc/webapps/my_django_app/myproject/media/uploads/','')
        ##temp_file = '/home/catubc/webapps/my_django_app/myproject/media/saved_uploads/'+work_dir+randomized_filename
        ##np.savetxt(temp_file, data_file)

        ##Load spiketimes from first column
        #spiketimes=data_file[:,0]
        #spiketimes= np.array(spiketimes*tsf.SampleFrequency, dtype=np.float32) #Convert spiketimes to timesteps

        ##Load unit ids from second column
        #spike_id=np.array(data_file[:,1],dtype=np.int32)

        ##NB: Extra steps needed for assigning map KK discountinuos unit ids onto 0 based sequential ids
        #unique_ids = np.unique(spike_id)

        ##print "Parsing Dan .csv file -> assigning to units"
        #self.units=[]
        #self.maxchan=[]
        #for i in range(len(unique_ids)):
            #indexes = np.where(spike_id==unique_ids[i])[0]
            #self.units.append(spiketimes[indexes])
            #self.maxchan.append(int(np.rint(np.average(maxchan[indexes]))))


        #for i in range(len(spiketimes)):
            #self.units[np.where(unique_ids==spike_id[i])[0][0]].append(spiketimes[i])
        #self.n_units=len(self.units)

        #self.size=np.zeros((self.n_units), dtype=np.float32)
        #for i in range(len(self.units)):
            #self.size[i]=len(self.units[i])
        ##print "No. of units in .csv: ", self.n_units

        #if False:  #Set this to false for now; do not load .tsf file
            #if (os.path.exists(fname[0:-4]+'_ptps.csv')==False):
                ##Compute PTP values from original .tsf file; needed for .csv sorted files
                #print "Manually recomputing each unit ptp and max channel values"
                #self.ptp=np.zeros((self.n_units), dtype=np.float32)
                #self.size=np.zeros((self.n_units), dtype=np.float32)
                #for i in range(self.n_units):
                    #self.size[i]=len(self.units[i])

                ##Use fortran to compute PTP very quickly
                #tsf.compute_ptp(self.units, tsf.n_electrodes, tsf.ec_traces, tsf.vscale_HP)
                #self.ptp=tsf.ptp
                #self.maxchan=tsf.maxchan
            #else:
                #self.ptp = np.loadtxt(fname[:-4]+'_ptps.csv', dtype=np.float32, delimiter=",")
                #self.maxchan = np.loadtxt(fname[:-4]+'_maxch.csv', dtype=np.int32, delimiter=",")
                #self.size = np.loadtxt(fname[:-4]+'_size.csv', dtype=np.int32, delimiter=",")

        #self.n_sorted_spikes = [None]*self.n_units
        #for k in range(self.n_units):
            #self.n_sorted_spikes[k] = len(self.units[k])
        

class Loadptcs(object):
    """Polytrode clustered spikes file neuron record"""
    def __init__(self, sorted_file, work_dir, ptcs_flag, save_timestamps):
        
        f = open(work_dir+sorted_file+'.ptcs', "rb")
        self.sorted_file = sorted_file
        self.name = sorted_file
        self.full_path = work_dir + sorted_file
        # call the appropriate method:
        self.VER2FUNC = {1: self.readHeader, 2: self.readHeader, 3: self.readHeader}

        self.readHeader(f, ptcs_flag, save_timestamps)
        
        self.nid = []  #Make unique unit id list for loading later.
        
        self.loadData(self.nsamplebytes, f, work_dir, ptcs_flag, save_timestamps)
        
        f.close()

    def __getstate__(self):
        """Instance methods must be excluded when pickling"""
        d = self.__dict__.copy()
        try: del d['VER2FUNC']
        except KeyError: pass
        return d

    def readHeader(self, f, ptcs_flag, save_timestamps):
        """Read in neuron record of .ptcs file version 3. 'zpos' field was replaced
        by 'sigma' field.
        nid: int64 (signed neuron id, could be -ve, could be non-contiguous with previous)
        ndescrbytes: uint64 (nbytes, keep as multiple of 8 for nice alignment, defaults to 0)
        descr: ndescrbytes of ASCII text
        (padded with null bytes if needed for 8 byte alignment)
        clusterscore: float64
        xpos: float64 (um)
        ypos: float64 (um)
        sigma: float64 (um) (Gaussian spatial sigma)
        nchans: uint64 (num chans in template waveforms)
        chanids: nchans * uint64 (0 based IDs of channels in template waveforms)
        maxchanid: uint64 (0 based ID of max channel in template waveforms)
        nt: uint64 (num timepoints per template waveform channel)
        nwavedatabytes: uint64 (nbytes, keep as multiple of 8 for nice alignment)
        wavedata: nwavedatabytes of nsamplebytes sized floats
        (template waveform data, laid out as nchans * nt, in uV,
        padded with null bytes if needed for 8 byte alignment)
        nwavestdbytes: uint64 (nbytes, keep as multiple of 8 for nice alignment)
        wavestd: nwavestdbytes of nsamplebytes sized floats
        (template waveform standard deviation, laid out as nchans * nt, in uV,
        padded with null bytes if needed for 8 byte alignment)
        nspikes: uint64 (number of spikes in this neuron)
        spike timestamps: nspikes * uint64 (us, should be sorted)
        """

        self.nid = int(np.fromfile(f, dtype=np.int64, count=1)) # nid
        self.ndescrbytes = int(np.fromfile(f, dtype=np.uint64, count=1)) # ndescrbytes
        self.descr = f.read(self.ndescrbytes).rstrip('\0 ') # descr

        if self.descr:
            try:
                self.descr = eval(self.descr) # might be a dict
            except: pass

        self.nneurons = int(np.fromfile(f, dtype=np.uint64, count=1)) # nneurons
        self.nspikes = int(np.fromfile(f, dtype=np.uint64, count=1)) # nspikes
        self.nsamplebytes = int(np.fromfile(f, dtype=np.uint64, count=1)) # nsamplebytes
        self.samplerate = int(np.fromfile(f, dtype=np.uint64, count=1)) # samplerate
        self.npttypebytes = int(np.fromfile(f, dtype=np.uint64, count=1)) # npttypebytes

        self.pttype = f.read(self.npttypebytes).rstrip('\0 ') # pttype

        self.nptchans = int(np.fromfile(f, dtype=np.uint64, count=1)) # nptchans
        self.chanpos = np.fromfile(f, dtype=np.float64, count=self.nptchans*2) # chanpos
        self.chanpos.shape = self.nptchans, 2 # reshape into rows of (x, y) coords
        self.nsrcfnamebytes = int(np.fromfile(f, dtype=np.uint64, count=1)) # nsrcfnamebytes
        self.srcfname = f.read(self.nsrcfnamebytes).rstrip('\0 ') # srcfname
        # maybe convert this to a proper Python datetime object in the Neuron:
        self.datetime = float(np.fromfile(f, dtype=np.float64, count=1)) # datetime (days)
        self.ndatetimestrbytes = int(np.fromfile(f, dtype=np.uint64, count=1)) # ndatetimestrbytes
        self.datetimestr = f.read(self.ndatetimestrbytes).rstrip('\0 ') # datetimestr

    def loadData(self, n_bytes, f, work_dir, ptcs_flag, save_timestamps):
        #call the appropriate method:
        #self.VER2FUNC = {1: self.read_ver_1, 2:self.read_ver_2, 3:self.read_ver_3}
        self.nsamplebytes = n_bytes
        self.wavedtype = {2: np.float16, 4: np.float32, 8: np.float64}[self.nsamplebytes]

        self.n_units=self.nneurons
        self.units=[None]*self.n_units
        self.uid = [None]*self.n_units  #Unique id for full track sorts
        self.n_sorted_spikes = [None]*self.n_units
        self.ptp=np.zeros((self.n_units), dtype=np.float32)
        self.size = []
        self.maxchan = []

        for k in range(self.n_units):
            self.readUnit(f,work_dir, ptcs_flag)
            self.units[k]= self.spikes
            if 'nick' in self.full_path:
                self.uid[k]= self.nid-1
            elif 'martin' in self.full_path:
                self.uid[k]= self.nid

            if ptcs_flag: #Martin's data has wrong flag for saves
                self.units[k]=[x*self.samplerate/1E+6 for x in self.units[k]] #Converts spiketimes from usec to timesteps
            else:
                self.units[k]=[x*self.samplerate/2/1E+6 for x in self.units[k]] #Converts spiketimes from usec to timesteps

            self.n_sorted_spikes[k] = len(self.units[k])
            self.size.append(self.nspikes)
            self.maxchan.append(self.maxchanu)
            self.ptp[k]=max(self.wavedata[np.where(self.chans==self.maxchanu)[0][0]]) - \
                        min(self.wavedata[np.where(self.chans==self.maxchanu)[0][0]]) #compute PTP of template;
        f.close()

        if (os.path.exists(work_dir+self.sorted_file+'_ptps.csv')==False):
            np.savetxt(work_dir+self.sorted_file+'_ptps.csv', self.ptp, delimiter=",")
            np.savetxt(work_dir+self.sorted_file+'_size.csv', self.size, delimiter=",")

    def readUnit(self,f, work_dir, ptcs_flag):
        self.nid = int(np.fromfile(f, dtype=np.int64, count=1)) # nid
        self.ndescrbytes = int(np.fromfile(f, dtype=np.uint64, count=1)) # ndescrbytes
        self.descr = f.read(self.ndescrbytes).rstrip('\0 ') # descr

        if self.descr:
            try:
                self.descr = eval(self.descr) # might be a dict
            except: pass

        self.clusterscore = float(np.fromfile(f, dtype=np.float64, count=1)) # clusterscore
        self.xpos = float(np.fromfile(f, dtype=np.float64, count=1)) # xpos (um)
        self.ypos = float(np.fromfile(f, dtype=np.float64, count=1)) # ypos (um)
        self.zpos = float(np.fromfile(f, dtype=np.float64, count=1)) # zpos (um)
        self.nchans = int(np.fromfile(f, dtype=np.uint64, count=1)) # nchans
        self.chans = np.fromfile(f, dtype=np.uint64, count=self.nchans) #NB: Some errors here from older .ptcs formats
        self.maxchanu = int(np.fromfile(f, dtype=np.uint64, count=1)) # maxchanid

        self.nt = int(np.fromfile(f, dtype=np.uint64, count=1)) # nt: number of time points in template

        self.nwavedatabytes, self.wavedata = self.read_wave(f) #TEMPLATE

        self.nwavestdbytes, self.wavestd = self.read_wave(f) #STANDARD DEVIATION
        self.nspikes = int(np.fromfile(f, dtype=np.uint64, count=1)) # nspikes

        # spike timestamps (us):
        self.spikes = np.fromfile(f, dtype=np.uint64, count=self.nspikes)

        # convert from unsigned to signed int for calculating intervals:
        self.spikes = np.asarray(self.spikes, dtype=np.float64)

    def read_wave(self, f):
        """Read wavedata/wavestd bytes"""
        # nwavedata/nwavestd bytes, padded:
        nbytes = int(np.fromfile(f, dtype=np.uint64, count=1))
        fp = f.tell()
        count = nbytes // self.nsamplebytes # trunc to ignore any pad bytes
        X = np.fromfile(f, dtype=self.wavedtype, count=count) # wavedata/wavestd (uV)
        if nbytes != 0:
            X.shape = self.nchans, self.nt # reshape
        f.seek(fp + nbytes) # skip any pad bytes
        return nbytes, X

    def rstrip(s, strip):
        """What I think str.rstrip should really do"""
        if s.endswith(strip):
            return s[:-len(strip)] # strip it
        else:
            return s

    def read(self):
        self.nid = self.parse_id()
        with open(self.fname, 'rb') as f:
            self.spikes = np.fromfile(f, dtype=np.int64) # spike timestamps (us)
        self.nspikes = len(self.spikes)

    def read_tsf(self,f):
        pass

class Loadcsv(object):

    def __init__(self, fname, Sort1, work_dir):

        self.fname = fname
        self.work_dir = work_dir
        self.loadSpikes(Sort1)

    def loadSpikes(self, Sort1):

        #Permanent save data file in "saved_uploads" media directory
        import shutil
        #NB: THIS FILE NAME SHOULD BE ASSIGNED DYNAMICALLY; Just removing all path except "randomized filename" part
        randomized_filename = self.fname.replace('/home/catubc/webapps/my_django_app/myproject/media/uploads/','')
        randomized_filename = randomized_filename.replace('../static/'+self.work_dir,'')
        temp_file = '/home/catubc/webapps/my_django_app/myproject/media/saved_uploads/'+self.work_dir +randomized_filename
        
        print "Saving files to saved_uploads: ", self.fname, temp_file
        shutil.copy(self.fname, temp_file)

        #NB: Ensure that 3 decimal places spiketimes (i.e. ms precision) is sufficient
        f = open(self.fname, "r")
        data_file = np.genfromtxt(f, delimiter=',', dtype=np.float32) #, usecols=(0,))
        f.close()
        
        #randomized_filename = fname.replace('/home/catubc/webapps/my_django_app/myproject/media/uploads/','')
        #temp_file = '/home/catubc/webapps/my_django_app/myproject/media/saved_uploads/'+work_dir+randomized_filename
        #np.savetxt(temp_file, data_file)

        #SPIKETIMES - 1ST COLUMN
        spiketimes=data_file[:,0]
        spiketimes= np.array(spiketimes*Sort1.samplerate, dtype=np.float32) #Convert spiketimes to timesteps

        #UNIT IDS - 2ND COLUMN
        spike_id=np.array(data_file[:,1],dtype=np.int32)

        #MAX CHS - 3RD COLUMN
        maxchan=np.array(data_file[:,2],dtype=np.int32)-1

        #NB: EXTRA STEP needed for assigning map KK discountinuos unit ids onto 0 based sequential ids
        unique_ids = np.unique(spike_id)
        self.n_units=len(unique_ids)

        #PARSING
        self.units=[]
        self.maxchan=[]
        self.ptp=[]
        for i in range(len(unique_ids)):
            print "Unit: ", i
            indexes = np.where(spike_id==unique_ids[i])[0]
            temp_spikes = spiketimes[indexes]; self.units.append(temp_spikes)
            #temp_maxchan = max(set(maxchan[indexes]), key=list(maxchan[indexes]).count)-1; 
            numbers =  np.unique(maxchan[indexes], return_index=True, return_counts=True)       #Return array with indexes, counts, number of counts
            #print maxchan[indexes]
            #print numbers
            temp_maxchan = numbers[0][np.argmax(numbers[2])]
            if i in [0,1,2,3,74]: 
                print numbers
            print "temp_maxchan: ", temp_maxchan,
            self.maxchan.append(temp_maxchan)
            
            temp_ptp = self.compute_ptp_maxchan(i, temp_spikes, temp_maxchan); self.ptp.append(temp_ptp)
            print "temp_ptp: ", temp_ptp
            print "\n\n"



        self.size=np.zeros((self.n_units), dtype=np.float32)
        for i in range(len(self.units)):
            self.size[i]=len(self.units[i])
        #print "No. of units in .csv: ", self.n_units

        #if False:  #Set this to false for now; do not load .tsf file
            #if (os.path.exists(fname[0:-4]+'_ptps.csv')==False):
                ##Compute PTP values from original .tsf file; needed for .csv sorted files
                #print "Manually recomputing each unit ptp and max channel values"
                #self.ptp=np.zeros((self.n_units), dtype=np.float32)
                #self.size=np.zeros((self.n_units), dtype=np.float32)
                #for i in range(self.n_units):
                    #self.size[i]=len(self.units[i])

                ##Use fortran to compute PTP very quickly
                #tsf.compute_ptp(self.units, tsf.n_electrodes, tsf.ec_traces, tsf.vscale_HP)
                #self.ptp=tsf.ptp
                #self.maxchan=tsf.maxchan
            #else:
                #self.ptp = np.loadtxt(fname[:-4]+'_ptps.csv', dtype=np.float32, delimiter=",")
                #self.maxchan = np.loadtxt(fname[:-4]+'_maxch.csv', dtype=np.int32, delimiter=",")
                #self.size = np.loadtxt(fname[:-4]+'_size.csv', dtype=np.int32, delimiter=",")

        self.n_sorted_spikes = [None]*self.n_units
        for k in range(self.n_units):
            self.n_sorted_spikes[k] = len(self.units[k])


    def compute_ptp_maxchan(self, unit, spikes, maxchan):
        #Load maxchan only from tsf file and compute PTP 
        
        fin = open('/home/catubc/webapps/my_django_app/myproject/static/' + self.work_dir + self.work_dir[:-1]+'_truth_filtered.tsf', "rb")
        self.header = fin.read(16)
        self.iformat = struct.unpack('i',fin.read(4))[0]
        self.SampleFrequency = struct.unpack('i',fin.read(4))[0]
        self.n_electrodes = struct.unpack('i',fin.read(4))[0]
        self.n_vd_samples = struct.unpack('i',fin.read(4))[0]
        self.vscale_HP = struct.unpack('f',fin.read(4))[0]

        self.Siteloc = np.zeros((2*self.n_electrodes), dtype=np.int16)
        self.Readloc = np.zeros((self.n_electrodes), dtype=np.int32)
        for i in range(self.n_electrodes):
            self.Siteloc[i*2] = struct.unpack('h', fin.read(2))[0]
            self.Siteloc[i*2+1] = struct.unpack('h', fin.read(2))[0]
            self.Readloc[i] = struct.unpack('i', fin.read(4))[0]
        
        #print self.Siteloc
        #print self.Readloc
        
        indent = 16+20+self.n_electrodes*8
        #import matplotlib.pyplot as plt
        #print maxchan
        for p in range(1):
            chan_indent =  maxchan # IMPORTANT TO USE 1 BASED INDEXES;

            #ax=plt.subplot(4,4,p+1)
            #ax.yaxis.set_ticks([])
            #ax.xaxis.set_ticks([])
            #plt.ylim(-250,250)
            #chan_indent = maxchan-1]+p
            #plt.title("CH: "+str(chan_indent))

            fin.seek(indent+chan_indent*2*self.n_vd_samples, os.SEEK_SET)         #Not 100% sure this indent is correct.
            ec_traces =  np.fromfile(fin, dtype=np.int16, count=self.n_vd_samples)

        
            #Compute ptp for a single channel
            ptp=0; counter=0.
            t = np.arange(-10,10,1)
            for k in spikes:#[0:100]:
                if k>10: 
                    temp_trace = ec_traces[int(k)-10:int(k)+10]
                    ptp += self.vscale_HP*(max(temp_trace)-min(temp_trace))
                    #plt.plot(t,temp_trace)
                    counter+=1.
            ptp = ptp/counter
            print "Ch: ", str(chan_indent), " ptp: ", ptp, " Readloc: ", self.Readloc[chan_indent]

        #plt.show()
        
        fin.close()
        return ptp
        
            
def Compare_sorts(Sort1, Sort2):

    import compare_sort
    #Hacked way of making a numpy array out of uneven list; using trailing zeros
    max_len = 0
    for i in range(len(Sort1.units)):
        if max_len< len(Sort1.units[i]): max_len=len(Sort1.units[i])
    spikes1 = np.zeros((len(Sort1.units),max_len), dtype=np.int32)
    for i in range(len(Sort1.units)):
        spikes1[i][0:len(Sort1.units[i])]=Sort1.units[i]

    spikes1=np.asfortranarray(spikes1)

    max_len = 0
    for i in range(len(Sort2.units)):
        if max_len< len(Sort2.units[i]): max_len=len(Sort2.units[i])
    spikes2 = np.zeros((len(Sort2.units),max_len), dtype=np.int32)
    for i in range(len(Sort2.units)):
        spikes2[i][0:len(Sort2.units[i])]=Sort2.units[i]

    spikes2=np.asfortranarray(spikes2)

    #Siteloc = np.asfortranarray(tsf.Siteloc)
    Siteloc = np.asfortranarray(Sort1.chanpos)

    print "... compare sort... into fortran..."
    out_array =  compare_sort.compare_sort(Sort1.samplerate, spikes1, spikes2, Siteloc, Sort1.maxchan, Sort2.maxchan)
    print "... out of fortran"
    
    Sort2.purity = out_array[0][0:len(Sort2.units)]
    Sort2.completeness= out_array[1][0:len(Sort2.units)]
    Sort2.match = out_array[2][0:len(Sort2.units)]

    return

def Plot_Composite_metrics(data_dir, sorted_file, sorted_dir, Sort2):
    #Plot_Composite_metrics(sim_dir, sorted_dirs, sorted_file, anonymous)
    #return Plot_Composite_metrics(data_dir, sorted_file, fname_save, Sort2)

    '''Compute metrics from sorted data and generate sorter .jpg file
    also compute and update .jpg files to be used by global metrics webpage 
    '''
    colors=['green','blue','violet','cyan','violet','mediumvioletred','indianred','lightsalmon','pink','darkolivegreen']

    purity=np.array(Sort2.purity)
    completeness = np.array(Sort2.completeness)
    ptp = np.array(Sort2.ptp)

    #Ordered data
    indexes = np.argsort(ptp)
    ptp_ordered = ptp[indexes]
    purity_ordered = purity[indexes]
    completeness_ordered = completeness[indexes]

    #************************************** MATPLOTLIB PLOTS ****************************************
    from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
    from matplotlib.figure import Figure
    import django
    import matplotlib.gridspec as gridspec
    
    fig=Figure()
    fig.set_size_inches(15, 15)
    font_size = 8
    
    gs = gridspec.GridSpec(21, 19)
    
    #******************* PLOT PURITY & COMPLETENESS *******************
    if True:
        ax = fig.add_subplot(gs[5:8, :])
        ax.set_ylim(0, 100)
        #weights=['bold','bold','bold','normal','normal','normal','normal','normal','normal','normal']
        text_out = "Individual Sorting Metrics \n\nPurity & Completeness (Units in Upload Order)"
        ax.set_title("Individual Sorting Metrics \n\nPurity & Completeness (Units in Upload Order)", fontweight='bold', fontsize=font_size+10)
        
        ax.set_ylabel('Percent (%)', fontsize=font_size+10)
        svm = []
        n_units = []

        for i in range(len(purity)):
            ax.bar(1+i*15, purity[i], 5, color=colors[0], alpha=1)
            ax.bar(6+i*15, completeness[i], 5, color=colors[1], alpha=1)
            svm.append(purity[i]*completeness[i]*1E-2)
        ax.set_xlim(0, (i+1)*15)

        x = np.arange(10,(i+1)*15+6,15)
        labels = np.arange(1,len(purity)+1,1)
        ax.set_xticks(x)
        ax.set_xticklabels(labels, fontsize=8)

        #plot 80% and 90% lines
        x = (0,(i+1)*15)
        a = (80,80)
        ax.plot(x,a, 'r--', color='black', linewidth=1)
        x = (0,(i+1)*15)
        a = (90,90)
        ax.plot(x,a, 'r--', color='black', linewidth=1)
      

    #******************* PLOT SVM - UPLOAD ORDER *******************
    if True:
        ax = fig.add_subplot(gs[9:12, :])
        ax.set_ylim(0, 100)
        ax.set_title('Single Value Metric (Purity * Completeness; Units in Upload Order)', fontsize=font_size+10)
        ax.set_ylabel('Percent (%)', fontsize=font_size+10)

        for i in range(len(svm)):
            ax.bar(1+i*15, svm[i], 10, color=colors[3], alpha=1)
        ax.set_xlim(0, (i+1)*15)

        x = np.arange(10,(i+1)*15+6,15)
        labels = np.arange(1,len(purity)+1,1)
        ax.set_xticks(x)
        ax.set_xticklabels(labels, fontsize=8)

        #plot 80% and 90% lines
        x = (0,(i+1)*15)
        a = (80,80)
        ax.plot(x,a, 'r--', color='black', linewidth=1)
        x = (0,(i+1)*15)
        a = (90,90)
        ax.plot(x,a, 'r--', color='black', linewidth=1)
        x = (0,(i+1)*15)
        a = (50,50)
        ax.plot(x,a, 'r--', color='black', linewidth=1)


    #******************* PLOT PURITY VS PTP *******************
    if True:
        ax = fig.add_subplot(gs[13:16, :])

        ax.set_ylim(0, 100)
        ax.set_xlim(0,ptp_ordered[-1]+10)
        ax.set_title('Purity vs PTP(uV) (Units in PTP Order)', fontsize=font_size+10)
        ax.set_ylabel('Percent (%)', fontsize=font_size+10)

        x = np.arange(0,max(ptp_ordered))
        #labels = np.arange(1,len(purity)+1,1)
        #ax.set_xticks(x)
        #ax.set_xticklabels(labels)
        
        #Plot 80% and 90% dashed lines
        x = (0,ptp_ordered[-1]+10)
        a = (80,80)
        ax.plot(x,a, 'r--', color='black', linewidth=1)
        ax.get_xaxis().set_visible(True)
        x = (0,ptp_ordered[-1]+10)
        a = (90,90)
        ax.plot(x,a, 'r--', color='black', linewidth=1)

        for i in range(len(purity_ordered)):
            ax.bar(ptp_ordered[i], purity_ordered[i], 2, color=colors[0], alpha=1)

    if True:
        ax = fig.add_subplot(gs[17:20, :])
        ax.set_xlim(0,ptp_ordered[-1]+10)

        ax.set_ylim(0, 100)
        ax.set_title('Completeness vs PTP(uV) (Units in PTP Order)', fontsize=font_size+10)
        ax.set_ylabel('Percent (%)', fontsize=font_size+10)
        
        #Plot 80% and 90% lines
        x = (0,ptp_ordered[-1]+10)
        a = (80,80)
        ax.plot(x,a, 'r--', color='black', linewidth=1)
        x = (0,ptp_ordered[-1]+10)
        a = (90,90)
        ax.plot(x,a, 'r--', color='black', linewidth=1)

        for i in range(len(purity)):
            ax.bar(ptp_ordered[i], completeness_ordered[i], 2, color=colors[1], alpha=1)


    #********************************** GLOBAL METRICS ************************************
    random_filename = sorted_dir.replace('/home/catubc/webapps/my_django_app/myproject/media/uploads/','')
    save_dir = '/home/catubc/webapps/my_django_app/myproject/media/saved_metrics/'+sorted_file+'/'

    #Load existing purity and completeness files: 
    import glob
    #MAKE SURE ORDERING DATA FILES OTHERWISE THEY MAY NOT MATCH UP LATER
    purity_files = sorted(glob.glob(save_dir +"*_purityx"))
    completeness_files = sorted(glob.glob(save_dir+"*_completenessx"))
    n_spikes_files = sorted(glob.glob(save_dir+"*_spikesx"))

    purity_array = []
    completeness_array = []
    for file_1, file_2 in zip(purity_files, completeness_files):
        purity_array.append(np.loadtxt(file_1))
        completeness_array.append(np.loadtxt(file_2))

    n_spikes_array = []

    for file_ in n_spikes_files:
        n_spikes_array.append(np.loadtxt(file_))
        
    print "Len purity array: ", len(purity_array)
    print "Len completeness array: ", len(completeness_array)
    
    #Save purity, completeness, n_spikes, into .../saved_metrics/... #DO IT after loading existing metrics
    #NB: Use "purityx" as filename in rare circumstances that someone uploads files using these terms
    np.savetxt(save_dir+random_filename+"_purityx", purity)
    np.savetxt(save_dir+random_filename+"_completenessx", completeness)
    n_spikes=0
    for k in range(len(Sort2.units)):
        n_spikes+=len(Sort2.units[i])
    f = open(save_dir+random_filename+"_spikesx", 'w')
    f.write(str(n_spikes))
    f.close()

    
    #******************* SINGLE VALUE METRIC 
    ax = fig.add_subplot(gs[0:3, 0:3])
    #Compute SVM for current sort
    total_svm = np.sum(np.array(svm))
    
    #Compute SVM for previous sorts
    svm_plot=[]
    svm_plot.append(total_svm) #Must append here in case value is off the plot
    for i in range(len(purity_array)):
        #For some reason some of the files have been corrupt and len purity array not eq completeness
        if len(purity_array[i])==len(completeness_array[i]):
            svm_plot.append(np.sum(np.array(purity_array[i])*np.array(completeness_array[i]))*1E-2)
    
    #Histogram all svm data
    svm_plot = np.array(svm_plot)
    bins_x = np.linspace(0,np.max(svm_plot)*1.1,25)
    hist, bins = np.histogram(svm_plot, bins = bins_x)
    width = 0.7 * (bins[1] - bins[0])
    center = (bins[:-1] + bins[1:]) / 2
    
    ax.bar(center, hist, align='center', width=width)

    #Plot current svm marker - red bar
    ax.axvspan(xmin=total_svm,xmax=total_svm+10, color='red')

    x_ticks=np.linspace(0,max(svm_plot)*1.1,4)
    
    ax.set_xticks(x_ticks)
    ax.tick_params(labelsize=8)
    ax.set_yticks([])
    ax.set_xlim(left=0)
    ax.set_ylim(bottom=0)
    ax.set_title('SVM: '+str(int(total_svm)), fontsize=font_size+5)

    #Save/Update SVM - global metrics
    if True: 
        font_size_jpg = 25
        fig_save=Figure()
        fig_save.set_size_inches(5, 5)
        
        ax_save = fig_save.add_subplot(111)
        #fig.set_canvas(ax_save.gcf().canvas)
        ax_save.bar(center, hist, align='center', width=width)

        ax_save.set_xticks(x_ticks)
        ax_save.tick_params(labelsize=font_size_jpg)
        ax_save.set_yticks([])
        ax_save.set_xlim(left=0)
        ax_save.set_ylim(bottom=0)
        ax_save.set_title('SVM', fontsize=font_size_jpg+10)

        canvas=FigureCanvas(fig_save)
        fig_save.savefig('/home/catubc/webapps/static_media/images/metrics/'+sorted_file+"_svm.png", format='png')
        #fig_save.close()
        
    #******************* # UNITS SORTED **************************
    ax = fig.add_subplot(gs[0:3, 4:7])

    #Histogram n_unit data
    n_units_array = []
    n_units_array.append(len(Sort2.units))
    for i in range(len(purity_array)):
        n_units_array.append(len(purity_array[i]))
        
    n_units_array = np.array(n_units_array)
    bins_x = np.linspace(0,np.max(n_units_array)*1.1,25)
    hist, bins = np.histogram(n_units_array, bins = bins_x)
    width = 0.7 * (bins[1] - bins[0])
    center = (bins[:-1] + bins[1:]) / 2
    ax.bar(center, hist, align='center', width=width)

    #Plot current svm marker - red bar
    ax.axvspan(xmin=len(Sort2.units),xmax=len(Sort2.units)+1, color='red')

    import math
    x_ticks=np.linspace(0,int(math.ceil(max(n_units_array) / 10.0)) * 10,4, dtype=np.int16)
    ax.set_xticks(x_ticks)
    ax.tick_params(labelsize=8)
    #ax.set_tick_params(axis='both', which='both', labelsize=8)
    
    ax.set_xlim(left=0)
    ax.set_ylim(bottom=0)
    ax.set_yticks([])
    ax.set_title('#Units: '+str(len(Sort2.units)), fontsize=font_size+5)

    #Save/Update SVM - global metrics
    if True: 
        fig_save=Figure()
        fig_save.set_size_inches(5, 5)
        ax_save = fig_save.add_subplot(111)

        ax_save.bar(center, hist, align='center', width=width)
        ax_save.set_xticks(x_ticks)
        ax_save.tick_params(labelsize=font_size_jpg)
        ax_save.set_xlim(left=0)
        ax_save.set_ylim(bottom=0)
        ax_save.set_yticks([])
        ax_save.set_title('#Units', fontsize=font_size_jpg+10)

        canvas=FigureCanvas(fig_save)
        fig_save.savefig('/home/catubc/webapps/static_media/images/metrics/'+sorted_file+"_n_units.png", format='png')

    #******************* # SPIKES SORTED 
    ax = fig.add_subplot(gs[0:3, 8:11])
    
    #Histogram n_spikes data
    n_spikes_array.append(n_spikes)
    n_spikes_array = np.array(n_spikes_array)
    bins_x = np.linspace(0,100000,25)
    hist, bins = np.histogram(n_spikes_array, bins = bins_x)
    width = 0.7 * (bins[1] - bins[0])
    center = (bins[:-1] + bins[1:]) / 2
    ax.bar(center, hist, align='center', width=width)

    #Plot current svm marker - red bar
    ax.axvspan(xmin=n_spikes,xmax=n_spikes+10, color='red')
    
    #Hardwire xticks here... good idea!?
    x_ticks=np.linspace(0,100000,4)

    #x_ticks=np.linspace(0,int(math.ceil(n_spikes / 10000.0)) * 10000,4)
    ax.set_xticks(x_ticks)
    ax.set_xlim(0,100000)
    ax.tick_params(labelsize=8)
    #ax.set_tick_params(axis='both', which='both', labelsize=8)
    
    ax.set_yticks([])
    #ax.set_xlim(left=0)
    ax.set_ylim(bottom=0)
    ax.set_title('#Spikes: '+str(n_spikes), fontsize=font_size+5)

    #Save/Update SVM - global metrics
    if True: 
        fig_save=Figure()
        fig_save.set_size_inches(5, 5)
        ax_save = fig_save.add_subplot(111)

        ax_save.bar(center, hist, align='center', width=width)
        #Hard wire 100,000 spikes as max # spikes
        #x_ticks=np.linspace(0,100000,4)
        ax_save.set_xticks(x_ticks)
        ax_save.tick_params(labelsize=font_size_jpg-7)
        ax_save.set_yticks([])
        ax_save.set_xlim(0,100000)
        ax_save.set_ylim(bottom=0)
        ax_save.set_title('#Spikes', fontsize=font_size_jpg+10)

        canvas=FigureCanvas(fig_save)
        fig_save.savefig('/home/catubc/webapps/static_media/images/metrics/'+sorted_file+"_n_spikes.png", format='png')

    #******************* AVE. UNIT PURITY 
    ax = fig.add_subplot(gs[0:3, 12:15])
    
    #Histogram n_spikes data
    purity_array.append(Sort2.purity)
    purity_ave = []
    for i in range(len(purity_array)):
        purity_ave.append(np.average(purity_array[i]))

    purity_ave = np.array(purity_ave)
    bins_x = np.linspace(0,np.max(purity_ave)*1.1,25)
    hist, bins = np.histogram(purity_ave, bins = bins_x)
    width = 0.7 * (bins[1] - bins[0])
    center = (bins[:-1] + bins[1:]) / 2
    ax.bar(center, hist, align='center', width=width)

    #Plot current svm marker - red bar
    ax.axvspan(xmin=np.average(Sort2.purity),xmax=np.average(Sort2.purity)+2, color='red')
        
    x_ticks=np.linspace(0,100,4, dtype=np.int16)
    ax.set_xticks(x_ticks)
    ax.tick_params(labelsize=8)
    ax.set_xlim(0,100)
    
    ax.set_yticks([])
    ax.set_title('Ave Purity: '+str(int(np.average(Sort2.purity)))+"%", fontsize=font_size+5)

    #Save/Update SVM - global metrics
    if True: 
        fig_save=Figure()
        fig_save.set_size_inches(5, 5)
        ax_save = fig_save.add_subplot(111)

        ax_save.bar(center, hist, align='center', width=width)
        ax_save.set_xticks(x_ticks)
        ax_save.tick_params(labelsize=font_size_jpg)
        ax_save.set_yticks([])
        ax_save.set_title('Ave Purity (%)', fontsize=font_size_jpg+10)
        ax_save.set_xlim(0,100)

        canvas=FigureCanvas(fig_save)
        fig_save.savefig('/home/catubc/webapps/static_media/images/metrics/'+sorted_file+"_purity.png", format='png')

    #******************* AVE UNIT COMPLETENESS 
    ax = fig.add_subplot(gs[0:3, 16:19])
    
    #Histogram n_spikes data
    completeness_array.append(Sort2.completeness)
    completeness_ave = []
    for i in range(len(completeness_array)):
        completeness_ave.append(np.average(completeness_array[i]))

    completeness_ave = np.array(completeness_ave)
    bins_x = np.linspace(0,np.max(completeness_ave)*1.1,25)
    hist, bins = np.histogram(completeness_ave, bins = bins_x)
    width = 0.7 * (bins[1] - bins[0])
    center = (bins[:-1] + bins[1:]) / 2
    ax.bar(center, hist, align='center', width=width)

    #Plot current svm marker - red bar
    ax.axvspan(xmin=np.average(Sort2.completeness),xmax=np.average(Sort2.completeness)+2, color='red')
        
    x_ticks=np.linspace(0,100,4, dtype=np.int16)
    ax.set_xticks(x_ticks)
    ax.tick_params(labelsize=8)
    ax.set_xlim(0,100)

    ax.set_yticks([])
    ax.set_title('Ave Completeness: '+str(int(np.average(Sort2.completeness)))+"%", fontsize=font_size+5)
    

    #Save/Update SVM - global metrics
    if True: 
        fig_save=Figure()
        fig_save.set_size_inches(5, 5)
        ax_save = fig_save.add_subplot(111)

        ax_save.bar(center, hist, align='center', width=width)
        ax_save.set_xticks(x_ticks)
        ax_save.tick_params(labelsize=font_size_jpg)
        ax_save.set_yticks([])
        ax_save.set_title('Ave Completeness (%)', fontsize=font_size_jpg+5)
        ax_save.set_xlim(0,100)

        canvas=FigureCanvas(fig_save)
        fig_save.savefig('/home/catubc/webapps/static_media/images/metrics/'+sorted_file+"_completeness.png", format='png')

    #***************** SHIP IT!!
    
    fig.suptitle("Sorted File: " + random_filename+"\n\nGlobal Metrics & Distributions", fontweight = 'bold', fontsize = font_size +10)
    canvas=FigureCanvas(fig)
    #response=django.http.HttpResponse(content_type='image/png')
    #fig.savefig(response, format='png')
    #return response

    #Save output .jpg to file and send name of location
    saved_output_jpg = '/home/catubc/webapps/static_media/saved_sort_images/' + random_filename+'.jpg'
    fig.savefig(saved_output_jpg)
    
    #Save purity and completeness .txt files
    np.savetxt('/home/catubc/webapps/static_media/saved_sort_images/purity_' + random_filename+'.txt', np.array(Sort2.purity))
    np.savetxt('/home/catubc/webapps/static_media/saved_sort_images/completeness_' + random_filename+'.txt', np.array(Sort2.completeness))

    return random_filename
    

    #******************************* BOKEH PLOTS ************************************
    if False:
        
        from bokeh.plotting import figure
        from bokeh.resources import CDN
        from bokeh.embed import file_html
        from bokeh.io import gridplot, output_file, show
        from bokeh.charts import Bar, output_file, show, color, defaults
        from pandas import DataFrame

        #*********** Purity
        dict = {'values':purity, "color": colors[0]}
        df = DataFrame(dict)
        ax1 = Bar(df, values='values', color='color',title="Purity (Upload Order)", xlabel='Sorted Units', ylabel='%', plot_width=700, plot_height=250)

        #*********** Purity
        #Order data by ptp
        indexes = np.argsort(ptp)
        ptp = np.round_(ptp[indexes],2)
        purity = np.round_(purity[indexes],2)
        
        dict = {'values':purity, 'label':ptp, "color": colors[1]}
        df = DataFrame(dict)
        ax2 = Bar(df, label='label', values='values', color='color',title="Purity (PTP Order)", xlabel='Unit PTP (uV)', ylabel='%', plot_width=700, plot_height=250)
        
        #********** Completeness
        dict = {'values':completeness, "color": colors[2]}
        df = DataFrame(dict)
        ax3 = Bar(df, values='values', color='color',title="Completeness", xlabel='Sorted units', ylabel='%', plot_width=700, plot_height=250)

        #********** Completeness
        completeness = np.round_(completeness[indexes],2)
        dict = {'values':completeness,'label':ptp, "color": colors[3]}
        df = DataFrame(dict)
        ax4 = Bar(df, label='label', values='values', color='color',title="Completeness (PTP Order)", xlabel='Sorted units', ylabel='%', plot_width=700, plot_height=250)

        #Send bokeh data to gridplot
        p = gridplot([[ax1], [ax2], [ax3], [ax4]]) #Can display these plots side-by-side e.g. [ax1,ax2]


        html = file_html(p, CDN, "my plot")
    
        return HttpResponse(html)

    ##SAVED WORKING FORMAT OF PLOTTING
    ##from bokeh.plotting import figure
    ##from bokeh.resources import CDN
    ##from bokeh.embed import file_html
    ##plot = figure()
    ##plot.circle([1,2], [3,4])
    ##html = file_html(fig, CDN, "my plot")
    ##return HttpResponse(html)

def Compute_global_metrics():
    pass

def find_nearest(array,value,dist):
    #print "Looking for: ", value
    #print "Diff to closest value: ", array[np.abs(array-value).argmin()]-value
    if abs(array[np.abs(array-value).argmin()]-value) < dist:
        idx = (np.abs(array-value)).argmin()
        return idx
    else:
        return None

