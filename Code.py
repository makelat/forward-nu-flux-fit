# Load Libraries and Define General Functions

import sys
import os
import numpy as np
import matplotlib
from matplotlib import pyplot as plt
import uproot
import json
import math
from scipy import interpolate
from scipy.stats import poisson
import datetime

### ---------------------------------------------------------
### input/output

def readfile_txt(filename):
    list_of_lists = []
    with open(filename) as f:
        for line in f:
            if len(line)==0: continue
            if line[0]=="#":continue
            inner_list = [float(elt.strip()) for elt in line.split( )]
            if len(inner_list)==0: continue
            list_of_lists.append(inner_list)
    return np.array(list_of_lists)

def read_and_interpolate(filename,mass):
    br=readfile_txt(filename)
    return np.interp(mass, br.T[0],br.T[1])
    
def save_json(filename, file):
    jason = json.dumps(file)
    f = open(filename+".json","w")
    f.write(jason)
    f.close()
    
def load_json(filename):
    f = open(filename+".json")
    return json.load(f)
    
### ---------------------------------------------------------
### Setup Experiments

# length [cm], rho [g/cm3], posz [m], x/y [m]
experiments = {
    "FASERv"   : {"length": 100, "rho": 19.3, "material": "W", "posz": 480, "lumi":150,
                  "xmin": -.125, "xmax": .125, "ymin": -.125, "ymax": .125 },
    "FASERv2"  : {"length": 660., "rho": 19.3, "material": "W", "posz": 620, "lumi":3000,
                  "xmin": -.2, "xmax": .2, "ymin": -.2, "ymax": .2 },
    "NuTeV@FPF-Pb" : {"length": 400., "rho": 11.29, "material": "Pb", "posz": 620, "lumi":3000,
                      "xmin": -.5, "xmax": .5, "ymin": -1.5, "ymax": 1.5 },
    "NuTeV@FPF-Fe" : {"length": 400., "rho": 7.874, "material": "Fe", "posz": 620, "lumi":3000,
                      "xmin": -.5, "xmax": .5, "ymin": -1.5, "ymax": 1.5 },
    "SND"      : {"length": 27.25, "rho": 19.3, "material": "W", "posz": 480, "lumi":150,
                  "xmin": 0.080, "xmax": 0.470, "ymin": 0.155, "ymax": 0.545 },
    "AndSND"   : {"length": 30., "rho": 19.3, "material": "W", "posz": 620, "lumi":3000,
                  "xmin": 0.10, "xmax": 0.70, "ymin": -0.30, "ymax": 0.30 },
    "FLARE"    : {"length": 700, "rho": 1.3954, "material": "Ar", "posz": 620, "lumi":3000,
                  "xmin": -0.50, "xmax": 0.50, "ymin": -0.50, "ymax": 0.50 },
    "FLARE100" : {"length": 3000, "rho": 1.3954, "material": "Ar", "posz": 620, "lumi":3000,
                  "xmin": -0.80, "xmax": 0.80, "ymin": -0.80, "ymax": 0.80 },
    "FCCFLARE" : {"length": 700, "rho": 1.3954, "material": "Ar", "posz": 1500, "lumi":30000,
                  "xmin": -0.50, "xmax": 0.50, "ymin": -0.50, "ymax": 0.50 },
}

### ---------------------------------------------------------
### Interaction Probability

def get_interaction_probability(experiment,energy,pid,procin=""):

    # length and densitity
    rho = experiments[experiment]["rho"] # [g/cm3]
    lengthcm = experiments[experiment]["length"] # [cm]
    material = experiments[experiment]["material"]
    mN = 1.6605*10.**(-24) # Atomic mass unit in [g]
        
    proc = procin
    fscltag=''  #Final state charged lepton tag
    for fscl in ['_ee','_emu','_etau','_mutau','_tautau']:
        if fscl in procin:
            fscltag = fscl
            proc = proc.replace(fscl,'')
            break
   
    # path to cross section per nucleon, default/main notebook case
    xsecfilename = "vxsecs/xs_"+material+"_"+pid+".txt"
    
    # path to cross section per nucleon, trident/dimuon cases
    if len(procin)!=0:
        procdir = proc
        proctag = proc
        if '4F' in proc: 
            procdir = 'tridents'
            proctag = '4F_dGV-0_dGA-0'  #Baseline delta g_V, delta g_A values, scaled later in model
        xsecfilename = 'files/trident_mc/'+procdir+'/xsec/cm2_nucleon/'+proctag+'_'+material+'_'+pid+fscltag+'.txt'
        # if material-specific xsec not found, try proton xsec
        # n.b. this does not need to be scaled by atomic mass, 
        # since the x-sec is per nucleon not nucleus
        if not os.path.exists(xsecfilename):
            xsecfilenamep = 'files/trident_mc/'+procdir+'/xsec/cm2_nucleon/'+proctag+'_p_'+pid+fscltag+'.txt'
            if os.path.exists(xsecfilenamep): xsecfilename = xsecfilenamep
            else: print('ERROR: neither '+xsecfilename+' nor '+xsecfilenamep+' found')
    
    # interpolate cross section value
    sigma = read_and_interpolate(xsecfilename,energy)
       
    # interaction probability
    prob_interact = sigma*rho*lengthcm/mN
    return prob_interact

### ---------------------------------------------------------
### Function to remove i:th col &row from array for projections
def dimDel(arr, i):
    ret = np.delete(arr, i, 0)
    ret = np.delete(ret, i, 1)
    return ret

### ---------------------------------------------------------
### Function To Convert Text File into ROOT file

def convert_to_root(inputfiles, filename="example.root"):
    
    # create root file
    file = uproot.recreate(filename)
    
    # initialize flux tree
    file["flux"] = {"wgt": [], "vtxx": [], "vtxy": [], "vtxz": [],
        "dist": [], "px": [],  "py": [], "pz": [], "E": [], "pdg": []}
    
    # initialze meta data
    maxEnergy, minWgt, maxWgt = 0, 10.**10, 0
    
    # loop through files and fill flux tree
    for inputfile, scale in inputfiles:
        # load file
        data = readfile_txt(inputfile)
        # write to file
        file["flux"].extend({
            "wgt"    : data.T[8]*scale,
            "vtxx"   : data.T[2]*1000.,
            "vtxy"   : data.T[3]*1000.,
            "vtxz"   : np.zeros(len(data)),
            "dist"   : data.T[4]*1000.,
            "px"     : data.T[5]*data.T[7],
            "py"     : data.T[6]*data.T[7],
            "pz"     : data.T[7]*np.sqrt(1-data.T[5]**2-data.T[6]**2),
            "E"      : data.T[7],
            "pdg"    : data.T[0],
        })
        emax, minw, maxw = max(data.T[7]), min(data.T[8]), max(data.T[8])
        if emax>maxEnergy: maxEnergy=emax
        if minw<minWgt:    minWgt=minw
        if maxw>maxWgt:    maxWgt=maxw
    
    # write meta tree
    file["meta"] = {
        "pdglist"   : [[12,-12,14,-14,16,-16]],
        "maxEnergy" : [maxEnergy],
        "minWgt"    : [minWgt],
        "maxWgt"    : [maxWgt],
        "protons"   : [0.001],
        "windowBase": [[-0.5, -0.5, 0.0]],
        "windowDir1": [[ 1.0,  0.0, 0.0]],
        "windowDir2": [[ 0.0,  1.0, 0.0]],
        "auxintname": [[]],
        "auxdblname": [[]],
        "infiles"   : [[]],
        "seed"      : [0],
        "metakey"   : [0],
    }

### ---------------------------------------------------------
### Inits required for decay length calculations

clight = 2.9979e8 #m/s

Dmesons = ['D+','D_s+','D0']
Dch     = ['D+','D_s+']
DchLtX  = [r'$D^{+}$', r'$D_{s}^{+}$']
Dnt     = ['D0']
DntLtX  = [r'$D^{0}$']

#Average lifetimes and uncertainties [s] (PDG)
tau={}
tau['D+'          ] = [1040.0e-15, 7.0e-15]
tau['D-'          ] = tau['D+']
tau['D0'          ] = [410.0e-15, 1.5e-15]
tau['Dbar0'       ] = tau['D0']
tau['D_s+'        ] = [500.0e-15, 7.0e-15]
tau['D_s-'        ] = tau['D_s+']
tau['omega'       ] = [7.58e-23, 0.11e-23]
tau['Lambda0'     ] = [2.631e-10, 0.020e-10]
tau['Lambdabar0'  ] = tau['Lambda0']
tau['Lambda_bbar0'] = [1409e-15, 55e-15]
tau['Lambda_c+'   ] = [200e-15, 6e-15]
tau['Lambda_cbar-'] = tau['Lambda_c+']
tau['eta'         ] = [5.0e-19, 0.3e-19]
tau["eta'"        ] = [3.2e-21, 0.2e-21]
tau['rho0'        ] = [4.5e-24, 0.0]                                     
tau['mu-'         ] = [21969811e-13, 22e-13]
tau['mu+'         ] = tau['mu-']
tau['tau-'        ] = [2.903e-13, 0.005e-13]
tau['tau+'        ] = tau['tau-']
tau['Sigma-'      ] = [1.479e-10, 0.011e-10]
tau['Sigmabar+'   ] = tau['Sigma-']
tau['phi'         ] = [1.55e-22, 0.01e-22]
tau['Omega_c0'    ] = [268e-15, 24e-15]
tau['B-'          ] = [1638e-15, 4e-15]
tau['B+'          ] = tau['B-']
tau['B0'          ] = [1519e-15, 4e-15]
tau['Bbar0'       ] = tau['B0']
tau['B_sbar0'     ] = [1515e-15, 4e-15]
tau['Xi-'         ] = [1.639e-10,0.015e-10]
tau['Xibar+'      ] = tau['Xi-']
tau['Xi_c+'       ] = [456e-15, 5e-15]
tau['Xi_c0'       ] = [153e-15, 6e-15]
tau['Xi_cbar0'    ] = tau['Xi_c0']
tau['Xi_bbar0'    ] = [1480e-15, 30e-15]


#Masses [GeV]
mass={}
mass['D+'          ] = 1.86958
mass['D-'          ] = mass['D+']
mass['D0'          ] = 1.86483
mass['Dbar0'       ] = mass['D0']
mass['D_s+'        ] = 1.96827
mass['D_s-'        ] = mass['D_s+']
mass['omega'       ] = 0.78266
mass['Lambda0'     ] = 1.115683
mass['Lambdabar0'  ] = mass['Lambda0']
mass['Lambda_bbar0'] = 5.61960
mass['Lambda_c+'   ] = 2.28646
mass['Lambda_cbar-'] = mass['Lambda_c+']
mass['eta'         ] = 0.547862
mass["eta'"        ] = 0.95778
mass['rho0'        ] = 0.77
mass['mu-'         ] = 0.1056583755
mass['mu+'         ] = mass['mu-']
mass['tau-'        ] = 1.77686
mass['tau+'        ] = mass['tau-']
mass['Sigma-'      ] = 1.197449
mass['Sigmabar+'   ] = mass['Sigma-']
mass['phi'         ] = 1.019461
mass['Omega_c0'    ] = 2.6975
mass['B-'          ] = 5.27934
mass['B+'          ] = mass['B-']
mass['B0'          ] = 5.27965
mass['Bbar0'       ] = mass['B0']
mass['B_sbar0'     ] = 5.36688
mass['Xi-'         ] = 1.32171
mass['Xibar+'      ] = mass['Xi-']
mass['Xi_c+'       ] = 2.46794
mass['Xi_c0'       ] = 2.47090
mass['Xi_cbar0'    ] = mass['Xi_c0']
mass['Xi_bbar0'    ] = 5.7919

def decayLength(name,px,py,pz):
    if name not in list(tau.keys()):
        print('WARNING: Mean lifetime unspecified for '+str(name))
        return -1
    if name not in list(mass.keys()):
        print('WARNING: mass unspecified for '+str(name))
        return -1
    #Smear mean lifetime by the given uncertainties
    rndm = np.random.normal()
    lifetime = tau[name][0] + rndm*tau[name][1]
    momentum = np.sqrt(px**2 + py**2 + pz**2)
    return lifetime*clight*momentum/mass[name]

### ---------------------------------------------------------
### Compute N bin uncertainty
def uncbands(axes,ebins,rbins,ecenters,detector,baseline,infomat,cid):
    
    global vpidstrs
    global Nptypes
    global Ngp
    
    #Init
    uncertainty = [0.0 for _ in vpidstrs]
    Nrbins = len(rbins)-1
    values, vectors = linalg.eig(infomat)
    
    #Add all r-bins's contributions to the uncertainty for each vpid
    for rbin in range(Nrbins):           
        rmin, rmax = rbins[rbin], rbins[rbin+1] 
        for ivpidsgn,vpid in enumerate(vpidstrs):
            base_entry = np.array(baseline[vpid]['n_int']).T[rbin]
            for value,vector in zip(values, vectors.T):
                if value<=0.0: continue  #Skip flat directions
                point=vector/np.sqrt(value)
                lambdavars = [point[sum([Ngp[m]-1 for m in range(n)]):sum([Ngp[m]-1 for m in range(n+1)])]\
                              for n in range(Nptypes)]
                lambdavars = [np.add(lambdavars[i],lambdaBL[i]) for i in range(len(lambdavars))]
                varied = model(detector=detector,radN='_rad'+str(Nrbins), lambdamat=lambdavars,cid=cid)
                entry = np.array(varied[vpid]['n_int']).T[rbin]
                uncertainty[ivpidsgn] += (entry - base_entry)**2            
    
    #Combine uncertainties for N radial bins
    for ivpidsgn,vpid in enumerate(vpidstrs):  #N.B. possibly 6 vpids here
        base_entry_sum = np.array(baseline[vpid]['n_int']).T[0]
        for rbin in range(1,Nrbins):
            base_entry_sum = np.add(base_entry_sum, np.array(baseline[vpid]['n_int']).T[rbin])
        uncertainty[ivpidsgn] = np.sqrt(uncertainty[ivpidsgn])
        wgtp = 1+uncertainty[ivpidsgn]/base_entry_sum   #Uncertainties for fractional...
        wgtm = 1-uncertainty[ivpidsgn]/base_entry_sum   #...uncs in lower panel
        errhi = np.multiply(wgtp,base_entry_sum)  #Uncertainties for spectra...
        errlo = np.multiply(wgtm,base_entry_sum)  #...histos (upper panel)

        #Plot N radial bin combined uncertainty
        ivpid = int(np.floor(0.5*int(ivpidsgn))) if cid else ivpidsgn
        lstyle = 'dotted' if '-' in vpid else 'solid'
        lwidth = 1.4 if '-' in vpid else 0.7
        #TODO plotTest
        axes[1,ivpid].hist(x=ecenters,weights=wgtp,bins=ebins,\
                           histtype='step',color='black',ls=lstyle,lw=lwidth)
        #TODO plotTest
        axes[1,ivpid].hist(x=ecenters,weights=wgtm,bins=ebins,\
                           histtype='step',color='black',ls=lstyle,lw=lwidth)
        axes[0,ivpid].bar(x=ebins[:-1],\
                          height=np.add(errhi,np.multiply(-1,errlo)),\
                          bottom=errlo,\
                          width=np.diff(ebins),\
                          align='edge', linewidth=0, color='gray',\
                          alpha=0.25, zorder=-1)
        axes[1,ivpid].bar(x=ebins[:-1],\
                          height=np.add(wgtp,np.multiply(-1,wgtm)),\
                          bottom=wgtm,\
                          width=np.diff(ebins),\
                          align='edge', linewidth=0, color='gray',\
                          alpha=0.25, zorder=-1)    

### ---------------------------------------------------------
### AUX function to format error messages in get_llr and terminate
def llrERROR(key,errtype,pardict):
    print('ERROR llr '+key+' is '+errtype+' for:')
    for key in list(pardict.keys()):
        print(key+'='+str(pardict[key]))
    sys.exit()

### ---------------------------------------------------------
### Fetch interpolation function for given acceptance table (Emu,efficiency)
def acctableIP(acctable):
    f = open(acctable, "r")
    Emu,eff=[],[]
    lines=[]
    lines=f.readlines()
    lines.pop(0)    #Remove header line
    for l in lines:
        lsplit = l.split()
        Emu.append(float(lsplit[0]))
        eff.append(float(lsplit[1]))
    f.close()
    #Return an interpolation function giving efficiency for input Emu
    fip = interpolate.interp1d(Emu,eff,fill_value='extrapolate')
    return fip

### ---------------------------------------------------------
### Fetch low-E FS muon and IS neutrino energies per event from P8 output
def EnuEmuDimuons():

    #Init
    Enu=[]
    Emu=[]
    inEvt = False       #True when currently reading an event  
    
    #Read Pythia output
    f = open('files/trident_mc/dimuons/pythia_output_20230406/0000_all_particles.txt', "r")
    lines=[]
    lines=f.readlines()
    f.close()
    
    for i in range(len(lines)):
            
        if 'PYTHIA Event Listing  (complete event)' in lines[i]:
            EnuEvt = 0.
            EmuEvt = []
            inEvt = True  #Event-reading phase entered
            mudtrs=[]     #Store mu dtrs' number in P8 output, to not double-count mus emitting radiation
            continue        
        if inEvt:
            if 'Charge sum:' not in lines[i]:
                lsplit = lines[i].split()
                if not len(lsplit)==0 and 'id  name' not in lines[i]:
                    pdgid  = int(lsplit[1])
                    name   = lsplit[2]
                    status = int(lsplit[3])
                    #Incoming particles have status +-12, check incident E_nu 
                    if abs(status)==12 and 'nu' in name: 
                        EnuEvt = float(lsplit[13])
                    #Count muons: ensure considered events have >=2 muons to be dimuon
                    if '(' not in name and abs(pdgid)==13 and lsplit[0] not in mudtrs: 
                        EmuEvt.append(float(lsplit[13]))  #Store muon energies
                        mudtrs.append(lsplit[6])       #Store indices of muons' daughter...
                        mudtrs.append(lsplit[7])       #...particles in the event
            #Finished reading this event, store or discard before moving to next one
            else:
                inEvt = False
                #Skip lines between 'Charge sum' line and possible 'iEvent' id-assignment line
                while i<len(lines) and 'End PYTHIA Event Listing' not in lines[i-1]: i+=1
                #Ensure this event was assigned a non-unique id in the file, i.e. labeled as an N-muon event
                if 'iEvent' in lines[i] and len(EmuEvt) >= 2:
                    #Make event ids unique by combining them with E-bin id
                    Enu.append(EnuEvt)
                    Emu.append(min(EmuEvt))  #Store smallest E mu
    return [Enu,Emu]

### ---------------------------------------------------------
### Fetch low-E FS muon and IS neutrino energies per event from trident generated events
def EnuEmuTridents(nE,detector='FASERv2',nuid='14'):
    material = experiments[detector]['material']
    nustr='nu'
    if int(nuid)<0: nustr += 'bar'
    if   abs(int(nuid))==12: nustr += 'el'
    elif abs(int(nuid))==16: nustr += 'tau'
    else:                    nustr += 'mu'
    Enu,Emu=[],[]
    
    #Generated trident events are given in separate files per incoming Enu, loop over them
    for iE in range(1,nE+1):
        striE = str(iE)
        if iE<10: striE = '0' + striE
        fname = 'files/trident_mc/tridents/generated_events/'+material+'_SM_'+ nustr + '_' + striE + '.txt'
        #Ensure file is not empty
        if (os.stat(fname).st_size == 0):
            print('Trident file empty at iE = '+str(iE))
            return [Enu,Emu]
        
        #Here we only consider the production of a mu-mubar pair
        lepton1 = int(nuid)/abs(int(nuid))*13
        lepton2 = -1*lepton1
        #If we wouldn't only consider dimuon production, could deduce lepton 1 & 2 PDG IDs (left4ref):
        #lepton1 = int(nuid)/abs(int(nuid))*(abs(int(nuid))-1)
        #lepton2 = -1*lepton1
        
        f = open(fname, "r")
        lines=[]
        lines=f.readlines()
        f.close()
        inEvt = False  #Flag if currently reading event lines
        mu1found = False
        mu2found = False
        Enufound = False
        p4_1 = [0.,0.,0.,0.]
        p4_2 = [0.,0.,0.,0.]
        for line in lines:
            #Check where event listings start and end
            if  '<event>' in line: 
                inEvt = True
                mu1found = False
                mu2found = False
                Enufound = False
                p4_1 = [0.,0.,0.,0.]  #1st muon 4-momentum
                p4_2 = [0.,0.,0.,0.]  #2nd muon 4-momentum
                continue
            if '</event>' in line: 
                inEvt = False
                continue
            if inEvt:
                lsplit = line.split()
                if lsplit[0]==nuid and lsplit[1]=='-1' and not Enufound: 
                    Enu.append(float(lsplit[5]))
                    Enufound = True
                elif int(lsplit[0])==lepton2 and lsplit[1]=='1' and not mu2found: 
                    mu2found = True
                    p4_2 = [float(lsplit[2]), float(lsplit[3]), float(lsplit[4]), float(lsplit[5])]
                elif int(lsplit[0])==lepton1 and lsplit[1]=='1' and not mu1found: 
                    mu1found = True
                    p4_1 = [float(lsplit[2]), float(lsplit[3]), float(lsplit[4]), float(lsplit[5])]
                else: continue
                if mu1found and mu2found:
                    Emu.append(min(p4_1[3],p4_2[3]))  #Store lower-E mu only
    return [Enu,Emu]

### ---------------------------------------------------------
### Fetch incoming nu E, deduce bin boundaries based on given central values
def fetchEbins(proc):
    lines=[]
    Evals = []
    if 'dimuons' in proc:
        f = open('files/trident_mc/dimuons/pythia_output_20230406/0000_Ebins_vpTOlj_for_dimu_v2.txt', "r")
        lines=f.readlines()
        for line in lines:
            if 'nevents' in line: continue
            Evals.append(float(line.split()[2]))
        f.close()
    else: Evals = np.loadtxt('files/trident_mc/tridents/xsec/cm2_nucleon/tridents_W_14.txt',usecols=0)
    Evalsbd = [0.5*(Evals[i]+Evals[i+1]) for i in range(len(Evals)-1)]
    Evalsbd = [2.0*Evals[0]-Evalsbd[0]] + Evalsbd      #Add lowest bin bd.
    Evalsbd = Evalsbd + [2.0*Evals[-1] - Evalsbd[-1]]  #Add highest bin bd.
    return [Evals,Evalsbd]

def efficiency(proc,acctable,detector,nuid='14'):
    print('Finding efficiency interpolation function for ',proc)
    #Step 1: weighted vs unweighted histos
    Evals,Ebins = fetchEbins(proc)
    Enu,Emu=[],[]
    if 'dimuons' in proc: Enu,Emu = EnuEmuDimuons()
    else: Enu,Emu = EnuEmuTridents(nE=len(Evals),detector=detector,nuid=nuid)
    accIP = acctableIP(acctable)
    effs = np.array([accIP(E) for E in Emu])
    hwgtd   = plt.hist(x=Enu,weights=effs,bins=Ebins,histtype='step',color='black',ls='solid',lw=1.2,label='weighted')
    hunwgtd = plt.hist(x=Enu,             bins=Ebins,histtype='step',color='red',  ls='solid',lw=1.2,label='unweighted')
    plt.title(proc)
    plt.xlabel(r'$E_\nu$ [GeV]')
    plt.xscale('log')
    plt.legend(frameon=False)
    plt.savefig('plots/efficiency_1_wgt_vs_unwgt_counts_unnormalized_'+detector+'_'+proc+'_'+nuid+'.pdf')
    plt.clf()
    
    #Step 2 obtain ratio from the histos
    ratio = [hwgtd[0][i]/hunwgtd[0][i] if hunwgtd[0][i]>0 else 0 for i in range(len(hwgtd[0]))]
    plt.plot(Evals,ratio)
    plt.title(proc)
    plt.xscale('log')
    plt.ylabel('Acceptance')
    plt.xlabel(r'$E_\nu$ [GeV]')
    plt.savefig('plots/efficiency_2_ratio_'+detector+'_'+proc+'_'+nuid+'.pdf')
    ip_eff = interpolate.interp1d(Evals,ratio,fill_value='extrapolate')
    return ip_eff

### ---------------------------------------------------------
### Find the required number of events for a model to be distinct from the
### SM according Poisson distribution/difference in the number of events
#Param  detectorlist  Considered detector( combination)s: e.g. ['det1'] or ['det1','det2']...
#       SM_nints      The number of events according to SM
#       CL            Desired confidence level, e.g. 0.95 for 95% CL
def findPoissonConstraint(detectorlist,SM_nints,CL,dmu_i=1.0):
    #Find the constraint for the number of events (muconstr)
    ktmp = max(round(SM_nints),0)
    dmu = dmu_i
    muconstr = ktmp
    converged = False
    approxInfty = 10000
    nrange = range(ktmp,approxInfty) if dmu>0 else range(ktmp)
    while not converged:
        muconstrlast = muconstr
        muconstr += dmu
        P = sum([poisson.pmf(mu=muconstr,k=n) for n in nrange])
        #When requested CL reached, go back and reduce step size to achieve higher accuracy
        if P > CL:
            if abs(dmu)<0.01:
                converged = True
            else:
                dmu *= 0.1
                muconstr = muconstrlast
    print('&'.join(detectorlist)+' SM expected n.o. events='+str(ktmp)\
          +', dmu='+str(dmu)+', constraint='+str(muconstr)+', final P = '+str(P))
    return muconstr

### ---------------------------------------------------------
### Test automation for checking that revisions reproduce previous plots
#Param  plotdata    an object containing the data to be printed
#       testtag     the printout filename
#       testsubdir  printouts produced into tests/testsubdir; e.g. notebook name
#       mode        'plot'/'histo', assumes plotdata is plt.plot/hist return object.
#                   If unspecified, assumes a list of vectors/lists, output as is
def plotTest(plotdata,testtag,testsubdir='',mode=''):
    
    #Ensure relevant test directory structure exists
    testdir = 'tests/'+testsubdir
    if testdir[-1]!='/': testdir+='/'
    testdirsplit = [t for t in testdir.split('/') if len(t)!=0]
    for i in range(1,len(testdirsplit)+1):
        currentpath = '/'.join(testdirsplit[:i])
        if not os.path.exists(currentpath): 
            os.mkdir(currentpath)
            
    #The first output is named "*_REFERENCE"; compare later-produced outputs with that
    suffix = '.txt'
    initmode = True
    ref = '_REFERENCE'
    outputpath = testdir + testtag
    if os.path.exists(outputpath+ref+suffix): 
        initmode = False
        ref=''
    
    #Write x and y axes's data to output file, one column for each
    f = open(outputpath+ref+suffix,"w")
    plotdataT = plotdata  #Will be transposed/formatted if need be
    if mode == 'plot': plotdataT = np.array(plotdata[-1].get_data()).T
    elif mode == 'histo':
        binslo = plotdata[1][:-1]
        binshi = plotdata[1][1:]
        vals   = plotdata[0]
        plotdataT = [binslo, binshi, vals]
    for vec in plotdataT:
        for x in vec:
            f.write(str(x)+"    ")
        f.write('\n')
    f.close()

    #If a reference file already exists, compare current output to it
    if not initmode:
        ref = np.loadtxt(outputpath + '_REFERENCE' + suffix).flatten()
        dat = np.loadtxt(outputpath                + suffix).flatten()
        ref = np.array([r for r in ref if not (np.isinf(r) or np.isnan(r))])
        dat = np.array([d for d in dat if not (np.isinf(d) or np.isnan(d))])
        f = open('tests.log','a')
        if False in np.isclose(ref,dat).flatten():
            f.write(str(datetime.datetime.now())+' WARNING! test '+outputpath+' failed!\n')
        else: f.write(str(datetime.datetime.now())+' test '+outputpath+' passed\n')
        f.close()
