# ## Load Libraries and Define General Functions

import sys
import os
import numpy as np
import matplotlib
from matplotlib import pyplot as plt
import uproot
import json

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

def get_interaction_probability(experiment,energy, pid):

    # length and densitity
    rho = experiments[experiment]["rho"] # [g/cm3]
    lengthcm = experiments[experiment]["length"] # [cm]
    material = experiments[experiment]["material"]
    mN = 1.6605*10.**(-24) # [g]
        
    # cross section
    sigma = read_and_interpolate("vxsecs/xs_"+material+"_"+pid+".txt",energy)
       
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
