#A Python3 script to parse the run commands via which the relevant xsec and event files
#can be produced using the code neutrino_tridents.cxx.
#To get started, you need an executable for the C++ code:
#
#  g++ --std=gnu++0x -o tridents.exe neutrino_tridents.cpp
#
#Then, have a look at the bottom of this script, containing a list of printCmds calls.
#Since a lot of runs will be required, consider commenting out a few of these and only
#running the code for some at a time, especially if all features are not required.
#Once confident with the selection, call
#
#  python3 genRunCmds.py
#
#This generates files cmds_0.txt - cmds_ncpu.txt. 
#Then, open ncpu amount of terminals, copy-paste the contents of each file into one, 
#and [use your imagination] while running. After everything succesful runs, move the
#resulting xsec files with names of the form
#
#  process_material_IncomingNuID_FSleptons.txt
#
#under
#
#  xsec/cm2_nucleon/
#
#and any generated event files with names of the form
#
#  Material_model_EnergyBinIndex.txt
#
#under 
#
#  generated_events/
#
#The files are not moved automatically in order to avoid risk of destroying and
#beginning to overwrite any existing files.

ncpu=8 #Specify max number of terminals to run in

def printCmds(mode,material,model):
    #Incoming nu_mu resulting in mumu final states. Relevant for SM and 4F studies
    procs=[['8','numu'], ['11','nubarmu']]
    
    #Cases only relevant for SM studies
    if model=='SM':
        
        #Incoming nu_e, nu_tau resulting in mumu final states
        procs+=[['2','nuel'], ['5','nubarel'], ['23','nutau'], ['28','nubartau']]
        
        #All other non-mumu FS, relevant mainly for producing the table in the paper
        if mode=='xsec':
            procs=[['1',  '12_ee'],    ['7',  '14_ee'],    ['22', '16_ee'],
                   ['4', '-12_ee'],    ['10','-14_ee'],    ['27','-16_ee'],
                   ['3',  '12_emu'],   ['9',  '14_emu'],
                   ['6', '-12_emu'],   ['12','-14_emu'],
                   ['14', '12_etau'],  ['24', '16_etau'],
                   ['16','-12_etau'],  ['29','-16_etau'],
                   ['18', '14_mutau'], ['25', '16_mutau'],
                   ['20','-14_mutau'], ['30','-16_mutau'],
                   ['13', '12_tautau'],['17', '14_tautau'],['21', '16_tautau'],
                   ['15','-12_tautau'],['19','-14_tautau'],['26','-16_tautau']]
    
    Enu = ['1.00000000e+00',\
           '1.25892541e+00',\
           '1.58489319e+00',\
           '1.99526231e+00',\
           '2.51188643e+00',\
           '3.16227766e+00',\
           '3.98107171e+00',\
           '5.01187234e+00',\
           '6.30957344e+00',\
           '7.94328235e+00',\
           '1.00000000e+01',\
           '1.25892541e+01',\
           '1.58489319e+01',\
           '1.99526231e+01',\
           '2.51188643e+01',\
           '3.16227766e+01',\
           '3.98107171e+01',\
           '5.01187234e+01',\
           '6.30957344e+01',\
           '7.94328235e+01',\
           '1.00000000e+02',\
           '1.25892541e+02',\
           '1.58489319e+02',\
           '1.99526231e+02',\
           '2.51188643e+02',\
           '3.16227766e+02',\
           '3.98107171e+02',\
           '5.01187234e+02',\
           '6.30957344e+02',\
           '7.94328235e+02',\
           '1.00000000e+03',\
           '1.25892541e+03',\
           '1.58489319e+03',\
           '1.99526231e+03',\
           '2.51188643e+03',\
           '3.16227766e+03',\
           '3.98107171e+03',\
           '5.01187234e+03',\
           '6.30957344e+03',\
           '7.94328235e+03',\
           '1.00000000e+04']
    
    #Interpolation grid for 4F studies
    if model=='4F':
        BSMstrs=['dGV','dGA']
        BSMs = []
        dGAvals = ['-2.00', '-1.75', '-1.50', '-1.25', '-1.00',\
                   '-0.75', '-0.50', '-0.25',  '0.00',  '0.25',\
                    '0.50',  '0.75',  '1.00',  '1.25',  '1.50',\
                    '1.75',  '2.00',  '2.25',  '2.50',  '2.75',\
                    '3.00',  '3.25',  '3.50',  '3.75',  '4.00']
        dGVvals = ['-5.00', '-4.75', '-4.50', '-4.25', '-4.00',\
                   '-3.75', '-3.50', '-3.25', '-3.00', '-2.75',\
                   '-2.50', '-2.25', '-2.00', '-1.75', '-1.50',\
                   '-1.25', '-1.00', '-0.75', '-0.50', '-0.25',\
                    '0.00',  '0.25',  '0.50',  '0.75',  '1.00']
        BSMs = [[dGV, dGA] for dGV in dGVvals for dGA in dGAvals]
    elif model=='SM': BSMs=[[]]  #No need for g', m_Z' in SM computations
    
    retcmds=[]
    for proc in procs:
        for BSM in BSMs:
            proccmds=[]
            for i in range(len(Enu)):
                cmd = './tridents.exe ' + proc[0] + ' ' + material + ' 1 ' + Enu[i]
                cmd += ' ' + model
                if model!='SM': cmd += ' ' + BSM[0] + ' ' + BSM[1]
                if   mode=='xsec':
                    cmd += ' CrossSection'
                elif mode=='evts': 
                    cmd += ' GenerateEvents 100000 '  #Modify number of events here
                    cmd += material
                    if model=='SM':
                        cmd += '_SM'                
                    else:
                        cmd += '_' + BSMstrs[0] + '-' + BSM[0].replace('.','d')
                        cmd += '_' + BSMstrs[1] + '-' + BSM[1].replace('.','d')
                    cmd += '_' + proc[1]        #numu, nubarmu, ...
                    Enustr = str(i+1)
                    if len(Enustr)<2: Enustr = '0'+Enustr
                    cmd += '_' + Enustr + '.txt'
                proccmds.append(cmd)
            retcmds.append(proccmds)
    return retcmds

#############################################################
# Print all commands into NCPU files (quasi-multithreading) #
#############################################################

import numpy as np

def divideCmds(nproc,cmds):
    nproc = min(nproc,len(cmds))
    step = int(np.ceil(len(cmds)/nproc))
    cmds_divided = [cmds[n*step:(n+1)*step] for n in range(nproc)]
    for i,cmds_cpu in enumerate(cmds_divided):
        f = open('cmds_'+str(i)+'.txt.','w')
        for cmd in np.array(cmds_cpu).flatten(): f.write(cmd+'\n')
        f.close()

#######################
# MAIN FUNCTION CALLS #
#######################

cmds=[]

#For SM results, both cross sections and event listings are needed
#N.B. in case you consider also detectors other than FASERv(2) in detail,
#     make sure to add the 'evts' calls for the corresponding materials here
cmds += printCmds(mode='xsec',material='W', model='SM')
cmds += printCmds(mode='evts',material='W', model='SM')

#Additional SM xsec relevant for producing the event count table in the paper
cmds += printCmds(mode='xsec',material='Ar',model='SM')
cmds += printCmds(mode='xsec',material='Fe',model='SM')
cmds += printCmds(mode='xsec',material='Pb',model='SM')

#Relevant cross sections for four-fermi (4F) studies
cmds += printCmds(mode='xsec',material='W', model='4F')

#Produce files
divideCmds(nproc=ncpu,cmds=cmds)
