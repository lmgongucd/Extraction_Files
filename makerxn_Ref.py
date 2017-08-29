#!/usr/bin/env python

import warnings
# Suppress warnings from Molecule class.
warnings.simplefilter("ignore")
import os, sys, re
#import networkx as nx
import numpy as np
#import copy
#import ast
#from collections import OrderedDict, defaultdict, Counter
import argparse
from nanoreactor import Nanoreactor
from nanoreactor import Molecule
from nanoreactor.nifty import commadash 
#import itertools
#import time

# Arguments:
# 1) --in .xyz file for trajectory.xyz
# 2) --qs .txt file for charge-spin.txt
# 3) --frm starting frame (indexed from zero)
# 4) --stride length of stride, default as 1 (no lengthened stride)
# 5) --end ending frame (indexed from zero)
# 6) --f(list of index numbers of interest throughout the reaction)

parser=argparse.ArgumentParser()
parser.add_argument('xyzin', help='XYZ input file')
parser.add_argument('--qs',type=str)
parser.add_argument('--frm',type=int)
parser.add_argument('--end',type=int)
parser.add_argument('--stride', type=int, default=1)
parser.add_argument('--f',type=str)
parser.add_argument('--N', help='if atoms selected are not yet neutralized', type=str)
#self, xyzin=None, qsin=None, ftype=None, stride=1, frames=0, xyzout='out.xyz'
args=parser.parse_args()
aindex=args.f.split(',')
atomindex=[]
for i in range(len(aindex)):
  atomindex.append(int(aindex[i]))
atomindex=np.array(atomindex)
stride = int(args.stride)
frameInit = int(args.frm)
frameFinal = int(args.end)
frames=np.arange(frameInit,frameFinal+1,stride)

M = Molecule(args.xyzin, ftype='xyz', build_topology=False)
P = Molecule(args.qs, ftype='xyz', build_topology=False)

Mcut=M[frameInit:frameFinal:stride]
Pcut=P[frameInit:frameFinal:stride]

#here I think to extract neutralizing charges
# want to neutralize the system following the Nanoreactor class in nanoreactor.py
# --input of coord, charge, atom index, frames
# ---- neutralize and chop frames with newly selected atoms
# output of frames charges, ready for Refine.py!!

PcutK1=Pcut.atom_select(8, build_topology=False)
PcutK2=Pcut.atom_select(13, build_topology=False)
PcutK3=Pcut.atom_select(34, build_topology=False)
PcutK4=Pcut.atom_select(42, build_topology=False)
Mcut=Mcut.atom_select(atomindex, build_topology=False)
Pcut=Pcut.atom_select(atomindex, build_topology=False)

#checking what the selection is
Mmid='Rxn_mid.xyz'
Pmid='Rxn_mid.pop'
Mcut.write(Mmid)
Pcut.write(Pmid,ftype='xyz')

#M2 = Nanoreactor(xyzin=Mmid, qsin=Pmid, ftype='xyz')
# is there a way to select which frames to not include?
framesel=np.arange(0,frameFinal-frameInit)
#print 'starting neutralizing'
#Reactant_Atoms = atomindex
#Neutral_Atoms, Neutral_Isos=M2.GetNeutralizing(Reactant_Atoms,framesel)
#RxnNum=0
#RxnList=[]
#Slices = []        # A list of lists of Molecule objects to be saved to disk
#Firsts = []        # A list of lists of first frames
#Lasts  = []        # A list of lists of last frames
#SliceIndices = []  # A list of atom indices corresponding to each saved trajectory
#SliceFrames = []   # A list of frames corresponding to each saved trajectory
#BufferTime = 0
##PadTime = self.PadTime
#
##currently no Neutral_Atoms are being selected by the Get.Neutralizing function
##currently also Neutral_Atoms is a different type, and so cannot "cast" with Reactant_Atoms
#Reactant_Atoms += Neutral_Atoms
##iso0 += Neutral_Isos
##iso1 += Neutral_Isos
#Reactant_Atoms = sorted(list(set(Reactant_Atoms)))
##Slice = self.atom_select(Reactant_Atoms)[FrameSel]
##rx0 = '+'.join([("%i" % iso0.count(i) if iso0.count(i) > 1 else "") + self.Isomers[i].ef() for i in sorted(set(iso0))])
##rx1 = '+'.join([("%i" % iso1.count(i) if iso1.count(i) > 1 else "") + self.Isomers[i].ef() for i in sorted(set(iso1))])
##evector="%s -> %s" % (rx0 if EndPoint else rx1, rx1 if EndPoint else rx0)
##tx0 = str([commadash(sorted([Reactant_Atoms.index(i) for i in self.TimeSeries[g]['graph'].L()])) for g in gid0]).replace(" ","")
##tx1 = str([commadash(sorted([Reactant_Atoms.index(i) for i in self.TimeSeries[g]['graph'].L()])) for g in gid1]).replace(" ","")
##tvector="%s -> %s" % (tx0 if EndPoint else tx1, tx1 if EndPoint else tx0)
##Slice.comms = ["Reaction: formula %s atoms %s frame %s charge %+.3f sz %+.3f sz^2 %.3f" % (evector, tvector, str(frame), sum(self.Charges[frame][Reactant_Atoms]),
##
##        sum(self.Spins[frame][Reactant_Atoms]), sum([j**2 for j in self.Spins[frame][Reactant_Atoms]])) for frame in FrameSel]
#print "finished neutralizing"

Pdata = np.array(Pcut.xyzs)
Pdata1 = np.array(PcutK1.xyzs)
Pdata2 = np.array(PcutK2.xyzs)
Pdata3 = np.array(PcutK3.xyzs)
Pdata4 = np.array(PcutK4.xyzs)
print "K1 mean,stdev axis=0"
print Pdata1[:,:,0].mean(axis=0), Pdata1[:,:,0].std(axis=0)
print "K2 mean,stdev axis=0"
print Pdata2[:,:,0].mean(axis=0), Pdata2[:,:,0].std(axis=0)
print "K3 mean,stdev axis=0"
print Pdata3[:,:,0].mean(axis=0), Pdata3[:,:,0].std(axis=0)
print "K4 mean,stdev axis=0"
print Pdata4[:,:,0].mean(axis=0), Pdata4[:,:,0].std(axis=0)

charges = Pdata[:,:,0].sum(axis=1)
spins = Pdata[:,:,1].sum(axis=1)
spins2 = spins**2
#M=M.atom_select(atomindex, build_topology=False)
#P=P.atom_select(atomindex, build_topology=False)

# Slice.comms = ["Product: formula %s atoms %s frame %s charge %+.3f sz %+.3f sz^2 %.3f" % (I.ef(), commadash(Indices[i]), str(frame), sum(self.Charges[frame][Indices[i]]), sum(self.Spins[frame][Indices[i]]), sum([j**2 for j in self.Spins[frame][Indices[i]]])) for frame in FrameSel]
Mcut.comms = ['Product: atoms %s frame %d charge %+.3f sz %+.3f sz^2 %.3f' %(commadash(atomindex), frames[i], charges[i], spins[i], spins2[i]) for i in range(Mcut.ns)]
Pcut.comms = Mcut.comms
Mout = 'Reaction.xyz'
Pout = 'Reaction.pop'
Mcut.write(Mout)
Pcut.write(Pout, ftype='xyz')
print "Atom extraction complete! see Reaction.xyz and Reaction.pop"
