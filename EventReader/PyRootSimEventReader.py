#!/usr/bin/env python
"""
------------------------------------------------------------------------

A script to read in simulated MEGALib events & store infomation as pandas dataframes.

Author: Henrike Fleischhack (fleischhack@cua.edu)
Date: August 13th, 2020

------------------------------------------------------------------------
"""

import os
import time
import sys
import fileinput
import numpy
from itertools import product, combinations, repeat
from collections import OrderedDict
import glob
import pandas
import math
import pickle

import gzip

from astropy.coordinates import cartesian_to_spherical, spherical_to_cartesian

from matplotlib import pyplot

import ROOT

ROOT.gSystem.Load("$(MEGALIB)/lib/libMEGAlib.so")

# Initialize MEGAlib
G = ROOT.MGlobal()
G.Initialize()

#We can use these to make printing of MEGAlib (and other root classes) easier.
def MegaPrint(self):
    try:
        return self.Data()
    except:
        return ":("
    
def MegaToStringPrint(self):
    try:
        return self.ToString().Data()
    except:
        return ":("

def MegaNamePrint(self):
    try:
        return self.GetName().Data()
    except:
        return ":("

setattr(ROOT.MString, '__str__', MegaPrint)
setattr(ROOT.MVector, '__str__', MegaToStringPrint)
setattr(ROOT.MSimHT, '__str__', MegaToStringPrint)
setattr(ROOT.MSimIA, '__str__', MegaToStringPrint)
setattr(ROOT.MComptonEvent, '__str__', MegaToStringPrint)
setattr(ROOT.MPairEvent, '__str__', MegaToStringPrint)
setattr(ROOT.MPhysicalEvent, '__str__', MegaToStringPrint)
setattr(ROOT.MDDetector, '__str__', MegaNamePrint)
setattr(ROOT.MDVoxel3D,  '__str__', MegaNamePrint)
setattr(ROOT.MDStrip2D,  '__str__', MegaNamePrint)


"""Project 3D hit position on one of the ACD panels to a common 2D coordinate frame for easier plotting later. Hardcoded panel positions for now."""
def ProjectSideACD( hit ):

    try:
        pos = hit.GetPosition()
    except:
        pos = hit
    
    
    #side panels
    
    if pos.GetX() == 53.25:
        return pos.GetY(), pos.GetZ(), 1
        
    if pos.GetY() == 53.25:
        return pos.GetX(), pos.GetZ(), 2
            
    if pos.GetX() == -53.25:
        return pos.GetY(), pos.GetZ(), 3

    if pos.GetY() == -53.25:
        return pos.GetX(), pos.GetZ(), 4
     
    #top panel
    if pos.GetZ() == 67.25:
        return pos.GetX(), pos.GetY(), 0
    
    
    print(u"😬", pos )
    assert False
 

"""Distance between a point and the line defined by a reference point and a direction, given as MEGALIB MVector """
def ImpactDistance( Point, ReferencePoint, Direction ):

    d = ReferencePoint - Point
    n = Direction.Unit()
    
    return (d - (d.Dot(n))*n ).Mag()



def MegaSimEventToDict(Event, Geometry):

    theDict = {}
       
    theACDHits = []
       
    theDict["ID"] = Event.GetID()
    
    pos1 = Event.GetIAAt(1).GetPosition()
    vol = Geometry.GetVolume( pos1 )
    
    theDict["TrueFirstIAVolume"] = vol.GetName().Data()
    
    if vol.IsSensitive():
        theDict["TrueFirstIADetector"] = vol.GetDetector().GetName().Data()
    else:
        theDict["TrueFirstIADetector"] = "passive"
    
    theDict["TrueFirstIAX"] = pos1.GetX()
    theDict["TrueFirstIAY"] = pos1.GetY()
    theDict["TrueFirstIAZ"] = pos1.GetZ()
       
    #First interaction type
    theDict["TrueType"] = Event.GetIAAt(1).GetProcess().Data()

    #The ED keyword (Total energy deposit in active material)
    ED = 0.0;
    
    meanHitPosition = ROOT.MVector(0,0,0)
    
    for iH in range(0, Event.GetNHTs()):
        theHit = Event.GetHTAt(iH)
        Ehit = theHit.GetEnergy()
        ED += Ehit
        meanHitPosition += Ehit*theHit.GetPosition()
        
        detectorType = theHit.GetDetector()
           
        theKey = "NHit_{}".format( detectorType )
        theDict[ theKey ] = theDict.get(theKey, 0) + 1
               
        theKey = "EHit_{}".format( detectorType )
        theDict[ theKey ] = theDict.get(theKey, 0) + Ehit/1e3
        
        hitDict = {}
        hitPos = theHit.GetPosition()
        hitDict["ID"] = theDict["ID"]
               
        hitDict["HitX"] = hitPos.GetX()
        hitDict["HitY"] = hitPos.GetY()
        hitDict["HitZ"] = hitPos.GetZ()
        hitDict["HitEnergy"] = Ehit/1e3
            
        #NOTE: The following lines are an example from when I was looking at the
        #hit distribution within the ACD panels, in this case simulated as voxel detectors.
        #if detectorType == 8: #voxel detector
        #    pX, pZ, type = ProjectSideACD( theHit )
        #    hitDict["HitProjectedX"] = pX
        #    hitDict["HitProjectedZ"] = pZ
        #    hitDict["HitPanel"] = type
                
        hitDict["HitDetector"] = Geometry.GetDetector( theHit.GetPosition() ).GetName().Data()
        theHits.append( hitDict )
           
    for iH in range(0, Event.GetNGRs()):
        ED += Event.GetGRAt(iH).GetEnergy()
    theDict["ED"] = ED
    if ED > 0:
        meanHitPosition = meanHitPosition * (1.0 / ED)
    theDict["TrueMeanHitX"] = meanHitPosition.GetX()
    theDict["TrueMeanHitY"] = meanHitPosition.GetY()
    theDict["TrueMeanHitZ"] = meanHitPosition.GetZ()
         
    #The EC keyword (Escapes)
    EC = 0.0;
    for iA in range(0, Event.GetNIAs()):
        if (Event.GetIAAt(iA).GetProcess().Data() == "ESCP"):
            EC += Event.GetIAAt(iA).GetMotherEnergy()
    theDict["EC"] = EC
         
    #Deposits in non-sensitive material
    theDict["NS"] = Event.GetEnergyDepositNotSensitiveMaterial()
       
    #incoming gamma-ray energy
    theDict["TrueEnergy"] = Event.GetICEnergy()/1e3
       
    #First interaction type
    theDict["TrueType"] = Event.GetIAAt(1).GetProcess().Data()
       
    initialPhotonDirection = -1*Event.GetICOrigin()
       
    theDict["TrueTheta"] = initialPhotonDirection.Theta()
    theDict["TruePhi"] =  initialPhotonDirection.Phi()

    if theDict["TrueType"] == "COMP":
        theDict["TrueThetaC"] = Event.GetICScatterAngle()       #True Compton angle
        theDict["TrueThetaE"] = Event.GetICElectronD().Angle(Event.GetICOrigin()) #True electron angle
        
    if theDict["TrueType"] == "PAIR":
        theDict["TruePairAngle"] = Event.GetIPPositronDir().Angle(Event.GetIPElectronDir() )
    
    
    return theDict, theHits


def readSimFile( FileName, GeometryName ):
    
    fileinfo = {}
    fileinfo["filename"] = FileName

    theEvents = []
    theHits = []

    #parse some info by hand
    
    if FileName[-2:] == "gz":
        import gzip
        op = gzip.open
    else:
        op = open
    
    with op(FileName, mode="rt") as fp:
        for line in fp:
            fields = line.replace(";", " ").split()
        
            if len(fields) < 1:
                continue
            
            keyword = fields[0]

            if keyword == "SimulationStartAreaFarField":
                fileinfo["Area"] = float(fields[1])
            elif keyword == "BeamType" and fields[1] == "FarFieldPointSource":
                fileinfo["SourceTheta"] = numpy.deg2rad(float(fields[2]))
                fileinfo["SourcePhi"] = numpy.deg2rad(float(fields[3]))
            elif keyword == "SpectralType" and fields[1] == "Mono":
                fileinfo["SourceEnergy"] = float(fields[2])/1000 #keV to MeV
            
            #first event! We'll use MEGAlib for the rest.
            if keyword == "SE":
                break

    # Use MEGAlib reader for the rest
    Geometry = ROOT.MDGeometryQuest()
    
    if Geometry.ScanSetupFile(ROOT.MString(GeometryName)) == True:
        print("Geometry " + GeometryName + " loaded!")
    else:
        print("Unable to load geometry " + GeometryName + " - Aborting!")
        quit()

    Reader = ROOT.MFileEventsSim(Geometry)
    if Reader.Open(ROOT.MString(FileName)) == False:
        print("Unable to open file " + FileName + ". Aborting!")
        quit()

    fileinfo["Area"] = Reader.GetSimulationStartAreaFarField()

    while True:
    
        Event = Reader.GetNextEvent()
        if not Event:
            break
        ROOT.SetOwnership(Event, True)
        
        eventDict, hitList = MegaSimEventToDict(Event, Geometry)
        theEvents.append( eventDict )
        theHits.extend( hitList )
               
    fileinfo["thrownEvents"] = Reader.GetSimulatedEvents()

    df = pandas.DataFrame(theEvents)
    
    #cleanup
    df.fillna(0, inplace=True)
    for b in ["NHit_1", "NHit_2", "NHit_4", "NHit_8"]:
        if b in df.columns:
            df = df.astype( {b:'int32' } )
        
    df.set_index("ID", inplace = True)


    df2 = pandas.DataFrame(theACDHits)
    df2.fillna(0, inplace = True)
    
    if "HitPanel" in df2.columns:
        df2 = df2.astype({'HitPanel': 'int32' })

    try:
        df2.set_index("ID", inplace = True)
    except:
        pass

    return df, df2, fileinfo
    
    
def MegaTraEventToDict( Event):

    theDict = {}
    
    theDict["ID"] = Event.GetId()
           
    theType = Event.GetTypeString().Data()
    
    theDict["RecoType"] = theType
        
    if theType == "Unidentifiable":
        theDict["RecoType"] = "bad"

        theDict["Subtype"] = " ({})".format(Event.GetBadString().Data() )
            

    if theType == "Compton":
        
        comptonPhotonEnergy = Event.Eg()/1e3
        comptonElectronEnergy = Event.Ee()/1e3
        
        totalEnergy = comptonPhotonEnergy+comptonElectronEnergy
        theDict["RecoEnergy"] = totalEnergy

        theDict["RecoThetaC_Cal"] = Event.Phi()
        theDict["RecoThetaE_Cal"] = Event.Epsilon()

        theDict["RecoDeltaTheta"] = Event.DeltaTheta()
                
        if Event.HasTrack():
            theDict["Subtype"] = " (tracked)"
        else:
            theDict["Subtype"] = " (untracked)"
                
        originalGammaDirection = -Event.GetOIDirection()
        scatteredGammaDirection = Event.C2() - Event.C1()
        theDict["RecoThetaC_Track"] = originalGammaDirection.Angle( scatteredGammaDirection )
    
    elif theType == "Pair":
    
        pairEnergy = Event.GetEnergy()/1e3
        theDict["RecoEnergy"] = pairEnergy
        theDict["RecoTheta"] = Event.GetOrigin().Theta()
        theDict["RecoPhi"] = Event.GetOrigin().Phi()
        theDict["RecoPairAngle"] = Event.GetOpeningAngle()
        theDict["RecoDelAngle"] = Event.GetOrigin().Angle(-Event.GetOIDirection())
        
    return theDict



def readTraFile( FileName, df = None ):
    
    theEvents = []

    Reader = ROOT.MFileEventsTra()
    if Reader.Open(ROOT.MString(FileName)) == False:
        print("Unable to open file " + FileName + ". Aborting!")
        quit()

    while True:
    
        Event = Reader.GetNextEvent()
        if not Event:
            break
        ROOT.SetOwnership(Event, True)
        
        eventDict = MegaTraEventToDict(Event)
        theEvents.append( eventDict )

    if len(theEvents) > 0:
        df = pandas.DataFrame(theEvents)
       
        #cleanup
        df["Subtype"] = df["Subtype"].fillna( "" )
        df.set_index("ID", inplace = True)

        return df
    
    else:
        return None
        

def readOneSetOfSims( geometry, simFileName, traFileName = None, **kwargs ):
   
    #print (geometry, simFileName, traFileName)
   
    if traFileName is None:
        traFileName = simFileName.replace(".sim", ".tra")
    
    simDfName = simFileName + ".csv"
    hitDfName = simFileName.replace(".sim", ".hits.csv" ).replace(".gz", "")
    simInfoName = simFileName + ".pkl"
    dfName = traFileName + ".csv"
    
    summaryFileName = traFileName.replace(".tra", ".eventsByType.csv").replace(".gz", "")
    summaryFileName_good = traFileName.replace(".tra", ".GoodEventsByType.csv").replace(".gz", "")

    if os.path.exists( simDfName ) and os.path.exists( hitDfName):
        simData = pandas.read_csv( simDfName, index_col="ID" )
        try:
            hitData = pandas.read_csv( hitDfName, index_col="ID" )
        except:
            hitData = pandas.DataFrame()
        with open(simInfoName,"rb") as f:
            simInfo = pickle.load(f)
    
    else:
        simData, hitData, simInfo = readSimFile(simFileName, geometry, impactDistributionFile)
        simData.to_csv( simDfName )
        hitData.to_csv( hitDfName )
        with open(simInfoName,"wb") as f:
            pickle.dump(simInfo,f)


    if os.path.exists(dfName):
        df = pandas.read_csv(dfName, index_col = "ID")

    else:
        recData = readTraFile(traFileName)
        
        #0 events in file
        if recData is None:
            return simData, hitData, simInfo, None, None, None
        
        df = simData.join(recData)

        df["RecoType"] = df["RecoType"].fillna( "no trigger" )
        df["Subtype"] = df["Subtype"].fillna( "" )
                
        df["LongType"] = df.RecoType + df.Subtype

        #Select compton events with reco energy within 10% of true energy, pair events within 50%
        good_compton = ( (df.RecoType == "Compton") & ( numpy.abs( 1.0 - df.RecoEnergy/df.TrueEnergy ) < 0.1  ) )
        good_pair = ( (df.RecoType == "Pair") & ( numpy.abs( 1.0 - df.RecoEnergy/df.TrueEnergy ) < 0.5  ) )

        df["RecoEnergyCut"] = good_compton | good_pair

        df = df[~ (df.TrueType.isin( ["RAYL", "PHOT"]))]
        df.to_csv(dfName)
 
 
    if os.path.exists(summaryFileName):
        summary = pandas.read_csv(summaryFileName)

    else:
        summary = df.groupby(["TrueType", "LongType"]).size().unstack(level=-1).fillna(0).reset_index()

        summary["ThrownEvents"] = simInfo["thrownEvents"]
        summary["Area"] = simInfo["Area"]
        summary["SourceEnergy"] = simInfo["SourceEnergy" ]
        summary["SourceTheta"] = simInfo["SourceTheta"]
        
        for key, value in kwargs.items():
            summary[key] = value
    
        summary.to_csv(summaryFileName)

 
    if os.path.exists(summaryFileName_good):
        summary_good = pandas.read_csv(summaryFileName_good)

    else:
        summary_good = df[df.RecoEnergyCut].groupby(["TrueType", "LongType"]).size().unstack(level=-1).fillna(0).reset_index()

        summary_good["ThrownEvents"] = simInfo["thrownEvents"]
        summary_good["Area"] = simInfo["Area"]
        summary_good["SourceEnergy"] = simInfo["SourceEnergy" ]
        summary_good["SourceTheta"] = simInfo["SourceTheta"]
        
        for key, value in kwargs.items():
            summary_good[key] = value
    
        summary_good.to_csv(summaryFileName_good)

    s2s = {}

    try:
        for n, d in df.groupby("TrueType"):
    
            summary2FileGlob = traFileName.replace(".tra", ".eventsByTypeAndFirstIA.{}.N*.csv".format(n))
        
            if len( glob.glob(summary2FileGlob) ) == 1 and os.path.exists(glob.glob(summary2FileGlob)[0] ):
                summary2 = pandas.read_csv(glob.glob(summary2FileGlob)[0], index_col = "LongType" )

            else:
                summary2 = d.groupby(["TrueFirstIADetector", "LongType"]).size().unstack(level=0).fillna(0)
                nEvents = summary2.values.sum()
        
                summary2["Total"] = summary2.sum(axis=1)
                summary2.loc["Total"] = summary2.sum(axis=0)
        
                summary2 = summary2/nEvents
                print (n, nEvents)
                print(summary2)
                s2s[(n, nEvents)] = summary2
            
                summary2FileName = traFileName.replace(".tra", ".eventsByTypeAndFirstIA.{}.N{}.csv".format(n, int(nEvents))).replace(".gz", "")
                summary2.to_csv( summary2FileName )
    except:
        pass
        

    return simData, hitData, simInfo, df, summary, s2s, summary_good
    
    

def readAllTheFiles(basedir, theSims, simDir, recoDir, geometry, **kwargs ):

    retDF = pandas.DataFrame()
    
    retDF2 = pandas.DataFrame()

    theSimFiles = glob.glob( os.path.join(basedir, simDir, theSims ) )
        
    for simFileName in theSimFiles:

        traFileName = os.path.join( basedir, recoDir, os.path.basename(simFileName).replace(".sim", ".tra") )
        if not os.path.exists(traFileName):
            print("No", traFileName, u"😭")
            continue
        simData, hitData, simInfo, df, summary, s2s, summary_good = readOneSetOfSims( geometry, simFileName, traFileName, **kwargs)

        if summary is not None:
            retDF = retDF.append(summary, ignore_index = True, sort = True )
    
        if summary_good is not None:
            retDF2 = retDF2.append(summary_good, ignore_index = True, sort = True )
    
    return retDF.fillna(0), retDF2.fillna(0)
    
