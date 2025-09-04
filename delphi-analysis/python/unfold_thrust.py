import ROOT
import numpy as np
import sys
import os
import argparse
from array import array

from ROOT import RooUnfoldResponse
from ROOT import RooUnfold
from ROOT import RooUnfoldBayes
from common_functions import *

def Proj2D_Y(h,xmin,xmax,hname="XXX"):
    imin=h.GetXaxis().FindBin(xmin)
    imax=h.GetXaxis().FindBin(xmax)-1
    
    proj_y=h.ProjectionY(hname, imin, imax)
    ROOT.SetOwnership(proj_y,True)
    return proj_y

def Proj2D_X(h,ymin,ymax,hname="XXX",Debug=False):
    imin=h.GetYaxis().FindBin(ymin)
    imax=h.GetYaxis().FindBin(ymax)-1

    proj_x=h.ProjectionX(hname, imin, imax)
    ROOT.SetOwnership(proj_x,True)
    return proj_x

def convert_to_roounfold_response(reco_hist, gen_hist, response_hist, name="response"):
    """
    Much simpler approach using direct constructor!
    """
    
    print(f"Creating RooUnfoldResponse '{name}' directly from histograms")
    
    # Just use the direct constructor - that's it!
    response = ROOT.RooUnfoldResponse(
        reco_hist,      # measured TH1
        gen_hist,       # truth TH1
        response_hist,   # response TH2
        name,           # name
        name,           # title
        False           # do_overflow
    )
    
    return response

# Your binning arrays
eijbins = [0.0, 0.0001, 0.0002, 0.0005, 0.00075, 0.001, 0.00125, 0.0015, 0.00175, 0.002, 0.00225, 0.0025, 0.00275, 0.003, 0.0035, 0.004, 0.005, 0.007, 0.01, 0.02, 0.03, 0.04, 0.05, 0.07, 0.10, 0.15, 0.20, 0.3, 1]

rbins = calcBinEdge(0.002, np.pi/2, 100)
zbins = calcBinEdge(0.000001, 0.5, 100)

tbins = np.linspace(0.5, 1.1, 61)
logtbins = np.linspace(-10, 0, 101)

eijbins = np.array(eijbins)
rbins = np.array(rbins)
zbins = np.array(zbins)
tbins = np.array(tbins)
logtbins = np.array(logtbins)

if __name__ == '__main__':

    #filenamein = "response_kk2f4146_qqpy_91.25_v14.root"
    #filenamein = "response_kk2f4146_qqardcy_91.25_v1.root"
    #filenamein = "response_apacic105_91.25_v1.root"
    filenamein = "response_pythia8_91.25_v2.root"
    #filenamein = "response_pythia8_dire_91.25_v1.root"

    datafile = 'h_94c_v5.root'
    dataname = 'ThrustMissPNCLog2'
    #dataname = 'ThrustMissPNCDelphi'

    #datafile = 'h_kk2f4146_qqpy_91.25_v5.root'
    #dataname = 'ThrustMissPNCDelphi'
    #dataname = 'ThrustMissPNCLog2'

    #filenameout = "unfolded_thrustDelphi_kk2f4146_qqpy_91.25_v2.root"
    #filenameout = "unfolded_data_thrustDelphi_kk2f4146_qqpy_91.25_v2.root"
    #filenameout = "unfolded_data_thrustDelphi_kk2f4146_qqardcy_91.25_v2.root"
    #filenameout = "unfolded_data_thrustDelphi_apacic_91.25_v2.root"
    #filenameout = "unfolded_data_thrustDelphi_pythia8_91.25_v2.root"
    #filenameout = "unfolded_data_thrustDelphi_pythia8_dire_91.25_v2.root"
    
    #filenameout = "unfolded_thrust_kk2f4146_qqpy_91.25_v2.root"
    #filenameout = "unfolded_data_thrust_kk2f4146_qqpy_91.25_v2.root"
    #filenameout = "unfolded_data_thrust_kk2f4146_qqardcy_91.25_v2.root"
    #filenameout = "unfolded_data_thrust_apacic_91.25_v2.root"
    filenameout = "unfolded_data_thrust_pythia8_91.25_v2.root"
    #filenameout = "unfolded_data_thrust_pythia8_dire_91.25_v2.root"

    #filenamein = "response_kk2f4146_qqpy_91.25_v14.root"
    #datafile = 'response_kk2f4146_qqpy_91.25_v14.root'
    #filenameout = "unfolded_thrust_kk2f4146_qqpy_91.25_v2.root"
    #dataname='reco_thrust_log2'

    fin = ROOT.TFile.Open(filenamein,'r')

    # Get original histograms
    counter = fin.Get("counter").Clone("N")
    _response = fin.Get("response_thrust_log2")
    reco = fin.Get('reco_thrust_log2')
    gen= fin.Get("gen_thrust_log2")

    #_response = fin.Get("response_thrustDelphi")
    #reco = fin.Get('reco_thrustDelphi')
    #gen= fin.Get("gen_thrustDelphi")

    fdata = ROOT.TFile.Open(datafile,'r')
    data = fdata.Get(dataname)

    normalization = fin.Get("counter").GetBinContent(2)
    try:
        n = fdata.Get('counter').GetBinContent(2)
    except:
        n = fdata.Get('N').GetBinContent(2)
    
    data.Scale(float(normalization)/n)

    
    # Save results
    fout = ROOT.TFile(filenameout, 'recreate')
    fout.cd()
    
    _response.Write("response")
    response=convert_to_roounfold_response(reco, gen, _response)

    # Check response matrix properties
    RESPONSE = response.Mresponse()
    singular = ROOT.TDecompSVD(RESPONSE)
    print("Response matrix singular values:")
    singular.GetSig().Print()

    # Perform unfolding
    unfold = ROOT.RooUnfoldBayes(response, data, 4)

    hUnf = unfold.Hunfold().Clone("unfolded")
    hErr = unfold.Eunfold()

    counter.Write("N")
    
    # Save restricted histograms used for unfolding
    reco.Write("reco")
    gen.Write("gen") 
    data.Write("data")
    
    # Save unfolding results
    hUnf.Write("unfolded")
    hErr.Write("unfolding_errors")

    fout.Close()

    print("Unfolding completed successfully!")
#    print(f"Original histogram bins: {_reco_orig.GetNbinsX()}")
#    print(f"Restricted histogram bins: {_reco.GetNbinsX()}")
#    print(f"Response matrix size: {_response.GetNbinsX()}x{_response.GetNbinsY()}")
