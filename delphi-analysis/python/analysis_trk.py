import ROOT
import numpy as np, dataclasses as dc
import sys
import os
import argparse
import math
import itertools
from array import array
from common_functions import *
from binning_and_selections import *

# Global histogram containers
histograms_1d = {}
histograms_2d = {}
histograms_nd = {}

def bookHistograms():
    global histograms_1d, histograms_2d, histograms_nd
    
    # Define 1D histogram configurations using numpy arrays for bin edges
    # Format: "name": bin_edges_array
    hist_1d_configs = {
        # Event counter
        "N": np.array(np.linspace(0, 10, 11)),
        
        # Charged track histograms  
        "trk_eff": np.array(np.linspace(0, 100, 101)),
        "trk_eff_theta": np.array(np.linspace(0, 180, 181)),
        "trk_eff_phi": np.array(np.linspace(0, 180, 181)),  # phi also 0-180
        "trk_fake": np.array(np.linspace(0, 100, 101)),
        "trk_fake_theta": np.array(np.linspace(0, 180, 181)),
        "trk_fake_phi": np.array(np.linspace(0, 180, 181)),  # phi also 0-180
        
        # Neutral track histograms
        "ntr_eff": np.array(np.linspace(0, 100, 101)),
        "ntr_eff_theta": np.array(np.linspace(0, 180, 181)),
        "ntr_eff_phi": np.array(np.linspace(0, 180, 181)),  # phi also 0-180
        "ntr_fake": np.array(np.linspace(0, 100, 101)),
        "ntr_fake_theta": np.array(np.linspace(0, 180, 181)),
        "ntr_fake_phi": np.array(np.linspace(0, 180, 181)),  # phi also 0-180
        
        # Photon histograms
        "pho_eff": np.array(np.linspace(0, 100, 101)),
        "pho_eff_theta": np.array(np.linspace(0, 180, 181)),
        "pho_eff_phi": np.array(np.linspace(0, 180, 181)),  # phi also 0-180
        "pho_fake": np.array(np.linspace(0, 100, 101)),
        "pho_fake_theta": np.array(np.linspace(0, 180, 181)),
        "pho_fake_phi": np.array(np.linspace(0, 180, 181)),  # phi also 0-180
    }
    
    # Define detector regions for theta binning
    theta_regions = np.array([0, 40, 140, 180])  # Forward endcap, barrel, backward endcap
    
    # Fine binning for angular response (theta_reco - theta_gen) from -pi to pi
    angular_response_bins = np.array(np.linspace(-np.pi, np.pi, 20001))  # 20000 bins with fine resolution
    
    # Define 2D histogram configurations using numpy arrays
    # resp_*_* histograms: x=response_value, y=theta_detector_regions  
    hist_2d_configs = {
        # Response histograms: X=angular_response (theta/phi_reco - theta/phi_gen), Y=detector_regions
        "resp_trk_theta": [angular_response_bins, theta_regions],
        "resp_trk_phi": [angular_response_bins, theta_regions],  # phi response vs detector regions
        "resp_ntr_theta": [angular_response_bins, theta_regions],
        "resp_ntr_phi": [angular_response_bins, theta_regions],
        "resp_pho_theta": [angular_response_bins, theta_regions],
        "resp_pho_phi": [angular_response_bins, theta_regions],
        
        # Matching histograms (keep original structure)
        "match_trk": [np.array(np.linspace(0, 100, 101)), np.array(np.linspace(0, 100, 101))],
        "match_trk_theta": [np.array(np.linspace(0, 180, 181)), np.array(np.linspace(0, 180, 181))],
        "match_trk_phi": [np.array(np.linspace(0, 180, 181)), np.array(np.linspace(0, 180, 181))],
        "match_ntr": [np.array(np.linspace(0, 100, 101)), np.array(np.linspace(0, 100, 101))],
        "match_ntr_theta": [np.array(np.linspace(0, 180, 181)), np.array(np.linspace(0, 180, 181))],
        "match_ntr_phi": [np.array(np.linspace(0, 180, 181)), np.array(np.linspace(0, 180, 181))],
        "match_pho": [np.array(np.linspace(0, 100, 101)), np.array(np.linspace(0, 100, 101))],
        "match_pho_theta": [np.array(np.linspace(0, 180, 181)), np.array(np.linspace(0, 180, 181))],
        "match_pho_phi": [np.array(np.linspace(0, 180, 181)), np.array(np.linspace(0, 180, 181))],
    }
    
    # Define pT/E bins for different particle types
    charged_pt_bins = np.array([0.4, 0.6, 1., 2., 4., 7., 10., 15., 20., 25., 30., 35., 45., 200])
    neutral_pt_bins = np.array([0.5, 1., 2., 4., 8., 16., 32., 200.])
    
    # Define response resolution binning (pt_reco/pt_gen) with fine binning
    response_resolution_bins = np.array(np.linspace(0, 2, 2001))  # 2000 bins for resolution studies
    
    # Define ND sparse histogram configurations  
    # 4D structure: axis(0)=response_resolution, axis(1)=pt_bins, axis(2)=detector_regions, axis(3)=dummy
    dummy_axis = np.array([0, 1])  # Single bin for unused axis(3)
    
    hist_nd_configs = {
        "resp_trk": [
            response_resolution_bins,  # axis(0): response value (pt_reco/pt_gen) 0-2
            charged_pt_bins,          # axis(1): pT bins for charged particles
            theta_regions,            # axis(2): detector regions [0,40], [40,140], [140,180]  
            dummy_axis                # axis(3): not useful, single bin
        ],
        "resp_ntr": [
            response_resolution_bins,  # axis(0): response value (pt_reco/pt_gen) 0-2
            neutral_pt_bins,          # axis(1): pT bins for neutral particles  
            theta_regions,            # axis(2): detector regions
            dummy_axis                # axis(3): not useful, single bin
        ],
        "resp_pho": [
            response_resolution_bins,  # axis(0): response value (E_reco/E_gen) 0-2
            neutral_pt_bins,          # axis(1): Energy bins for photons (same as neutral)
            theta_regions,            # axis(2): detector regions
            dummy_axis                # axis(3): not useful, single bin
        ],
    }
    
    # Book 1D histograms from numpy arrays - CORRECTED VERSION
    for name, bin_edges in hist_1d_configs.items():
        nbins = len(bin_edges) - 1
        # Convert numpy array to ROOT array format
        bin_array = ROOT.std.vector('double')(bin_edges)
        
        hist = ROOT.TH1D(name, name, nbins, bin_array.data())
        hist.SetDirectory(0)
        histograms_1d[name] = hist
    
    # Book 2D histograms from numpy arrays - CORRECTED VERSION
    for name, [x_edges, y_edges] in hist_2d_configs.items():
        nbinsx = len(x_edges) - 1
        nbinsy = len(y_edges) - 1

        x_edges = np.asarray(x_edges, dtype=np.double)
        y_edges = np.asarray(y_edges, dtype=np.double)  
        
        hist = ROOT.TH2D(name, name, nbinsx, x_edges, nbinsy, y_edges)
        hist.SetDirectory(0)
        histograms_2d[name] = hist
    
    # Book ND histograms from numpy arrays - WORKING VERSION
    for name, bin_edges_list in hist_nd_configs.items():
        ndim = len(bin_edges_list)
        
        # Create individual TAxis objects with variable binning
        axes = ROOT.std.vector('TAxis')()
        
        for dim, edges in enumerate(bin_edges_list):
            nbins = len(edges) - 1
            
            # Create TAxis with variable binning
            edges = np.asarray(edges, dtype=np.double)
            axis = ROOT.TAxis(nbins, edges)
            axis.SetName(f"axis_{dim}")
            axes.push_back(axis)
        
        # Create THnSparseD using the axes vector constructor
        hist = ROOT.THnSparseD(name, name, axes)
        histograms_nd[name] = hist

if __name__ == "__main__":

    filename = '/eos/user/z/zhangj/DELPHI/simulation/v94c/91.25/kk2f4146_qqpy/nanoaod_kk2f4146_qqpy_91.25_40001.sdst.root'
    filenameout = "trk_test.root"

    parser = argparse.ArgumentParser()
    parser.add_argument("infiles", nargs='?', default=filename, help="name of input files")
    parser.add_argument("outfile", nargs='?', default=filenameout, help="name of output files")
    args = parser.parse_args()

    fin = ROOT.TFile.Open(args.infiles, 'r')
    treco = fin.Get("t")
    tgen = fin.Get("tgenBefore")

    nevt = treco.GetEntries()

    bookHistograms()

    for ievt in range(nevt):

        treco.GetEntry(ievt)
        tgen .GetEntry(ievt)

        if ievt % 1000 == 0:
            print(f"Processing event {ievt}/{nevt}")

        E_reco   = treco.Energy
        E_gen   = tgen.Energy
        if abs(E_gen - 91.25) > 1: continue  ## this should always be False in fixed energy MC

        get      = lambda *names: (np.asarray(getattr(treco, n)) for n in names)
        get_gen  = lambda *names: (np.asarray(getattr(tgen , n)) for n in names)

        px,py,pz,m,q,th,pt,eta,phi,d0,z0,pwflag,hp,idx,cspidx = get(
            'px','py','pz','mass','charge','theta','pt','eta','phi','d0','z0','pwflag','highPurity','index','correspondenceIndex')
        
        px_gen,py_gen,pz_gen,m_gen,q_gen,th_gen,pt_gen,eta_gen,phi_gen,pwflag_gen,hp_gen,idx_gen,cspidx_gen = get_gen(
            'px','py','pz','mass','charge','theta','pt','eta','phi','pwflag','highPurity','index','correspondenceIndex')

        # Gen level selection neccessary for DELPHI for V0 particle complication
        # For ALEPH, this is always true anyway
        e = np.sqrt(px_gen**2 + py_gen**2 + pz_gen**2 + m_gen**2)
        if (np.sum(e) - E_gen) > 0.1: continue

        histograms_1d['N'].Fill(0.5)

        all_results = apply_track_selection_delphi(px=px, py=py, pz=pz, m=m, q=q, th=th,
                                                    pt=pt, eta=eta, phi=phi, d0=d0, z0=z0, pwflag=pwflag, hp=hp,
                                                    px_gen=px_gen, py_gen=py_gen, pz_gen=pz_gen, m_gen=m_gen, q_gen=q_gen,
                                                    th_gen=th_gen, pt_gen=pt_gen, eta_gen=eta_gen, phi_gen=phi_gen,
                                                    pwflag_gen=pwflag_gen, hp_gen=hp_gen)

        sel_c = all_results['sel_c']          # Reco charged particles
        sel = all_results['sel']              # All reco particles
        sel_c_gen = all_results['sel_c_gen']  # Gen charged particles
        sel_gen = all_results['sel_gen']      # All gen particles

        px_c, py_c, pz_c, m_c, q_c, pt_c, th_c, phi_c, idx_c, cspidx_c = (
            v1[sel_c] for v1 in (px, py, pz, m, q, pt, th, phi, idx, cspidx)
        )

        px_gen_c, py_gen_c, pz_gen_c, m_gen_c, q_gen_c, pt_gen_c, th_gen_c, phi_gen_c, idx_gen_c, cspidx_gen_c = (
            v2[sel_c_gen] for v2 in (px_gen, py_gen, pz_gen, m_gen, q_gen, pt_gen, th_gen, phi_gen, idx_gen, cspidx_gen)
        )

        px_n, py_n, pz_n, m_n, q_n, pt_n, th_n, phi_n = (
            v3[sel] for v3 in (px, py, pz, m, q, pt, th, phi)
        )

        px_gen_n, py_gen_n, pz_gen_n, m_gen_n, q_gen_n, pt_gen_n, th_gen_n, phi_gen_n = (
            v4[sel_gen] for v4 in (px_gen, py_gen, pz_gen, m_gen, q_gen, pt_gen, th_gen, phi_gen)
        )

        # --- Event selection and thrust ---
        e_c = np.sqrt(px_c**2 + py_c**2 + pz_c**2 + m_c**2)
        e_n = np.sqrt(px_n**2 + py_n**2 + pz_n**2 + m_n**2)

        e_gen_c = np.sqrt(px_gen_c**2 + py_gen_c**2 + pz_gen_c**2 + m_gen_c**2)
        e_gen_n = np.sqrt(px_gen_n**2 + py_gen_n**2 + pz_gen_n**2 + m_gen_n**2)

        rec_c = P4Block.build(px_c ,py_c ,pz_c ,q_c ,
                                pt_c ,th_c ,phi_c ,e_c ,m_c ,idx_c, cspidx_c)
        gen_c = P4Block.build(px_gen_c ,py_gen_c ,pz_gen_c ,q_gen_c ,
                                pt_gen_c ,th_gen_c ,phi_gen_c ,e_gen_c, m_gen_c, idx_gen_c, cspidx_gen_c)

        rec = P4Block.build(px_n ,py_n ,pz_n ,q_n ,
                                pt_n ,th_n ,phi_n ,e_n)
        gen = P4Block.build(px_gen_n ,py_gen_n ,pz_gen_n ,q_gen_n ,
                                pt_gen_n ,th_gen_n ,phi_gen_n ,e_gen_n)

        axis_n_met, T_n_met = thrust_axis_fast(rec.p, include_met=True)
        axis_gen_n_met, T_gen_n_met = thrust_axis_fast(gen.p, include_met=True)

        theta_Tu = thrust_theta(axis_n_met, T_n_met, fold=False)
        theta_gen_Tu = thrust_theta(axis_gen_n_met, T_n_met, fold=False)

        results = apply_event_selection_delphi(
            e_c, e_n, e_gen_c, e_gen_n, theta_Tu, theta_gen_Tu, E_reco, E_gen
        )

        pass_reco = results['pass_reco']

        #if not (pass_reco and pass_gen):
        if not pass_reco:
            continue
            
        histograms_1d['N'].Fill(1.5)

        ireco, igen, imiss, ifake = match_angular(rec_c, gen_c, 0.05)

        
        # Convert angles to degrees for histogram filling
        th_c_deg = np.degrees(th_c)
        phi_c_deg = np.degrees(phi_c)
        th_gen_c_deg = np.degrees(th_gen_c)
        phi_gen_c_deg = np.degrees(phi_gen_c)
        
        th_n_deg = np.degrees(th_n)
        phi_n_deg = np.degrees(phi_n)
        th_gen_n_deg = np.degrees(th_gen_n)
        phi_gen_n_deg = np.degrees(phi_gen_n)
        
        # =================================================================
        # CHARGED TRACK EFFICIENCY AND FAKE RATE FILLING
        # =================================================================
        
        # Fill efficiency histograms with ALL gen charged particles (denominator for efficiency)
        for i in range(len(pt_gen_c)):
            histograms_1d['trk_eff'].Fill(pt_gen_c[i])
            histograms_1d['trk_eff_theta'].Fill(th_gen_c_deg[i])
            histograms_1d['trk_eff_phi'].Fill(phi_gen_c_deg[i])
        print("DEBUG1", histograms_1d['trk_eff'].Integral())
        
        # Fill fake rate histograms with ALL reco charged particles
        for i in range(len(pt_c)):
            histograms_1d['trk_fake'].Fill(pt_c[i])
            histograms_1d['trk_fake_theta'].Fill(th_c_deg[i])
            histograms_1d['trk_fake_phi'].Fill(phi_c_deg[i])
        
        # =================================================================
        # CHARGED TRACK RESPONSE HISTOGRAMS (2D)
        # =================================================================
        
        # Fill 2D response histograms for matched particles only
        for reco_idx, gen_idx in zip(ireco, igen):
            # Angular response: reco - gen (in radians)
            theta_response = th_c[reco_idx] - th_gen_c[gen_idx]
            phi_response = phi_c[reco_idx] - phi_gen_c[gen_idx]
            
            # Fill 2D response histograms: X=angular_response, Y=gen_theta_deg
            histograms_2d['resp_trk_theta'].Fill(theta_response, th_gen_c_deg[gen_idx])
            histograms_2d['resp_trk_phi'].Fill(phi_response, th_gen_c_deg[gen_idx])
        
        # =================================================================
        # CHARGED TRACK 4D RESPONSE HISTOGRAMS
        # =================================================================
        
        # Fill 4D sparse histograms for resolution studies
        for reco_idx, gen_idx in zip(ireco, igen):
            # Resolution: pt_reco/pt_gen
            resolution = pt_c[reco_idx] / pt_gen_c[gen_idx] 
            
            # Fill 4D histogram: [resolution, pt_gen, theta_gen_deg, 0]
            histograms_nd['resp_trk'].Fill(np.array([resolution, pt_gen_c[gen_idx], th_gen_c_deg[gen_idx], 0]))
        
        # =================================================================
        # CHARGED TRACK MATCHING HISTOGRAMS (2D)
        # =================================================================
        
        # Fill matching correlation histograms for matched particles only
        print(zip(ireco, igen))
        for reco_idx, gen_idx in zip(ireco, igen):
            print(reco_idx, gen_idx)
            # True vs reco pT correlation
            histograms_2d['match_trk'].Fill(pt_c[reco_idx], pt_gen_c[gen_idx])
            
            histograms_2d['match_trk_theta'].Fill(th_c_deg[reco_idx], th_gen_c_deg[gen_idx])
            histograms_2d['match_trk_phi'].Fill(phi_c_deg[reco_idx], phi_gen_c_deg[gen_idx])

        print("DEBUG2", histograms_2d['match_trk'].Integral())
        
        # =================================================================
        # NEUTRAL TRACK HISTOGRAMS (using same matching function)
        # =================================================================
        
        # Match neutral particles using the same function as charged tracks
        ireco_n, igen_n, imiss_n, ifake_n = match_angular(rec, gen, 0.05)
        
        # Fill neutral track efficiency histograms with ALL gen neutral particles (using energy)
        for i in range(len(e_gen_n)):
            histograms_1d['ntr_eff'].Fill(e_gen_n[i])
            histograms_1d['ntr_eff_theta'].Fill(th_gen_n_deg[i])
            histograms_1d['ntr_eff_phi'].Fill(phi_gen_n_deg[i])
        
        # Fill neutral fake rate histograms with ALL reco neutral particles (using energy)
        for i in range(len(e_n)):
            histograms_1d['ntr_fake'].Fill(e_n[i])
            histograms_1d['ntr_fake_theta'].Fill(th_n_deg[i])
            histograms_1d['ntr_fake_phi'].Fill(phi_n_deg[i])
        
        # Fill neutral 2D response histograms for matched particles only
        for reco_idx, gen_idx in zip(ireco_n, igen_n):
            # Angular response
            theta_response = th_n[reco_idx] - th_gen_n[gen_idx]
            phi_response = phi_n[reco_idx] - phi_gen_n[gen_idx]
            
            # Fill 2D response histograms: X=angular_response, Y=gen_theta_deg  
            histograms_2d['resp_ntr_theta'].Fill(theta_response, th_gen_n_deg[gen_idx])
            histograms_2d['resp_ntr_phi'].Fill(phi_response, th_gen_n_deg[gen_idx])
        
        # Fill neutral 4D response histograms for matched particles only
        for reco_idx, gen_idx in zip(ireco_n, igen_n):
            # Resolution: E_reco/E_gen (using energy for neutrals)
            resolution = e_n[reco_idx] / e_gen_n[gen_idx] 

            # Fill 4D histogram: [resolution, energy_gen, theta_gen_deg, 0]
            histograms_nd['resp_ntr'].Fill(np.array([resolution, e_gen_n[gen_idx], th_gen_n_deg[gen_idx], 0]))
        
        # Fill neutral matching histograms for matched particles only
        for reco_idx, gen_idx in zip(ireco_n, igen_n):
            # Use energy for neutral particle matching plots
            histograms_2d['match_ntr'].Fill(e_n[reco_idx], e_gen_n[gen_idx])
            
            histograms_2d['match_ntr_theta'].Fill(th_n_deg[reco_idx], th_gen_n_deg[gen_idx])
            histograms_2d['match_ntr_phi'].Fill(phi_n_deg[reco_idx], phi_gen_n_deg[gen_idx])

    fout = ROOT.TFile.Open(args.outfile, 'recreate')

    for k in histograms_1d.keys():
        histograms_1d[k].Write()
    for k in histograms_2d.keys():
        histograms_2d[k].Write()
    for k in histograms_nd.keys():
        histograms_nd[k].Write()

    fout.Close()