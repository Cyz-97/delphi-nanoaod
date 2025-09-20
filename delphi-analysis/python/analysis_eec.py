import ROOT
import math
import numpy as np
from array import array
import itertools
import argparse
import sys
from common_functions import *
from binning_and_selections import *

# 0) A one‐bin "counter" histogram
h0 = ROOT.TH1D("N", "", 2, 0, 2)

# 1) Define all your 1D histos and their bin edges
h1d_defs = {
    "ETotal"        : np.linspace(0, 300, 301),     

    "EEC_r"         : rbins,
    "EEC_z"         : zbins,
    "EEC_r_trkSele" : rbins,
    "EEC_z_trkSele" : zbins,
    "EEC_r_before"  : rbins,
    "EEC_z_before"  : zbins,
    "EEC_r_pos"     : rbins,
    "EEC_z_pos"     : rbins,
    "EEC_r_neg"     : rbins,
    "EEC_z_neg"     : rbins,
    "EEC_r_cross"   : rbins,
    "EEC_z_cross"   : rbins,
    "theta_conv_ele": np.linspace(0, 4, 41),
    
    "Thrust_before"   : tbins,
    "Thrust_before2"  : tbins2,
    "Thrust_beforeDelphi" : tbinsDelphi,

    "Thrust_before_log" : logtbins, 
    "Thrust_before_log2": logtbins2,
    
    "ThrustC"       : tbins,
    "ThrustMissPC"  : tbins,
    "ThrustNC"      : tbins,
    "ThrustNCLog"   : logtbins,
    
    "ThrustMissPNC" : tbins,
    "ThrustMissPNC2" : tbins2,
    "ThrustMissPNCDelphi" : tbinsDelphi,
    
    "ThrustMissPNCLog" : logtbins,
    "ThrustMissPNCLog2" : logtbins2,

    "ThrustMissPNC_ALEPH" : tbins,
    "ThrustMissPNC2_ALEPH" : tbins2,
    "ThrustMissPNCLog_ALEPH" : logtbins,
    "ThrustMissPNCLog2_ALEPH" : logtbins2,

    "ThrustThetaMissPNC" : np.linspace(0, 180, 181),
    
    "MissPC"        : np.linspace(0, 100, 1001),
    "MissPNC"       : np.linspace(0, 100, 1001),

    "TrackE"        : ebins,
    "TrackPT"       : np.linspace(0, 100, 1001),
    "NuetralE"      : np.linspace(0, 100, 1001),
    "SumE"          : np.linspace(0, 200, 2001),
    "SumTrackE"     : np.linspace(0, 200, 2001),
    "TrackTheta"    : np.linspace(0, 180, 181),
    "TrackPhi"      : np.linspace(0, 180, 181),

    "TrackESele"    : ebins,
    "TrackPTSele"   : np.linspace(0, 100, 1001),
    "NuetralESele"  : np.linspace(0, 100, 1001),
    "SumESele"      : np.linspace(0, 200, 2001),
    "HJM"           : np.linspace(0, 200, 2001),
    "SumTrackESele" : np.linspace(0, 200, 2001),
    "TrackThetaSele": np.linspace(0, 180, 181),
    "TrackPhiSele"  : np.linspace(0, 180, 181),
    "NuetralThetaSele": np.linspace(0, 180, 181),
    "NuetralPhiSele": np.linspace(0, 180, 181),
    
    "NCharged"      : np.linspace(0, 100, 101),
    "NChargedTrkSele": np.linspace(0, 100, 101),
    "NChargedSele"  : np.linspace(0, 100, 101),
    "NPart"         : np.linspace(0, 200, 201),
}

# Build them in a single comprehension
h1d = {
    name: ROOT.TH1D(name, "", len(edges)-1, np.array(edges))
    for name, edges in h1d_defs.items()
}

# 2) Define your 2D histos: x–bins come from bins_theta/z, y–bins from eijbins2
h2d_defs = {
    "EEC2d_r"       : (rbins, eijbins),
    "EEC2d_z"       : (zbins,     eijbins),
    "EEC2d_r_pos"   : (rbins, eijbins),
    "EEC2d_z_pos"   : (zbins,     eijbins),
    "EEC2d_r_neg"   : (rbins, eijbins),
    "EEC2d_z_neg"   : (zbins,     eijbins),
    "EEC2d_r_cross" : (rbins, eijbins),
    "EEC2d_z_cross" : (zbins,     eijbins),
}

h2d = {
    name: ROOT.TH2D(
        name, "",
        len(xedges)-1, np.array(xedges),
        len(yedges)-1, np.array(yedges)
    )
    for name, (xedges, yedges) in h2d_defs.items()
}

HISTS = {
    'all':    (h1d['EEC_r'], h1d['EEC_z'],
               h2d['EEC2d_r'], h2d['EEC2d_z']),
    'pos':    (h1d['EEC_r_pos'], h1d['EEC_z_pos'],
               h2d['EEC2d_r_pos'], h2d['EEC2d_z_pos']),
    'neg':    (h1d['EEC_r_neg'], h1d['EEC_z_neg'],
               h2d['EEC2d_r_neg'], h2d['EEC2d_z_neg']),
    'cross':  (h1d['EEC_r_cross'], h1d['EEC_z_cross'],
               h2d['EEC2d_r_cross'], h2d['EEC2d_z_cross']),
}

# ---------------------------------------------------------------------------
# 3) Un-modified fill_pair – we just add a *lookup* argument
# ---------------------------------------------------------------------------
def fill_pair(i, j, r, z, eec, charges, table, weight=1):
    cat = 'neg' if charges[i]<0 and charges[j]<0 else \
          'pos' if charges[i]>0 and charges[j]>0 else 'cross'
    for tag in ('all', cat):
        h1,h2,h3,h4 = table[tag]
        h1.Fill(r, eec*weight)
        h2.Fill(z, eec*weight)
        h3.Fill(r, eec, weight)
        h4.Fill(z, eec, weight)

def calc_weight(pt):
    pt = np.asarray(pt)  # make sure it works for scalars & arrays
    return np.where(
        pt < 30, 
        1,
        np.where(pt <= 50, pt/30., 1./0.6)
    )

def calculate_event_eec_histogram(pairs_data, temp_hist, n):
    """
    Calculate single-event EEC histogram eec^(k)
    
    pairs_data: list of tuples [(jacobian_val, eij_val, weight), ...]
    template_hist: ROOT 2D histogram to get binning from
    nx, ny: total bins including overflow/underflow
    
    Returns: event histogram vector eec^(k) (flattened), including overflow bins
    """
    total_bins = n + 2
    
    # Create event histogram vector eec^(k)
    event_eec = np.zeros(total_bins)
    
    # Fill the event histogram - include ALL bins (no overflow check)
    for jacobian_val, weight in pairs_data:
        i_bin = template_hist.GetXaxis().FindBin(jacobian_val)  # Can be 0 (underflow) or nx+1 (overflow)
        event_eec[i_bin] += weight
    
    return event_eec

if __name__ == "__main__":

    filename = '/eos/user/z/zhangj/DELPHI/simulation/v94c/91.25/kk2f4146_qqpy/nanoaod_kk2f4146_qqpy_91.25_40001.sdst.root'
    #filename = "/afs/cern.ch/user/z/zhangj/private/DELPHI/delphi-nanoaod/nanoaod_pythia8_1.root"
    filenameout = 'h_test_with_covariance.root'
    isGen = False
    doChargedThrust = False
    
    parser = argparse.ArgumentParser()
    parser.add_argument("infiles", nargs='?', default=filename, help="name of input files")
    parser.add_argument("outfile", nargs='?', default=filenameout, help="name of output files")
    parser.add_argument("--jacobian", type=str, default="r", help="which variable to calculate covariance for (r or z)")
    parser.add_argument("--use_weights", action='store_true', default=False, help="use pt-dependent weights (default: no weights)")
    args = parser.parse_args()

    treename = "t"
    isALEPH  = False           # always defined

    if "ALEPH" in args.infiles:
      isALEPH = False # True / False
    if isGen:
        treename = "tgenBefore" 
    
    fin = ROOT.TFile.Open(args.infiles, 'r')
    t_hadrons = fin.Get(treename)

    fout = ROOT.TFile(args.outfile, 'RECREATE')

    # Set up covariance calculation
    jacobian_var = args.jacobian
    if jacobian_var == "r":
        template_hist = h1d["EEC_r"]
        print("Calculating covariance for r distribution")
        total_bins = len(rbins)-1
    elif jacobian_var == "z":
        template_hist = h1d["EEC_z"]
        print("Calculating covariance for z distribution")
        total_bins = len(zbins)-1
    else:
        print(f"ERROR: Unknown jacobian {jacobian_var}")
        sys.exit(1)

    # Initialize covariance calculation following your algorithm
    sum_of_eecs = np.zeros(total_bins+2)                       # Sum of eec^(k) vectors
    sum_of_eec_products = np.zeros((total_bins+2, total_bins+2)) # Sum of eec^(k)[i] * eec^(k)[j]

    N=0

    for iEvt in range(t_hadrons.GetEntries()):
        if iEvt % 1000 == 0:
            print(f"Processing event {iEvt}/{t_hadrons.GetEntries()}")
            
        t_hadrons.GetEntry(iEvt)

        E   = t_hadrons.Energy
        if abs(E - 91.25) > 1: continue
            
        get = lambda *names: (np.array(getattr(t_hadrons,n)) for n in names)

        if isGen:
            if isALEPH:
                px,py,pz,m,q,th,pt,eta,phi,pwflag = get(
                  'px','py','pz','mass','charge','theta','pt','eta','phi','pwflag')
            else:
                px,py,pz,m,q,th,pt,eta,phi,pwflag,hp = get(
                    'px','py','pz','mass','charge','theta','pt','eta','phi','pid','highPurity')
        else:
            px,py,pz,m,q,th,pt,eta,phi,d0,z0,pwflag = get(
                'px','py','pz','mass','charge','theta','pt','eta','phi','d0','z0','pwflag')

        sphericity = calculate_sphericity_with_fallback(px, py, pz)

        if isGen:
            if isALEPH:
                conversion = conversion_veto_mask(eta, phi, q, pwflag, 0.05, 0.05)
                mask_conv = ~conversion
                theta_conv = th
                theta_conv = theta_conv[mask_conv]

                for t in theta_conv:
                    h1d["theta_conv_ele"].Fill(t)
                
                ## For ALEPH sample, always veto conversion electron first
                px,py,pz,m,q,th,pt,eta,phi= (
                    arr[conversion] for arr in [px,py,pz,m,q,th,pt,eta,phi]
                )
            else:
                ## Pythia sample ISR photon tagged as highPurity == False
                ## DELPHI sample highPurity always True

                isr = (hp > -1)
                px,py,pz,m,q,th,pt,eta,phi= (
                    arr[isr] for arr in [px,py,pz,m,q,th,pt,eta,phi]
                )

        # --- Before any selection
        e  = np.sqrt(px**2 + py**2 + pz**2 + m**2)
        if isGen and np.sum(e)>E+0.1: continue
        h0.Fill(0.5)
        h1d["ETotal"].Fill(E)
        h1d["SumE"].Fill(np.sum(e))
        h1d["NPart"].Fill(len(px))

        sumtrke = 0
        sumtrk = 0
        for k in range(len(pt)):
            if abs(q[k]) > 0.1:
                h1d["TrackE"].Fill(e[k])
                h1d["TrackPT"].Fill(pt[k])
                sumtrke+=e[k]
                sumtrk+=1
                h1d["TrackTheta"].Fill(math.degrees(th[k]))
                h1d["TrackPhi"].Fill(math.degrees(phi[k]))
        h1d["SumTrackE"].Fill(sumtrke)
        h1d["NCharged"].Fill(sumtrk)
        
        p3 = np.stack((px, py, pz), axis=1)

        axis, T = thrust_axis_fast(p3, include_met=True)
        h1d["Thrust_before"].Fill(T)
        h1d["Thrust_before2"].Fill(T)
        h1d["Thrust_beforeDelphi"].Fill(T)
        
        h1d["Thrust_before_log"].Fill(np.log(1-T))
        h1d["Thrust_before_log2"].Fill(np.log(1-T))

        sel_before = (abs(q) > 0.5)
        p3 = p3[sel_before]
        e = e[sel_before]
        
        # |p|   (N,1) and (1,N) norms
        norm  = np.linalg.norm(p3, axis=1, keepdims=True)       # shape (N,1)
        denom = norm * norm.T                                  

        # cos(theta_ij)  =  (p_i · p_j) / (|p_i||p_j|)
        dot   = p3 @ p3.T                                     # shape (N,N)
        cos   = np.divide(dot, denom, out=np.ones_like(dot), where=denom>0)

        r     = np.arccos(np.clip(cos, -1.0, 1.0))                # angles (N,N)
        z     = 0.5 * (1 - cos)                                   # same shape

        # outer product of energies
        ee    = np.outer(e, e) / (E * E)                      # (N,N)

        # upper-triangle indices (i < j)
        iu, ju = np.triu_indices(len(p3), k=1)

        for i,j in zip(iu,ju):       # iu_sel from the reduced list
            h1d["EEC_r_before"].Fill(r[i, j], ee[i,j])
            h1d["EEC_z_before"].Fill(z[i, j], ee[i,j])
        
        # --- Track selection
        if not isGen:
            reco_results = apply_track_selection_delphi(px=px, py=py, pz=pz, m=m, q=q, th=th, pt=pt, eta=eta, phi=phi, d0=d0, z0=z0, pwflag=pwflag)
            sel_c = reco_results['sel_c']        # Charged particles only
            sel_n = reco_results['sel_n']        # Neutral particles only
            sel = reco_results['sel']            # All particles (charged + neutral)
        else:
            gen_results = apply_track_selection_delphi(px_gen=px, py_gen=py, pz_gen=pz, m_gen=m, q_gen=q, th_gen=th, pt_gen=pt, eta_gen=eta, phi_gen=phi, pwflag_gen=pwflag, hp_gen=hp)
            sel_c = gen_results['sel_c_gen']     # Generated charged particles only
            sel_n = gen_results['sel_n_gen']     # Generated neutral particles only
            sel = gen_results['sel_gen']         # All generated particles

        # --- Calc thrust and EEC
        px_c, py_c, pz_c, m_c, q_c, pt_c, th_c, phi_c = (
            v1[sel_c] for v1 in (px, py, pz, m, q, pt, th, phi)
        )
        
        e_c  = np.sqrt(px_c**2 + py_c**2 + pz_c**2 + m_c**2)
        p3_c = np.stack((px_c, py_c, pz_c), axis=1)

        px_nc, py_nc, pz_nc, m_nc, th_nc, phi_nc = (
            v2[sel] for v2 in (px, py, pz, m, th, phi)
        )
        e_nc  = np.sqrt(px_nc**2 + py_nc**2 + pz_nc**2 + m_nc**2)
        p3_nc = np.stack((px_nc, py_nc, pz_nc), axis=1)

        p4_nc = np.stack((e_nc, px_nc, py_nc, pz_nc), axis=1)

        px_n, py_n, pz_n, m_n, th_n, phi_n = (
            v3[sel_n] for v3 in (px, py, pz, m, th, phi)
        )
        e_n  = np.sqrt(px_n**2 + py_n**2 + pz_n**2 + m_n**2)
        p3_n = np.stack((px_n, py_n, pz_n), axis=1)
        
        axis_nc, T_nc = thrust_axis_fast(p3_nc)
        axis_nc_met, T_nc_met = thrust_axis_fast(p3_nc, include_met=True)

        theta_Tu = thrust_theta(axis_nc_met, T_nc_met, fold=False)

        h1d["NChargedTrkSele"].Fill(len(p3_c))

        # |p|   (N,1) and (1,N) norms
        norm  = np.linalg.norm(p3_c, axis=1, keepdims=True)       # shape (N,1)
        denom = norm * norm.T                                      

        # cos(theta_ij)  =  (p_i · p_j) / (|p_i||p_j|)
        dot   = p3_c @ p3_c.T                                     # shape (N,N)
        cos   = np.divide(dot, denom, out=np.ones_like(dot), where=denom>0)

        r     = np.arccos(np.clip(cos, -1.0, 1.0))                # angles (N,N)
        z     = 0.5 * (1 - cos)                                   # same shape

        # outer product of energies
        ee    = np.outer(e_c, e_c) / (E * E)                      # (N,N)

        # upper-triangle indices (i < j)
        iu, ju = np.triu_indices(len(p3_c), k=1)

        for i,j in zip(iu,ju):
            h1d["EEC_r_trkSele"].Fill(r[i, j], ee[i,j])
            h1d["EEC_z_trkSele"].Fill(z[i, j], ee[i,j])

        if len(p3_c) >= 7 and np.sum(e_nc) > 0.5*E:
            h1d["ThrustThetaMissPNC"].Fill(np.degrees(theta_Tu))

        pass_delphi = apply_event_selection_delphi(e_c=e_c, e_n=e_nc, theta_Tu=theta_Tu, E_reco=E)['pass_reco']
        pass_aleph = apply_event_selection_aleph(e_c=e_c, e_nc=e_nc, sphericity=sphericity["cos_theta_v1"])['pass_reco']

        if pass_aleph:
            h1d["ThrustMissPNC_ALEPH"].Fill(T_nc_met)
            h1d["ThrustMissPNC2_ALEPH"].Fill(T_nc_met)
            h1d["ThrustMissPNCLog_ALEPH"].Fill(np.log(1-T_nc_met))
            h1d["ThrustMissPNCLog2_ALEPH"].Fill(np.log(1-T_nc_met))
        
        if pass_delphi:
            h0.Fill(1.5)
            
            h1d["NChargedSele"].Fill(len(p3_c))
            h1d["SumESele"].Fill(np.sum(e_nc))
            h1d["SumTrackESele"].Fill(np.sum(e_c))
            for k in range(len(pt_c)):
                weight = calc_weight(pt_c[k]) if args.use_weights else 1.0
                h1d["TrackPTSele"].Fill(pt_c[k], weight)
                h1d["TrackESele"].Fill(e_c[k], weight)
                h1d["TrackThetaSele"].Fill(math.degrees(th_c[k]), weight)
                h1d["TrackPhiSele"].Fill(math.degrees(phi_c[k]), weight)

            for l in range(len(e_n)):
                h1d["NuetralESele"].Fill(e_n[l])
                h1d["NuetralThetaSele"].Fill(math.degrees(th_n[l]))
                h1d["NuetralPhiSele"].Fill(math.degrees(phi_n[l]))

            p_miss_c = missing_p(p3_c)
            p_miss_nc = missing_p(p3_nc)

            h1d["MissPC"].Fill(np.linalg.norm(p_miss_c))
            h1d["MissPNC"].Fill(np.linalg.norm(p_miss_nc))

            if doChargedThrust:
                axis_c, T_c = thrust_axis(p3_c)
                axis_c_met, T_c_met = thrust_axis(p3_c, include_met=True)
                h1d["ThrustC"].Fill(T_c)
                h1d["ThrustMissPC"].Fill(T_c_met)

            h1d["ThrustNC"].Fill(T_nc)
            h1d["ThrustMissPNC"].Fill(T_nc_met)
            h1d["ThrustMissPNC2"].Fill(T_nc_met)
            h1d["ThrustMissPNCDelphi"].Fill(T_nc_met)
            h1d["ThrustNCLog"].Fill(np.log(1-T_nc))
            h1d["ThrustMissPNCLog"].Fill(np.log(1-T_nc_met))
            h1d["ThrustMissPNCLog2"].Fill(np.log(1-T_nc_met))

            M_h = heavy_jet_mass(p4_nc, axis_nc_met)
            h1d["HJM"].Fill(M_h)

            # ===== COVARIANCE CALCULATION =====
            # Step 1: Calculate single-event EEC histogram eec^(k)
            event_pairs = []
            
            for idx_i, idx_j in zip(iu, ju):
                # Get the values we need
                r_val = r[idx_i, idx_j]
                z_val = z[idx_i, idx_j]
                eij_val = ee[idx_i, idx_j]
                
                # Choose which jacobian to use
                if jacobian_var == "r":
                    jacobian_val = r_val
                elif jacobian_var == "z":
                    jacobian_val = z_val
                
                # Apply weight if requested
                if args.use_weights:
                    weight = calc_weight(pt_c[idx_i]) * calc_weight(pt_c[idx_j])
                else:
                    weight = 1.0
                
                # Store for covariance calculation (jacobian_val for binning, eij_val as weight)
                event_pairs.append((jacobian_val, eij_val))
                
                # Fill normal histograms
                fill_pair(idx_i, idx_j, r_val, z_val, eij_val, q_c, HISTS, weight)
            
            # Step 2: Calculate this event's EEC histogram eec^(k)
            event_eec = calculate_event_eec_histogram(event_pairs, template_hist, total_bins)
            
            # Step 3: Add to running sums following your algorithm
            sum_of_eecs += event_eec
            # Calculate outer product eec^(k)[i] * eec^(k)[j] and add to running sum
            sum_of_eec_products += np.outer(event_eec, event_eec)

        N += 1

    print("Processed", N, "events")
    
    # Write results - only the raw sums for later combination
    fout.cd()
    h0.Write()
    for key in h1d.keys():
        h1d[key].Write()
    for key in h2d.keys():
        h2d[key].Write()


    h_sum_of_eec=ROOT.TH1D("sum_of_eec", "", 200, 0, 200)
    h_sum_of_eec_products=ROOT.TH2D("sum_of_eec_products", "", 200, 0, 200, 200, 0, 200)

    for i in range(1, 201):
        h_sum_of_eec.SetBinContent(i, sum_of_eecs[i])
        for j in range(1, 201):
            h_sum_of_eec_products.SetBinContent(i, j, sum_of_eec_products[i,j])

    
    h_sum_of_eec.Write()
    h_sum_of_eec_products.Write()
