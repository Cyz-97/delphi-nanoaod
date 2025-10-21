import ROOT
import numpy as np
import argparse
from array import array
from common_functions import *
from binning_and_selections import *

class ThrustResponse:

    def __init__(self, reco_tree, gen_tree, weight=False):
        self._hists = {}
        self._treco = reco_tree
        self._tgen = gen_tree
        self._evt_counter = 0
        self._isALEPH = False
        self._weight = weight

    def setALEPH(self):
        self._isALEPH = True

    def writeToFile(self, output):
        fout = ROOT.TFile(output, 'recreate')
        fout.cd()
        for key in self._hists.keys():
            self._hists[key].Write()
        fout.Close()

    def bookHistograms(self):
        """Book only thrust-related histograms"""
        h1d_defs = {
            "counter"             : np.array([0., 1., 2., 3.], dtype="f"),
            "gen_thrust"          : tbins,
            "reco_thrust"         : tbins,
            "gen_thrust_log"      : logtbins,
            "reco_thrust_log"     : logtbins,
            "gen_thrust2"         : tbins2,
            "reco_thrust2"        : tbins2,
            "gen_thrustDelphi"    : tbinsDelphi,
            "reco_thrustDelphi"   : tbinsDelphi,
            "gen_thrustDelphi_c"    : tbinsDelphi,
            "reco_thrustDelphi_c"   : tbinsDelphi,
            "gen_thrust_log2"     : logtbins2,
            "reco_thrust_log2"    : logtbins2,
            "reco_thrust_theta"   : np.array(np.linspace(0, 180, 181)),
            "gen_thrust_theta"    : np.array(np.linspace(0, 180, 181)),
        }

        h2d_defs = {
            "response_thrust"         : (tbins, tbins),
            "response_thrust_log"     : (logtbins, logtbins),
            "response_thrust2"        : (tbins2, tbins2),
            "response_thrustDelphi"   : (tbinsDelphi, tbinsDelphi),
            "response_thrustDelphi_c"   : (tbinsDelphi, tbinsDelphi),
            "response_thrust_log2"    : (logtbins2, logtbins2),
        }

        self._hists = {}

        # Book 1-D histograms
        self._hists.update({
            name: ROOT.TH1F(name, "", len(edges) - 1, edges)
            for name, edges in h1d_defs.items()
        })

        # Book 2-D response matrices
        self._hists.update({
            name: ROOT.TH2D(name, "", len(xe) - 1, xe, len(ye) - 1, ye)
            for name, (xe, ye) in h2d_defs.items()
        })

    def loop(self):
        nevt = self._treco.GetEntries()
        for ievt in range(nevt):
            if ievt % 1000 == 0:
                print(f"Processing event {ievt}/{nevt}")

            treco, tgen = self._treco, self._tgen
            treco.GetEntry(ievt)
            tgen.GetEntry(ievt)

            try:
                E_reco = treco.Energy
            except:
                E_reco = 91.25

            try:
                E_gen = tgen.Energy
            except:
                E_gen = 91.25

            get = lambda *names: (np.asarray(getattr(treco, n)) for n in names)
            get_gen = lambda *names: (np.asarray(getattr(tgen, n)) for n in names)

            if self._isALEPH:
                px, py, pz, m, q, th, pt, eta, phi, d0, z0, pwflag, hp = get(
                    'px', 'py', 'pz', 'mass', 'charge', 'theta', 'pt', 'eta', 'phi', 'd0', 'z0', 'pwflag', 'highPurity')
                px_gen, py_gen, pz_gen, m_gen, q_gen, th_gen, pt_gen, eta_gen, phi_gen, pwflag_gen, hp_gen = get_gen(
                    'px', 'py', 'pz', 'mass', 'charge', 'theta', 'pt', 'eta', 'phi', 'pwflag', 'highPurity')
            else:
                px, py, pz, m, q, th, pt, eta, phi, d0, z0, pwflag, hp, idx, cspidx = get(
                    'px', 'py', 'pz', 'mass', 'charge', 'theta', 'pt', 'eta', 'phi', 'd0', 'z0', 'pwflag', 'highPurity', 'index', 'correspondenceIndex')
                px_gen, py_gen, pz_gen, m_gen, q_gen, th_gen, pt_gen, eta_gen, phi_gen, pwflag_gen, hp_gen, idx_gen, cspidx_gen = get_gen(
                    'px', 'py', 'pz', 'mass', 'charge', 'theta', 'pt', 'eta', 'phi', 'pwflag', 'highPurity', 'index', 'correspondenceIndex')

            sphericity = calculate_sphericity_with_fallback(px, py, pz)
            sphericity_gen = calculate_sphericity_with_fallback(px_gen, py_gen, pz_gen)

            # Gen level selection --- this might create selection bias, be careful!
            e = np.sqrt(px_gen**2 + py_gen**2 + pz_gen**2 + m_gen**2)
            #if (np.sum(e) - E_gen) > 0.1:
            #    continue

            self._hists['counter'].Fill(0.5)

            # Apply track selections
            if self._isALEPH:
                all_results = apply_track_selection_delphi(
                    px=px, py=py, pz=pz, m=m, q=q, th=th, pt=pt, eta=eta, phi=phi, d0=d0, z0=z0, pwflag=pwflag, hp=hp,
                    px_gen=px_gen, py_gen=py_gen, pz_gen=pz_gen, m_gen=m_gen, q_gen=q_gen,
                    th_gen=th_gen, pt_gen=pt_gen, eta_gen=eta_gen, phi_gen=phi_gen,
                    pwflag_gen=pwflag_gen, hp_gen=hp_gen)
            else:
                all_results = apply_track_selection_delphi(
                    px=px, py=py, pz=pz, m=m, q=q, th=th, pt=pt, eta=eta, phi=phi, d0=d0, z0=z0, pwflag=pwflag, hp=hp,
                    px_gen=px_gen, py_gen=py_gen, pz_gen=pz_gen, m_gen=m_gen, q_gen=q_gen,
                    th_gen=th_gen, pt_gen=pt_gen, eta_gen=eta_gen, phi_gen=phi_gen,
                    pwflag_gen=pwflag_gen, hp_gen=hp_gen)


            sel_c = all_results['sel_c']  
            sel_gen_c = all_results['sel_c_gen']  

            px_c, py_c, pz_c, m_c, q_c, pt_c, th_c, phi_c = (
                v[sel_c] for v in (px, py, pz, m, q, pt, th, phi)
            )

            px_gen_c, py_gen_c, pz_gen_c, m_gen_c, q_gen_c, pt_gen_c, th_gen_c, phi_gen_c = (
                v[sel_gen_c] for v in (px_gen, py_gen, pz_gen, m_gen, q_gen, pt_gen, th_gen, phi_gen)
            )

            e_c = np.sqrt(px_c**2 + py_c**2 + pz_c**2 + m_c**2)  # Charged only
            e_gen_c = np.sqrt(px_gen_c**2 + py_gen_c**2 + pz_gen_c**2 + m_gen_c**2)  # Charged only

            rec_c = P4Block.build(px_c, py_c, pz_c, q_c, pt_c, th_c, phi_c, e_c)
            gen_c = P4Block.build(px_gen_c, py_gen_c, pz_gen_c, q_gen_c, pt_gen_c, th_gen_c, phi_gen_c, e_gen_c)

            sel = all_results['sel']
            sel_gen = all_results['sel_gen']

            px_n, py_n, pz_n, m_n, q_n, pt_n, th_n, phi_n = (
                v[sel] for v in (px, py, pz, m, q, pt, th, phi)
            )

            px_gen_n, py_gen_n, pz_gen_n, m_gen_n, q_gen_n, pt_gen_n, th_gen_n, phi_gen_n = (
                v[sel_gen] for v in (px_gen, py_gen, pz_gen, m_gen, q_gen, pt_gen, th_gen, phi_gen)
            )

            e_n = np.sqrt(px_n**2 + py_n**2 + pz_n**2 + m_n**2)
            e_gen_n = np.sqrt(px_gen_n**2 + py_gen_n**2 + pz_gen_n**2 + m_gen_n**2)

            rec = P4Block.build(px_n, py_n, pz_n, q_n, pt_n, th_n, phi_n, e_n)
            gen = P4Block.build(px_gen_n, py_gen_n, pz_gen_n, q_gen_n, pt_gen_n, th_gen_n, phi_gen_n, e_gen_n)

            # Calculate thrust
            axis_n_met, T_n_met = thrust_axis_fast(rec.p, include_met=True)
            axis_gen_n_met, T_gen_n_met = thrust_axis_fast(gen.p, include_met=False)

            axis_c, T_c = thrust_axis_fast(rec_c.p, include_met=False)
            axis_gen_c, T_gen_c = thrust_axis_fast(gen_c.p, include_met=False)

            theta_Tu = thrust_theta(axis_n_met, T_n_met, fold=False)
            theta_gen_Tu = thrust_theta(axis_gen_n_met, T_gen_n_met, fold=False)

            self._hists["reco_thrust_theta"].Fill(np.degrees(theta_Tu))
            self._hists["gen_thrust_theta"].Fill(np.degrees(theta_gen_Tu))

            # Event selection
            if self._isALEPH:
                results = apply_event_selection_aleph(
                    e_c, e_n, e_gen_c, e_gen_n, sphericity["cos_theta_v1"], sphericity_gen["cos_theta_v1"]
                )
            else:
                results = apply_event_selection_delphi(
                    e_c, e_n, e_gen_c, e_gen_n, theta_Tu, theta_gen_Tu, E_reco, E_gen
                )

            pass_reco = results['pass_reco']

            if not pass_reco:
                continue

            if self._weight:
                evt_weight = calc_multiplicity_weight_linear(len(px_n))
            else:
                evt_weight = 1.0

            self._hists['counter'].Fill(1.5, evt_weight)
            self._evt_counter += 1

            # Fill thrust histograms
            self._hists["reco_thrust"].Fill(1 - T_n_met, evt_weight)
            self._hists["reco_thrust_log"].Fill(np.log(1 - T_n_met), evt_weight)
            self._hists["gen_thrust"].Fill(1 - T_gen_n_met, evt_weight)
            self._hists["gen_thrust_log"].Fill(np.log(1 - T_gen_n_met), evt_weight)
            self._hists["response_thrust"].Fill(1 - T_n_met, 1 - T_gen_n_met, evt_weight)
            self._hists["response_thrust_log"].Fill(np.log(1 - T_n_met), np.log(1 - T_gen_n_met), evt_weight)

            self._hists["reco_thrust2"].Fill(1 - T_n_met, evt_weight)
            self._hists["reco_thrustDelphi"].Fill(1 - T_n_met, evt_weight)
            self._hists["reco_thrustDelphi_c"].Fill(1 - T_c, evt_weight)
            self._hists["reco_thrust_log2"].Fill(np.log(1 - T_n_met), evt_weight)
            self._hists["gen_thrust2"].Fill(1 - T_gen_n_met, evt_weight)
            self._hists["gen_thrustDelphi"].Fill(1 - T_gen_n_met, evt_weight)
            self._hists["gen_thrustDelphi_c"].Fill(1 - T_gen_c, evt_weight)
            self._hists["gen_thrust_log2"].Fill(np.log(1 - T_gen_n_met), evt_weight)
            self._hists["response_thrust2"].Fill(1 - T_n_met, 1 - T_gen_n_met, evt_weight)
            self._hists["response_thrustDelphi"].Fill(1 - T_n_met, 1 - T_gen_n_met, evt_weight)
            self._hists["response_thrustDelphi_c"].Fill(1 - T_c, 1 - T_gen_c, evt_weight)
            self._hists["response_thrust_log2"].Fill(np.log(1 - T_n_met), np.log(1 - T_gen_n_met), evt_weight)


if __name__ == "__main__":
    filename = '/eos/user/z/zhangj/DELPHI/simulation/v94c/91.25/kk2f4146_qqpy/nanoaod_kk2f4146_qqpy_91.25_40001.sdst.root'
    filenameout = "response_thrust_only.root"

    parser = argparse.ArgumentParser()
    parser.add_argument("infiles", nargs='?', default=filename, help="name of input files")
    parser.add_argument("outfile", nargs='?', default=filenameout, help="name of output file")
    parser.add_argument("--use_evt_weights", action='store_true', default=False, help="use multiplicity-dependent weights (default: no weights)")
    args = parser.parse_args()

    treco = 't'
    if "ALEPH" in args.infiles:
        tgen = 'tgen'
    else:
        tgen = 'tgenBefore'

    t_reco = ROOT.TChain(treco)
    t_gen = ROOT.TChain(tgen)

    print("Reading input from:", args.infiles)
    InputRootFiles = []
    if args.infiles.find(".root") > -1:
        InputRootFiles.append(args.infiles)
    else:
        InputRootFiles = ReadFilesFromList(args.infiles)

    for f in InputRootFiles:
        t_reco.Add(f)
        t_gen.Add(f)

    fnameout = args.outfile

    weight = args.use_evt_weights
    response = ThrustResponse(t_reco, t_gen, weight)
    if "ALEPH" in args.infiles:
        response.setALEPH()
        print("ALEPH mode enabled")
    
    response.bookHistograms()
    response.loop()
    response.writeToFile(fnameout)
    
    print(f"Thrust response matrices written to {fnameout}")
    print(f"Total events processed: {response._evt_counter}")