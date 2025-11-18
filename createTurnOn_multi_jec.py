#!/usr/bin/env python

import argparse
from array import array
import math
import numpy as np
import os
import re
import sys
import ROOT

parser = argparse.ArgumentParser(description='Create turn on curves.')
parser.add_argument('--input-data', required=True, type=str, help="skimmed data input")
parser.add_argument('--input-dy-mc', required=True, type=str, help="skimmed DY MC input")
parser.add_argument('--output', required=True, type=str, help="output file prefix")
parser.add_argument('--channels', required=False, type=str, default='ditaujet_jetleg60', help="channels to process")
# parser.add_argument('--channels', required=False, type=str, default='ditau,ditau_withptiso_nobitcut,ditau_withptiso_bit1,ditau_withptiso_bit1_bit17,ditau_withptiso_bit1_bit17_0bit18', help="channels to process")
# parser.add_argument('--decay-modes', required=False, type=str, default='all,0,1,10,11,1011', help="decay modes to process")
parser.add_argument('--working-points', required=False, type=str,
                    default='VVVLoose,VVLoose,VLoose,Loose,Medium,Tight,VTight,VVTight',
                    help="working points to process")
args = parser.parse_args()

path_prefix = '' if 'TAU-Trigger-NANO' in os.getcwd() else 'TAU-Trigger-NANO/'
sys.path.insert(0, path_prefix + 'Common/python')
from AnalysisTypes import *
from AnalysisTools import *
import RootPlotting
ROOT.ROOT.EnableImplicitMT(4)
ROOT.gROOT.SetBatch(True)
ROOT.TH1.SetDefaultSumw2()
RootPlotting.ApplyDefaultGlobalStyle()

bin_scans = {
    # 2:  [ 0.01 ],
    5:  [ 0.01, 0.05 ],
    10: [ 0.05, 0.1 ],
    20: [ 0.1, 0.2 ],
    30: [0.2, 0.4],
    50: [ 0.4, 0.5 ],
    100: [ 0.5 ],
    150: [ 1 ]
}
bin_scan_pairs = []
for max_bin_delta_pt, max_rel_err_vec in bin_scans.items():
    for max_rel_err in max_rel_err_vec:
        bin_scan_pairs.append([max_bin_delta_pt, max_rel_err])


def CreateBins(max_pt, for_fitting):
    if for_fitting:
        step=1
        return np.arange(40, 1000+step, step=step), False
        #high_pt_bins = np.arange(100, 501, step=5)
    else:
        bins = np.arange(20, 30, step=5)
        bins = np.append(bins, np.arange(30, 70, step=4))
        #bins = np.append(bins, np.arange(60, 100, step=10))
        high_pt_bins = [ 70, 100, 150, 200]
        use_logx = max_pt > 200
        return np.append(bins, high_pt_bins), use_logx

class TurnOnData:
    def __init__(self):
        self.hist_total = None
        self.hist_passed = None
        self.eff = None

def CreateHistograms(input_file, channels, discr_name, working_points, hist_models, label, var,
                     output_file):
    df = ROOT.RDataFrame('Events', input_file)
    turnOn_data = {}
    dm_labels = {}

    # for wp in working_points:
    #     wp_bit = ParseEnum(DiscriminatorWP, wp)
        # df_wp = df_dm.Filter('({} & (1 << {})) != 0'.format(discr_name, wp_bit))
        # df = df.Filter('{0} >= {1}'.format(discr_name, wp_bit))
    # for var_mc in ["_nom", "_corrRegrouped_BBEC1_down","_corrRegrouped_BBEC1_up","_corrRegrouped_HF_down","_corrRegrouped_HF_up","_corrRegrouped_EC2_2024_down","_corrRegrouped_EC2_2024_up","_corrRegrouped_EC2_down","_corrRegrouped_EC2_up","_corrRegrouped_RelativeSample_2024_down","_corrRegrouped_RelativeSample_2024_up","_corrRegrouped_BBEC1_2024_down","_corrRegrouped_BBEC1_2024_up","_corrRegrouped_Absolute_down","_corrRegrouped_Absolute_up","_corrRegrouped_FlavorQCD_down","_corrRegrouped_FlavorQCD_up","_corrRegrouped_HF_2024_down","_corrRegrouped_HF_2024_up","_corrRegrouped_RelativeBal_down","_corrRegrouped_RelativeBal_up","_corrRegrouped_Total_down","_corrRegrouped_Total_up","_corrRegrouped_Absolute_2024_down","_corrRegrouped_Absolute_2024_up"]:
    for var_mc in ["_nom"]:
        var_used = var
        if label == 'mc':
            var_used = var + var_mc
        turnOn_data[var_mc] = {}
        for title, eta_cut in abs_jet_eta_cuts:
            eta_cut_label = title.strip().replace("<", "_").replace(">", "_")
            turnOn_data[var_mc][eta_cut_label] = {}
            if label == 'mc':
                eta_cut = eta_cut.replace("eta", "eta" + var_mc)
            for channel in channels:
                turnOn_data[var_mc][eta_cut_label][channel] = {}
                if label == 'mc':
                    channel_used = channel + var_mc
                else:
                    channel_used = channel
                df_tot = df.Filter(eta_cut).Filter('pass_ditaujet > 0.5 && leading_jet_pt > 20')
                print(eta_cut)
                df_ch = df_tot.Filter('pass_{} > 0.5'.format(channel_used))
                for model_name, hist_model in hist_models.items():
                    turn_on = TurnOnData()
                    turn_on.hist_total = df_tot.Histo1D(hist_model, var_used, 'weight')
                    turn_on.hist_passed = df_ch.Histo1D(hist_model, var_used, 'weight')
                    turnOn_data[var_mc][eta_cut_label][channel][model_name] = turn_on

            for channel in channels:
                for model_name, hist_model in hist_models.items():
                    turn_on = turnOn_data[var_mc][eta_cut_label][channel][model_name]
                    name_pattern = '{}_{}_{}_{}_{}_{{}}'.format(var_mc, eta_cut_label, label, channel, model_name)
                    turn_on.name_pattern = name_pattern
                    if 'fit' in model_name:
                        passed, total, eff = AutoRebinAndEfficiency(turn_on.hist_passed.GetPtr(),
                                                                    turn_on.hist_total.GetPtr(), bin_scan_pairs)
                    else:
                        passed, total = turn_on.hist_passed.GetPtr(), turn_on.hist_total.GetPtr()
                        # # new add by botao
                        # if (passed.Integral() < 0) or (total.Integral() < 0):
                        #     continue
                        # # end add
                        try:
                            FixEfficiencyBins(passed, total)
                        except:
                            continue
                        turn_on.eff = ROOT.TEfficiency(passed, total)
                        eff = turn_on.eff
                    output_file.WriteTObject(total, name_pattern.format('total'), 'Overwrite')
                    output_file.WriteTObject(passed, name_pattern.format('passed'), 'Overwrite')
                    output_file.WriteTObject(eff, name_pattern.format('eff'), 'Overwrite')
                    # print(name_pattern)
                    # print('hist_total {}'.format(turn_on.hist_total.GetPtr().GetNbinsX()))
                    # for n in range(turn_on.hist_total.GetPtr().GetNbinsX() + 1):
                    #     print('\t{} {} +/- {} {} +/- {}'.format(turn_on.hist_total.GetPtr().GetBinLowEdge(n+1), turn_on.hist_total.GetPtr().GetBinContent(n+1), turn_on.hist_total.GetPtr().GetBinError(n+1), turn_on.hist_passed.GetPtr().GetBinContent(n+1), turn_on.hist_passed.GetPtr().GetBinError(n+1)))
                    # print('hist_passed {}'.format(turn_on.hist_passed.GetPtr().GetNbinsX()))
                    # for n in range(turn_on.hist_passed.GetPtr().GetNbinsX() + 1):
                    #     print('\t{} {}'.format(turn_on.hist_passed.GetPtr().GetBinLowEdge(n+1), ))
                    #if 'fit' not in model_name:
                    #    turn_on.eff = ROOT.TEfficiency(turn_on.hist_passed.GetPtr(), turn_on.hist_total.GetPtr())
                    #    output_file.WriteTObject(turn_on.eff, name_pattern.format('passed'), 'Overwrite')

    return turnOn_data

output_file = ROOT.TFile(args.output + '.root', 'RECREATE')
input_files = [ args.input_data, args.input_dy_mc ]
n_inputs = len(input_files)
labels = [ 'data', 'mc' ]
var = 'leading_jet_pt'
title, x_title = 'jet p_{T}', 'jet p_{T} (GeV)'
channels = args.channels.split(',')
abs_jet_eta_cuts = [('|jet #eta| < 2.5', 'abs(leading_jet_eta) < 2.5',),('2.5 < |jet #eta| < 3', 'abs(leading_jet_eta) < 3 && abs(leading_jet_eta) > 2.5',),('3 < |jet #eta| < 5', 'abs(leading_jet_eta) < 5 && abs(leading_jet_eta) > 3',)]
working_points = args.working_points.split(',')
bins, use_logx = CreateBins(200, False)
bins_fit, _ = CreateBins(200, True)
hist_models = {
    'plot': ROOT.RDF.TH1DModel(var, var, len(bins) - 1, array('d', bins)),
    'fit': ROOT.RDF.TH1DModel(var, var, len(bins_fit) - 1, array('d', bins_fit))
}
turnOn_data = [None] * n_inputs
for input_id in range(n_inputs):
    print("Creating {} histograms...".format(labels[input_id]))
    turnOn_data[input_id] = CreateHistograms(input_files[input_id], channels, 'leading_tau_idDeepTauVSjet', # tau_idDeepTau2017v2p1VSjet,
                                             ["Medium"], hist_models, labels[input_id], var, output_file)

colors = [ ROOT.kRed, ROOT.kBlack ]
canvas = RootPlotting.CreateCanvas()

n_plots = len(abs_jet_eta_cuts) * len(channels) * 1
plot_id = 0
# for var_mc in ["_nom", "_corrRegrouped_BBEC1_down","_corrRegrouped_BBEC1_up","_corrRegrouped_HF_down","_corrRegrouped_HF_up","_corrRegrouped_EC2_2024_down","_corrRegrouped_EC2_2024_up","_corrRegrouped_EC2_down","_corrRegrouped_EC2_up","_corrRegrouped_RelativeSample_2024_down","_corrRegrouped_RelativeSample_2024_up","_corrRegrouped_BBEC1_2024_down","_corrRegrouped_BBEC1_2024_up","_corrRegrouped_Absolute_down","_corrRegrouped_Absolute_up","_corrRegrouped_FlavorQCD_down","_corrRegrouped_FlavorQCD_up","_corrRegrouped_HF_2024_down","_corrRegrouped_HF_2024_up","_corrRegrouped_RelativeBal_down","_corrRegrouped_RelativeBal_up","_corrRegrouped_Total_down","_corrRegrouped_Total_up","_corrRegrouped_Absolute_2024_down","_corrRegrouped_Absolute_2024_up"]:
for var_mc in ["_nom"]:
    var_used = var
    for eta_title, eta_cut in abs_jet_eta_cuts:
        eta_cut_label = eta_title.strip().replace("<", "_").replace(">", "_")
        for channel in channels:
            ratio_graph = None
            ref_hist = hist_models['plot'].GetHistogram()
            ratio_ref_hist = ref_hist.Clone()
            turnOns = [None] * n_inputs
            curves = [None] * n_inputs
            for input_id in range(n_inputs):
                turnOns[input_id] = turnOn_data[input_id][var_mc][eta_cut_label][channel]['plot']
                curves[input_id] = turnOns[input_id].eff
            y_min, y_max = (0, 1)
            y_title = 'Efficiency'
            title = '{}'.format(channel)
            plain_title = '{}'.format(channel)
            main_pad, ratio_pad, title_controls = RootPlotting.CreateTwoPadLayout(canvas, ref_hist, ratio_ref_hist,
                                                                                    log_x=use_logx, title=title)
            RootPlotting.ApplyAxisSetup(ref_hist, ratio_ref_hist, x_title=x_title, y_title=y_title,
                                        ratio_y_title='Ratio', y_range=(y_min, y_max * 1.1), max_ratio=1.5)
            legend = RootPlotting.CreateLegend(pos=(0.78, 0.28), size=(0.2, 0.15))
            for input_id in range(n_inputs):
                curve = curves[input_id]
                try:
                    curve.Draw('SAME')
                except:
                    continue
                RootPlotting.ApplyDefaultLineStyle(curve, colors[input_id])
                legend.AddEntry(curve, labels[input_id], 'PLE')

                if input_id < n_inputs - 1:
                    ratio_graph = RootPlotting.CreateEfficiencyRatioGraph(turnOns[input_id].hist_passed,
                                                                            turnOns[input_id].hist_total,
                                                                            turnOns[-1].hist_passed,
                                                                            turnOns[-1].hist_total)
                    if ratio_graph:
                        output_file.WriteTObject(ratio_graph, 'ratio_{}'.format(plain_title), 'Overwrite')
                        ratio_pad.cd()
                        ratio_color = colors[input_id] if n_inputs > 2 else ROOT.kBlack
                        RootPlotting.ApplyDefaultLineStyle(ratio_graph, ratio_color)
                        ratio_graph.Draw("0PE SAME")
                        main_pad.cd()
            legend.Draw()

            canvas.Update()
            output_file.WriteTObject(canvas, 'canvas_{}'.format(plain_title), 'Overwrite')
            RootPlotting.PrintAndClear(canvas, args.output + '.pdf', plain_title, plot_id, n_plots,
                                        [ main_pad, ratio_pad ])
            plot_id += 1
output_file.Close()
