import ROOT

#TH1F *hnew = (TH1F*)gDirectory->Get("hnew");

f = ROOT.TFile("out_mc_full_matched_byReference.root")
dirr = f.GetDirectory("demo");
f.cd("demo")
# f.ls()
# dirr.ls()

myTree = dirr.Get("kaon_tree")
myTree.ls()
#gPad->Update();
myTree.Draw("v0_Ks_inv_m_pi:hps_pt_1>>hname", "v0_Ks_inv_m_pi>0 && hps_pt_1>0")

ROOT.gDirectory.ls()

hname = ROOT.gDirectory.Get("hname")
hname.SetDirectory(0)
hname.Print("all")
ROOT.gDirectory.ls()

hname.FitSlicesY()
hname_1 = ROOT.gDirectory.Get("hname_1");
hname_1.Draw()
ROOT.gDirectory.ls()
# wait for you to close all open canvases before exiting
# wait() will have no effect if ROOT is in batch mode:
#ROOT.gROOT.SetBatch(True)
raw_input('Press Enter to exit')