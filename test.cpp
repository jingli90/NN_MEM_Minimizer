#include "TROOT.h"
#include "TMath.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TInterpreter.h"
#include "TString.h"
#include "TLorentzVector.h"

#include "Math/Integrator.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "Math/GSLMCIntegrator.h"
#include "Math/IntegratorMultiDim.h"
#include "Math/IOptions.h"
#include "Math/IntegratorOptions.h"
#include "Math/AllIntegrationTypes.h"
#include "Math/AdaptiveIntegratorMultiDim.h"

#include "/afs/cern.ch/work/j/jing/public/IPHCNtuple/2017_04_19_80X-branch/MyAnalysis/CMSSW_8_0_9/src/ttH/MEM/src/MEPhaseSpace.h"
#include "/afs/cern.ch/work/j/jing/public/IPHCNtuple/2017_04_19_80X-branch/MyAnalysis/CMSSW_8_0_9/src/ttH/MEM/src/HypIntegrator.h"
#include "/afs/cern.ch/work/j/jing/public/IPHCNtuple/2017_04_19_80X-branch/MyAnalysis/CMSSW_8_0_9/src/ttH/MEM/src/MultiLepton.h"
#include "/afs/cern.ch/work/j/jing/public/IPHCNtuple/2017_04_19_80X-branch/MyAnalysis/CMSSW_8_0_9/src/ttH/MEM/src/ReadGenFlatTree.h"
#include "/afs/cern.ch/work/j/jing/public/IPHCNtuple/2017_04_19_80X-branch/MyAnalysis/CMSSW_8_0_9/src/ttH/MEM/src/ConfigParser.h"
#include "/afs/cern.ch/work/j/jing/public/IPHCNtuple/2017_04_19_80X-branch/MyAnalysis/CMSSW_8_0_9/src/ttH/MEM/src/Permutations.h"

#include "/afs/cern.ch/work/j/jing/public/IPHCNtuple/2017_04_19_80X-branch/MyAnalysis/CMSSW_8_0_9/src/ttH/MEM/src/TransferFunctions.h"

#include "NN_MEM.h"
#include "simpleNN.h"
#include "NN_MEM_Eval.h"

int main(int argc, char *argv[]){
	TString Filename_ttZ = "/afs/cern.ch/work/j/jing/public/IPHCNtuple/fromXavier/trees_from_Xavier_20170405/ttZJets_madgraphMLM.root";
	TString Filename_ttW = "/afs/cern.ch/work/j/jing/public/IPHCNtuple/fromXavier/trees_from_Xavier_20170405/ttWJets_madgraphMLM.root";
	TString Filename_ttH = "/afs/cern.ch/work/j/jing/public/IPHCNtuple/fromXavier/trees_from_Xavier_20170405/ttHToNonbb.root";
	TString treeName = "Tree";
	TString outFilename = "output";
	TString configFile="/afs/cern.ch/work/j/jing/private/FCPPL/learning_NN/simpleBP_MEM/NN_MEM_moreInput/NN_MEM_git/NN_MEM_Minimizer/test/config.cfg";

	NN_MEM m(Filename_ttZ, outFilename, "ttZ", configFile);
	m.SetNLayer(1);
	m.SetNNodes(1,10);
	m.SetNum_outputNN(5);
	//m.SetNEpochs(1000);
	m.SetNEpochs(-1);
	m.SetNeuronType("tanh");

	m.SetCut("is_3l_TTH_SR==1 && mc_ttZhypAllowed==1 && catJets==0");
	m.SetWeightExpression("weight");
	//m.AddVariable("nJet25_Recl");
	//m.AddVariable("max_Lep_eta");
	//m.AddVariable("MT_met_lep1");
	//m.AddVariable("mindr_lep1_jet");
	//m.AddVariable("mindr_lep2_jet");
	//m.AddVariable("LepGood_conePt0");
	//m.AddVariable("LepGood_conePt1");
	m.AddVariable("multilepton_Bjet1_P4->E()");
	m.AddVariable("multilepton_Bjet1_P4->Theta()");
	m.AddVariable("multilepton_Bjet1_P4->Phi()");
	m.AddVariable("multilepton_Bjet2_P4->E()");
	m.AddVariable("multilepton_Bjet2_P4->Theta()");
	m.AddVariable("multilepton_Bjet2_P4->Phi()");
	m.AddVariable("multilepton_JetClosestMw1_P4->E()");
	m.AddVariable("multilepton_JetClosestMw1_P4->Theta()");
	m.AddVariable("multilepton_JetClosestMw1_P4->Phi()");
	m.AddVariable("multilepton_JetClosestMw2_P4->E()");
	m.AddVariable("multilepton_JetClosestMw2_P4->Theta()");
	m.AddVariable("multilepton_JetClosestMw2_P4->Phi()");
	m.AddVariable("multilepton_Lepton1_P4->E()");
	m.AddVariable("multilepton_Lepton1_P4->Theta()");
	m.AddVariable("multilepton_Lepton1_P4->Phi()");
	m.AddVariable("multilepton_Lepton2_P4->E()");
	m.AddVariable("multilepton_Lepton2_P4->Theta()");
	m.AddVariable("multilepton_Lepton2_P4->Phi()");
	m.AddVariable("multilepton_Lepton3_P4->E()");
	m.AddVariable("multilepton_Lepton3_P4->Theta()");
	m.AddVariable("multilepton_Lepton3_P4->Phi()");
	m.AddVariable("multilepton_mET->Pt()");
	m.AddVariable("multilepton_mET->Phi()");
	m.AddSpectator("is_3l_TTH_SR");
	m.AddSpectator("mc_ttZhypAllowed");
	m.AddSpectator("catJets");

	vector<double>* xL;
	vector<double>* xU;
	xU = new vector<double>;
	xL = new vector<double>;
	xU->push_back(1500);
	xU->push_back(1.563771);
	xU->push_back(1500);
	xU->push_back(6.283185);
	xU->push_back(1.563771);
	xL->push_back(4.7);
	xL->push_back(-1.545340);
	xL->push_back(4.7);
	xL->push_back(0);
	xL->push_back(-1.545340);

	m.SetXLXU(xL, xU);

	vector<double>* var_max_int;
	vector<double>* var_min_int;
	var_max_int = new vector<double>;
	var_min_int = new vector<double>;
	var_min_int->push_back(25); // multilepton_Bjet1_P4->E()
	var_max_int->push_back(2603);
	var_min_int->push_back(0); // multilepton_Bjet1_P4->Theta()
	var_max_int->push_back(3);
	var_min_int->push_back(-4); // multilepton_Bjet1_P4->Phi()
	var_max_int->push_back(4);
	var_min_int->push_back(0); // multilepton_Bjet2_P4->E()
	var_max_int->push_back(2203);
	var_min_int->push_back(0); // multilepton_Bjet2_P4->Theta()
	var_max_int->push_back(3);
	var_min_int->push_back(-4); // multilepton_Bjet2_P4->Phi()
	var_max_int->push_back(4);
	var_min_int->push_back(0); // multilepton_JetClosestMw1_P4->E()
	var_max_int->push_back(3170);
	var_min_int->push_back(0); // multilepton_JetClosestMw1_P4->Theta()
	var_max_int->push_back(3);
	var_min_int->push_back(-4); // multilepton_JetClosestMw1_P4->Phi()
	var_max_int->push_back(4);
	var_min_int->push_back(0); // multilepton_JetClosestMw2_P4->E()
	var_max_int->push_back(2879);
	var_min_int->push_back(0); // multilepton_JetClosestMw2_P4->Theta()
	var_max_int->push_back(3);
	var_min_int->push_back(-4); // multilepton_JetClosestMw2_P4->Phi()
	var_max_int->push_back(4);
	var_min_int->push_back(21); // multilepton_Lepton1_P4->E()
	var_max_int->push_back(2437);
	var_min_int->push_back(0); // multilepton_Lepton1_P4->Theta()
	var_max_int->push_back(3);
	var_min_int->push_back(-4); // multilepton_Lepton1_P4->Phi()
	var_max_int->push_back(4);
	var_min_int->push_back(10); // multilepton_Lepton2_P4->E()
	var_max_int->push_back(1343);
	var_min_int->push_back(0); // multilepton_Lepton2_P4->Theta()
	var_max_int->push_back(3);
	var_min_int->push_back(-4); // multilepton_Lepton2_P4->Phi()
	var_max_int->push_back(4);
	var_min_int->push_back(0); // multilepton_Lepton3_P4->E()
	var_max_int->push_back(755);
	var_min_int->push_back(0); // multilepton_Lepton3_P4->Theta()
	var_max_int->push_back(3);
	var_min_int->push_back(-4); // multilepton_Lepton3_P4->Phi()
	var_max_int->push_back(4);
	var_min_int->push_back(0); // multilepton_mET->Pt()
	var_max_int->push_back(860);
	var_min_int->push_back(-4); // multilepton_mET->Phi() 
	var_max_int->push_back(4);
	m.SetMinMax(var_min_int, var_max_int);

	int kMode=kMEM_TTLL_TopAntitopDecay;
	int kCatJets=kCat_3l_2b_2j;
	m.SetModeCat(kMode, kCatJets);

	//m.SetMinimizationType("Minuit2Minimizer");
	m.SetMinimizationType("GSLSimAnMinimizer");

	m.SetFileInit("/afs/cern.ch/work/j/jing/private/FCPPL/learning_NN/simpleBP_MEM/NN_MEM_moreInput/NN_MEM_git/NN_MEM_Minimizer/test/weightInit.txt");

	m.AddOtherFiles(Filename_ttW, "ttW");
	m.AddOtherFiles(Filename_ttH, "ttH");

	m.myana();

}
