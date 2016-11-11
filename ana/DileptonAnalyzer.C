#include "TROOT.h"
#include "TSystem.h"
#include "TRint.h"
#include "TStyle.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TH2D.h"
#include "THStack.h"
#include "TTree.h"
#include "TLegend.h"
#include "TMath.h"
#include "TRandom3.h"

#include "MuonPogTree.h"
#include "CMSAnalysis.h"

using namespace muon_pog;

class WmunuAnalysis : public CMSAnalysis {
public:
  WmunuAnalysis(){};
  virtual ~WmunuAnalysis(){};
  
  // Preselection
  bool Preselect();
  
  // Trigger filter
  bool isTriggered();

  //leptonSF
  void Apply_lepton_SF(Muon* mu, double& weight, bool notrig=false);
  void Muon_triggerSF(Muon* mu, double& weight);
  void Apply_dimuon_triggerSF(Muon* mu1, Muon* mu2, double& weight);

  //Checks all the input parameters
  Double_t deltaPhi(Double_t Phi0, Double_t Phi1);

  // Selection of a Wprime event: NULL if rejected, link to muon if accepted
  std::vector<muon_pog::Muon*> GoodZprimeMuons();

private:
  
};

//Input Arguments
TString mode = "DUMMY";
TString InjectedScaleNAME = "DUMMY";
float InjectedScale = 0.;
TString assignment = "DUMMY";
TString folder = "DUMMY";

///Other stuff
FILE *Monitor_file;

int main(int argc, char** argv){

  mode = argv[1];
  InjectedScaleNAME = argv[2];
  InjectedScale = float(InjectedScaleNAME.Atof());
  assignment = argv[3];
  folder = argv[4];

  if (mode == "Scale" || mode == "Validation"){
    printf("Selected Mode: %s \n", mode.Data());
  }else{
    printf("USE a valid Mode: ./cmsTreeRun DileptonAnalyzer.C [Mode] [Bias] [Assignment] [OutputFolder] \n");
    return 0;
  }
  if (assignment == "TUNEP" || assignment == "INNER" || assignment == "PICKY_HIGHPT" || assignment == "DYT" || assignment == "PICKY" || assignment == "TPFMS" || assignment == "GLOBAL"){
    printf("Selected assignment: %s \n", assignment.Data());
    printf("Selected BIAS: %s \n", InjectedScaleNAME.Data());
  }else{
    printf("USE a valid Assignment: ./cmsTreeRun DileptonAnalyzer.C [Mode] [Bias] [Assignment] [OutputFolder] \n");
    return 0; 
  }
  if (folder == "DUMMY"){
    printf("Select a folder, USE: ./cmsTreeRun DileptonAnalyzer.C [Mode] [Bias] [Assignment] [OutputFolder] \n");
    return 0;
  }else{
    printf("Selected histogram folder : %s \n", folder.Data());
  }


  // Initialize analysis structure
  WmunuAnalysis cms;

  // FAST TESTS => process "1/stepData" fraction of real data events; 
  //               set stepData=1 for final running
  unsigned int stepData = 1; //Do not change.
  
  // Add data sample. Input is: (titleId, file, luminosity)
  // One should add this sample first, to define the luminosity for normalizations
  TString lumi_invfb = ""; 
  double lumi_invpb = 0.; 
  
  if (mode == "Validation" || mode == "Scale"){
    lumi_invfb = "1.0"; // golden json preliminary (2.07+10%)
    lumi_invpb = lumi_invfb.Atof()*1000./stepData; 
//    cms.AddDataSample("Data",  "/afs/cern.ch/work/r/rradogna/Scale/CMSSW_8_0_8/src/GeneralizedEndpoint_80/ntuples/skim_mumu_50_15_CMSTrees_SingleMuon_2016_Prompt_v2_MINIAOD.root", lumi_invpb);
      //cms.AddDataSample("Data",  "root://cms-xrd-global.cern.ch///store/group/phys_muon/Commissioning/Ntuples/Commissioning2016/data25ns/skimmed/output_data10June/muonPOGNtuple_SingleMu2016B-v2-Golden-10-June_skimmed_2_muons_20_GeV_nohits_and_matches_no_HLT_2.root", lumi_invpb);
      printf("reading data");
      cms.AddDataSample("Data",  "/afs/cern.ch/work/r/rradogna/Scale/ScalePOG/CMSSW_8_0_8/src/ntuples/muonPOGNtuple_8_0_8_ZToMuMu_powheg_M_1400_2300_13.root", lumi_invpb);
  }

  TString chtit_data = TString::Format("2016 Data, L=%d fb^{-1}", int(lumi_invpb/1.e3+0.5));
   
  int prescale = -1; // use all MC events, final plots DO NOT CHANGE.
  prescale = 1; // use all MC events, final plots DO NOT CHANGE.

  // MC Cross sections
  double dy_k_factor_nnlo=1.006;
  double z_50_120_powheg_xsec = 1975.*dy_k_factor_nnlo;
  double z_120_200_powheg_xsec = 19.32*dy_k_factor_nnlo;
  double z_200_400_powheg_xsec = 2.731*dy_k_factor_nnlo;
  double z_400_800_powheg_xsec = 0.241*dy_k_factor_nnlo;
  double z_800_1400_powheg_xsec = 0.01678*dy_k_factor_nnlo;
  double z_1400_2300_powheg_xsec = 0.00139*dy_k_factor_nnlo;
  double z_2300_3500_powheg_xsec = 8.948e-05*dy_k_factor_nnlo;
  double z_3500_4500_powheg_xsec = 4.135e-06*dy_k_factor_nnlo;
  //
  double z_pt_150_madgraph_xsec =  18.39*1.3;//9.53 18.39 pb

  //
  double w_xsec = 20508.9*3; //Factor 3 due possible LNu decays. NNLO
  //

  double tt_xsec = 831.76;//730.0;//831.76; // Inclusive top cross section, NNLO+NNLL, Filter Eff for 700 to 1000 wrt inclusive = 0.0921, 1000 to inf... 0.02474
  double singlet_xsec = 35.6; // 
  double singletbar_xsec = 35.6; // 
  double tt_700_1000_xsec = 831.76*0.0921; // Inclusive top cross section, NNLO+NNLL.
  double tt_1000_xsec = 831.76*0.02474; // Inclusive top cross section, NNLO+NNLL.
  //

    double VV_k_factor_nlo=1;//1.3;
  double ww_xsec = 118.7*VV_k_factor_nlo; //NNLO QCD
    double wz_xsec = 47.13;//(40.2+25.)*VV_k_factor_nlo; //NLO
  double zz_xsec = 16.523*VV_k_factor_nlo; //NLO

  // Initialize the yields
  double w_in_total = 24184766.;   double w_in_total_signed = 16541248.;//25ns
  
  double z_50_120_in_powheg_total = 2967200;
  double z_120_200_in_powheg_total = 99200;
  double z_200_400_in_powheg_total = 100000;
  double z_400_800_in_powheg_total = 100000;
  double z_800_1400_in_powheg_total = 100000;
  double z_1400_2300_in_powheg_total = 100000;
  double z_2300_3500_in_powheg_total = 99200;
  double z_3500_4500_in_powheg_total = 95200;

  double z_pt_150_in_madgraph_total =  2355383;//2355383

  double tt_in_total = 182123200;//25ns
  double singlet_in_total = 998400; //
  double singletbar_in_total = 985000; //
  double tt_700_1000_in_total = 3868007;//25ns
  double tt_1000_in_total = 2360497;//25ns

  double ww_in_total = 993214;//25ns
  double wz_in_total = 1000000;//25ns
  double zz_in_total = 989312;//25ns

//  TString dir = "/afs/cern.ch/work/r/rradogna/Scale/CMSSW_8_0_8/src/GeneralizedEndpoint_80/ntuples/";//25ns
  TString dir = "/afs/cern.ch/work/r/rradogna/Scale/CMSSW_8_0_8/src/GeneralizedEndpoint_80/ntuples/SKIM_Zprime/";//25ns
  TString SkimSuffix = "_ZPRIME_dimuon_skimmed"; // 25ns

///cms.AddMCSample("SM_Z", dir+"/test/ZToMuMu_NNPDF30_13TeV-powheg_M_50_120"+SkimSuffix+"_WPRIME_skimmed_PtZ_0_250.root", z_50_120_in_powheg_total, z_50_120_powheg_xsec, z_50_120_in_powheg_total, prescale);
///  cms.AddMCSample("SM_Z", dir+"/test/DYJetsToLL_M-50_Zpt-150toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8"+SkimSuffix+"_WPRIME_skimmed_PtZ_250_3000.root", z_pt_150_in_madgraph_total, z_pt_150_madgraph_xsec, z_pt_150_in_madgraph_total, prescale);
//  cms.AddMCSample("SM_Z", "/afs/cern.ch/work/r/rradogna/Scale/CMSSW_8_0_8/src/GeneralizedEndpoint_80/ntuples/CMSTrees_DYJetsToLL_M-50_Zpt-150toInf_TuneCUETP8M1_13TeV_80X_asymptotic_MINIAODv2_WPRIME_skimmed_PtZ_150_Inf.root", z_pt_150_in_madgraph_total, z_pt_150_madgraph_xsec, z_pt_150_in_madgraph_total, prescale);
//    
//  cms.AddMCSample("SM_Z", dir+"ZToMuMu_NNPDF30_13TeV-powheg_M_50_120"+SkimSuffix+"_WPRIME_skimmed_PtZ_0_150.root", z_50_120_in_powheg_total, z_50_120_powheg_xsec, z_50_120_in_powheg_total, prescale);
//  cms.AddMCSample("SM_Z", dir+"ZToMuMu_NNPDF30_13TeV-powheg_M_120_200"+SkimSuffix+"_WPRIME_skimmed_PtZ_0_150.root", z_120_200_in_powheg_total, z_120_200_powheg_xsec, z_120_200_in_powheg_total, prescale);
//  cms.AddMCSample("SM_Z", dir+"ZToMuMu_NNPDF30_13TeV-powheg_M_200_400"+SkimSuffix+"_WPRIME_skimmed_PtZ_0_150.root", z_200_400_in_powheg_total, z_200_400_powheg_xsec, z_200_400_in_powheg_total, prescale);
//  cms.AddMCSample("SM_Z", dir+"ZToMuMu_NNPDF30_13TeV-powheg_M_400_800"+SkimSuffix+"_WPRIME_skimmed_PtZ_0_150.root", z_400_800_in_powheg_total, z_400_800_powheg_xsec, z_400_800_in_powheg_total, prescale);
//  cms.AddMCSample("SM_Z", dir+"ZToMuMu_NNPDF30_13TeV-powheg_M_800_1400"+SkimSuffix+".root", z_800_1400_in_powheg_total, z_800_1400_powheg_xsec, z_800_1400_in_powheg_total, prescale);
//  cms.AddMCSample("SM_Z", dir+"ZToMuMu_NNPDF30_13TeV-powheg_M_1400_2300"+SkimSuffix+".root", z_1400_2300_in_powheg_total, z_1400_2300_powheg_xsec, z_1400_2300_in_powheg_total, prescale);
//  cms.AddMCSample("SM_Z", dir+"ZToMuMu_NNPDF30_13TeV-powheg_M_2300_3500"+SkimSuffix+".root", z_2300_3500_in_powheg_total, z_2300_3500_powheg_xsec, z_2300_3500_in_powheg_total, prescale);
//  cms.AddMCSample("SM_Z", dir+"ZToMuMu_NNPDF30_13TeV-powheg_M_3500_4500"+SkimSuffix+".root", z_3500_4500_in_powheg_total, z_3500_4500_powheg_xsec, z_3500_4500_in_powheg_total, prescale);
//
////  cms.AddMCSample("SM_W", dir+"WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8"+SkimSuffix+".root", w_in_total, w_xsec, w_in_total_signed, prescale);
//  cms.AddMCSample("t#bar{t}", "/afs/cern.ch/work/r/rradogna/Scale/CMSSW_8_0_8/src/GeneralizedEndpoint_80/ntuples/skim_mumu_50_15__CMSTrees_TT_TuneCUETP8M1_13TeV-powheg-pythia8_80X_asymptotic_MINIAODv2.root", tt_in_total, tt_xsec, tt_in_total, prescale);
//  cms.AddMCSample("t#bar{t}", dir+"CMSTrees_ST_tW_top_5f_13TeV-powheg-pythia8_TuneCUETP8M1_80X_asymptotic_MINIAODv2"+SkimSuffix+".root", singlet_in_total, singlet_xsec, singlet_in_total, prescale);
//  cms.AddMCSample("t#bar{t}", dir+"CMSTrees_ST_tW_antitop_5f_13TeV-powheg-pythia8_TuneCUETP8M1_80X_asymptotic_MINIAODv2"+SkimSuffix+".root", singletbar_in_total, singletbar_xsec, singletbar_in_total, prescale);
//
//  cms.AddMCSample("VV", dir+"CMSTrees_WW_TuneCUETP8M1_13TeV-pythia8_80X_asymptotic_MINIAODv2"+SkimSuffix+".root", ww_in_total, ww_xsec, ww_in_total, prescale);
//  cms.AddMCSample("VV", dir+"CMSTrees_WZ_TuneCUETP8M1_13TeV-pythia8_80X_asymptotic_MINIAODv2"+SkimSuffix+".root", wz_in_total, wz_xsec, wz_in_total, prescale);
//  cms.AddMCSample("VV", dir+"CMSTrees_ZZ_TuneCUETP8M1_13TeV-pythia8_80X_asymptotic_MINIAODv2"+SkimSuffix+".root", zz_in_total, zz_xsec, zz_in_total, prescale);

//  cms.AddMCSample("VV", "root://cms-xrd-global.cern.ch///store/caf/user/oglez/NTUPLES/80X_CMSTrees_2016_06_20/CMSTrees_ZZ_TuneCUETP8M1_13TeV-pythia8_80X_asymptotic_MINIAODv2.root", zz_in_total, zz_xsec, zz_in_total, prescale);
  cms.AddMCSample("SM_Z", "/afs/cern.ch/work/r/rradogna/Scale/ScalePOG/CMSSW_8_0_8/src/ntuples/muonPOGNtuple_8_0_8_ZToMuMu_powheg_M_1400_2300_13.root", 1000, z_1400_2300_powheg_xsec, 1000, prescale);


  unsigned int REBIN = 1; // Use REBIN>1 to get less bins in most histograms

  TString lepton_name_for_axis = "#mu";
  TString lepton_id = "HighPt";
  float iso_lepton=50.;
  float pt_lepton=25.;
  float z_peak_low=65.;
  float z_peak_high=120.;
  float eta_lepton=2.5;
  
  if (mode=="Validation"){
    //Pt
    cms.AddPlot1D("hPtmuPeak_mu1",  lepton_name_for_axis+" p_{T} [GeV]", 50/REBIN, 53., 300.);
    cms.AddPlot1D("hPtmuPeak_mu2",  lepton_name_for_axis+" p_{T} [GeV]", 50/REBIN, 25., 300.);
    cms.AddPlot1D("hPtmuTail",  lepton_name_for_axis+" p_{T} [GeV]", 50/REBIN, pt_lepton, 1000.);

    //Boson Pt
    cms.AddPlot1D("hPtZPeak", "p_{T} dimuom pair [GeV]", 50./REBIN, 50, 200.);
    cms.AddPlot1D("hPtZTail", "p_{T} dimuom pair [GeV]", 50./REBIN, 50, 1500.);

    //Other
    cms.AddPlot1D("hPhimu", lepton_name_for_axis+" #phi", 50/REBIN, -M_PI, M_PI);
    cms.AddPlot1D("hEtamu", lepton_name_for_axis+" #eta", 50/REBIN, -1*eta_lepton, eta_lepton);
    cms.AddPlot1D("hDxymu", lepton_name_for_axis+" d_{xy} [cm]", 200/REBIN, -0.02, 0.02);
    cms.AddPlot1D("hDzmu",  lepton_name_for_axis+" d_{z} [cm]", 200/REBIN, -0.05, 0.05);  
    cms.AddPlot1D("hTrkIso", lepton_name_for_axis+" isolation", 50/REBIN, 0.0, iso_lepton);
    cms.AddPlot1D("hTrkIsoOverPt", lepton_name_for_axis+" isolation over pt", 50/REBIN, 0.0, 0.1);

    //Negative
    cms.AddPlot1D("hPtmu_N",  lepton_name_for_axis+" p_{T} [GeV], q=-1", 50/REBIN, pt_lepton, 750.);
    cms.AddPlot1D("hPtmu_NBarrel",  lepton_name_for_axis+" p_{T} [GeV], q=-1, |#eta|<1.2", 50/REBIN, pt_lepton, 750.);
    cms.AddPlot1D("hPtmu_NEndcap",  lepton_name_for_axis+" p_{T} [GeV], q=-1, 1.2>|#eta|>2.1", 50/REBIN, pt_lepton, 750.);
    cms.AddPlot1D("hPtmu_NVeryEndcap",  lepton_name_for_axis+" p_{T} [GeV], q=-1, |#eta|>2.1", 50/REBIN, pt_lepton, 750.);
    cms.AddPlot1D("hPtmu_NLargeTail",  lepton_name_for_axis+" p_{T} [GeV], q=-1", 50/REBIN, pt_lepton, 1500.);
    cms.AddPlot1D("hPtmu_NLargeTailBarrel",  lepton_name_for_axis+" p_{T} [GeV], q=-1, |#eta|<1.2", 50/REBIN, pt_lepton, 1500.);
    cms.AddPlot1D("hPtmu_NLargeTailEndcap",  lepton_name_for_axis+" p_{T} [GeV], q=-1, 1.2>|#eta|>2.1", 50/REBIN, pt_lepton, 1500.);
    cms.AddPlot1D("hPtmu_NLargeTailForward",  lepton_name_for_axis+" p_{T} [GeV], q=-1, |#eta|>2.1", 50/REBIN, pt_lepton, 1500.);
    cms.AddPlot1D("hPhimu_N", lepton_name_for_axis+" #phi, q=-1", 50/REBIN, -M_PI, M_PI);
    cms.AddPlot1D("hEtamu_N", lepton_name_for_axis+" #eta, q=-1", 50/REBIN, -1*eta_lepton, eta_lepton);

    //Positive
    cms.AddPlot1D("hPtmu_P",  lepton_name_for_axis+" p_{T} [GeV], q=+1", 50/REBIN, pt_lepton, 750.);
    cms.AddPlot1D("hPtmu_PBarrel",  lepton_name_for_axis+" p_{T} [GeV], q=+1, |#eta|<1.2", 50/REBIN, pt_lepton, 750.);
    cms.AddPlot1D("hPtmu_PEndcap",  lepton_name_for_axis+" p_{T} [GeV], q=+1, 1.2>|#eta|>2.1", 50/REBIN, pt_lepton, 750.);
    cms.AddPlot1D("hPtmu_PVeryEndcap",  lepton_name_for_axis+" p_{T} [GeV], q=+1, |#eta|>2.1", 50/REBIN, pt_lepton, 750.);
    cms.AddPlot1D("hPtmu_PLargeTail",  lepton_name_for_axis+" p_{T} [GeV], q=+1", 50/REBIN, pt_lepton, 1500.);
    cms.AddPlot1D("hPtmu_PLargeTailBarrel",  lepton_name_for_axis+" p_{T} [GeV], q=+1, |#eta|<1.2", 50/REBIN, pt_lepton, 1500.);
    cms.AddPlot1D("hPtmu_PLargeTailEndcap",  lepton_name_for_axis+" p_{T} [GeV], q=+1, 1.2>|#eta|>2.1", 50/REBIN, pt_lepton, 1500.);
    cms.AddPlot1D("hPtmu_PLargeTailForward",  lepton_name_for_axis+" p_{T} [GeV], q=+1, |#eta|>2.1", 50/REBIN, pt_lepton, 1500.);
    cms.AddPlot1D("hKMuon",  "q/p_{T} [c/TeV] ", 40, -5, 5.);
    cms.AddPlot1D("hKMuonLowMass",  "q/p_{T} [c/TeV] ", 80, -10, 10.);
    cms.AddPlot1D("hKMuonLowMassBarrel",  "q/p_{T} [c/TeV] ", 80, -10, 10.);
    cms.AddPlot1D("hKMuonPeak", "p_{T} [c/TeV] ", 160, -20, 20.);

    cms.AddPlot1D("hPhimu_P", lepton_name_for_axis+" #phi, q=+1", 50/REBIN, -M_PI, M_PI);
    cms.AddPlot1D("hEtamu_P", lepton_name_for_axis+" #eta, q=+1", 50/REBIN, -1*eta_lepton, eta_lepton);

    //Inv masses
    cms.AddPlot1D("hInvMassPeak", "M_{#mu #mu} [GeV]", 50./REBIN, z_peak_low, z_peak_high);
    cms.AddPlot1D("hInvMassPeakEndcap", "M_{#mu #mu} [GeV]", 50./REBIN, z_peak_low, z_peak_high);
    cms.AddPlot1D("hInvMassPeakVeryEndcap", "M_{#mu #mu} [GeV]", 50./REBIN, z_peak_low, z_peak_high);
    cms.AddPlot1D("hInvMassTail", "M_{#mu #mu} [GeV]", 50./REBIN, z_peak_low, 1500.);
    cms.AddPlot1D("hInvMassTailv2", "M_{#mu #mu} [GeV]", 110./REBIN, z_peak_low, 1500.);
    cms.AddPlot1D("hInvMassLargeTail", "M_{#mu #mu} [GeV]", 50./REBIN, z_peak_low, 2500.);
  }


  bool scale=false;
  bool ScaleMap = false;
  if (mode == "Scale") {
    scale = true;
    ScaleMap = true;
  }
  
  // Scan range to be synchornized with the production.
  double IniBias = -1.0; 
  double EndBias = 1.0; 
  double StepBias = 0.01; 
  int nbin = 0.;  
  float GoodBias = 0.;  

  if (scale == true){
    for (float kbias =IniBias; kbias<=EndBias; kbias = kbias + StepBias){    
      if (fabs(kbias-InjectedScale) > 0.00001) {
	//	printf("nbin %d kbias %f \n", nbin, kbias-InjectedScale);
	nbin++;
      }else{	
	GoodBias = kbias;       	
	printf("Identifier %d InjectedBias %f \n", nbin, GoodBias);
	cms.AddPlot1D("hKMuon_"+TString(Form("%d",nbin)),  "q/p_{T} [c/TeV] if #kappa_{b}="+TString(Form("%.2f",GoodBias)), 40, -5, 5.);
	cms.AddPlot1D("hKMuonBarrel_"+TString(Form("%d",nbin)),  "q/p_{T} [c/TeV] , |#eta|<1.2, if #kappa_{b}="+TString(Form("%.2f",GoodBias)), 40, -5., 5.);
	cms.AddPlot1D("hKMuonEndcap_"+TString(Form("%d",nbin)),  "q/p_{T} [c/TeV] |#eta|>1.2, if #kappa_{b}="+TString(Form("%.2f",GoodBias)), 40, -5., 5.);
	cms.AddPlot1D("hKMuonPosEndcap_"+TString(Form("%d",nbin)),  "q/p_{T} [c/TeV] #eta>1.2, if #kappa_{b}="+TString(Form("%.2f",GoodBias)), 40, -5., 5.);
	cms.AddPlot1D("hKMuonNegEndcap_"+TString(Form("%d",nbin)),  "q/p_{T} [c/TeV] #eta<-1.2, if #kappa_{b}="+TString(Form("%.2f",GoodBias)), 40, -5., 5.);
	cms.AddPlot1D("hKMuonPeak_"+TString(Form("%d",nbin)),  "q/p_{T} [c/TeV] in ZPeak if #kappa_{b}="+TString(Form("%.2f",GoodBias)), 160, -20, 20.);
	cms.AddPlot1D("hKMuonLowMass_"+TString(Form("%d",nbin)),  "q/p_{T} [c/TeV] in ZPeak if #kappa_{b}="+TString(Form("%.2f",GoodBias)), 160, -10, 10.);
    cms.AddPlot1D("hKMuonLowMassBarrel_"+TString(Form("%d",nbin)),  "q/p_{T} [c/TeV] in ZPeak if #kappa_{b}="+TString(Form("%.2f",GoodBias)), 160, -10, 10.);
	cms.AddPlot1D("hKMuonLowMass_P_"+TString(Form("%d",nbin)),  "q/p_{T} [c/TeV], q=+1, if #kappa_{b}="+TString(Form("%.2f",GoodBias)), 43, 0.70, 10.);
	cms.AddPlot1D("hKMuonLowMass_N_"+TString(Form("%d",nbin)),  "q/p_{T} [c/TeV], q=-1, if #kappa_{b}="+TString(Form("%.2f",GoodBias)), 43, 0.70, 10.);
	if (ScaleMap == true){
	  float Range0 = 9.; //This values settle the curvature cut for the MAP eta >2.1
	  float Range1 = 9.; //This values settle the curvature cut for the MAP 1.2 > eta > 2.1
	  float Range2 = 9.; //This values settle the curvature cut for the MAP eta < 1.2
	  if (assignment == "PICKY_HIGHPT"){
	    Range0 = 9.; 
	    Range1 = 6.6;
	    Range2 = 5.5;	    
	  }
//	  // -2.4,-2.1 	  
	  cms.AddPlot1D("hKMuon_0_1_"+TString(Form("%d",nbin)),  "q/p_{T} [c/TeV] -2.4<#eta<-2.1, -180< #phi <-60  if #kappa_{b}="+TString(Form("%.2f",GoodBias)), 160, -1*Range0, Range0);
	  cms.AddPlot1D("hKMuon_0_2_"+TString(Form("%d",nbin)),  "q/p_{T} [c/TeV] -2.4<#eta<-2.1, -60< #phi <60   if #kappa_{b}="+TString(Form("%.2f",GoodBias)), 160, -1*Range0, Range0);
	  cms.AddPlot1D("hKMuon_0_3_"+TString(Form("%d",nbin)),  "q/p_{T} [c/TeV] -2.4<#eta<-2.1,  60< #phi <180  if #kappa_{b}="+TString(Form("%.2f",GoodBias)), 160, -1*Range0, Range0);
//	  // -2.1,-1.2
//	  cms.AddPlot1D("hKMuon_1_1_"+TString(Form("%d",nbin)),  "q/p_{T} [c/TeV] -2.1<#eta<-1.2, -180< #phi <-60  if #kappa_{b}="+TString(Form("%.2f",GoodBias)), 160, -1*Range1, Range1);
//	  cms.AddPlot1D("hKMuon_1_2_"+TString(Form("%d",nbin)),  "q/p_{T} [c/TeV] -2.1<#eta<-1.2, -60< #phi <60   if #kappa_{b}="+TString(Form("%.2f",GoodBias)), 160, -1*Range1, Range1);
//	  cms.AddPlot1D("hKMuon_1_3_"+TString(Form("%d",nbin)),  "q/p_{T} [c/TeV] -2.1<#eta<-1.2,  60< #phi <180  if #kappa_{b}="+TString(Form("%.2f",GoodBias)), 160, -1*Range1, Range1);
//	  // -1.2,0 
//	  cms.AddPlot1D("hKMuon_2_1_"+TString(Form("%d",nbin)),  "q/p_{T} [c/TeV] -1.2<#eta<0, -180< #phi <-60  if #kappa_{b}="+TString(Form("%.2f",GoodBias)), 160, -1*Range2, Range2);
//	  cms.AddPlot1D("hKMuon_2_2_"+TString(Form("%d",nbin)),  "q/p_{T} [c/TeV] -1.2<#eta<0, -60< #phi <60   if #kappa_{b}="+TString(Form("%.2f",GoodBias)), 160, -1*Range2, Range2);
//	  cms.AddPlot1D("hKMuon_2_3_"+TString(Form("%d",nbin)),  "q/p_{T} [c/TeV] -1.2<#eta<0,  60< #phi <180  if #kappa_{b}="+TString(Form("%.2f",GoodBias)), 160, -1*Range2, Range2);
//	  // 0,1.2 
//	  cms.AddPlot1D("hKMuon_3_1_"+TString(Form("%d",nbin)),  "q/p_{T} [c/TeV] 0<#eta<1.2, -180< #phi <-60  if #kappa_{b}="+TString(Form("%.2f",GoodBias)), 160, -1*Range2, Range2);
//	  cms.AddPlot1D("hKMuon_3_2_"+TString(Form("%d",nbin)),  "q/p_{T} [c/TeV] 0<#eta<1.2, -60< #phi <60   if #kappa_{b}="+TString(Form("%.2f",GoodBias)), 160, -1*Range2, Range2);
//	  cms.AddPlot1D("hKMuon_3_3_"+TString(Form("%d",nbin)),  "q/p_{T} [c/TeV] 0<#eta<1.2,  60< #phi <180  if #kappa_{b}="+TString(Form("%.2f",GoodBias)), 160, -1*Range2, Range2);
//	  // 1.2,2.1 
//	  cms.AddPlot1D("hKMuon_4_1_"+TString(Form("%d",nbin)),  "q/p_{T} [c/TeV] 1.2<#eta<2.1, -180< #phi <-60  if #kappa_{b}="+TString(Form("%.2f",GoodBias)), 160, -1*Range1, Range1);
//	  cms.AddPlot1D("hKMuon_4_2_"+TString(Form("%d",nbin)),  "q/p_{T} [c/TeV] 1.2<#eta<2.1, -60< #phi <60   if #kappa_{b}="+TString(Form("%.2f",GoodBias)), 160, -1*Range1, Range1);
//	  cms.AddPlot1D("hKMuon_4_3_"+TString(Form("%d",nbin)),  "q/p_{T} [c/TeV] 1.2<#eta<2.1,  60< #phi <180  if #kappa_{b}="+TString(Form("%.2f",GoodBias)), 160, -1*Range1, Range1);
//	  // 2.1,2.4 
	  cms.AddPlot1D("hKMuon_5_1_"+TString(Form("%d",nbin)),  "q/p_{T} [c/TeV] 2.1<#eta<2.4, -180< #phi <-60  if #kappa_{b}="+TString(Form("%.2f",GoodBias)), 160, -1*Range0, Range0);
	  cms.AddPlot1D("hKMuon_5_2_"+TString(Form("%d",nbin)),  "q/p_{T} [c/TeV] 2.1<#eta<2.4, -60< #phi <60   if #kappa_{b}="+TString(Form("%.2f",GoodBias)), 160, -1*Range0, Range0);
	  cms.AddPlot1D("hKMuon_5_3_"+TString(Form("%d",nbin)),  "q/p_{T} [c/TeV] 2.1<#eta<2.4,  60< #phi <180  if #kappa_{b}="+TString(Form("%.2f",GoodBias)), 160, -1*Range0, Range0);
        // -2.4,-2.1
//        cms.AddPlot1D("hKMuon_0_1_"+TString(Form("%d",nbin)),  "q/p_{T} [c/TeV] -2.4<#eta<-2.1, -180< #phi <-60  if #kappa_{b}="+TString(Form("%.2f",GoodBias)), 160, -5,5);
//        cms.AddPlot1D("hKMuon_0_2_"+TString(Form("%d",nbin)),  "q/p_{T} [c/TeV] -2.4<#eta<-2.1, -60< #phi <60   if #kappa_{b}="+TString(Form("%.2f",GoodBias)), 160, -5,5);
//        cms.AddPlot1D("hKMuon_0_3_"+TString(Form("%d",nbin)),  "q/p_{T} [c/TeV] -2.4<#eta<-2.1,  60< #phi <180  if #kappa_{b}="+TString(Form("%.2f",GoodBias)), 160, -5,5);
        // -2.1,-1.2
        cms.AddPlot1D("hKMuon_1_1_"+TString(Form("%d",nbin)),  "q/p_{T} [c/TeV] -2.1<#eta<-1.2, -180< #phi <-60  if #kappa_{b}="+TString(Form("%.2f",GoodBias)), 160, -5,5);
        cms.AddPlot1D("hKMuon_1_2_"+TString(Form("%d",nbin)),  "q/p_{T} [c/TeV] -2.1<#eta<-1.2, -60< #phi <60   if #kappa_{b}="+TString(Form("%.2f",GoodBias)), 160, -5,5);
        cms.AddPlot1D("hKMuon_1_3_"+TString(Form("%d",nbin)),  "q/p_{T} [c/TeV] -2.1<#eta<-1.2,  60< #phi <180  if #kappa_{b}="+TString(Form("%.2f",GoodBias)), 160, -5,5);
        // -1.2,0
        cms.AddPlot1D("hKMuon_2_1_"+TString(Form("%d",nbin)),  "q/p_{T} [c/TeV] -1.2<#eta<0, -180< #phi <-60  if #kappa_{b}="+TString(Form("%.2f",GoodBias)), 160, -5,5);
        cms.AddPlot1D("hKMuon_2_2_"+TString(Form("%d",nbin)),  "q/p_{T} [c/TeV] -1.2<#eta<0, -60< #phi <60   if #kappa_{b}="+TString(Form("%.2f",GoodBias)), 160, -5,5);
        cms.AddPlot1D("hKMuon_2_3_"+TString(Form("%d",nbin)),  "q/p_{T} [c/TeV] -1.2<#eta<0,  60< #phi <180  if #kappa_{b}="+TString(Form("%.2f",GoodBias)), 160, -5,5);
        // 0,1.2
        cms.AddPlot1D("hKMuon_3_1_"+TString(Form("%d",nbin)),  "q/p_{T} [c/TeV] 0<#eta<1.2, -180< #phi <-60  if #kappa_{b}="+TString(Form("%.2f",GoodBias)), 160, -5,5);
        cms.AddPlot1D("hKMuon_3_2_"+TString(Form("%d",nbin)),  "q/p_{T} [c/TeV] 0<#eta<1.2, -60< #phi <60   if #kappa_{b}="+TString(Form("%.2f",GoodBias)), 160, -5,5);
        cms.AddPlot1D("hKMuon_3_3_"+TString(Form("%d",nbin)),  "q/p_{T} [c/TeV] 0<#eta<1.2,  60< #phi <180  if #kappa_{b}="+TString(Form("%.2f",GoodBias)), 160, -5,5);
        // 1.2,2.1
        cms.AddPlot1D("hKMuon_4_1_"+TString(Form("%d",nbin)),  "q/p_{T} [c/TeV] 1.2<#eta<2.1, -180< #phi <-60  if #kappa_{b}="+TString(Form("%.2f",GoodBias)), 160, -5,5);
        cms.AddPlot1D("hKMuon_4_2_"+TString(Form("%d",nbin)),  "q/p_{T} [c/TeV] 1.2<#eta<2.1, -60< #phi <60   if #kappa_{b}="+TString(Form("%.2f",GoodBias)), 160, -5,5);
        cms.AddPlot1D("hKMuon_4_3_"+TString(Form("%d",nbin)),  "q/p_{T} [c/TeV] 1.2<#eta<2.1,  60< #phi <180  if #kappa_{b}="+TString(Form("%.2f",GoodBias)), 160, -5,5);
        // 2.1,2.4
//        cms.AddPlot1D("hKMuon_5_1_"+TString(Form("%d",nbin)),  "q/p_{T} [c/TeV] 2.1<#eta<2.4, -180< #phi <-60  if #kappa_{b}="+TString(Form("%.2f",GoodBias)), 160, -5,5);
//        cms.AddPlot1D("hKMuon_5_2_"+TString(Form("%d",nbin)),  "q/p_{T} [c/TeV] 2.1<#eta<2.4, -60< #phi <60   if #kappa_{b}="+TString(Form("%.2f",GoodBias)), 160, -5,5);
//        cms.AddPlot1D("hKMuon_5_3_"+TString(Form("%d",nbin)),  "q/p_{T} [c/TeV] 2.1<#eta<2.4,  60< #phi <180  if #kappa_{b}="+TString(Form("%.2f",GoodBias)), 160, -5,5);
	  break;
	}
	//	nbin++;
      }
    }
    //    nbin=0;
  }
    
  // Loop on samples
  int nsamples = cms.GetNumberOfSamples();
  //Monitor Intersting events. 
  Monitor_file = fopen("GeneralizedEndpointEvents_ZMuMu.txt","w+");

  for (unsigned int iSample=0; iSample<nsamples; ++iSample) {
    if (iSample>0) printf("iSample %i Equivalent Luminosity %f fb, NEvents %i, Xsec %f, Weight %f \n", 
			  iSample,  cms.GetSampleEquivalentLuminosity(iSample)/1000., cms.GetNumberOfEvents(iSample), cms.GetSampleXsection(iSample), cms.GetSampleWeight(iSample));    
    
    // Set tree
    cms.SetTree(iSample);

    // Loop on events for the current sample
    unsigned int step = 1;     int nevents = cms.GetNumberOfEvents(iSample);
    if (cms.isDataFile()) step = stepData;
    //if (!cms.isModelFile() ) continue; // if SignalMC only running over signal MC samples
    
    for (int iEvent=0; iEvent<nevents; iEvent+=step) {
      // Get next entry
      if (cms.GetEntry(iEvent)<0) break;
        if (cms.isDataFile()){
            if (!cms.isTriggered()) continue;
                printf("data Triggered\n");
        }

      // Fast preselection using the summary information
      //if (!cms.Preselect()) continue;

      // Fill histograms
      double weight = 1.;
      // This is needed only for aMC@NLO
//      if (!cms.isDataFile() && cms.genInfo().PDFweight < 0) {
//	weight = -1*weight;
//      }//TO BE FIXED
      
      //SF for Id
      // 2 good muons, with the right Z mass                                            
      std::vector<muon_pog::Muon*> goodMuons = cms.GoodZprimeMuons();
        printf("data good Zprime\n");

      if (goodMuons.size()!=2) continue;
      muon_pog::Muon* mu1 = goodMuons[0];
      muon_pog::Muon* mu2 = goodMuons[1];
      if (!cms.isDataFile()){
          // Scale Facctors!!!
	cms.Apply_lepton_SF(mu1, weight);
	cms.Apply_lepton_SF(mu2, weight, true);
    cms.Apply_dimuon_triggerSF(mu1, mu2, weight);
    weight = weight *0.95;
          // FIX ME.... check the SF for new data and re-reco
	// weight = weight *0.987*0.987;
	// if (fabs(mu2->pt) < 52.) weight = weight *0.97;
      }
      double zMass2 = sqrt(2*mu1->pt_tuneP*mu2->pt_tuneP*(cosh(mu1->eta-mu2->eta)- cos(cms.deltaPhi(mu1->phi,mu2->phi))));

      //Double check Invariant mass
      double px_mu1 = mu1->pt_tuneP*cos(mu1->phi);
      double py_mu1 = mu1->pt_tuneP*sin(mu1->phi);
      double pz_mu1 = mu1->pt_tuneP*sinh(mu1->eta);
      double px_mu2 = mu2->pt_tuneP*cos(mu2->phi);
      double py_mu2 = mu2->pt_tuneP*sin(mu2->phi);
      double pz_mu2 = mu2->pt_tuneP*sinh(mu2->eta);
      double px_Z= px_mu1+px_mu2;
      double py_Z= py_mu1+py_mu2;
      double pz_Z= pz_mu1+pz_mu2;
      double en_mu1= sqrt(0.1057*0.1057+px_mu1*px_mu1+py_mu1*py_mu1+pz_mu1*pz_mu1);
      double en_mu2= sqrt(0.1057*0.1057+px_mu2*px_mu2+py_mu2*py_mu2+pz_mu2*pz_mu2);
      double InvMass = sqrt((en_mu1+en_mu2)*(en_mu1+en_mu2)-(px_Z)*(px_Z)-(py_Z)*(py_Z)-(pz_Z)*(pz_Z));
      ///////

      const Event& ev = cms.event();
      if (zMass2 > 800 && cms.isDataFile()) printf("zMass2 %f, InvMass %f , ratio %f, mu1pt %f, mu2pt %f EVENT: %i:%i:%i \n", zMass2, InvMass, InvMass/zMass2, mu1->pt_tuneP, mu2->pt_tuneP, ev.runNumber, ev.luminosityBlockNumber, ev.eventNumber);


      double pZ =  sqrt(px_Z*px_Z+py_Z*py_Z+pz_Z*pz_Z);
      double pt_Z = sqrt(px_Z*px_Z+py_Z*py_Z);
      if ((fabs(mu1->pt_tuneP) > 200 || fabs(mu2->pt_tuneP) > 200 ) && fabs(mu1->pt_tuneP) > 53 && fabs(mu2->pt_tuneP) > 25 && cms.isDataFile() == true) fprintf(Monitor_file,"%i:%i:%i, %f %f %f ::: pZ= %f and ptZ %f, Delta Phi %f, Delta Eta %f \n",  ev.runNumber,  ev.luminosityBlockNumber, ev.eventNumber, zMass2, mu1->pt_tuneP, mu2->pt_tuneP, pZ, pt_Z, cms.deltaPhi(mu1->phi,mu2->phi), mu1->eta-mu2->eta);
      

      if (mode=="Validation"){
	//printf("weight %f \n", weight);
	cms.FillPlot1D("hPtmuPeak_mu1", iSample, fabs(mu1->pt), weight);
	cms.FillPlot1D("hPtmuTail", iSample, fabs(mu1->pt), weight);
	cms.FillPlot1D("hPhimu", iSample, mu1->phi, weight);
	cms.FillPlot1D("hEtamu", iSample, mu1->eta, weight);
	cms.FillPlot1D("hDxymu", iSample, mu1->dxy, weight);
	cms.FillPlot1D("hDzmu", iSample, mu1->dz, weight);
	cms.FillPlot1D("hTrkIso", iSample, mu1->trackerIso, weight);
	cms.FillPlot1D("hTrkIsoOverPt", iSample, mu1->trackerIso/mu1->pt, weight);
	
	cms.FillPlot1D("hPtmuPeak_mu2", iSample, fabs(mu2->pt), weight);
	cms.FillPlot1D("hPtmuTail", iSample, fabs(mu2->pt), weight);
	cms.FillPlot1D("hPhimu", iSample, mu2->phi, weight);
	cms.FillPlot1D("hEtamu", iSample, mu2->eta, weight);
	cms.FillPlot1D("hDxymu", iSample, mu2->dxy, weight);
	cms.FillPlot1D("hDzmu", iSample, mu2->dz, weight);
	cms.FillPlot1D("hTrkIso", iSample, mu2->trackerIso, weight);
	cms.FillPlot1D("hTrkIsoOverPt", iSample, mu2->trackerIso/mu2->pt, weight);
	
	cms.FillPlot1D("hKMuon", iSample, (1000*mu1->charge)/fabs(mu1->pt), weight);
	cms.FillPlot1D("hKMuon", iSample, (1000*mu2->charge)/fabs(mu2->pt), weight);
	cms.FillPlot1D("hKMuonLowMass", iSample, (1000*mu1->charge)/fabs(mu1->pt), weight);
	cms.FillPlot1D("hKMuonLowMass", iSample, (1000*mu2->charge)/fabs(mu2->pt), weight);
	cms.FillPlot1D("hKMuonPeak", iSample, (1000*mu1->charge)/fabs(mu1->pt), weight);
	cms.FillPlot1D("hKMuonPeak", iSample, (1000*mu2->charge)/fabs(mu2->pt), weight);

	if (mu1->charge ==1){
	  cms.FillPlot1D("hPtmu_P", iSample, fabs(mu1->pt), weight);
      if (fabs(mu1->eta) < 1.2 )  cms.FillPlot1D("hKMuonLowMassBarrel", iSample, (1000*mu1->charge)/fabs(mu1->pt), weight);
	  if (fabs(mu1->eta) < 1.2 ) cms.FillPlot1D("hPtmu_PBarrel", iSample, mu1->pt, weight);
	  if (fabs(mu1->eta) > 1.2 ) cms.FillPlot1D("hPtmu_PEndcap", iSample, mu1->pt, weight);
	  if (fabs(mu1->eta) > 2.1 ) cms.FillPlot1D("hPtmu_PVeryEndcap", iSample, mu1->pt, weight);
	  cms.FillPlot1D("hPtmu_PLargeTail", iSample, fabs(mu1->pt), weight);
	  if (fabs(mu1->eta) < 1.2 ) cms.FillPlot1D("hPtmu_PLargeTailBarrel", iSample, mu1->pt, weight);
	  if (fabs(mu1->eta) > 1.2 && fabs(mu1->eta) < 2.0) cms.FillPlot1D("hPtmu_PLargeTailEndcap", iSample, mu1->pt, weight);
	  if (fabs(mu1->eta) > 2.1 ) cms.FillPlot1D("hPtmu_PLargeTailForward", iSample, mu1->pt, weight);
	  cms.FillPlot1D("hPhimu_P", iSample, mu1->phi, weight);
	  cms.FillPlot1D("hEtamu_P", iSample, mu1->eta, weight);
	  cms.FillPlot1D("hPtmu_N", iSample, mu2->pt, weight);
	  if (fabs(mu2->eta) < 1.2 ) cms.FillPlot1D("hPtmu_NBarrel", iSample, mu2->pt, weight);
	  if (fabs(mu2->eta) > 1.2 ) cms.FillPlot1D("hPtmu_NEndcap", iSample, mu2->pt, weight);
	  if (fabs(mu2->eta) > 2.1 ) cms.FillPlot1D("hPtmu_NVeryEndcap", iSample, mu2->pt, weight);
	  cms.FillPlot1D("hPtmu_NLargeTail", iSample, mu2->pt, weight);
	  if (fabs(mu2->eta) < 1.2 ) cms.FillPlot1D("hPtmu_NLargeTailBarrel", iSample, mu2->pt, weight);
	  if (fabs(mu2->eta) > 1.2 && fabs(mu2->eta) < 2.1) cms.FillPlot1D("hPtmu_NLargeTailEndcap", iSample, mu2->pt, weight);
	  if (fabs(mu2->eta) > 2.1 ) cms.FillPlot1D("hPtmu_NLargeTailForward", iSample, mu2->pt, weight);
	  cms.FillPlot1D("hPhimu_N", iSample, mu2->phi, weight);
	  cms.FillPlot1D("hEtamu_N", iSample, mu2->eta, weight);
	}else{
	  cms.FillPlot1D("hPtmu_P", iSample, fabs(mu2->pt), weight);
      if (fabs(mu2->eta) < 1.2 )  cms.FillPlot1D("hKMuonLowMassBarrel", iSample, (1000*mu2->charge)/fabs(mu2->pt), weight);
	  if (fabs(mu2->eta) < 1.2 ) cms.FillPlot1D("hPtmu_PBarrel", iSample, mu2->pt, weight);
	  if (fabs(mu2->eta) > 1.2 ) cms.FillPlot1D("hPtmu_PEndcap", iSample, mu2->pt, weight);
	  if (fabs(mu2->eta) > 2.1 ) cms.FillPlot1D("hPtmu_PVeryEndcap", iSample, mu2->pt, weight);
	  cms.FillPlot1D("hPtmu_PLargeTail", iSample, fabs(mu2->pt), weight);
	  if (fabs(mu2->eta) < 1.2 ) cms.FillPlot1D("hPtmu_PLargeTailBarrel", iSample, mu2->pt, weight);
	  if (fabs(mu2->eta) > 1.2 && fabs(mu2->eta) < 2.1) cms.FillPlot1D("hPtmu_PLargeTailEndcap", iSample, mu2->pt, weight);
	  if (fabs(mu2->eta) > 2.1 ) cms.FillPlot1D("hPtmu_PLargeTailForward", iSample, mu2->pt, weight);
	  cms.FillPlot1D("hPhimu_P", iSample, mu2->phi, weight);
	  cms.FillPlot1D("hEtamu_P", iSample, mu2->eta, weight);
	  cms.FillPlot1D("hPtmu_N", iSample, mu1->pt, weight);
	  if (fabs(mu1->eta) < 1.2 ) cms.FillPlot1D("hPtmu_NBarrel", iSample, mu1->pt, weight);
	  if (fabs(mu1->eta) > 1.2 ) cms.FillPlot1D("hPtmu_NEndcap", iSample, mu1->pt, weight);
	  if (fabs(mu1->eta) > 2.1 ) cms.FillPlot1D("hPtmu_NVeryEndcap", iSample, mu1->pt, weight);
	  cms.FillPlot1D("hPtmu_NLargeTail", iSample, mu1->pt, weight);
	  if (fabs(mu1->eta) < 1.2 ) cms.FillPlot1D("hPtmu_NLargeTailBarrel", iSample, mu1->pt, weight);
	  if (fabs(mu1->eta) > 1.2 && fabs(mu1->eta) < 2.0) cms.FillPlot1D("hPtmu_NLargeTailEndcap", iSample, mu1->pt, weight);
	  if (fabs(mu1->eta) > 2.1 ) cms.FillPlot1D("hPtmu_NLargeTailForward", iSample, mu1->pt, weight);
	  cms.FillPlot1D("hPhimu_N", iSample, mu1->phi, weight);
	  cms.FillPlot1D("hEtamu_N", iSample, mu1->eta, weight);	
	}
	cms.FillPlot1D("hInvMassPeak", iSample, zMass2, weight);	
	if (fabs(mu1->eta) > 1.2 || fabs(mu2->eta) > 1.2) cms.FillPlot1D("hInvMassPeakEndcap", iSample, zMass2, weight);	
	if (fabs(mu1->eta) > 2.1 || fabs(mu2->eta) > 2.1) cms.FillPlot1D("hInvMassPeakVeryEndcap", iSample, zMass2, weight);	
	cms.FillPlot1D("hInvMassTail", iSample, zMass2, weight);	
	cms.FillPlot1D("hInvMassTailv2", iSample, InvMass, weight);	
	cms.FillPlot1D("hInvMassLargeTail", iSample, zMass2, weight);	
	
	cms.FillPlot1D("hPtZPeak", iSample, pt_Z, weight);	
	cms.FillPlot1D("hPtZTail", iSample, pt_Z, weight);	
      }      


      if (scale == true){
	// if (!cms.isDataFile()){
	//   //Apply bias in TuneP assigment
	//   kappa_mu1 = kappa_mu1 + GoodBias;
	//   kappa_mu2 = kappa_mu2 + GoodBias;
	// }	    
	// float kappa_mu1_refit = 1./(mu1->ptPICKY/1000.); //InTeV 
	// float kappa_mu2_refit = 1./(mu2->ptPICKY/1000.); //InTeV
	float kappa_mu1_refit = mu1->charge/(mu1->pt_tuneP/1000.); //InTeV and TuneP
	float kappa_mu2_refit = mu2->charge/(mu2->pt_tuneP/1000.); //InTeV
	
	if (assignment == "INNER"){
	  kappa_mu1_refit = 1./(mu1->pt_tracker/1000.); //InTeV
	  kappa_mu2_refit = 1./(mu2->pt_tracker/1000.); //InTeV
	}
//	if (assignment == "TPFMS"){
//	  kappa_mu1_refit = 1./(mu1->ptTPFMS/1000.); //InTeV
//	  kappa_mu2_refit = 1./(mu2->ptTPFMS/1000.); //InTeV
//	}/// NEW
	if (assignment == "GLOBAL"){
	  kappa_mu1_refit = 1./(mu1->pt_global/1000.); //InTeV
	  kappa_mu2_refit = 1./(mu2->pt_global/1000.); //InTeV
	}
//	if (assignment == "DYT"){
//	  kappa_mu1_refit = 1./(mu1->ptDYT/1000.); //InTeV
//	  kappa_mu2_refit = 1./(mu2->ptDYT/1000.); //InTeV
//	}
//	if (assignment == "PICKY" || assignment == "PICKY_HIGHPT"){
//	  kappa_mu1_refit = 1./(mu1->ptPICKY/1000.); //InTeV
//	  kappa_mu2_refit = 1./(mu2->ptPICKY/1000.); //InTeV	      
//	}/// NEW
	if (assignment == "TUNEP"){
	  kappa_mu1_refit = mu1->charge/(mu1->pt_tuneP/1000.); //InTeV mu1->pt assignment is not signed in CIEMAT tuples.
	  kappa_mu2_refit = mu2->charge/(mu2->pt_tuneP/1000.); //InTeV
	}
	if (!cms.isDataFile()){
	  //Apply bias in MC refit assignment
	  kappa_mu1_refit = kappa_mu1_refit + GoodBias;
	  kappa_mu2_refit = kappa_mu2_refit + GoodBias;
	}
	
	// General Curvature
	//Above 200 GeV
	cms.FillPlot1D("hKMuon_"+TString(Form("%d",nbin)), iSample, kappa_mu1_refit, weight);
	cms.FillPlot1D("hKMuon_"+TString(Form("%d",nbin)), iSample, kappa_mu2_refit, weight);		
	//Above 100 GeV
	cms.FillPlot1D("hKMuonPeak_"+TString(Form("%d",nbin)), iSample, kappa_mu1_refit, weight);
	cms.FillPlot1D("hKMuonPeak_"+TString(Form("%d",nbin)), iSample, kappa_mu2_refit, weight);
	//Above 50 GeV
	cms.FillPlot1D("hKMuonLowMass_"+TString(Form("%d",nbin)), iSample, kappa_mu1_refit, weight);
	cms.FillPlot1D("hKMuonLowMass_"+TString(Form("%d",nbin)), iSample, kappa_mu2_refit, weight);

	// Inclusive Results
	if (fabs(mu1->eta) < 1.2 ){
        cms.FillPlot1D("hKMuonLowMassBarrel_"+TString(Form("%d",nbin)), iSample, kappa_mu1_refit, weight);
	  cms.FillPlot1D("hKMuonBarrel_"+TString(Form("%d",nbin)), iSample, kappa_mu1_refit, weight);
	}else{
	  cms.FillPlot1D("hKMuonEndcap_"+TString(Form("%d",nbin)), iSample, kappa_mu1_refit, weight);
	  if (mu1->eta > 1.2 ) cms.FillPlot1D("hKMuonPosEndcap_"+TString(Form("%d",nbin)), iSample, kappa_mu1_refit, weight);
	  if (mu1->eta < -1.2 ) cms.FillPlot1D("hKMuonNegEndcap_"+TString(Form("%d",nbin)), iSample, kappa_mu1_refit, weight);
	}
	if (fabs(mu2->eta) < 1.2 ){
      cms.FillPlot1D("hKMuonLowMassBarrel_"+TString(Form("%d",nbin)), iSample, kappa_mu2_refit, weight);
	  cms.FillPlot1D("hKMuonBarrel_"+TString(Form("%d",nbin)), iSample, kappa_mu2_refit, weight);
	}else{
	  cms.FillPlot1D("hKMuonEndcap_"+TString(Form("%d",nbin)), iSample, kappa_mu2_refit, weight);
	  if (mu2->eta > 1.2 ) cms.FillPlot1D("hKMuonPosEndcap_"+TString(Form("%d",nbin)), iSample, kappa_mu2_refit, weight);
	  if (mu2->eta < -1.2 ) cms.FillPlot1D("hKMuonNegEndcap_"+TString(Form("%d",nbin)), iSample, kappa_mu2_refit, weight);	      
	}
	
	
	//Differential Analysis
	if (ScaleMap == true){
	  //First Muon
	  // -2.4,-2.1 
	  if (( mu1->eta > -2.4 && mu1->eta <= -2.1) && (mu1->phi > -M_PI && mu1->phi <= -M_PI*60/180.)) cms.FillPlot1D("hKMuon_0_1_"+TString(Form("%d",nbin)), iSample, kappa_mu1_refit, weight);
	  if (( mu1->eta > -2.4 && mu1->eta <= -2.1) && (mu1->phi > -M_PI*60/180. && mu1->phi <= M_PI*60/180.)) cms.FillPlot1D("hKMuon_0_2_"+TString(Form("%d",nbin)), iSample, kappa_mu1_refit, weight);
	  if (( mu1->eta > -2.4 && mu1->eta <= -2.1) && (mu1->phi > M_PI*60/180. && mu1->phi <= M_PI)) cms.FillPlot1D("hKMuon_0_3_"+TString(Form("%d",nbin)), iSample, kappa_mu1_refit, weight);
	  // -2.1, -1.2
	  if (( mu1->eta > -2.1 && mu1->eta <= -1.2) && (mu1->phi > -M_PI && mu1->phi <= -M_PI*60/180.)) cms.FillPlot1D("hKMuon_1_1_"+TString(Form("%d",nbin)), iSample, kappa_mu1_refit, weight);
	  if (( mu1->eta > -2.1 && mu1->eta <= -1.2) && (mu1->phi > -M_PI*60/180. && mu1->phi <= M_PI*60/180.)) cms.FillPlot1D("hKMuon_1_2_"+TString(Form("%d",nbin)), iSample, kappa_mu1_refit, weight);
	  if (( mu1->eta > -2.1 && mu1->eta <= -1.2) && (mu1->phi > M_PI*60/180. && mu1->phi <= M_PI)) cms.FillPlot1D("hKMuon_1_3_"+TString(Form("%d",nbin)), iSample, kappa_mu1_refit, weight);
	  // -1.2, 0
	  if (( mu1->eta > -1.2 && mu1->eta <= 0.) && (mu1->phi > -M_PI && mu1->phi <= -M_PI*60/180.)) cms.FillPlot1D("hKMuon_2_1_"+TString(Form("%d",nbin)), iSample, kappa_mu1_refit, weight);
	  if (( mu1->eta > -1.2 && mu1->eta <= 0.) && (mu1->phi > -M_PI*60/180. && mu1->phi <= M_PI*60/180.)) cms.FillPlot1D("hKMuon_2_2_"+TString(Form("%d",nbin)), iSample, kappa_mu1_refit, weight);
	  if (( mu1->eta > -1.2 && mu1->eta <= 0.) && (mu1->phi > M_PI*60/180. && mu1->phi <= M_PI)) cms.FillPlot1D("hKMuon_2_3_"+TString(Form("%d",nbin)), iSample, kappa_mu1_refit, weight);
	  // 0, 1.2
	  if (( mu1->eta > 0. && mu1->eta <= 1.2) && (mu1->phi > -M_PI && mu1->phi <= -M_PI*60/180.)) cms.FillPlot1D("hKMuon_3_1_"+TString(Form("%d",nbin)), iSample, kappa_mu1_refit, weight);
	  if (( mu1->eta > 0. && mu1->eta <= 1.2) && (mu1->phi > -M_PI*60/180. && mu1->phi <= M_PI*60/180.)) cms.FillPlot1D("hKMuon_3_2_"+TString(Form("%d",nbin)), iSample, kappa_mu1_refit, weight);
	  if (( mu1->eta > 0. && mu1->eta <= 1.2) && (mu1->phi > M_PI*60/180. && mu1->phi <= M_PI)) cms.FillPlot1D("hKMuon_3_3_"+TString(Form("%d",nbin)), iSample, kappa_mu1_refit, weight);
	  // 1.2, 2.1
	  if (( mu1->eta > 1.2 && mu1->eta <= 2.1) && (mu1->phi > -M_PI && mu1->phi <= -M_PI*60/180.)) cms.FillPlot1D("hKMuon_4_1_"+TString(Form("%d",nbin)), iSample, kappa_mu1_refit, weight);
	  if (( mu1->eta > 1.2 && mu1->eta <= 2.1) && (mu1->phi > -M_PI*60/180. && mu1->phi <= M_PI*60/180.)) cms.FillPlot1D("hKMuon_4_2_"+TString(Form("%d",nbin)), iSample, kappa_mu1_refit, weight);
	  if (( mu1->eta > 1.2 && mu1->eta <= 2.1) && (mu1->phi > M_PI*60/180. && mu1->phi <= M_PI)) cms.FillPlot1D("hKMuon_4_3_"+TString(Form("%d",nbin)), iSample, kappa_mu1_refit, weight);
	  // 2.1, 2.4
	  if (( mu1->eta > 2.1 && mu1->eta <= 2.4) && (mu1->phi > -M_PI && mu1->phi <= -M_PI*60/180.)) cms.FillPlot1D("hKMuon_5_1_"+TString(Form("%d",nbin)), iSample, kappa_mu1_refit, weight);
	  if (( mu1->eta > 2.1 && mu1->eta <= 2.4) && (mu1->phi > -M_PI*60/180. && mu1->phi <= M_PI*60/180.)) cms.FillPlot1D("hKMuon_5_2_"+TString(Form("%d",nbin)), iSample, kappa_mu1_refit, weight);
	  if (( mu1->eta > 2.1 && mu1->eta <= 2.4) && (mu1->phi > M_PI*60/180. && mu1->phi <= M_PI)) cms.FillPlot1D("hKMuon_5_3_"+TString(Form("%d",nbin)), iSample, kappa_mu1_refit, weight);
	  //Second Muon
	  // -2.4,-2.1 
	  if (( mu2->eta > -2.4 && mu2->eta <= -2.1) && (mu2->phi > -M_PI && mu2->phi <= -M_PI*60/180.)) cms.FillPlot1D("hKMuon_0_1_"+TString(Form("%d",nbin)), iSample, kappa_mu2_refit, weight);
	  if (( mu2->eta > -2.4 && mu2->eta <= -2.1) && (mu2->phi > -M_PI*60/180. && mu2->phi <= M_PI*60/180.)) cms.FillPlot1D("hKMuon_0_2_"+TString(Form("%d",nbin)), iSample, kappa_mu2_refit, weight);
	  if (( mu2->eta > -2.4 && mu2->eta <= -2.1) && (mu2->phi > M_PI*60/180. && mu2->phi <= M_PI)) cms.FillPlot1D("hKMuon_0_3_"+TString(Form("%d",nbin)), iSample, kappa_mu2_refit, weight);
	  // -2.1, -1.2
	  if (( mu2->eta > -2.1 && mu2->eta <= -1.2) && (mu2->phi > -M_PI && mu2->phi <= -M_PI*60/180.)) cms.FillPlot1D("hKMuon_1_1_"+TString(Form("%d",nbin)), iSample, kappa_mu2_refit, weight);
	  if (( mu2->eta > -2.1 && mu2->eta <= -1.2) && (mu2->phi > -M_PI*60/180. && mu2->phi <= M_PI*60/180.)) cms.FillPlot1D("hKMuon_1_2_"+TString(Form("%d",nbin)), iSample, kappa_mu2_refit, weight);
	  if (( mu2->eta > -2.1 && mu2->eta <= -1.2) && (mu2->phi > M_PI*60/180. && mu2->phi <= M_PI)) cms.FillPlot1D("hKMuon_1_3_"+TString(Form("%d",nbin)), iSample, kappa_mu2_refit, weight);
	  // -1.2, 0
	  if (( mu2->eta > -1.2 && mu2->eta <= 0.) && (mu2->phi > -M_PI && mu2->phi <= -M_PI*60/180.)) cms.FillPlot1D("hKMuon_2_1_"+TString(Form("%d",nbin)), iSample, kappa_mu2_refit, weight);
	  if (( mu2->eta > -1.2 && mu2->eta <= 0.) && (mu2->phi > -M_PI*60/180. && mu2->phi <= M_PI*60/180.)) cms.FillPlot1D("hKMuon_2_2_"+TString(Form("%d",nbin)), iSample, kappa_mu2_refit, weight);
	  if (( mu2->eta > -1.2 && mu2->eta <= 0.) && (mu2->phi > M_PI*60/180. && mu2->phi <= M_PI)) cms.FillPlot1D("hKMuon_2_3_"+TString(Form("%d",nbin)), iSample, kappa_mu2_refit, weight);
	  // 0, 1.2
	  if (( mu2->eta > 0. && mu2->eta <= 1.2) && (mu2->phi > -M_PI && mu2->phi <= -M_PI*60/180.)) cms.FillPlot1D("hKMuon_3_1_"+TString(Form("%d",nbin)), iSample, kappa_mu2_refit, weight);
	  if (( mu2->eta > 0. && mu2->eta <= 1.2) && (mu2->phi > -M_PI*60/180. && mu2->phi <= M_PI*60/180.)) cms.FillPlot1D("hKMuon_3_2_"+TString(Form("%d",nbin)), iSample, kappa_mu2_refit, weight);
	  if (( mu2->eta > 0. && mu2->eta <= 1.2) && (mu2->phi > M_PI*60/180. && mu2->phi <= M_PI)) cms.FillPlot1D("hKMuon_3_3_"+TString(Form("%d",nbin)), iSample, kappa_mu2_refit, weight);
	  // 1.2, 2.1
	  if (( mu2->eta > 1.2 && mu2->eta <= 2.1) && (mu2->phi > -M_PI && mu2->phi <= -M_PI*60/180.)) cms.FillPlot1D("hKMuon_4_1_"+TString(Form("%d",nbin)), iSample, kappa_mu2_refit, weight);
	  if (( mu2->eta > 1.2 && mu2->eta <= 2.1) && (mu2->phi > -M_PI*60/180. && mu2->phi <= M_PI*60/180.)) cms.FillPlot1D("hKMuon_4_2_"+TString(Form("%d",nbin)), iSample, kappa_mu2_refit, weight);
	  if (( mu2->eta > 1.2 && mu2->eta <= 2.1) && (mu2->phi > M_PI*60/180. && mu2->phi <= M_PI)) cms.FillPlot1D("hKMuon_4_3_"+TString(Form("%d",nbin)), iSample, kappa_mu2_refit, weight);
	  // 2.1, 2.4
	  if (( mu2->eta > 2.1 && mu2->eta <= 2.4) && (mu2->phi > -M_PI && mu2->phi <= -M_PI*60/180.)) cms.FillPlot1D("hKMuon_5_1_"+TString(Form("%d",nbin)), iSample, kappa_mu2_refit, weight);
	  if (( mu2->eta > 2.1 && mu2->eta <= 2.4) && (mu2->phi > -M_PI*60/180. && mu2->phi <= M_PI*60/180.)) cms.FillPlot1D("hKMuon_5_2_"+TString(Form("%d",nbin)), iSample, kappa_mu2_refit, weight);
	  if (( mu2->eta > 2.1 && mu2->eta <= 2.4) && (mu2->phi > M_PI*60/180. && mu2->phi <= M_PI)) cms.FillPlot1D("hKMuon_5_3_"+TString(Form("%d",nbin)), iSample, kappa_mu2_refit, weight);	  
	}
      }     
    }
  }
  

  TString suffix = "ZprimeMuMu_13TeV";
  TString plot_title = "CMS RunII preliminary "+lumi_invfb+" fb^{-1}, #mu #mu, #sqrt{s} = 13 TeV ";
  
  bool mc_study = false; //set this to true if your analysis only has MC (Used only in preparations for data taking)
  // Histograms before kin selection
  if (mode == "Validation"){
      printf("Validation*");
    cms.DrawPlot1D(folder,"hPtmuPeak_mu1", suffix, mc_study, plot_title, "log");
    cms.DrawPlot1D(folder,"hPtmuPeak_mu2", suffix, mc_study, plot_title, "log"); 
    cms.DrawPlot1D(folder,"hPtmuTail", suffix, mc_study, plot_title, "log"); 
    cms.DrawPlot1D(folder,"hPhimu", suffix, mc_study, plot_title); 
    cms.DrawPlot1D(folder,"hEtamu", suffix, mc_study, plot_title); 
    cms.DrawPlot1D(folder,"hDxymu", suffix, mc_study, plot_title);
    cms.DrawPlot1D(folder,"hDzmu", suffix, mc_study, plot_title); 
    cms.DrawPlot1D(folder,"hTrkIso", suffix, mc_study, plot_title, "log"); 
    cms.DrawPlot1D(folder,"hTrkIsoOverPt", suffix, mc_study, plot_title, "log"); 
    
    cms.DrawPlot1D(folder,"hPtmu_P", suffix, mc_study, plot_title, "log"); 
    cms.DrawPlot1D(folder,"hPtmu_PBarrel", suffix, mc_study, plot_title, "log"); 
    cms.DrawPlot1D(folder,"hPtmu_PEndcap", suffix, mc_study, plot_title, "log"); 
    cms.DrawPlot1D(folder,"hPtmu_PVeryEndcap", suffix, mc_study, plot_title, "log"); 
    cms.DrawPlot1D(folder,"hPtmu_PLargeTail", suffix, mc_study, plot_title, "log"); 
    cms.DrawPlot1D(folder,"hPtmu_PLargeTailBarrel", suffix, mc_study, plot_title, "log"); 
    cms.DrawPlot1D(folder,"hPtmu_PLargeTailEndcap", suffix, mc_study, plot_title, "log"); 
    cms.DrawPlot1D(folder,"hPtmu_PLargeTailForward", suffix, mc_study, plot_title, "log"); 
    cms.DrawPlot1D(folder,"hPhimu_P", suffix, mc_study, plot_title);
    cms.DrawPlot1D(folder,"hEtamu_P", suffix, mc_study, plot_title);
    
    cms.DrawPlot1D(folder,"hPtmu_N", suffix, mc_study, plot_title, "log"); 
    cms.DrawPlot1D(folder,"hPtmu_NBarrel", suffix, mc_study, plot_title, "log"); 
    cms.DrawPlot1D(folder,"hPtmu_NEndcap", suffix, mc_study, plot_title, "log"); 
    cms.DrawPlot1D(folder,"hPtmu_NVeryEndcap", suffix, mc_study, plot_title, "log"); 
    cms.DrawPlot1D(folder,"hPtmu_NLargeTail", suffix, mc_study, plot_title, "log"); 
    cms.DrawPlot1D(folder,"hPtmu_NLargeTailBarrel", suffix, mc_study, plot_title, "log"); 
    cms.DrawPlot1D(folder,"hPtmu_NLargeTailEndcap", suffix, mc_study, plot_title, "log"); 
    cms.DrawPlot1D(folder,"hPtmu_NLargeTailForward", suffix, mc_study, plot_title, "log"); 
    cms.DrawPlot1D(folder,"hPhimu_N", suffix, mc_study, plot_title);
    cms.DrawPlot1D(folder,"hEtamu_N", suffix, mc_study, plot_title);
  
    cms.DrawPlot1D(folder,"hInvMassPeak", suffix, mc_study, plot_title); 
    cms.DrawPlot1D(folder,"hInvMassPeakEndcap", suffix, mc_study, plot_title); 
    cms.DrawPlot1D(folder,"hInvMassPeakVeryEndcap", suffix, mc_study, plot_title); 
    cms.DrawPlot1D(folder,"hInvMassTail", suffix, mc_study, plot_title, "log"); 
    cms.DrawPlot1D(folder,"hInvMassTailv2", suffix, mc_study, plot_title, "log"); 
    cms.DrawPlot1D(folder,"hInvMassLargeTail", suffix, mc_study, plot_title, "log"); 
    
    cms.DrawPlot1D(folder,"hPtZPeak", suffix, mc_study, plot_title, "log"); 
    cms.DrawPlot1D(folder,"hPtZTail", suffix, mc_study, plot_title, "log"); 

    cms.DrawPlot1D(folder,"hKMuon", suffix, mc_study, plot_title, "log");
    cms.DrawPlot1D(folder,"hKMuonLowMass", suffix, mc_study, plot_title, "log");
    cms.DrawPlot1D(folder,"hKMuonLowMassBarrel", suffix, mc_study, plot_title, "log");
    cms.DrawPlot1D(folder,"hKMuonPeak", suffix, mc_study, plot_title, "log");
  }
  if (scale == true){
    cms.DrawPlot1D(folder,"hKMuon_"+TString(Form("%d",nbin)), suffix, mc_study, plot_title, "log");
    cms.DrawPlot1D(folder,"hKMuonBarrel_"+TString(Form("%d",nbin)), suffix, mc_study, plot_title, "log");
    cms.DrawPlot1D(folder,"hKMuonEndcap_"+TString(Form("%d",nbin)), suffix, mc_study, plot_title, "log");
    cms.DrawPlot1D(folder,"hKMuonPosEndcap_"+TString(Form("%d",nbin)), suffix, mc_study, plot_title, "log");
    cms.DrawPlot1D(folder,"hKMuonNegEndcap_"+TString(Form("%d",nbin)), suffix, mc_study, plot_title, "log");
    cms.DrawPlot1D(folder,"hKMuonPeak_"+TString(Form("%d",nbin)), suffix, mc_study, plot_title, "log");
    cms.DrawPlot1D(folder,"hKMuonLowMass_"+TString(Form("%d",nbin)), suffix, mc_study, plot_title, "log");
    cms.DrawPlot1D(folder,"hKMuonLowMassBarrel_"+TString(Form("%d",nbin)), suffix, mc_study, plot_title, "log");
	
    if (ScaleMap == true){      
      // -2.4,-2.1 
      cms.DrawPlot1D(folder,"hKMuon_0_1_"+TString(Form("%d",nbin)), suffix, mc_study, plot_title, "log");
      cms.DrawPlot1D(folder,"hKMuon_0_2_"+TString(Form("%d",nbin)), suffix, mc_study, plot_title, "log");
      cms.DrawPlot1D(folder,"hKMuon_0_3_"+TString(Form("%d",nbin)), suffix, mc_study, plot_title, "log");
      // -2.1,-1.2 
      cms.DrawPlot1D(folder,"hKMuon_1_1_"+TString(Form("%d",nbin)), suffix, mc_study, plot_title, "log");
      cms.DrawPlot1D(folder,"hKMuon_1_2_"+TString(Form("%d",nbin)), suffix, mc_study, plot_title, "log");
      cms.DrawPlot1D(folder,"hKMuon_1_3_"+TString(Form("%d",nbin)), suffix, mc_study, plot_title, "log");
      // -1.2,0 
      cms.DrawPlot1D(folder,"hKMuon_2_1_"+TString(Form("%d",nbin)), suffix, mc_study, plot_title, "log");
      cms.DrawPlot1D(folder,"hKMuon_2_2_"+TString(Form("%d",nbin)), suffix, mc_study, plot_title, "log");
      cms.DrawPlot1D(folder,"hKMuon_2_3_"+TString(Form("%d",nbin)), suffix, mc_study, plot_title, "log");
      // 0,1.2 
      cms.DrawPlot1D(folder,"hKMuon_3_1_"+TString(Form("%d",nbin)), suffix, mc_study, plot_title, "log");
      cms.DrawPlot1D(folder,"hKMuon_3_2_"+TString(Form("%d",nbin)), suffix, mc_study, plot_title, "log");
      cms.DrawPlot1D(folder,"hKMuon_3_3_"+TString(Form("%d",nbin)), suffix, mc_study, plot_title, "log");
      // 1.2,2.1 
      cms.DrawPlot1D(folder,"hKMuon_4_1_"+TString(Form("%d",nbin)), suffix, mc_study, plot_title, "log");
      cms.DrawPlot1D(folder,"hKMuon_4_2_"+TString(Form("%d",nbin)), suffix, mc_study, plot_title, "log");
      cms.DrawPlot1D(folder,"hKMuon_4_3_"+TString(Form("%d",nbin)), suffix, mc_study, plot_title, "log");
      // 2.1,2.4 
      cms.DrawPlot1D(folder,"hKMuon_5_1_"+TString(Form("%d",nbin)), suffix, mc_study, plot_title, "log");
      cms.DrawPlot1D(folder,"hKMuon_5_2_"+TString(Form("%d",nbin)), suffix, mc_study, plot_title, "log");
      cms.DrawPlot1D(folder,"hKMuon_5_3_"+TString(Form("%d",nbin)), suffix, mc_study, plot_title, "log");      
    }
  }
    
  fclose(Monitor_file);
  return 0;
}
  
bool WmunuAnalysis::Preselect() {
  // Cut on global event information to save processing time
  const Event& ev = event();
//  if (ev.nMuons<1) return false; // already in skim
    if (2<1) return false; /// NEW
  return true;
}


bool WmunuAnalysis::isTriggered() {
  bool trigger = false;

  muon_pog::HLT& hlt = this->hlt();
  unsigned int nhlt = hlt.triggers.size();
  for (unsigned int ih=0; ih<nhlt; ++ih) {
    if (hlt.triggers[ih].find("HLT_Mu50_v")!=std::string::npos) trigger = true;
    if (hlt.triggers[ih].find("HLT_TkMu50_")!=std::string::npos) trigger = true;
    if (hlt.triggers[ih].find("HLT_IsoMu22_")!=std::string::npos) trigger = true;
    if (hlt.triggers[ih].find("HLT_IsoTkMu22_")!=std::string::npos) trigger = true;
    //if (hlt.triggers[ih].find("HLT_")!=std::string::npos) trigger = true;      
    //printf("Triggers fired? %d %s \n",ih, hlt.triggers[ih].data());    	  
    //if (hlt.triggers[ih].find("HLT_Mu40")!=std::string::npos) trigger = true;      
    //if (hlt.triggers[ih].find("HLT_MET80_")!=std::string::npos) trigger = true;
    if (trigger) break;
  }
  return trigger;
}


std::vector<muon_pog::Muon*> WmunuAnalysis::GoodZprimeMuons() {
  // Return already if there is nothing to do
  std::vector<muon_pog::Muon*> selectedMuons;
//    printf("inside good muons");

  bool verbose = true;//For debugging the seleciton
  // Reject if two or more muons with pt>25 GeV
  unsigned nmuons25 = 0;
  if (verbose == true )printf("Number of Muons %i \n", nMuons());

  for (unsigned int im=0; im<nMuons(); ++im) {
//      printf("inside good muons");

    if (muon(im).pt_tuneP>25.) nmuons25++;

  }
  if (nmuons25<2) return selectedMuons;
  if (verbose == true )printf("Number of Muons above 25 %i \n", nmuons25);

  for (unsigned int im1=0; im1<nMuons(); ++im1) {
    Muon* mu1 = pmuon(im1);
    if (verbose == true ) printf("Muon1 pt %f \n", mu1->pt_tuneP);
    if (fabs(mu1->pt_tuneP)<53.) continue;
    if (!mu1->isHighPt) continue; // high-pt Id
    if (fabs(mu1->eta)>2.4) continue; // eta cut
    if (mu1->bestMuPtErr>0.3*mu1->pt_tuneP) continue; // reasonable momentum uncertainty
    if (mu1->isTracker>0.1*mu1->pt) continue; // track-isolated
    if (fabs(mu1->dxy)>0.02) continue; // d0 cut

    for (unsigned int im2=im1+1; im2<nMuons(); ++im2) {
      Muon* mu2 = pmuon(im2);
      if (verbose == true ) printf("Muon2 pt %f \n", mu2->pt_tuneP);
      if (mu1->charge*mu2->charge == 1) continue;
      if (fabs(mu2->pt_tuneP)<25.) continue;
      //if (fabs(mu2->pt)<53.) continue;
      if (!mu2->isHighPt) continue; // high-pt Id
      if (fabs(mu2->eta)>2.4) continue; // eta cut
      if (mu2->bestMuPtErr>0.3*mu2->pt_tuneP) continue; // reasonable momentum uncertainty
      if (mu2->isTracker>0.1*mu2->pt) continue; // track-isolated
      if (fabs(mu2->dxy)>0.02) continue; // d0 cut

      double zMass2 = 2*mu1->pt_tuneP*mu2->pt_tuneP*(cosh(mu1->eta-mu2->eta)-cos(deltaPhi(mu1->phi, mu2->phi)));
      zMass2 = sqrt(zMass2);
      if (verbose == true )  printf("Invariant mass %f \n", zMass2);
      if (zMass2<50.) continue;

      if (mu1->pt_tuneP > mu2->pt_tuneP){
        selectedMuons.push_back(mu1);
        selectedMuons.push_back(mu2);
      }else{
        selectedMuons.push_back(mu2);
        selectedMuons.push_back(mu1);
      }
      return selectedMuons;
    }
  }
  return selectedMuons;      
}
void WmunuAnalysis::Apply_lepton_SF(Muon* mu, double& weight, bool notrig){
  unsigned int etaBINS = 4;

  //Reference: https://twiki.cern.ch/twiki/bin/view/CMS/MuonReferenceEffsRun2                                                                                                                                                  

//  //Here I set the SF for trigger                                                                                                                                                                                              
//  double SF_trigger[etaBINS];
//  if (fabs(mu->pt) < 60){
//    SF_trigger[0] = 0.911513; SF_trigger[1] = 0.882274; SF_trigger[2] = 0.839448; SF_trigger[3] = 0.723469;
//  }else if(fabs(mu->pt) < 140){
//    SF_trigger[0] = 0.914351; SF_trigger[1] = 0.883462; SF_trigger[2] = 0.838485; SF_trigger[3] = 0.7526;
//  }else{
//    SF_trigger[0] = 0.91233; SF_trigger[1] = 0.856823; SF_trigger[2] = 0.801867; SF_trigger[3] = 0.708108;
//  }
  //Here I set the SF for the HighPt muonID
  double SF_id[etaBINS];
  if (fabs(mu->pt) < 60){
    SF_id[0] = 0.97894; SF_id[1] = 0.97200; SF_id[2] = 0.99446; SF_id[3] = 0.97859;
  }else if(fabs(mu->pt) < 120){
    SF_id[0] = 0.97914; SF_id[1] = 0.97798; SF_id[2] = 0.996493; SF_id[3] = 0.97932;
  }else{
    SF_id[0] = 0.97914; SF_id[1] = 0.97798; SF_id[2] = 0.996493; SF_id[3] = 0.97932;
  }
  //Here I set the SF for the Isolation                                                                                                                                                                                        
  double SF_iso[etaBINS];
  if (fabs(mu->pt) < 60){
    SF_iso[0] = 0.999325; SF_iso[1] = 0.99888; SF_iso[2] = 0.998823; SF_iso[3] = 0.998515;
  }else if(fabs(mu->pt) < 120){
    SF_iso[0] = 0.999048; SF_iso[1] = 0.999334; SF_iso[2] = 0.999504; SF_iso[3] = 0.999914;
  }else{
    SF_iso[0] = 0.999048; SF_iso[1] = 0.999334; SF_iso[2] = 0.999504; SF_iso[3] = 0.999914;
  }

  //Eta bins definition.                                                                                                                                                                                                       
  double etabin[etaBINS+1];
  etabin[0]=0.; etabin[1]=0.9; etabin[2]=1.2; etabin[3]=2.1; etabin[4]=2.4;

  unsigned int ibin = etaBINS;
  for (unsigned int kbin=0; kbin<=etaBINS; ++kbin) {
    if (fabs(mu->eta)<etabin[kbin+1]) {
      ibin = kbin;
      break;
    }
  }
//  double sf = SF_trigger[ibin]*SF_id[ibin]*SF_iso[ibin];
//  if (notrig == true && fabs(mu->pt) > 51.5) sf = SF_id[ibin]*SF_iso[ibin];
  double sf = SF_id[ibin]*SF_iso[ibin];
  weight = weight*sf;
					
}
void WmunuAnalysis::Muon_triggerSF(Muon* mu, double& weight){
    unsigned int etaBINS = 4;

    //Here I set the SF for trigger
    double SF_trigger[etaBINS];
    if (fabs(mu->pt) < 60){
        SF_trigger[0] = 0.911513; SF_trigger[1] = 0.882274; SF_trigger[2] = 0.839448; SF_trigger[3] = 0.723469;
    }else if(fabs(mu->pt) < 140){
        SF_trigger[0] = 0.914351; SF_trigger[1] = 0.883462; SF_trigger[2] = 0.838485; SF_trigger[3] = 0.7526;
    }else{
        SF_trigger[0] = 0.91233; SF_trigger[1] = 0.856823; SF_trigger[2] = 0.801867; SF_trigger[3] = 0.708108;
    }
    
    //Eta bins definition.
    double etabin[etaBINS+1];
    etabin[0]=0.; etabin[1]=0.9; etabin[2]=1.2; etabin[3]=2.1; etabin[4]=2.4;
    
    unsigned int ibin = etaBINS;
    for (unsigned int kbin=0; kbin<=etaBINS; ++kbin) {
        if (fabs(mu->eta)<etabin[kbin+1]) {
            ibin = kbin;
            break;
        }
    }
    double sf = SF_trigger[ibin];
    weight = weight*sf;
}
void WmunuAnalysis::Apply_dimuon_triggerSF(Muon* mu1, Muon* mu2, double& weight){
    double eff1, eff2;
    eff1=1; eff2=1;
    Muon_triggerSF(mu1, eff1);
    Muon_triggerSF(mu2, eff2);
    weight = weight*(eff1+eff2-(eff1*eff2));
    
    
}

Double_t WmunuAnalysis::deltaPhi(Double_t Phi0, Double_t Phi1) {

   Double_t dphi = Phi0-Phi1;

   while (dphi >= TMath::Pi()) dphi -= TMath::TwoPi();
   while (dphi < -TMath::Pi()) dphi += TMath::TwoPi();

   return dphi;
    
    
}
