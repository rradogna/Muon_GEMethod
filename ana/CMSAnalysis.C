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
#include "TLatex.h"
#include "TColor.h"
#include "TMath.h"
#include "TF1.h"
#include "TLine.h"
#include "tdrstyle.C"

#include "MuonPogTree.h"
#include "CMSAnalysis.h"

CMSAnalysis::CMSAnalysis() {
      _file = 0;
      _tree = 0; 
      _currentIndex = -1; 
      _dataIndex = -1;

      _bEvent = 0;

      _iEvent = -1;
      _pevent = 0;

};

CMSAnalysis::~CMSAnalysis() { 
      if (_file) _file->Close();
};

void CMSAnalysis::AddDataSample(const TString& id, const TString& file, double luminosity) {
  _lumi = luminosity;
  _sampleId.push_back(id);
  _sampleFile.push_back(file);

  //TFile file_tmp(file,"READONLY");
  TFile* file_tmp=TFile::Open(file,"READONLY");
  TTree* tree_tmp = 0;
  file_tmp->GetObject("MUONPOGTREE",tree_tmp);
  if (!tree_tmp) file_tmp->GetObject("MuonPogTree",tree_tmp);
  if (!tree_tmp) file_tmp->GetObject("MuonPogTree/MUONPOGTREE",tree_tmp);
  if (!tree_tmp) file_tmp->GetObject("MUONPOGTREE/MuonPogTree",tree_tmp);
  if (!tree_tmp) printf("Error reading tree for file %s; crash expected!!\n", file.Data());
  _sampleNevents.push_back(tree_tmp->GetEntriesFast());
  int nDump = 1; while (_sampleNevents.back()/nDump>50) nDump *= 10;
  _nDump10.push_back(nDump);

  _sampleXsection.push_back(_sampleNevents.back()/_lumi);
  _sampleWeight.push_back(1.);
  _dataIndex = _sampleId.size()-1;
  _signalFlag.push_back(false); // set to false to start with
  _plotFlag.push_back(true); // set to true to start with
  _modelFlag.push_back(false); // set to false to start with
  printf("Added Data: %s lumi %.2f \n", 
	 file. Data(), _lumi );

}

void CMSAnalysis::AddMCSample(const TString& id, const TString& file, double maxevents, double xsec, double total_events_for_xsection, double prescale) {
  if (_dataIndex<0) {
      printf(">>> Warning: you should call AddDataSample first, to define the luminosity!!!\n");
      printf(">>> SAMPLE NOT ADDED!\n");
      return;
  }
  _sampleId.push_back(id);
  _sampleFile.push_back(file);

//  TFile file_tmp(file,"READONLY");
  TFile* file_tmp=TFile::Open(file,"READONLY");
  TTree* tree_tmp = 0;
  file_tmp->GetObject("MUONPOGTREE",tree_tmp);
  if (!tree_tmp) file_tmp->GetObject("MuonPogTree",tree_tmp);
  if (!tree_tmp) file_tmp->GetObject("MuonPogTree/MUONPOGTREE",tree_tmp);
  if (!tree_tmp) file_tmp->GetObject("MUONPOGTREE/MuonPogTree",tree_tmp);
  if (!tree_tmp) printf("Error reading tree for file %s; crash expected!!\n", file.Data());
  int maxeventsInTreeUnprescaled = tree_tmp->GetEntriesFast();
  int maxeventsInTree = maxeventsInTreeUnprescaled/prescale;
  _sampleNevents.push_back(maxeventsInTree);

  int nDump = 1; while (_sampleNevents.back()/nDump>50) nDump *= 10;
  _nDump10.push_back(nDump);

  _sampleXsection.push_back(xsec);
  _signalFlag.push_back(false); // set to false to start with
  _plotFlag.push_back(true); // set to true to start with
  _modelFlag.push_back(false); // set to false to start with

  double equivalent_lumi = maxeventsInTree/xsec;
  if (total_events_for_xsection>0 && xsec>0) equivalent_lumi = total_events_for_xsection/xsec;
  //printf("Sample %s with negative events? %f, TotalEventsForXsec %i,  MaxEvents %i \n ", id.Data(), total_events_for_xsection/maxevents, total_events_for_xsection, maxevents);
  if (maxevents>0 && xsec>0) equivalent_lumi = (total_events_for_xsection/prescale)/xsec;
  _sampleWeight.push_back(_lumi/equivalent_lumi);
  printf("Added MC: Name %s lumi %.2f, xsec %f, NTupleEvents %.0d, Processing %.0d, FullSampleEvents %.0f, NegFrac %.2f, EqLumi %f, Weight %f, prescale %.0f\n", 
	 id.Data(), _lumi, _sampleXsection.back(), maxeventsInTreeUnprescaled, _sampleNevents.back(), maxevents, total_events_for_xsection/maxevents, equivalent_lumi, _sampleWeight.back(), prescale);
}

void CMSAnalysis::AddMCSignalSample(const TString& id, const TString& file, double maxevents, double xsec, double total_events_for_xsection, double prescale) {

  if (_dataIndex<0) {
      printf(">>> Warning: you should call AddDataSample first, to define the luminosity!!!\n");
      printf(">>> SAMPLE NOT ADDED!\n");
      return;
  }
  _sampleId.push_back(id);
  _sampleFile.push_back(file);

  TFile file_tmp(file,"READONLY");
  TTree* tree_tmp = 0;
  file_tmp.GetObject("MUONPOGTREE",tree_tmp);
  if (!tree_tmp) file_tmp.GetObject("MuonPogTree",tree_tmp);
  if (!tree_tmp) file_tmp.GetObject("MuonPogTree/MUONPOGTREE",tree_tmp);
  if (!tree_tmp) file_tmp.GetObject("MUONPOGTREE/MuonPogTree",tree_tmp);
  if (!tree_tmp) printf("Error reading tree for file %s; crash expected!!\n", file.Data());
  int maxeventsInTreeUnprescaled = tree_tmp->GetEntriesFast();
  int maxeventsInTree = maxeventsInTreeUnprescaled/prescale;
  _sampleNevents.push_back(maxeventsInTree);

  int nDump = 1; while (_sampleNevents.back()/nDump>50) nDump *= 10;
  _nDump10.push_back(nDump);

  _sampleXsection.push_back(xsec);
  _signalFlag.push_back(true); // set to false to start with
  _plotFlag.push_back(true); // set to true to start with
  _modelFlag.push_back(false); // set to false to start with

  double equivalent_lumi = maxeventsInTree/xsec;
  //if (total_events_for_xsection>0 && xsec>0) equivalent_lumi = total_events_for_xsection/xsec;  
  //  if (total_events_for_xsection != maxeventsInTree ) printf("Sample %s with negative events? %f \n ", id.Data(), total_events_for_xsection/maxevents);
  if (total_events_for_xsection>0 && xsec>0) equivalent_lumi = maxevents*(total_events_for_xsection/maxevents)/xsec;
  //if (maxevents>0 && xsec>0) equivalent_lumi = maxevents/xsec;
  if (maxevents>0 && xsec>0) equivalent_lumi = (total_events_for_xsection/prescale)/xsec;

  _sampleWeight.push_back(_lumi/equivalent_lumi);
  printf("Added MC Signal: Name %s lumi %.2f, xsec %f, NTupleEvents %.0d, Processing %.0d, FullSampleEvents %.0f, NegFrac %.2f, EqLumi %f, Weight %f, prescale %.0f\n", 
	 id.Data(), _lumi, _sampleXsection.back(), maxeventsInTreeUnprescaled, _sampleNevents.back(), maxevents, total_events_for_xsection/maxevents, equivalent_lumi, _sampleWeight.back(), prescale);
}

void CMSAnalysis::AddModel(const TString& id, const TString& file, double maxevents, double xsec, double total_events_for_xsection, double prescale, bool plot) {
  if (_dataIndex<0) {
      printf(">>> Warning: you should call AddDataSample first, to define the luminosity!!!\n");
      printf(">>> SAMPLE NOT ADDED!\n");
      return;
  }
  _sampleId.push_back(id);
  _sampleFile.push_back(file);

  TFile file_tmp(file,"READONLY");
  TTree* tree_tmp = 0;
  file_tmp.GetObject("MUONPOGTREE",tree_tmp);
  if (!tree_tmp) file_tmp.GetObject("MuonPogTree",tree_tmp);
  if (!tree_tmp) file_tmp.GetObject("MuonPogTree/MUONPOGTREE",tree_tmp);
  if (!tree_tmp) file_tmp.GetObject("MUONPOGTREE/MuonPogTree",tree_tmp);
  if (!tree_tmp) printf("Error reading tree for file %s; crash expected!!\n", file.Data());
  int maxeventsInTreeUnprescaled = tree_tmp->GetEntriesFast();
  int maxeventsInTree = maxeventsInTreeUnprescaled/prescale;
  _sampleNevents.push_back(maxeventsInTree);

  int nDump = 1; while (_sampleNevents.back()/nDump>50) nDump *= 10;
  _nDump10.push_back(nDump);

  _sampleXsection.push_back(xsec);
  _signalFlag.push_back(false); // set to false to start with
  _plotFlag.push_back(plot); 
  _modelFlag.push_back(true); // set to false to start with

  double equivalent_lumi = maxeventsInTree/xsec;
  //if (total_events_for_xsection>0 && xsec>0) equivalent_lumi = total_events_for_xsection/xsec;
  if (maxevents>0 && xsec>0) equivalent_lumi = (total_events_for_xsection/prescale)/xsec;
  _sampleWeight.push_back(_lumi/equivalent_lumi);
  printf("Added MC Signal: Name %s lumi %.2f, xsec %f, NTupleEvents %.0d, Processing %.0d, FullSampleEvents %.0f, NegFrac %.2f, EqLumi %f, Weight %f, prescale %.0f\n", 
	 id.Data(), _lumi, _sampleXsection.back(), maxeventsInTreeUnprescaled, _sampleNevents.back(), maxevents, total_events_for_xsection/maxevents, equivalent_lumi, _sampleWeight.back(), prescale);

}

void CMSAnalysis::AddPlot1D(const TString& name, const TString& title, int nbins, double xmin, double xmax) {
      for (unsigned int i=0; i<_sampleId.size(); ++i) {
            bool existing = false;
            for (unsigned int j=0; j<hists_1D.size(); ++j) {
                  TString thisname = hists_1D[j]->GetName();
                  if (thisname == _sampleId[i]+"_"+name) {
                        existing = true; 
                        break;
                  }
            }
            if (existing) continue;

            hists_1D.push_back(new TH1D(_sampleId[i]+"_"+name, title, nbins, xmin, xmax));
            hists_1D[hists_1D.size()-1]->Sumw2();
      }
}
void CMSAnalysis::AddPlot1DRB(const TString& name, const TString& title, int nbins, double xbins[150]) {
      for (unsigned int i=0; i<_sampleId.size(); ++i) {
            bool existing = false;
            for (unsigned int j=0; j<hists_1D.size(); ++j) {
                  TString thisname = hists_1D[j]->GetName();
                  if (thisname == _sampleId[i]+"_"+name) {
                        existing = true; 
                        break;
                  }
            }
            if (existing) continue;

            hists_1D.push_back(new TH1D(_sampleId[i]+"_"+name, title, nbins , xbins));
            hists_1D[hists_1D.size()-1]->Sumw2();
      }
}
  
void CMSAnalysis::FillPlot1D(const TString& name, int isample, double value, double weight) {
      for (unsigned int j=0; j<hists_1D.size(); ++j) {
            if (hists_1D[j]->GetName()==_sampleId[isample]+"_"+name) {
                  hists_1D[j]->Fill(value,_sampleWeight[isample]*weight);
                  return;
            }
      }
}
void CMSAnalysis::NormPlot1D(const TString& name, int nbins, double xbins[150], TString ytitle) {
  bool Wsample= false;
  for (unsigned int j=0; j<hists_1D.size(); ++j) {
    for (unsigned int i=0; i<_sampleId.size(); ++i) {
      if (hists_1D[j]->GetName()==_sampleId[i]+"_"+name){
	if (Wsample  == true && _sampleId[i]=="W+jets") continue;
	if (_sampleId[i]=="W+jets" ) Wsample = true;	
	//	printf("j sample %i %s \n", j, _sampleId[i].Data());	
	for (unsigned int kbin = 1; kbin <=  nbins+1; kbin++) {
	  //	  printf("kbin %i/%i bincontent %f with %f \n", kbin,nbins,hists_1D[j]->GetBinContent(kbin), xbins[kbin-1]);
	  hists_1D[j]->SetBinContent(kbin, hists_1D[j]->GetBinContent(kbin)/xbins[kbin-1]);	  
	  hists_1D[j]->SetBinError(kbin, hists_1D[j]->GetBinError(kbin)/xbins[kbin-1]);	  
	  if (i==0 ) hists_1D[j]->SetYTitle(ytitle);
	}
      }
    }
  }
}

void CMSAnalysis::DrawPlot1D(const TString& folder, const TString& name, const TString& suffix, const bool mc_study, const TString& title, const TString& type) {
  setTDRStyle();
  TFile *f = new TFile(folder+name+"-"+suffix+".root","RECREATE");

  //gROOT->SetStyle("Pub"); 
  gROOT->SetStyle("Plain"); 
      
  gStyle->SetPadGridX(true); 
  gStyle->SetPadGridY(true);
  //gStyle->SetOptStat(1111111);
  //gStyle->SetOptStat(1111);
  gStyle->SetOptStat(0);

  //TLegend* leg = new TLegend(0.52,0.6,0.85,0.85);
  TLegend* leg = new TLegend(0.70,0.60,.90,.90); // to see the whole histogram
  leg->SetFillColor(0);

  //Histograms for templates
  TH1D* histData;
  TH1D* hists[90];
  //
  THStack* hMCStack = new THStack(name,name+" histograms");
  TH1D* hData;      
  TH1D* hDenMC;      
  // Get data; stack all MC except signal MC
  unsigned int nhists = hists_1D.size();
  //int colors[10] = {46, 47, 48, 49, 50, 51, 52, 53, 54, 55};
  //int colors[10] = {6, 7, 8, 9, 11, 44, 46, 2, 3, 4}; 
  //int colors[10] = {TColor::GetColor("#ff55ff"), TColor::GetColor("#00ffff"), TColor::GetColor("#ff5500"), 9, 11, 44, 46, 2, 3, 4}; 
  //      int colors[10] = {TColor::GetColor("#00ddff"), TColor::GetColor("#0055ff"), kMagenta, kRed, 11, 44, 46, 2, 3, 4}; 
  int colors[10] = {TColor::GetColor("#00ddff"), TColor::GetColor("#0055ff"), kRed, 8, 6, 44, 46, 5, 3, 4}; 
  int colors_model[10] = {6, 2, 4, 8, 6, 44, 46, 5, 3, 4}; 
  int mcindex = -1;
  int nmodels = 0;
  for (unsigned int j=0; j<nhists; ++j) {
    TString histname = hists_1D[j]->GetName();
    TString suffix = "_" + name;
    if (!histname.EndsWith(suffix)) continue;    
    for (unsigned int i=0; i<_sampleFile.size(); ++i) {
      TString prefix = _sampleId[i] + "_";
      if (!histname.BeginsWith(prefix)) continue;
      if (histname==_sampleId[_dataIndex]+"_"+name) {
	hData = hists_1D[j];
	hData->Sumw2();			
	hData->SetMarkerStyle(20);
	hData->SetMarkerSize(1.0);
	if (mc_study == false) leg->AddEntry(hData,_sampleId[i].Data(),"P");
	
	// Needed to save data histogram for templates.
	histData = (TH1D*)hists_1D[j]->Clone();
	histData->Sumw2();
	histData->SetMarkerStyle(20);
	histData->SetMarkerSize(1.0);
	TString tmp = hData->GetName();
	tmp.Append("_ratio");
	hData->SetName(tmp);
	//
	//			printf("histname ==  %s\n",histname.Data());
	break;			
      } else {
	if (i == 1 ){// First MC Sample
	  hDenMC = (TH1D*)hists_1D[j]->Clone();
	}else if (_modelFlag[i] == false) {
	  hDenMC->Add(hists_1D[j]);					      
	  hDenMC->Sumw2();				      		      
	}
	if (_modelFlag[i] == false) {
	  mcindex++;
	  int color = colors[mcindex%10];
	  hists_1D[j]->SetLineWidth(3);
	  hists_1D[j]->SetFillColor(color);
	  hists_1D[j]->SetLineColor(TColor::GetColorDark(color));
	  leg->AddEntry(hists_1D[j],_sampleId[i].Data(),"F");
	}
	if (_modelFlag[i] == true){
	  nmodels++;
	  if( _plotFlag[i]==true) {	    
	    int color_model = colors_model[nmodels%10];
	    hists_1D[j]->SetLineWidth(4);
	    hists_1D[j]->SetLineColor(TColor::GetColorDark(color_model));
	    leg->AddEntry(hists_1D[j],_sampleId[i].Data(),"l");
	  }
	}
	break;
      }      
    }
  }
  // Needed for MC templates.
  TString tmp = hDenMC->GetName();
  tmp.Append("_sumallMC");
  hDenMC->SetName(tmp);
  
  // Add the signal component to the stack now
  int firstSignalSample = 1;
  for (unsigned int j=nhists-1; j>=0; j--) {
    if (hists_1D[j]->GetName()==_sampleId[_sampleFile.size()-1]+"_"+name ) {		    
      firstSignalSample=j-nmodels;
      //	  printf("first signal sample = %i\n", j);
      for (unsigned int i=firstSignalSample; i>=firstSignalSample-mcindex; i--) {//mcindex = #MC histogrms -1
	hMCStack->Add(hists_1D[i]);
      }
      break;
    }
  }
  // Needed for MC templates.
  int count = 0;
  for (unsigned int j=1; j<nhists; ++j) {
    TString histname = hists_1D[j]->GetName();
    if (!histname.EndsWith(name)) continue;
    if ( count < 90 ) {
      hists[count] = (TH1D*)hists_1D[j]->Clone();
      hists[count]->SetFillColor(0);
    }
    count++;
  }
  
  
  TString c1name = "c1_" + name;
  TCanvas* c1 = new TCanvas(c1name.Data(),c1name.Data(),10,10,600,600);
  c1->cd();

  if (mc_study ==  false ){ // Do we need ratio plot for MC studies?
    TPad *pad1 = new TPad("pad1","top pad",0,0.3,1,1);
    pad1->SetBottomMargin(0);
    pad1->Draw();
    //
    TPad *pad2 = new TPad("pad2","bottom pad",0,0.0,1,0.3);
    pad2->SetTopMargin(0);
    pad2->SetBottomMargin(0.22);
    pad2->Draw(); 
    pad2->cd();	       
    TH1D *hNumData = (TH1D*)hData->Clone();    
    hNumData->Sumw2();
    hNumData->Divide(hDenMC);
    hNumData->Draw("ep");
    hNumData->SetMinimum(0.5);
    hNumData->SetMaximum(1.5);
    hNumData->SetTitle("");
    TAxis* xax = hNumData->GetXaxis();
    TAxis* yax = hNumData->GetYaxis();
    xax->SetTitleSize(0.10);
    xax->SetLabelSize(0.09);
    xax->SetTitleOffset(0.9);
    xax->SetTitle(hData->GetTitle());
    yax->SetTitleSize(0.10);
    yax->SetTitle("Data/MC");
    yax->SetLabelSize(0.09);
    yax->SetTitleOffset(0.4);
    yax->SetNdivisions(405);      
    Double_t minxaxis = xax->GetXmin();
    Double_t maxxaxis = xax->GetXmax();
    TLine *l=new TLine(minxaxis,1.0,maxxaxis,1.0);
    l->SetLineColor(2);
    l->SetLineWidth(2);
    l->Draw();
    pad1->cd();
  }

  // Drawing the data histogram, the stack of montecarlos and the models.
  hData->Draw("e");
  hMCStack->Draw("samehist");

  //Drawing the models
  for (unsigned int j=0; j<nhists; ++j) {
    TString histname = hists_1D[j]->GetName();
    TString suffix = "_" + name;
    if (!histname.EndsWith(suffix)) continue;
    for (unsigned int i=0; i<_sampleFile.size(); ++i) {
      TString prefix = _sampleId[i] + "_";
      if (!histname.BeginsWith(prefix)) continue;
      //      printf("Histoname %s maximum %f, total max %f \n", histname.Data(), hists_1D[j]->GetMaximum(), hData->GetMaximum());
      if (_modelFlag[i] ==  true && _plotFlag[i]==true) {
	if (hData->GetMaximum() < 1.6*hists_1D[j]->GetMaximum() ) hData->SetMaximum(1.6*hists_1D[j]->GetMaximum());
	hists_1D[j]->Draw("same");
      }
    }
  }
	
  // This draws, statistical error and data over all simulations.
  if (mc_study == false){
    hDenMC->SetFillColor(1);
    hDenMC->SetFillStyle(3002);
    hDenMC->Draw("e2 same");
    hData->Draw("esame");
    leg->AddEntry(hDenMC,"MC Uncertainty","F");
  }


  hData->SetXTitle(hData->GetTitle());
  hData->SetTitle("");
  hData->SetTitleOffset(1.2);

  // Setting the maximum of the histogram
  if (hData->GetMinimum()>0.) hData->SetMinimum(.1);
  if (type =="log" ) {
    gPad->SetLogy(1);
    hData->SetMaximum(1000*hData->GetMaximum());
    hData->SetMinimum(0.0001);
  }
  if (type =="normal" ) {
    hData->SetMaximum(1.2*hData->GetMaximum());
  }

  gPad->SetTicks(1,1);
  gPad->RedrawAxis();
  
  TLatex* preliminary = new TLatex(0.25,0.92,title);
  preliminary->SetTextSize(0.035);
  preliminary->SetNDC();
  preliminary->SetTextFont(42);
  preliminary->Draw();

  leg->Draw();

  //Ratio plot Data Over MC/Stack

  if (suffix!="") {
    c1->SaveAs(folder+name+"_"+suffix+".root");
    c1->SaveAs(folder+name+"_"+suffix+".jpg");
    c1->SaveAs(folder+name+"_"+suffix+".pdf");
    c1->SaveAs(folder+name+"_"+suffix+".eps");
  } else {
    c1->SaveAs(folder+name+".root");
    c1->SaveAs(folder+name+".jpg");
    c1->SaveAs(folder+name+".pdf");
    c1->SaveAs(folder+name+".eps");
  }
  f->Write();
  delete f;
}
  
void CMSAnalysis::SetTree(int i) {
      if (_file) {
            _file->Close();
            _file = NULL;
      }

      printf("Processing sample '%s'...\n", _sampleFile[i].Data());
      _currentIndex = i;
      TFile* _file = TFile::Open(_sampleFile[i],"READONLY");
      _tree = 0;
      _file->GetObject("MUONPOGTREE",_tree);
      if (!_tree) _file->GetObject("MuonPogTree",_tree);
      if (!_tree) _file->GetObject("MuonPogTree/MUONPOGTREE",_tree);
      if (!_tree) _file->GetObject("MUONPOGTREE/MuonPogTree",_tree);

      _bEvent = _tree->GetBranch("event");
      _bEvent->SetAddress(&_pevent);

      int nentriesInTree = _tree->GetEntriesFast();
      printf("\tReading %d entries from a total of %d\n", _sampleNevents[i], nentriesInTree);

}

Int_t CMSAnalysis::GetEntry(int iEvent) {
      if (_tree->LoadTree(iEvent)<0) return -1;
      if (iEvent%_nDump10[_currentIndex]==0) printf("... event index %d\n", iEvent);

      _iEvent = iEvent;
      if (_pevent) delete _pevent; _pevent = 0;
      return 0;
}

// EVENT
muon_pog::Event* CMSAnalysis::pevent() {
      if (_pevent==0) _bEvent->GetEntry(_iEvent); 
      return (_pevent);
}
muon_pog::Event& CMSAnalysis::event(){return *pevent();};

//// GENINFO
std::vector<muon_pog::GenInfo> CMSAnalysis::genInfos() {
    if (_pevent==0) _bEvent->GetEntry(_iEvent);
    //      if (!(_pevent->isMC)) return 0;
    if (_pevent->genInfos.size())
        return (_pevent->genInfos);
}
muon_pog::GenInfo* CMSAnalysis::pgenInfo(unsigned int i) {
    if (_pevent==0) _bEvent->GetEntry(_iEvent);
    if (!(_pevent->genInfos.size())) return 0;
    return (&(_pevent->genInfos.at(i)));
}
muon_pog::GenInfo& CMSAnalysis::genInfo(unsigned int i){return *pgenInfo(i);};

// HLT
muon_pog::HLT* CMSAnalysis::phlt() {
     if (_pevent==0) _bEvent->GetEntry(_iEvent);
      return (&(_pevent->hlt));
}
muon_pog::HLT& CMSAnalysis::hlt(){return *phlt();};

// MET
muon_pog::METs* CMSAnalysis::pmet() {
    if (_pevent==0) _bEvent->GetEntry(_iEvent);
      return (&(_pevent->mets));
}
muon_pog::METs& CMSAnalysis::met(){return *pmet();};

// GENPARTICLES
//unsigned int CMSAnalysis::nGenParticles() {
//      if (_pevent==0) _bEvent->GetEntry(_iEvent); 
//      if (!(_pevent->isMC)) return 0;
//      if (_pgenParticles==0 && _bGenParticles) _bGenParticles->GetEntry(_iEvent); 
//      return _pgenParticles->size();
//}
//std::vector<muon_pog::GenParticle> CMSAnalysis::genParticles() {
//      if (_pgenParticles==0) _bGenParticles->GetEntry(_iEvent); 
//      return (*_pgenParticles);
//}
//muon_pog::GenParticle* CMSAnalysis::pgenParticle(unsigned int i) {
//      if (nGenParticles()>i) {
//            if (_pgenParticles==0) _bGenParticles->GetEntry(_iEvent); 
//            return &(_pgenParticles->at(i));
//      } else {
//            return 0;
//      }
//}
//muon_pog::GenParticle& CMSAnalysis::genParticle(unsigned int i){return *pgenParticle(i);};

// MUONS
unsigned int CMSAnalysis::nMuons() {
      if (_pevent==0) _bEvent->GetEntry(_iEvent);
    return _pevent->muons.size();
//      return 1; ///NEW
}
std::vector<muon_pog::Muon> CMSAnalysis::muons() {
    if (_pevent==0) _bEvent->GetEntry(_iEvent);
    return (_pevent->muons);
}
muon_pog::Muon* CMSAnalysis::pmuon(unsigned int i) {
      if (nMuons()>i) {
          if (_pevent==0) _bEvent->GetEntry(_iEvent);
            return &(_pevent->muons.at(i));
      } else {
            return 0;
      }
}
muon_pog::Muon& CMSAnalysis::muon(unsigned int i){return *pmuon(i);};
