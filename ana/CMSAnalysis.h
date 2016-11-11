#include "TFile.h"
#include "TH1D.h"
#include "TTree.h"
#include "TString.h"

/// Complete analysis class to do CMS analysis using CMSTree2012 root files.
class CMSAnalysis {
  public:
      CMSAnalysis();
      virtual ~CMSAnalysis();

      /// Initialize and add data sample to the analysis (always first sample to add)
      void AddDataSample(const TString& id, const TString& file, double luminosity);
      /// Initialize and add the MC signal sample
      void AddMCSignalSample(const TString& id, const TString& file, double maxevents, double xsection, double total_events_for_xsection, double prescale);
      /// Initialize and add another MC background sample. All samples with the same "id" are added to the same histogram components at the end (stacked)
      void AddMCSample(const TString& id, const TString& file, double maxevents, double xsection, double total_events_for_xsection, double prescale);
      /// Initialize and add another MC background sample. All samples with the same "id" are added to the same histogram components at the end (stacked)
      void AddModel(const TString& id, const TString& file, double maxevents, double xsec, double total_events_for_xsection, double prescale, bool plot=false);
      /// Set tree to the sample with index i
      void SetTree(int i);

      /// Gets entry from the current tree; return 0 if OK, -1 if error or beyond limits
      Int_t GetEntry(int entry);

      /// True if the current file was declared to be real data
      bool isDataFile(){return _currentIndex==_dataIndex;};
      /// True if the current file was declared to be MC
      bool isMCFile(){return _currentIndex!=_dataIndex;};
      /// True if the current file was declared to be the signal MC
      bool isSignalFile(){return _signalFlag[_currentIndex]==true;};
      /// True if the current file was declared to be plotted
      bool isPlotFile(){return _plotFlag[_currentIndex]==true;};
      /// True if the current file was declared to be a model
      bool isModelFile(){return _modelFlag[_currentIndex]==true;};
      /// Get luminosity
      double GetLumi(){return _lumi;};
      /// Get the index to the data sample
      int GetDataIndex() {return _dataIndex;};
      /// Get pointer to current tree
      TTree* GetTree() {return _tree;};
      /// Get pointer to current file
      TFile* GetFile() {return _file;};
      /// Get the index of the sample currently in use
      int GetIndex() {return _currentIndex;};

      /// Get number of samples being processed
      int GetNumberOfSamples() {return _sampleId.size();};
      /// Get number of events to be read in sample with index i
      int GetNumberOfEvents(int i) {return _sampleNevents[i];};
      /// Get sampleId of sample with index i
      TString GetSampleId(int i) {return _sampleId[i];};
      /// Get file name of sample with index i
      TString GetSampleFile(int i) {return _sampleFile[i];};
      /// Get cross section for sample with index i (in data this gives number_of_events/luminosity)
      double GetSampleXsection(int i) {return _sampleXsection[i];};
      /// Get effective luminosity provided by sample with index i (in data this gives just the declared luminosity)
      double GetSampleEquivalentLuminosity(int i) {return _lumi/_sampleWeight[i];};
      /// Get sample weight for sample with index i (to normalize it to the total integrated luminosity)
      double GetSampleWeight(int i) {return _sampleWeight[i];};

      /// Book 1D Plot for data and all bkgd components. It follows TH1 conventions
      void AddPlot1D(const TString& name, const TString& title, int nbins, double xmin, double xmax);
      void AddPlot1DRB(const TString& name, const TString& title, int nbins, double xbins[31]);
      /// Fill histogram for the 1D Plot. It follows TH1 conventions.
      void FillPlot1D(const TString& name, int isample, double value, double weight=1.);
      void NormPlot1D(const TString& name, int nbins, double xbins[150], TString ytitle);
      /// Draw 1D plot with data and all bckg components. It follows TH1 conventions. It also produces .pdf, .png and .root versions of the histogram, with an optional suffix (to avoid overwriting other plots)
      //void DrawPlot1D(const TString& name, const TString& suffix="");
      void DrawPlot1D(const TString& folder, const TString& name, const TString& suffix, bool mc_study , const TString& title = "CMS preeliminary" , const TString& type = "normal");
      /// Book 2D Plot for data and all bkgd components. It follows TH1 conventions
      void AddPlot2D(const TString& name, const TString& title, int nbinsx, double xmin, double xmax, int nbinsy, double ymin, double ymax);
      /// Fill histogram for the 2D Plot. It follows TH1 conventions.
      void FillPlot2D(const TString& name, int isample, double xvalue, double yvalue, double weight=1.);
      /// Draw 2D plot with data and all bckg components. It follows TH1 conventions. It also produces .pdf, .png and .root versions of the histogram, with an optional suffix (to avoid overwriting other plots)
      void DrawPlot2D(const TString& name, const TString& suffix="");

      /// Get pointer to global information in current event
      muon_pog::Event* pevent();
      /// Get global information in current event
      muon_pog::Event& event();

      /// Get pointer to GenInfo information in current event
      muon_pog::GenInfo* pgenInfo();
      /// Get GenInfo information in current event
      muon_pog::GenInfo& genInfo();

      /// Get pointer to HLT information in current event
      muon_pog::HLT* phlt();
      /// Get HLT information in current event
      muon_pog::HLT& hlt();

      /// Get pointer to METFilter information in current event
//      muon_pog::HLT* pMETFilter();
      /// Get HLT information in current event
//      muon_pog::HLT& METFilter();

      /// Get pointer to MET information in current event
      muon_pog::METs* pmet();
      /// Get MET information in current event
      muon_pog::METs& met();

      /// Get number of genParticles in current event
      unsigned int nGenParticles();
      /// Get vector of genParticles in current event
      std::vector<muon_pog::GenParticle> genParticles();
      /// Get pointer to genParticle i in current event
      muon_pog::GenParticle* pgenParticle(unsigned int i);
      /// Get genParticle i in current event
      muon_pog::GenParticle& genParticle(unsigned int i);

      /// Get number of muons in current event
      unsigned int nMuons();
      /// Get vector of muons in current event
      std::vector<muon_pog::Muon> muons();
      /// Get pointer to muon i in current event
      muon_pog::Muon* pmuon(unsigned int i);
      /// Get muon i in current event
      muon_pog::Muon& muon(unsigned int i);

//      /// Get number of electrons in current event
//      unsigned int nElectrons();
//      /// Get vector of electrons in current event
//      std::vector<muon_pog::Electron> electrons();
//      /// Get pointer to electron i in current event
//      muon_pog::Electron* pelectron(unsigned int i);
//      /// Get electron i in current event
//      muon_pog::Electron& electron(unsigned int i);

//      /// Get number of jets in current event
//      unsigned int nJets();
//      /// Get vector of jets in current event
//      std::vector<muon_pog::Jet> jets();
//      /// Get pointer to jet i in current event
//      muon_pog::Jet* pjet(unsigned int i);
//      /// Get jet i in current event
//      muon_pog::Jet& jet(unsigned int i);

//      /// True if the current event is a real data event
//      bool isDataEvent(){return (!event().isMC);};
//      /// True if the current event is a MC event
//      bool isMCEvent(){return (event().isMC);};

  private:
      double _lumi; ///<
      std::vector<TString> _sampleFile; ///<
      std::vector<TString> _sampleId; ///<
      std::vector<int> _sampleNevents; ///<
      std::vector<int> _nDump10; ///<
      std::vector<double> _sampleXsection; ///<
      std::vector<double> _sampleWeight; ///<
      std::vector<bool> _signalFlag; ///<
      std::vector<bool> _plotFlag; ///<
      std::vector<bool> _modelFlag; ///<

      int _dataIndex; ///<

      int _currentIndex; ///<
      TFile* _file; ///<
      TTree* _tree; ///<

      TBranch* _bEvent; ///<
//      TBranch* _bGenInfo; ///<
//      TBranch* _bGenParticles; ///<
//      TBranch* _bMuons; ///<
//      TBranch* _bElectrons; ///<
//      TBranch* _bJets; ///<
//      TBranch* _bHLT; ///<
//      TBranch* _bMETFilter; ///<
//      TBranch* _bMET; ///<

      Int_t _iEvent; ///<
      muon_pog::Event* _pevent; ///<
      muon_pog::GenInfo* _pgenInfo; ///<
      std::vector<muon_pog::GenParticle>* _pgenParticles; ///<
      std::vector<muon_pog::Muon>* _pmuons; ///<
//      std::vector<muon_pog::Electron>* _pelectrons; ///<
//      std::vector<muon_pog::Jet>* _pjets; ///<
      muon_pog::HLT* _phlt; ///<
//      muon_pog::HLT* _pMETFilter; ///<
      muon_pog::METs* _pmet; ///<

      std::vector<TH1D*> hists_1D; ///<
};
