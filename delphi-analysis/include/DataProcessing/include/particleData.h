#ifndef PARTICLEDATA_H
#define PARTICLEDATA_H

//cpp dependencies
#include <string>
#include <vector>

//ROOT dependencies
#include "TTree.h"

//Local dependencies
#include "dataTools.h"

#include "eventData.h"
struct SOURCE{
  enum MCSOURCE{data=0, bbmc, qqmc, eeee=31, eemm, eeuu, eess, eecc, mumu=4, tautau, 
        lep1Data=60, lep1MC=61, lep2Data=70, lep2MC=71, 
        KQQ=8,
        GGBB=91, GGCC, GGSS, GGTT, GGUD, 
        TT=100, 
        KWW4F=110, KWENU, 
        PZZ=120, ZNN, PZEE,
        undefined=-1};
};

class particleData{
 public:
  static const int nMaxPart = 650;
  bool initMinimal;

  int nParticle;
  int EventNo;
  int RunNo;
  int year;
  int subDir;
  int process;
  int source;
  bool isMC;
  bool isOnres;
  unsigned long long uniqueID;
  float Energy;
  int bFlag;
  float particleWeight;
  float bx;
  float by;
  float ebx;
  float eby;
  float px[nMaxPart];
  float py[nMaxPart];
  float pz[nMaxPart];
  float pt[nMaxPart];
  float pmag[nMaxPart];
  float rap[nMaxPart];
  float eta[nMaxPart];
  float theta[nMaxPart];
  float phi[nMaxPart];
  float mass[nMaxPart];
  Short_t charge[nMaxPart];
  // Starting from 0, pwflag (via Marcello) - CHARGED_TRACK, CHARGED_LEPTONS1, CHARGED_LEPTONS2, V0, PHOTON, NEUTRAL_HADRON
  Short_t pwflag[nMaxPart];
  int pid[nMaxPart];
  float d0[nMaxPart];
  float z0[nMaxPart];
  Bool_t highPurity[nMaxPart];
  Short_t ntpc[nMaxPart];
  Short_t nitc[nMaxPart];
  Short_t nvdet[nMaxPart];
  float vx[nMaxPart];
  float vy[nMaxPart];
  float vz[nMaxPart];
  float weight[nMaxPart];

  //thrust axis variables
  float pt_wrtThrMissP[nMaxPart];
  float eta_wrtThrMissP[nMaxPart];
  float rap_wrtThrMissP[nMaxPart];
  float theta_wrtThrMissP[nMaxPart];
  float phi_wrtThrMissP[nMaxPart];

  float * pt_wrtThr;
  float * eta_wrtThr;
  float * rap_wrtThr;
  float * theta_wrtThr;
  float * phi_wrtThr;
  float * pt_wrtThrPerp;
  float * eta_wrtThrPerp;
  float * rap_wrtThrPerp;
  float * theta_wrtThrPerp;
  float * phi_wrtThrPerp;
  float * pt_wrtChThr;
  float * eta_wrtChThr;
  float * rap_wrtChThr;
  float * theta_wrtChThr;
  float * phi_wrtChThr;
  float * pt_wrtChThrPerp;
  float * eta_wrtChThrPerp;
  float * rap_wrtChThrPerp;
  float * theta_wrtChThrPerp;
  float * phi_wrtChThrPerp;
  float * pt_wrtNeuThr;
  float * eta_wrtNeuThr;
  float * rap_wrtNeuThr;
  float * theta_wrtNeuThr;
  float * phi_wrtNeuThr;
  float * pt_wrtNeuThrPerp;
  float * eta_wrtNeuThrPerp;
  float * rap_wrtNeuThrPerp;
  float * theta_wrtNeuThrPerp;
  float * phi_wrtNeuThrPerp;
  float * pt_wrtThrCorr;
  float * eta_wrtThrCorr;
  float * rap_wrtThrCorr;
  float * theta_wrtThrCorr;
  float * phi_wrtThrCorr;
  float * pt_wrtThrCorrPerp;
  float * eta_wrtThrCorrPerp;
  float * rap_wrtThrCorrPerp;
  float * theta_wrtThrCorrPerp;
  float * phi_wrtThrCorrPerp;
  float * pt_wrtThrCorrInverse;
  float * eta_wrtThrCorrInverse;
  float * rap_wrtThrCorrInverse;
  float * theta_wrtThrCorrInverse;
  float * phi_wrtThrCorrInverse;
  float * pt_wrtThrCorrInversePerp;
  float * eta_wrtThrCorrInversePerp;
  float * rap_wrtThrCorrInversePerp;
  float * theta_wrtThrCorrInversePerp;
  float * phi_wrtThrCorrInversePerp;
  float * pt_wrtThrMissPPerp;
  float * eta_wrtThrMissPPerp;
  float * rap_wrtThrMissPPerp;
  float * theta_wrtThrMissPPerp;
  float * phi_wrtThrMissPPerp;


  // artificial acceptance correction
  Bool_t passesArtificAccept[nMaxPart];
  float artificAcceptEffCorrection[nMaxPart];
  int anc[nMaxPart];
  int index[nMaxPart];
  int correspondenceIndex[nMaxPart];
    
  static const int nVar = 105;
  std::string varStr[nVar] = {"nParticle",
			      "EventNo",
			      "RunNo",
			      "year",
			      "subDir",
			      "process",
			      "source",
			      "isMC",
			      "isOnres",
			      "uniqueID",
			      "Energy",
			      "bFlag",
                              "particleWeight",
			      "bx",
			      "by",
			      "ebx",
			      "eby",
			      "px",
			      "py",
			      "pz",
			      "pt",
			      "pmag",
			      "rap",
			      "eta",
			      "theta",
			      "phi",
			      "mass",
			      "charge",
			      "pwflag",
			      "pid",
			      "d0",
			      "z0",
                              "highPurity",
			      "ntpc",
			      "nitc",
			      "nvdet",
			      "vx",
			      "vy",
			      "vz",
			      "weight",
			      "pt_wrtThr",
			      "eta_wrtThr",
			      "rap_wrtThr",
			      "theta_wrtThr",
			      "phi_wrtThr",
			      "pt_wrtThrPerp",
			      "eta_wrtThrPerp",
			      "rap_wrtThrPerp",
			      "theta_wrtThrPerp",
			      "phi_wrtThrPerp",
			      "pt_wrtChThr",
			      "eta_wrtChThr",
			      "rap_wrtChThr",
			      "theta_wrtChThr",
			      "phi_wrtChThr",
			      "pt_wrtChThrPerp",
			      "eta_wrtChThrPerp",
			      "rap_wrtChThrPerp",
			      "theta_wrtChThrPerp",
			      "phi_wrtChThrPerp",
			      "pt_wrtNeuThr",
			      "eta_wrtNeuThr",
			      "rap_wrtNeuThr",
			      "theta_wrtNeuThr",
			      "phi_wrtNeuThr",
			      "pt_wrtNeuThrPerp",
			      "eta_wrtNeuThrPerp",
			      "theta_wrtNeuThrPerp",
			      "phi_wrtNeuThrPerp",
			      "rap_wrtNeuThrPerp",
			      "pt_wrtThrCorr",
			      "eta_wrtThrCorr",
			      "rap_wrtThrCorr",
			      "theta_wrtThrCorr",
			      "phi_wrtThrCorr",
			      "pt_wrtThrCorrPerp",
			      "eta_wrtThrCorrPerp",
			      "rap_wrtThrCorrPerp",
			      "theta_wrtThrCorrPerp",
			      "phi_wrtThrCorrPerp",
			      "pt_wrtThrCorrInverse",
			      "eta_wrtThrCorrInverse",
			      "rap_wrtThrCorrInverse",
			      "theta_wrtThrCorrInverse",
			      "phi_wrtThrCorrInverse",
			      "pt_wrtThrCorrInversePerp",
			      "eta_wrtThrCorrInversePerp",
			      "rap_wrtThrCorrInversePerp",
			      "theta_wrtThrCorrInversePerp",
			      "phi_wrtThrCorrInversePerp",
			      "pt_wrtThrMissP",
			      "eta_wrtThrMissP",
			      "rap_wrtThrMissP",
			      "theta_wrtThrMissP",
			      "phi_wrtThrMissP",
			      "pt_wrtThrMissPPerp",
			      "eta_wrtThrMissPPerp",
			      "rap_wrtThrMissPPerp",
			      "theta_wrtThrMissPPerp",
			      "phi_wrtThrMissPPerp",
			      "passesArtificAccept",
			      "artificAcceptEffCorrection",
			      "anc",
			      "index",
			      "correspondenceIndex"};
  
  bool varIsGood[nVar];
  
  particleData(bool in_initMinimal=true);
  ~particleData();
  void SetStatusAndAddressRead(TTree* inTree_p, std::vector<std::string> inList);
  void SetBranchWrite(TTree* inTree_p, bool writeMinimal=0);
  void preFillClean();

  // quick particleData placeholder template 
  // If you want to do so, please change the definition of varStr first, and the do _init() -> _setRead() -> _setWrite() -> _setReducedPrecision(),
  // you can get the correspond codes instantaneously, but recommend the user should check it carefully, this is not seriously tested.
  void _init();
  void _setRead();
  void _setWrite();
  void _setReducedPrecision();
};

inline particleData::particleData(bool in_initMinimal) : initMinimal(in_initMinimal)
{
  //init values in event they are written to tree but are meaningless - treat -999 as danger
  nParticle = -999;
  EventNo = -999;
  RunNo = -999;
  year = -999;
  subDir = -999;
  process = -999;
  source = -999;
  isMC = false;
  isOnres = false;
  uniqueID = 0;
  Energy = 0;
  bFlag = -999;
  particleWeight = 1;//keep this default as 1 (only used for mixed events, but should be 1 otherwise)
  bx = -999;
  by = -999;
  ebx = -999;
  eby = -999;

  if (!initMinimal)
  {
    std::cout << "Init full branches ..." << std::endl;
    pt_wrtThr = new float[nMaxPart];
    eta_wrtThr = new float[nMaxPart];
    rap_wrtThr = new float[nMaxPart];
    theta_wrtThr = new float[nMaxPart];
    phi_wrtThr = new float[nMaxPart];
    pt_wrtThrPerp = new float[nMaxPart];
    eta_wrtThrPerp = new float[nMaxPart];
    rap_wrtThrPerp = new float[nMaxPart];
    theta_wrtThrPerp = new float[nMaxPart];
    phi_wrtThrPerp = new float[nMaxPart];
    pt_wrtChThr = new float[nMaxPart];
    eta_wrtChThr = new float[nMaxPart];
    rap_wrtChThr = new float[nMaxPart];
    theta_wrtChThr = new float[nMaxPart];
    phi_wrtChThr = new float[nMaxPart];
    pt_wrtChThrPerp = new float[nMaxPart];
    eta_wrtChThrPerp = new float[nMaxPart];
    rap_wrtChThrPerp = new float[nMaxPart];
    theta_wrtChThrPerp = new float[nMaxPart];
    phi_wrtChThrPerp = new float[nMaxPart];
    pt_wrtNeuThr = new float[nMaxPart];
    eta_wrtNeuThr = new float[nMaxPart];
    rap_wrtNeuThr = new float[nMaxPart];
    theta_wrtNeuThr = new float[nMaxPart];
    phi_wrtNeuThr = new float[nMaxPart];
    pt_wrtNeuThrPerp = new float[nMaxPart];
    eta_wrtNeuThrPerp = new float[nMaxPart];
    rap_wrtNeuThrPerp = new float[nMaxPart];
    theta_wrtNeuThrPerp = new float[nMaxPart];
    phi_wrtNeuThrPerp = new float[nMaxPart];
    pt_wrtThrCorr = new float[nMaxPart];
    eta_wrtThrCorr = new float[nMaxPart];
    rap_wrtThrCorr = new float[nMaxPart];
    theta_wrtThrCorr = new float[nMaxPart];
    phi_wrtThrCorr = new float[nMaxPart];
    pt_wrtThrCorrPerp = new float[nMaxPart];
    eta_wrtThrCorrPerp = new float[nMaxPart];
    rap_wrtThrCorrPerp = new float[nMaxPart];
    theta_wrtThrCorrPerp = new float[nMaxPart];
    phi_wrtThrCorrPerp = new float[nMaxPart];
    pt_wrtThrCorrInverse = new float[nMaxPart];
    eta_wrtThrCorrInverse = new float[nMaxPart];
    rap_wrtThrCorrInverse = new float[nMaxPart];
    theta_wrtThrCorrInverse = new float[nMaxPart];
    phi_wrtThrCorrInverse = new float[nMaxPart];
    pt_wrtThrCorrInversePerp = new float[nMaxPart];
    eta_wrtThrCorrInversePerp = new float[nMaxPart];
    rap_wrtThrCorrInversePerp = new float[nMaxPart];
    theta_wrtThrCorrInversePerp = new float[nMaxPart];
    phi_wrtThrCorrInversePerp = new float[nMaxPart];
    pt_wrtThrMissPPerp = new float[nMaxPart];
    eta_wrtThrMissPPerp = new float[nMaxPart];
    rap_wrtThrMissPPerp = new float[nMaxPart];
    theta_wrtThrMissPPerp = new float[nMaxPart];
    phi_wrtThrMissPPerp = new float[nMaxPart];
  }

  for(Int_t i = 0; i < nMaxPart; ++i){
    px[i] = -999;
    py[i] = -999;
    pz[i] = -999;
    pt[i] = -999;
    pmag[i] = -999;
    rap[i] = -999;
    eta[i] = -999;
    theta[i] = -999;
    phi[i] = -999;
    mass[i] = -999;
    charge[i] = -127;
    pwflag[i] = -127;
    pid[i] = -999;
    d0[i] = -999;
    z0[i] = -999;
    highPurity[i] = true;
    ntpc[i] = -127;
    nitc[i] = -127;
    nvdet[i] = -127;
    vx[i] = -999.;
    vy[i] = -999.;
    vz[i] = -999.;
    weight[i] = -999.;
    index[i] = -999;
    correspondenceIndex[i] = -999;
    
    if (!initMinimal)
    {
      pt_wrtThr[i] = -999;
      eta_wrtThr[i] = -999;
      rap_wrtThr[i] = -999;
      theta_wrtThr[i] = -999;
      phi_wrtThr[i] = -999;
      pt_wrtThrPerp[i] = -999;
      eta_wrtThrPerp[i] = -999;
      rap_wrtThrPerp[i] = -999;
      theta_wrtThrPerp[i] = -999;
      phi_wrtThrPerp[i] = -999;
      pt_wrtChThr[i] = -999;
      eta_wrtChThr[i] = -999;
      rap_wrtChThr[i] = -999;
      theta_wrtChThr[i] = -999;
      phi_wrtChThr[i] = -999;
      pt_wrtChThrPerp[i] = -999;
      eta_wrtChThrPerp[i] = -999;
      rap_wrtChThrPerp[i] = -999;
      theta_wrtChThrPerp[i] = -999;
      phi_wrtChThrPerp[i] = -999;
      pt_wrtNeuThr[i] = -999;
      eta_wrtNeuThr[i] = -999;
      rap_wrtNeuThr[i] = -999;
      theta_wrtNeuThr[i] = -999;
      phi_wrtNeuThr[i] = -999;
      pt_wrtNeuThrPerp[i] = -999;
      eta_wrtNeuThrPerp[i] = -999;
      theta_wrtNeuThrPerp[i] = -999;
      phi_wrtNeuThrPerp[i] = -999;
      rap_wrtNeuThrPerp[i] = -999;
      pt_wrtThrCorr[i] = -999;
      eta_wrtThrCorr[i] = -999;
      rap_wrtThrCorr[i] = -999;
      theta_wrtThrCorr[i] = -999;
      phi_wrtThrCorr[i] = -999;
      pt_wrtThrCorrPerp[i] = -999;
      eta_wrtThrCorrPerp[i] = -999;
      rap_wrtThrCorrPerp[i] = -999;
      theta_wrtThrCorrPerp[i] = -999;
      phi_wrtThrCorrPerp[i] = -999;
      pt_wrtThrCorrInverse[i] = -999;
      eta_wrtThrCorrInverse[i] = -999;
      rap_wrtThrCorrInverse[i] = -999;
      theta_wrtThrCorrInverse[i] = -999;
      phi_wrtThrCorrInverse[i] = -999;
      pt_wrtThrCorrInversePerp[i] = -999;
      eta_wrtThrCorrInversePerp[i] = -999;
      rap_wrtThrCorrInversePerp[i] = -999;
      theta_wrtThrCorrInversePerp[i] = -999;
      phi_wrtThrCorrInversePerp[i] = -999; 
      pt_wrtThrMissPPerp[i] = -999;
      eta_wrtThrMissPPerp[i] = -999;
      rap_wrtThrMissPPerp[i] = -999;
      theta_wrtThrMissPPerp[i] = -999;
      phi_wrtThrMissPPerp[i] = -999;
    }

    pt_wrtThrMissP[i] = -999;
    eta_wrtThrMissP[i] = -999;
    rap_wrtThrMissP[i] = -999;
    theta_wrtThrMissP[i] = -999;
    phi_wrtThrMissP[i] = -999;

      
    passesArtificAccept[i] = true;
    artificAcceptEffCorrection[i] = -999;
    anc[i] = -999;
  }
  
  for(int i = 0; i < nVar; ++i)
  {
    varIsGood[i] = (initMinimal && 
                    (  varStr[i].find("_wrt")==std::string::npos ||
                      (varStr[i].find("_wrtThrMissP")!=std::string::npos &&
                       varStr[i].find("_wrtThrMissPPerp")==std::string::npos )
                     ) )? true: false;
    // printf("[%d] %s: %o\n", i, varStr[i].c_str(), varIsGood[i]);
  }
  return;
}


inline particleData::~particleData()
{
  if (!initMinimal)
  {
    std::cout << "Delete branches ..." << std::endl;
    delete [] pt_wrtThr;
    delete [] eta_wrtThr;
    delete [] rap_wrtThr;
    delete [] theta_wrtThr;
    delete [] phi_wrtThr;
    delete [] pt_wrtThrPerp;
    delete [] eta_wrtThrPerp;
    delete [] rap_wrtThrPerp;
    delete [] theta_wrtThrPerp;
    delete [] phi_wrtThrPerp;
    delete [] pt_wrtChThr;
    delete [] eta_wrtChThr;
    delete [] rap_wrtChThr;
    delete [] theta_wrtChThr;
    delete [] phi_wrtChThr;
    delete [] pt_wrtChThrPerp;
    delete [] eta_wrtChThrPerp;
    delete [] rap_wrtChThrPerp;
    delete [] theta_wrtChThrPerp;
    delete [] phi_wrtChThrPerp;
    delete [] pt_wrtNeuThr;
    delete [] eta_wrtNeuThr;
    delete [] rap_wrtNeuThr;
    delete [] theta_wrtNeuThr;
    delete [] phi_wrtNeuThr;
    delete [] pt_wrtNeuThrPerp;
    delete [] eta_wrtNeuThrPerp;
    delete [] rap_wrtNeuThrPerp;
    delete [] theta_wrtNeuThrPerp;
    delete [] phi_wrtNeuThrPerp;
    delete [] pt_wrtThrCorr;
    delete [] eta_wrtThrCorr;
    delete [] rap_wrtThrCorr;
    delete [] theta_wrtThrCorr;
    delete [] phi_wrtThrCorr;
    delete [] pt_wrtThrCorrPerp;
    delete [] eta_wrtThrCorrPerp;
    delete [] rap_wrtThrCorrPerp;
    delete [] theta_wrtThrCorrPerp;
    delete [] phi_wrtThrCorrPerp;
    delete [] pt_wrtThrCorrInverse;
    delete [] eta_wrtThrCorrInverse;
    delete [] rap_wrtThrCorrInverse;
    delete [] theta_wrtThrCorrInverse;
    delete [] phi_wrtThrCorrInverse;
    delete [] pt_wrtThrCorrInversePerp;
    delete [] eta_wrtThrCorrInversePerp;
    delete [] rap_wrtThrCorrInversePerp;
    delete [] theta_wrtThrCorrInversePerp;
    delete [] phi_wrtThrCorrInversePerp;
    delete [] pt_wrtThrMissPPerp;
    delete [] eta_wrtThrMissPPerp;
    delete [] rap_wrtThrMissPPerp;
    delete [] theta_wrtThrMissPPerp;
    delete [] phi_wrtThrMissPPerp;
  }
}

inline void particleData::SetStatusAndAddressRead(TTree* inTree_p, std::vector<std::string> inList)
{
  if(inList.size() != 0)
  {
    for(int i = 0; i < nVar; ++i){varIsGood[i] = false;}

    for(unsigned int i = 0; i < inList.size(); ++i)
    {
      for(Int_t j = 0; j < nVar; ++j)
      {
      	if(inList.at(i).size() == varStr[j].size() && inList.at(i).find(varStr[j]) != std::string::npos)
        {
          varIsGood[j] = (initMinimal && 
                    (  varStr[j].find("_wrt")==std::string::npos ||
                      (varStr[j].find("_wrtThrMissP")!=std::string::npos &&
                       varStr[j].find("_wrtThrMissPPerp")==std::string::npos )
                     ) )? true: false;
          // printf("[%d] %s: %o\n", j, varStr[j].c_str(), varIsGood[j]);
      	  break;
      	}
      }
    }
  }

  for(Int_t i = 0; i < nVar; ++i){
    if(varIsGood[i]) {
      inTree_p->SetBranchStatus(varStr[i].c_str(), 1);
      // printf("SetBranchStatus %s\n", varStr[i].c_str());
    }
  }

  if(varIsGood[0]) inTree_p->SetBranchAddress("nParticle", &nParticle); 
  if(varIsGood[1]) inTree_p->SetBranchAddress("EventNo", &EventNo);
  if(varIsGood[2]) inTree_p->SetBranchAddress("RunNo", &RunNo);
  if(varIsGood[3]) inTree_p->SetBranchAddress("year", &year);
  if(varIsGood[4]) inTree_p->SetBranchAddress("subDir", &subDir);
  if(varIsGood[5]) inTree_p->SetBranchAddress("process", &process);
  if(varIsGood[6]) inTree_p->SetBranchAddress("source", &source);
  if(varIsGood[7]) inTree_p->SetBranchAddress("isMC", &isMC);
  if(varIsGood[8]) inTree_p->SetBranchAddress("isOnres", &isOnres);
  if(varIsGood[9]) inTree_p->SetBranchAddress("uniqueID", &uniqueID);
  if(varIsGood[10]) inTree_p->SetBranchAddress("Energy", &Energy);
  if(varIsGood[11]) inTree_p->SetBranchAddress("bFlag", &bFlag);
  if(varIsGood[12]) inTree_p->SetBranchAddress("particleWeight", &particleWeight);
  if(varIsGood[13]) inTree_p->SetBranchAddress("bx", &bx);
  if(varIsGood[14]) inTree_p->SetBranchAddress("by", &by);
  if(varIsGood[15]) inTree_p->SetBranchAddress("ebx", &ebx);
  if(varIsGood[16]) inTree_p->SetBranchAddress("eby", &eby);
  if(varIsGood[17]) inTree_p->SetBranchAddress("px", px);
  if(varIsGood[18]) inTree_p->SetBranchAddress("py", py);
  if(varIsGood[19]) inTree_p->SetBranchAddress("pz", pz);
  if(varIsGood[20]) inTree_p->SetBranchAddress("pt", pt);
  if(varIsGood[21]) inTree_p->SetBranchAddress("pmag", pmag);
  if(varIsGood[22]) inTree_p->SetBranchAddress("rap", rap);
  if(varIsGood[23]) inTree_p->SetBranchAddress("eta", eta);
  if(varIsGood[24]) inTree_p->SetBranchAddress("theta", theta);
  if(varIsGood[25]) inTree_p->SetBranchAddress("phi", phi);
  if(varIsGood[26]) inTree_p->SetBranchAddress("mass", mass);
  if(varIsGood[27]) inTree_p->SetBranchAddress("charge", charge);
  if(varIsGood[28]) inTree_p->SetBranchAddress("pwflag", pwflag);
  if(varIsGood[29]) inTree_p->SetBranchAddress("pid", pid);
  if(varIsGood[30]) inTree_p->SetBranchAddress("d0", d0);
  if(varIsGood[31]) inTree_p->SetBranchAddress("z0", z0);
  if(varIsGood[32]) inTree_p->SetBranchAddress("highPurity", highPurity);
  if(varIsGood[33]) inTree_p->SetBranchAddress("ntpc", ntpc);
  if(varIsGood[34]) inTree_p->SetBranchAddress("nitc", nitc);
  if(varIsGood[35]) inTree_p->SetBranchAddress("nvdet", nvdet);
  if(varIsGood[36]) inTree_p->SetBranchAddress("vx", vx);
  if(varIsGood[37]) inTree_p->SetBranchAddress("vy", vy);
  if(varIsGood[38]) inTree_p->SetBranchAddress("vz", vz);
  if(varIsGood[39]) inTree_p->SetBranchAddress("weight", weight);
  if(varIsGood[40]) inTree_p->SetBranchAddress("pt_wrtThr", pt_wrtThr);
  if(varIsGood[41]) inTree_p->SetBranchAddress("eta_wrtThr", eta_wrtThr);
  if(varIsGood[42]) inTree_p->SetBranchAddress("rap_wrtThr", rap_wrtThr);
  if(varIsGood[43]) inTree_p->SetBranchAddress("theta_wrtThr", theta_wrtThr);
  if(varIsGood[44]) inTree_p->SetBranchAddress("phi_wrtThr", phi_wrtThr);
  if(varIsGood[45]) inTree_p->SetBranchAddress("pt_wrtThrPerp", pt_wrtThrPerp);
  if(varIsGood[46]) inTree_p->SetBranchAddress("eta_wrtThrPerp", eta_wrtThrPerp);
  if(varIsGood[47]) inTree_p->SetBranchAddress("rap_wrtThrPerp", rap_wrtThrPerp);
  if(varIsGood[48]) inTree_p->SetBranchAddress("theta_wrtThrPerp", theta_wrtThrPerp);
  if(varIsGood[49]) inTree_p->SetBranchAddress("phi_wrtThrPerp", phi_wrtThrPerp);
  if(varIsGood[50]) inTree_p->SetBranchAddress("pt_wrtChThr", pt_wrtChThr);
  if(varIsGood[51]) inTree_p->SetBranchAddress("eta_wrtChThr", eta_wrtChThr);
  if(varIsGood[52]) inTree_p->SetBranchAddress("rap_wrtChThr", rap_wrtChThr);
  if(varIsGood[53]) inTree_p->SetBranchAddress("theta_wrtChThr", theta_wrtChThr);
  if(varIsGood[54]) inTree_p->SetBranchAddress("phi_wrtChThr", phi_wrtChThr);
  if(varIsGood[55]) inTree_p->SetBranchAddress("pt_wrtChThrPerp", pt_wrtChThrPerp);
  if(varIsGood[56]) inTree_p->SetBranchAddress("eta_wrtChThrPerp", eta_wrtChThrPerp);
  if(varIsGood[57]) inTree_p->SetBranchAddress("phi_wrtChThrPerp", phi_wrtChThrPerp);
  if(varIsGood[58]) inTree_p->SetBranchAddress("rap_wrtChThrPerp", rap_wrtChThrPerp);
  if(varIsGood[59]) inTree_p->SetBranchAddress("theta_wrtChThrPerp", theta_wrtChThrPerp);
  if(varIsGood[60]) inTree_p->SetBranchAddress("pt_wrtNeuThr", pt_wrtNeuThr);
  if(varIsGood[61]) inTree_p->SetBranchAddress("eta_wrtNeuThr", eta_wrtNeuThr);
  if(varIsGood[62]) inTree_p->SetBranchAddress("rap_wrtNeuThr", rap_wrtNeuThr);
  if(varIsGood[63]) inTree_p->SetBranchAddress("theta_wrtNeuThr", theta_wrtNeuThr);
  if(varIsGood[64]) inTree_p->SetBranchAddress("phi_wrtNeuThr", phi_wrtNeuThr);
  if(varIsGood[65]) inTree_p->SetBranchAddress("pt_wrtNeuThrPerp", pt_wrtNeuThrPerp);
  if(varIsGood[66]) inTree_p->SetBranchAddress("eta_wrtNeuThrPerp", eta_wrtNeuThrPerp);
  if(varIsGood[67]) inTree_p->SetBranchAddress("theta_wrtNeuThrPerp", theta_wrtNeuThrPerp);
  if(varIsGood[68]) inTree_p->SetBranchAddress("phi_wrtNeuThrPerp", phi_wrtNeuThrPerp);
  if(varIsGood[69]) inTree_p->SetBranchAddress("rap_wrtNeuThrPerp", rap_wrtNeuThrPerp);
  if(varIsGood[70]) inTree_p->SetBranchAddress("pt_wrtThrCorr", pt_wrtThrCorr);
  if(varIsGood[71]) inTree_p->SetBranchAddress("eta_wrtThrCorr", eta_wrtThrCorr);
  if(varIsGood[72]) inTree_p->SetBranchAddress("rap_wrtThrCorr", rap_wrtThrCorr);
  if(varIsGood[73]) inTree_p->SetBranchAddress("theta_wrtThrCorr", theta_wrtThrCorr);
  if(varIsGood[74]) inTree_p->SetBranchAddress("phi_wrtThrCorr", phi_wrtThrCorr);
  if(varIsGood[75]) inTree_p->SetBranchAddress("pt_wrtThrCorrPerp", pt_wrtThrCorrPerp);
  if(varIsGood[76]) inTree_p->SetBranchAddress("eta_wrtThrCorrPerp", eta_wrtThrCorrPerp);
  if(varIsGood[77]) inTree_p->SetBranchAddress("rap_wrtThrCorrPerp", rap_wrtThrCorrPerp);
  if(varIsGood[78]) inTree_p->SetBranchAddress("theta_wrtThrCorrPerp", theta_wrtThrCorrPerp);
  if(varIsGood[79]) inTree_p->SetBranchAddress("phi_wrtThrCorrPerp", phi_wrtThrCorrPerp);
  if(varIsGood[80]) inTree_p->SetBranchAddress("pt_wrtThrCorrInverse", pt_wrtThrCorrInverse);
  if(varIsGood[81]) inTree_p->SetBranchAddress("eta_wrtThrCorrInverse", eta_wrtThrCorrInverse);
  if(varIsGood[82]) inTree_p->SetBranchAddress("rap_wrtThrCorrInverse", rap_wrtThrCorrInverse);
  if(varIsGood[83]) inTree_p->SetBranchAddress("theta_wrtThrCorrInverse", theta_wrtThrCorrInverse);
  if(varIsGood[84]) inTree_p->SetBranchAddress("phi_wrtThrCorrInverse", phi_wrtThrCorrInverse);
  if(varIsGood[85]) inTree_p->SetBranchAddress("pt_wrtThrCorrInversePerp", pt_wrtThrCorrInversePerp);
  if(varIsGood[86]) inTree_p->SetBranchAddress("eta_wrtThrCorrInversePerp", eta_wrtThrCorrInversePerp);
  if(varIsGood[87]) inTree_p->SetBranchAddress("rap_wrtThrCorrInversePerp", rap_wrtThrCorrInversePerp);
  if(varIsGood[88]) inTree_p->SetBranchAddress("theta_wrtThrCorrInversePerp", theta_wrtThrCorrInversePerp);
  if(varIsGood[89]) inTree_p->SetBranchAddress("phi_wrtThrCorrInversePerp", phi_wrtThrCorrInversePerp);
  if(varIsGood[90]) inTree_p->SetBranchAddress("pt_wrtThrMissP", pt_wrtThrMissP);
  if(varIsGood[91]) inTree_p->SetBranchAddress("eta_wrtThrMissP", eta_wrtThrMissP);
  if(varIsGood[92]) inTree_p->SetBranchAddress("rap_wrtThrMissP", rap_wrtThrMissP);
  if(varIsGood[93]) inTree_p->SetBranchAddress("theta_wrtThrMissP", theta_wrtThrMissP);
  if(varIsGood[94]) inTree_p->SetBranchAddress("phi_wrtThrMissP", phi_wrtThrMissP);
  if(varIsGood[95]) inTree_p->SetBranchAddress("pt_wrtThrMissPPerp", pt_wrtThrMissPPerp);
  if(varIsGood[96]) inTree_p->SetBranchAddress("eta_wrtThrMissPPerp", eta_wrtThrMissPPerp);
  if(varIsGood[97]) inTree_p->SetBranchAddress("rap_wrtThrMissPPerp", rap_wrtThrMissPPerp);
  if(varIsGood[98]) inTree_p->SetBranchAddress("theta_wrtThrMissPPerp", theta_wrtThrMissPPerp);
  if(varIsGood[99]) inTree_p->SetBranchAddress("phi_wrtThrMissPPerp", phi_wrtThrMissPPerp);
  if(varIsGood[100]) inTree_p->SetBranchAddress("passesArtificAccept", passesArtificAccept);
  if(varIsGood[101]) inTree_p->SetBranchAddress("artificAcceptEffCorrection", artificAcceptEffCorrection);
  if(varIsGood[102]) inTree_p->SetBranchAddress("anc", anc);
  if(varIsGood[103]) inTree_p->SetBranchAddress("index", index);
  if(varIsGood[104]) inTree_p->SetBranchAddress("correspondenceIndex", correspondenceIndex);
  
    return;
}

inline void particleData::SetBranchWrite(TTree* inTree_p, bool writeMinimal)
{
  writeMinimal = writeMinimal && initMinimal;
  inTree_p->Branch("EventNo", &EventNo, "EventNo/I");
  inTree_p->Branch("RunNo", &RunNo, "RunNo/I");
  inTree_p->Branch("source", &source, "source/I");
  inTree_p->Branch("isMC", &isMC, "isMC/O");
  inTree_p->Branch("isOnres", &isOnres, "isOnres/O");
  inTree_p->Branch("Energy", &Energy, "Energy/F");
  inTree_p->Branch("subDir", &subDir, "subDir/I");
  inTree_p->Branch("particleWeight", &particleWeight, "particleWeight/F");
  inTree_p->Branch("nParticle", &nParticle, "nParticle/I"); 
  inTree_p->Branch("pt", pt, "pt[nParticle]/F");
  inTree_p->Branch("rap", rap, "rap[nParticle]/F");
  inTree_p->Branch("eta", eta, "eta[nParticle]/F");
  inTree_p->Branch("theta", theta, "theta[nParticle]/F");
  inTree_p->Branch("phi", phi, "phi[nParticle]/F");
  inTree_p->Branch("charge", charge, "charge[nParticle]/S");
  inTree_p->Branch("pwflag", pwflag, "pwflag[nParticle]/S");
  inTree_p->Branch("pid", pid, "pid[nParticle]/I");
  inTree_p->Branch("highPurity", highPurity, "highPurity[nParticle]/O");
  inTree_p->Branch("px", px, "px[nParticle]/F");
  inTree_p->Branch("py", py, "py[nParticle]/F");
  inTree_p->Branch("pz", pz, "pz[nParticle]/F");
  inTree_p->Branch("mass", mass, "mass[nParticle]/F");
  inTree_p->Branch("ntpc", ntpc, "ntpc[nParticle]/S");
  inTree_p->Branch("d0", d0, "d0[nParticle]/F");
  inTree_p->Branch("z0", z0, "z0[nParticle]/F");
  inTree_p->Branch("weight", weight, "weight[nParticle]/F");
  inTree_p->Branch("index", index, "index[nParticle]/I");
  inTree_p->Branch("correspondenceIndex", correspondenceIndex, "correspondenceIndex[nParticle]/I");
  
  if (!writeMinimal)
  {
    inTree_p->Branch("year", &year, "year/I");
    inTree_p->Branch("process", &process, "process/I");
    inTree_p->Branch("uniqueID", &uniqueID, "uniqueID/l");
    inTree_p->Branch("bFlag", &bFlag, "bFlag/I");
    inTree_p->Branch("bx", &bx, "bx/F");
    inTree_p->Branch("by", &by, "by/F");
    inTree_p->Branch("ebx", &ebx, "ebx/F");
    inTree_p->Branch("eby", &eby, "eby/F");


    inTree_p->Branch("pmag", pmag, "pmag[nParticle]/F");
    inTree_p->Branch("nitc", nitc, "nitc[nParticle]/S");
    inTree_p->Branch("nvdet", nvdet, "nvdet[nParticle]/S");
    inTree_p->Branch("vx", vx, "vx[nParticle]/F");
    inTree_p->Branch("vy", vy, "vy[nParticle]/F");
    inTree_p->Branch("vz", vz, "vz[nParticle]/F");

    inTree_p->Branch("pt_wrtThrPerp", pt_wrtThrPerp, "pt_wrtThrPerp[nParticle]/F");
    inTree_p->Branch("eta_wrtThrPerp", eta_wrtThrPerp, "eta_wrtThrPerp[nParticle]/F");
    inTree_p->Branch("rap_wrtThrPerp", rap_wrtThrPerp, "rap_wrtThrPerp[nParticle]/F");
    inTree_p->Branch("theta_wrtThrPerp", theta_wrtThrPerp, "theta_wrtThrPerp[nParticle]/F");
    inTree_p->Branch("phi_wrtThrPerp", phi_wrtThrPerp, "phi_wrtThrPerp[nParticle]/F");

    inTree_p->Branch("passesArtificAccept", passesArtificAccept, "passesArtificAccept[nParticle]/O");
    inTree_p->Branch("artificAcceptEffCorrection", artificAcceptEffCorrection, "artificAcceptEffCorrection[nParticle]/F");
    inTree_p->Branch("anc", anc, "anc[nParticle]/I");
  }
  
  if (!writeMinimal && do_chThrust )           {
    inTree_p->Branch("pt_wrtThr", pt_wrtThr, "pt_wrtThr[nParticle]/F");
    inTree_p->Branch("eta_wrtThr", eta_wrtThr, "eta_wrtThr[nParticle]/F");
    inTree_p->Branch("rap_wrtThr", rap_wrtThr, "rap_wrtThr[nParticle]/F");
    inTree_p->Branch("theta_wrtThr", theta_wrtThr, "theta_wrtThr[nParticle]/F");
    inTree_p->Branch("phi_wrtThr", phi_wrtThr, "phi_wrtThr[nParticle]/F");
    inTree_p->Branch("pt_wrtChThr", pt_wrtChThr, "pt_wrtChThr[nParticle]/F");
    inTree_p->Branch("eta_wrtChThr", eta_wrtChThr, "eta_wrtChThr[nParticle]/F");
    inTree_p->Branch("rap_wrtChThr", rap_wrtChThr, "rap_wrtChThr[nParticle]/F");
    inTree_p->Branch("theta_wrtChThr", theta_wrtChThr, "theta_wrtChThr[nParticle]/F");
    inTree_p->Branch("phi_wrtChThr", phi_wrtChThr, "phi_wrtChThr[nParticle]/F");
    inTree_p->Branch("pt_wrtChThrPerp", pt_wrtChThrPerp, "pt_wrtChThrPerp[nParticle]/F");
    inTree_p->Branch("eta_wrtChThrPerp", eta_wrtChThrPerp, "eta_wrtChThrPerp[nParticle]/F");
    inTree_p->Branch("rap_wrtChThrPerp", rap_wrtChThrPerp, "rap_wrtChThrPerp[nParticle]/F");
    inTree_p->Branch("theta_wrtChThrPerp", theta_wrtChThrPerp, "theta_wrtChThrPerp[nParticle]/F");
    inTree_p->Branch("phi_wrtChThrPerp", phi_wrtChThrPerp, "phi_wrtChThrPerp[nParticle]/F");
  }
  if (!writeMinimal && do_neuThrust )          {
    inTree_p->Branch("pt_wrtNeuThr", pt_wrtNeuThr, "pt_wrtNeuThr[nParticle]/F");
    inTree_p->Branch("eta_wrtNeuThr", eta_wrtNeuThr, "eta_wrtNeuThr[nParticle]/F");
    inTree_p->Branch("rap_wrtNeuThr", rap_wrtNeuThr, "rap_wrtNeuThr[nParticle]/F");
    inTree_p->Branch("theta_wrtNeuThr", theta_wrtNeuThr, "theta_wrtNeuThr[nParticle]/F");
    inTree_p->Branch("phi_wrtNeuThr", phi_wrtNeuThr, "phi_wrtNeuThr[nParticle]/F");
    inTree_p->Branch("pt_wrtNeuThrPerp", pt_wrtNeuThrPerp, "pt_wrtNeuThrPerp[nParticle]/F");
    inTree_p->Branch("eta_wrtNeuThrPerp", eta_wrtNeuThrPerp, "eta_wrtNeuThrPerp[nParticle]/F");
    inTree_p->Branch("theta_wrtNeuThrPerp", theta_wrtNeuThrPerp, "theta_wrtNeuThrPerp[nParticle]/F");
    inTree_p->Branch("phi_wrtNeuThrPerp", phi_wrtNeuThrPerp, "phi_wrtNeuThrPerp[nParticle]/F");
    inTree_p->Branch("rap_wrtNeuThrPerp", rap_wrtNeuThrPerp, "rap_wrtNeuThrPerp[nParticle]/F");
   } 
  if (!writeMinimal && do_thrustCorr )         {
    inTree_p->Branch("pt_wrtThrCorr", pt_wrtThrCorr, "pt_wrtThrCorr[nParticle]/F");
    inTree_p->Branch("eta_wrtThrCorr", eta_wrtThrCorr, "eta_wrtThrCorr[nParticle]/F");
    inTree_p->Branch("rap_wrtThrCorr", rap_wrtThrCorr, "rap_wrtThrCorr[nParticle]/F");
    inTree_p->Branch("theta_wrtThrCorr", theta_wrtThrCorr, "theta_wrtThrCorr[nParticle]/F");
    inTree_p->Branch("phi_wrtThrCorr", phi_wrtThrCorr, "phi_wrtThrCorr[nParticle]/F");
    inTree_p->Branch("pt_wrtThrCorrPerp", pt_wrtThrCorrPerp, "pt_wrtThrCorrPerp[nParticle]/F");
    inTree_p->Branch("eta_wrtThrCorrPerp", eta_wrtThrCorrPerp, "eta_wrtThrCorrPerp[nParticle]/F");
    inTree_p->Branch("rap_wrtThrCorrPerp", rap_wrtThrCorrPerp, "rap_wrtThrCorrPerp[nParticle]/F");
    inTree_p->Branch("theta_wrtThrCorrPerp", theta_wrtThrCorrPerp, "theta_wrtThrCorrPerp[nParticle]/F");
    inTree_p->Branch("phi_wrtThrCorrPerp", phi_wrtThrCorrPerp, "phi_wrtThrCorrPerp[nParticle]/F");
  }
  if (!writeMinimal && do_thrustCorrInverse )  {
    inTree_p->Branch("pt_wrtThrCorrInverse", pt_wrtThrCorrInverse, "pt_wrtThrCorrInverse[nParticle]/F");
    inTree_p->Branch("eta_wrtThrCorrInverse", eta_wrtThrCorrInverse, "eta_wrtThrCorrInverse[nParticle]/F");
    inTree_p->Branch("rap_wrtThrCorrInverse", rap_wrtThrCorrInverse, "rap_wrtThrCorrInverse[nParticle]/F");
    inTree_p->Branch("theta_wrtThrCorrInverse", theta_wrtThrCorrInverse, "theta_wrtThrCorrInverse[nParticle]/F");
    inTree_p->Branch("phi_wrtThrCorrInverse", phi_wrtThrCorrInverse, "phi_wrtThrCorrInverse[nParticle]/F");
    inTree_p->Branch("pt_wrtThrCorrInversePerp", pt_wrtThrCorrInversePerp, "pt_wrtThrCorrInversePerp[nParticle]/F");
    inTree_p->Branch("eta_wrtThrCorrInversePerp", eta_wrtThrCorrInversePerp, "eta_wrtThrCorrInversePerp[nParticle]/F");
    inTree_p->Branch("rap_wrtThrCorrInversePerp", rap_wrtThrCorrInversePerp, "rap_wrtThrCorrInversePerp[nParticle]/F");
    inTree_p->Branch("theta_wrtThrCorrInversePerp", theta_wrtThrCorrInversePerp, "theta_wrtThrCorrInversePerp[nParticle]/F");
    inTree_p->Branch("phi_wrtThrCorrInversePerp", phi_wrtThrCorrInversePerp, "phi_wrtThrCorrInversePerp[nParticle]/F");
  }
  if ( do_thrustMissP )        {
    inTree_p->Branch("pt_wrtThrMissP", pt_wrtThrMissP, "pt_wrtThrMissP[nParticle]/F");
    inTree_p->Branch("eta_wrtThrMissP", eta_wrtThrMissP, "eta_wrtThrMissP[nParticle]/F");
    inTree_p->Branch("rap_wrtThrMissP", rap_wrtThrMissP, "rap_wrtThrMissP[nParticle]/F");
    inTree_p->Branch("theta_wrtThrMissP", theta_wrtThrMissP, "theta_wrtThrMissP[nParticle]/F");
    inTree_p->Branch("phi_wrtThrMissP", phi_wrtThrMissP, "phi_wrtThrMissP[nParticle]/F");
    if (!writeMinimal){
      inTree_p->Branch("pt_wrtThrMissPPerp", pt_wrtThrMissPPerp, "pt_wrtThrMissPPerp[nParticle]/F");
      inTree_p->Branch("eta_wrtThrMissPPerp", eta_wrtThrMissPPerp, "eta_wrtThrMissPPerp[nParticle]/F");
      inTree_p->Branch("rap_wrtThrMissPPerp", rap_wrtThrMissPPerp, "rap_wrtThrMissPPerp[nParticle]/F");
      inTree_p->Branch("theta_wrtThrMissPPerp", theta_wrtThrMissPPerp, "theta_wrtThrMissPPerp[nParticle]/F");
      inTree_p->Branch("phi_wrtThrMissPPerp", phi_wrtThrMissPPerp, "phi_wrtThrMissPPerp[nParticle]/F");
    }
  }
  return;
}


inline void particleData::preFillClean()
{
  for(Int_t i = 0; i < nParticle; ++i){
    px[i] = reducedPrecision(px[i]);
    py[i] = reducedPrecision(py[i]);
    pz[i] = reducedPrecision(pz[i]);
    pt[i] = reducedPrecision(pt[i]);
    pmag[i] = reducedPrecision(pmag[i]);
    rap[i] = reducedPrecision(rap[i]);
    eta[i] = reducedPrecision(eta[i]);
    theta[i] = reducedPrecision(theta[i]);
    phi[i] = reducedPrecision(phi[i]);
    mass[i] = reducedPrecision(mass[i]);
    pid[i] = reducedPrecision(pid[i]);
    d0[i] = reducedPrecision(d0[i]);
    z0[i] = reducedPrecision(z0[i]);
    vx[i] = reducedPrecision(vx[i]);
    vy[i] = reducedPrecision(vy[i]);
    vz[i] = reducedPrecision(vz[i]);
    weight[i] = reducedPrecision(weight[i]);
    index[i] = reducedPrecision(index[i]);
    correspondenceIndex[i] = reducedPrecision(correspondenceIndex[i]);
    
    if (!initMinimal)
    {
      pt_wrtThr[i] = reducedPrecision(pt_wrtThr[i]);
      eta_wrtThr[i] = reducedPrecision(eta_wrtThr[i]);
      rap_wrtThr[i] = reducedPrecision(rap_wrtThr[i]);
      theta_wrtThr[i] = reducedPrecision(theta_wrtThr[i]);
      phi_wrtThr[i] = reducedPrecision(phi_wrtThr[i]);
      pt_wrtThrPerp[i] = reducedPrecision(pt_wrtThrPerp[i]);
      eta_wrtThrPerp[i] = reducedPrecision(eta_wrtThrPerp[i]);
      rap_wrtThrPerp[i] = reducedPrecision(rap_wrtThrPerp[i]);
      theta_wrtThrPerp[i] = reducedPrecision(theta_wrtThrPerp[i]);
      phi_wrtThrPerp[i] = reducedPrecision(phi_wrtThrPerp[i]);
      if ( do_chThrust )           {
        pt_wrtChThr[i] = reducedPrecision(pt_wrtChThr[i]);
        eta_wrtChThr[i] = reducedPrecision(eta_wrtChThr[i]);
        rap_wrtChThr[i] = reducedPrecision(rap_wrtChThr[i]);
        theta_wrtChThr[i] = reducedPrecision(theta_wrtChThr[i]);
        phi_wrtChThr[i] = reducedPrecision(phi_wrtChThr[i]);
        pt_wrtChThrPerp[i] = reducedPrecision(pt_wrtChThrPerp[i]);
        eta_wrtChThrPerp[i] = reducedPrecision(eta_wrtChThrPerp[i]);
        rap_wrtChThrPerp[i] = reducedPrecision(rap_wrtChThrPerp[i]);
        theta_wrtChThrPerp[i] = reducedPrecision(theta_wrtChThrPerp[i]);
        phi_wrtChThrPerp[i] = reducedPrecision(phi_wrtChThrPerp[i]);
      }
      if ( do_neuThrust )          {
        pt_wrtNeuThr[i] = reducedPrecision(pt_wrtNeuThr[i]);
        eta_wrtNeuThr[i] = reducedPrecision(eta_wrtNeuThr[i]);
        rap_wrtNeuThr[i] = reducedPrecision(rap_wrtNeuThr[i]);
        theta_wrtNeuThr[i] = reducedPrecision(theta_wrtNeuThr[i]);
        phi_wrtNeuThr[i] = reducedPrecision(phi_wrtNeuThr[i]);
        pt_wrtNeuThrPerp[i] = reducedPrecision(pt_wrtNeuThrPerp[i]);
        theta_wrtNeuThrPerp[i] = reducedPrecision(theta_wrtNeuThrPerp[i]);
        eta_wrtNeuThrPerp[i] = reducedPrecision(eta_wrtNeuThrPerp[i]);
        phi_wrtNeuThrPerp[i] = reducedPrecision(phi_wrtNeuThrPerp[i]);
        rap_wrtNeuThrPerp[i] = reducedPrecision(rap_wrtNeuThrPerp[i]);
      }
      if ( do_thrustCorr )         {
        pt_wrtThrCorr[i] = reducedPrecision(pt_wrtThrCorr[i]);
        eta_wrtThrCorr[i] = reducedPrecision(eta_wrtThrCorr[i]);
        rap_wrtThrCorr[i] = reducedPrecision(rap_wrtThrCorr[i]);
        theta_wrtThrCorr[i] = reducedPrecision(theta_wrtThrCorr[i]);
        phi_wrtThrCorr[i] = reducedPrecision(phi_wrtThrCorr[i]);
        pt_wrtThrCorrPerp[i] = reducedPrecision(pt_wrtThrCorrPerp[i]);
        eta_wrtThrCorrPerp[i] = reducedPrecision(eta_wrtThrCorrPerp[i]);
        rap_wrtThrCorrPerp[i] = reducedPrecision(rap_wrtThrCorrPerp[i]);
        theta_wrtThrCorrPerp[i] = reducedPrecision(theta_wrtThrCorrPerp[i]);
        phi_wrtThrCorrPerp[i] = reducedPrecision(phi_wrtThrCorrPerp[i]);
      }
      if ( do_thrustCorrInverse )  {
        pt_wrtThrCorrInverse[i] = reducedPrecision(pt_wrtThrCorrInverse[i]);
        eta_wrtThrCorrInverse[i] = reducedPrecision(eta_wrtThrCorrInverse[i]);
        rap_wrtThrCorrInverse[i] = reducedPrecision(rap_wrtThrCorrInverse[i]);
        theta_wrtThrCorrInverse[i] = reducedPrecision(theta_wrtThrCorrInverse[i]);
        phi_wrtThrCorrInverse[i] = reducedPrecision(phi_wrtThrCorrInverse[i]);
        pt_wrtThrCorrInversePerp[i] = reducedPrecision(pt_wrtThrCorrInversePerp[i]);
        eta_wrtThrCorrInversePerp[i] = reducedPrecision(eta_wrtThrCorrInversePerp[i]);
        rap_wrtThrCorrInversePerp[i] = reducedPrecision(rap_wrtThrCorrInversePerp[i]);
        theta_wrtThrCorrInversePerp[i] = reducedPrecision(theta_wrtThrCorrInversePerp[i]);
        phi_wrtThrCorrInversePerp[i] = reducedPrecision(phi_wrtThrCorrInversePerp[i]);
      }
      if ( do_thrustMissP )     {
        pt_wrtThrMissPPerp[i] = reducedPrecision(pt_wrtThrMissPPerp[i]);
        eta_wrtThrMissPPerp[i] = reducedPrecision(eta_wrtThrMissPPerp[i]);
        rap_wrtThrMissPPerp[i] = reducedPrecision(rap_wrtThrMissPPerp[i]);
        theta_wrtThrMissPPerp[i] = reducedPrecision(theta_wrtThrMissPPerp[i]);
        phi_wrtThrMissPPerp[i] = reducedPrecision(phi_wrtThrMissPPerp[i]);
      }
    }

    if ( do_thrustMissP )        {
      pt_wrtThrMissP[i] = reducedPrecision(pt_wrtThrMissP[i]);
      eta_wrtThrMissP[i] = reducedPrecision(eta_wrtThrMissP[i]);
      rap_wrtThrMissP[i] = reducedPrecision(rap_wrtThrMissP[i]);
      theta_wrtThrMissP[i] = reducedPrecision(theta_wrtThrMissP[i]);
      phi_wrtThrMissP[i] = reducedPrecision(phi_wrtThrMissP[i]);
    }
    artificAcceptEffCorrection[i] = reducedPrecision(artificAcceptEffCorrection[i]);
    anc[i] = reducedPrecision(anc[i]);
  }

  return;
}


inline void particleData::_init()
{
  printf("// Variable Declaration\n");
  for (int i = 0; i < nVar; ++i)
  {
    // bool
    if(varStr[i].find("isMC") != std::string::npos)                     printf("\tbool  %s;\n", varStr[i].c_str());
    else if(varStr[i].find("isOnres") != std::string::npos)             printf("\tbool  %s;\n", varStr[i].c_str());
    else if(varStr[i].find("highPurity") != std::string::npos)          printf("\tbool  %s[nMaxPart];\n", varStr[i].c_str());
    else if(varStr[i].find("passesArtificAccept") != std::string::npos) printf("\tbool  %s[nMaxPart];\n", varStr[i].c_str());

    // short
    else if(varStr[i].find("charge") != std::string::npos)     printf("\tShort_t  %s[nMaxPart];\n", varStr[i].c_str());
    else if(varStr[i].find("pwflag") != std::string::npos)     printf("\tShort_t  %s[nMaxPart];\n", varStr[i].c_str());
    else if(varStr[i].find("ntpc") != std::string::npos)       printf("\tShort_t  %s[nMaxPart];\n", varStr[i].c_str());
    else if(varStr[i].find("nitc") != std::string::npos)       printf("\tShort_t  %s[nMaxPart];\n", varStr[i].c_str());
    else if(varStr[i].find("nvdet") != std::string::npos)      printf("\tShort_t  %s[nMaxPart];\n", varStr[i].c_str());

    // unsigned long long
    else if(varStr[i].find("uniqueID") != std::string::npos)      printf("\tunsigned long long  %s;\n", varStr[i].c_str());

    // int    
    else if(varStr[i].find("EventNo") != std::string::npos)    printf("\tint  %s;\n", varStr[i].c_str());
    else if(varStr[i].find("RunNo") != std::string::npos)      printf("\tint  %s;\n", varStr[i].c_str());
    else if(varStr[i].find("year") != std::string::npos)       printf("\tint  %s;\n", varStr[i].c_str());
    else if(varStr[i].find("subDir") != std::string::npos)     printf("\tint  %s;\n", varStr[i].c_str());
    else if(varStr[i].find("process") != std::string::npos)    printf("\tint  %s;\n", varStr[i].c_str());
    else if(varStr[i].find("source") != std::string::npos)     printf("\tint  %s;\n", varStr[i].c_str());
    else if(varStr[i].find("bFlag") != std::string::npos)      printf("\tint  %s;\n", varStr[i].c_str());
    else if(varStr[i].find("nParticle") != std::string::npos)  printf("\tint  %s;\n", varStr[i].c_str());
    else if(varStr[i].find("pid") != std::string::npos)        printf("\tfloat  %s[nMaxPart];\n", varStr[i].c_str());

    // float
    else{
      if(varStr[i].find("Energy") != std::string::npos)                 printf("\tfloat  %s;\n", varStr[i].c_str());
      else if(varStr[i].find("particleWeight") != std::string::npos)    printf("\tfloat  %s;\n", varStr[i].c_str());
      else if(varStr[i].find("bx") != std::string::npos)                printf("\tfloat  %s;\n", varStr[i].c_str());
      else if(varStr[i].find("by") != std::string::npos)                printf("\tfloat  %s;\n", varStr[i].c_str());
      else if(varStr[i].find("ebx") != std::string::npos)               printf("\tfloat  %s;\n", varStr[i].c_str());
      else if(varStr[i].find("eby") != std::string::npos)               printf("\tfloat  %s;\n", varStr[i].c_str());

      else printf("\tfloat  %s[nMaxPart];\n", varStr[i].c_str());
    }
  }
  printf("// Variable Initialization (please rearrange the branches in for loop on your own)\n");
  for (int i = 0; i < nVar; ++i)
  {
    // bool
    if(varStr[i].find("isMC") != std::string::npos)                     printf("\t%s = false;\n", varStr[i].c_str());
    else if(varStr[i].find("isOnres") != std::string::npos)             printf("\t%s = false;\n", varStr[i].c_str());
    else if(varStr[i].find("highPurity") != std::string::npos)          printf("\t%s[i] = true;\n", varStr[i].c_str());
    else if(varStr[i].find("passesArtificAccept") != std::string::npos) printf("\t%s[i] = true;\n", varStr[i].c_str());

    // short
    else if(varStr[i].find("charge") != std::string::npos)     printf("\t%s[i] = -127;\n", varStr[i].c_str());
    else if(varStr[i].find("pwflag") != std::string::npos)     printf("\t%s[i] = -127;\n", varStr[i].c_str());
    else if(varStr[i].find("ntpc") != std::string::npos)       printf("\t%s[i] = -127;\n", varStr[i].c_str());
    else if(varStr[i].find("nitc") != std::string::npos)       printf("\t%s[i] = -127;\n", varStr[i].c_str());
    else if(varStr[i].find("nvdet") != std::string::npos)      printf("\t%s[i] = -127;\n", varStr[i].c_str());

    //unsigned long long
    else if(varStr[i].find("uniqueID") != std::string::npos)   printf("\t%s = 0;\n", varStr[i].c_str());

    // int    
    else if(varStr[i].find("EventNo") != std::string::npos)    printf("\t%s = -999;\n", varStr[i].c_str());
    else if(varStr[i].find("RunNo") != std::string::npos)      printf("\t%s = -999;\n", varStr[i].c_str());
    else if(varStr[i].find("year") != std::string::npos)       printf("\t%s = -999;\n", varStr[i].c_str());
    else if(varStr[i].find("subDir") != std::string::npos)     printf("\t%s = -999;\n", varStr[i].c_str());
    else if(varStr[i].find("process") != std::string::npos)    printf("\t%s = -999;\n", varStr[i].c_str());
    else if(varStr[i].find("source") != std::string::npos)     printf("\t%s = -999;\n", varStr[i].c_str());
    else if(varStr[i].find("bFlag") != std::string::npos)      printf("\t%s = -999;\n", varStr[i].c_str());
    else if(varStr[i].find("nParticle") != std::string::npos)  printf("\t%s = -999;\n", varStr[i].c_str());
    else if(varStr[i].find("pid") != std::string::npos)        printf("\t%s[i] = -999;\n", varStr[i].c_str());

    // float
    else{
      if(varStr[i].find("Energy") != std::string::npos)                 printf("\t%s = 0;\n", varStr[i].c_str());
      else if(varStr[i].find("particleWeight") != std::string::npos)    printf("\t%s = 1;\n", varStr[i].c_str());
      else if(varStr[i].find("bx") != std::string::npos)                printf("\t%s = -999;\n", varStr[i].c_str());
      else if(varStr[i].find("by") != std::string::npos)                printf("\t%s = -999;\n", varStr[i].c_str());
      else if(varStr[i].find("ebx") != std::string::npos)               printf("\t%s = -999;\n", varStr[i].c_str());
      else if(varStr[i].find("eby") != std::string::npos)               printf("\t%s = -999;\n", varStr[i].c_str());

      else printf("\t%s[i] = -999;\n", varStr[i].c_str());
    }
  }
}

inline void particleData::_setRead()
{
  for (int i = 0; i < nVar; ++i)
  {
    // bool
    if(varStr[i].find("isMC") != std::string::npos)                     printf("\tif(varIsGood[%d]) inTree_p->SetBranchAddress(\"isMC\", &isMC);\n", i);
    else if(varStr[i].find("isOnres") != std::string::npos)             printf("\tif(varIsGood[%d]) inTree_p->SetBranchAddress(\"isOnres\", &isOnres);\n", i);
    else if(varStr[i].find("highPurity") != std::string::npos)          printf("\tif(varIsGood[%d]) inTree_p->SetBranchAddress(\"highPurity\", highPurity);\n", i);
    else if(varStr[i].find("passesArtificAccept") != std::string::npos) printf("\tif(varIsGood[%d]) inTree_p->SetBranchAddress(\"passesArtificAccept\", passesArtificAccept);\n", i);

    // short
    else if(varStr[i].find("charge") != std::string::npos)     printf("\tif(varIsGood[%d]) inTree_p->SetBranchAddress(\"charge\", charge);\n", i);
    else if(varStr[i].find("pwflag") != std::string::npos)     printf("\tif(varIsGood[%d]) inTree_p->SetBranchAddress(\"pwflag\", pwflag);\n", i);
    else if(varStr[i].find("ntpc") != std::string::npos)       printf("\tif(varIsGood[%d]) inTree_p->SetBranchAddress(\"ntpc\", ntpc);\n", i);
    else if(varStr[i].find("nitc") != std::string::npos)       printf("\tif(varIsGood[%d]) inTree_p->SetBranchAddress(\"nitc\", nitc);\n", i);
    else if(varStr[i].find("nvdet") != std::string::npos)      printf("\tif(varIsGood[%d]) inTree_p->SetBranchAddress(\"nvdet\", nvdet);\n", i);
    
    // unsigned long long
    else if(varStr[i].find("uniqueID") != std::string::npos)   printf("\tif(varIsGood[%d]) inTree_p->SetBranchAddress(\"uniqueID\", &uniqueID);\n", i);

    // int    
    else if(varStr[i].find("EventNo") != std::string::npos)    printf("\tif(varIsGood[%d]) inTree_p->SetBranchAddress(\"EventNo\", &EventNo);\n", i);
    else if(varStr[i].find("RunNo") != std::string::npos)      printf("\tif(varIsGood[%d]) inTree_p->SetBranchAddress(\"RunNo\", &RunNo);\n", i);
    else if(varStr[i].find("year") != std::string::npos)       printf("\tif(varIsGood[%d]) inTree_p->SetBranchAddress(\"year\", &year);\n", i);
    else if(varStr[i].find("subDir") != std::string::npos)     printf("\tif(varIsGood[%d]) inTree_p->SetBranchAddress(\"subDir\", &subDir);\n", i);
    else if(varStr[i].find("process") != std::string::npos)    printf("\tif(varIsGood[%d]) inTree_p->SetBranchAddress(\"process\", &process);\n", i);
    else if(varStr[i].find("source") != std::string::npos)     printf("\tif(varIsGood[%d]) inTree_p->SetBranchAddress(\"source\", &source);\n", i);
    else if(varStr[i].find("bFlag") != std::string::npos)      printf("\tif(varIsGood[%d]) inTree_p->SetBranchAddress(\"bFlag\", &bFlag);\n", i);
    else if(varStr[i].find("nParticle") != std::string::npos)  printf("\tif(varIsGood[%d]) inTree_p->SetBranchAddress(\"nParticle\", &nParticle);\n", i);
    else if(varStr[i].find("pid") != std::string::npos)        printf("\tif(varIsGood[%d]) inTree_p->SetBranchAddress(\"pid\", pid);\n", i);

    // float
    else{
      if(varStr[i].find("Energy") != std::string::npos)                 printf("\tif(varIsGood[%d]) inTree_p->SetBranchAddress(\"Energy\", &Energy);\n", i);
      else if(varStr[i].find("particleWeight") != std::string::npos)    printf("\tif(varIsGood[%d]) inTree_p->SetBranchAddress(\"particleWeight\", &particleWeight);\n", i);
      else if(varStr[i].find("bx") != std::string::npos)                printf("\tif(varIsGood[%d]) inTree_p->SetBranchAddress(\"bx\", &bx);\n", i);
      else if(varStr[i].find("by") != std::string::npos)                printf("\tif(varIsGood[%d]) inTree_p->SetBranchAddress(\"by\", &by);\n", i);
      else if(varStr[i].find("ebx") != std::string::npos)               printf("\tif(varIsGood[%d]) inTree_p->SetBranchAddress(\"ebx\", &ebx);\n", i);
      else if(varStr[i].find("eby") != std::string::npos)               printf("\tif(varIsGood[%d]) inTree_p->SetBranchAddress(\"eby\", &eby);\n", i);

      else printf("\tif(varIsGood[%d]) inTree_p->SetBranchAddress(\"%s\", %s);\n", i, varStr[i].c_str(), varStr[i].c_str());
    }
  }
}

inline void particleData::_setWrite()
{
  for (int i = 0; i < nVar; ++i)
  {
    // bool
    if(varStr[i].find("isMC") != std::string::npos)                     printf("\tinTree_p->Branch(\"isMC\", &isMC, \"isMC/O\");\n");
    else if(varStr[i].find("isOnres") != std::string::npos)             printf("\tinTree_p->Branch(\"isOnres\", &isOnres, \"isOnres/O\");\n");
    else if(varStr[i].find("highPurity") != std::string::npos)          printf("\tinTree_p->Branch(\"highPurity\", highPurity, \"highPurity[nParticle]/O\");\n");
    else if(varStr[i].find("passesArtificAccept") != std::string::npos) printf("\tinTree_p->Branch(\"passesArtificAccept\", passesArtificAccept, \"passesArtificAccept[nParticle]/O\");\n");

    // short
    else if(varStr[i].find("charge") != std::string::npos)     printf("\tinTree_p->Branch(\"charge\", charge, \"charge[nParticle]/S\");\n");
    else if(varStr[i].find("pwflag") != std::string::npos)     printf("\tinTree_p->Branch(\"pwflag\", pwflag, \"pwflag[nParticle]/S\");\n");
    else if(varStr[i].find("ntpc") != std::string::npos)       printf("\tinTree_p->Branch(\"ntpc\", ntpc, \"ntpc[nParticle]/S\");\n");
    else if(varStr[i].find("nitc") != std::string::npos)       printf("\tinTree_p->Branch(\"nitc\", nitc, \"nitc[nParticle]/S\");\n");
    else if(varStr[i].find("nvdet") != std::string::npos)      printf("\tinTree_p->Branch(\"nvdet\", nvdet, \"nvdet[nParticle]/S\");\n");

    // unsigned long long
    else if(varStr[i].find("uniqueID") != std::string::npos)      printf("\tinTree_p->Branch(\"uniqueID\", &uniqueID, \"uniqueID/l\");\n");

    // int    
    else if(varStr[i].find("EventNo") != std::string::npos)    printf("\tinTree_p->Branch(\"EventNo\", &EventNo, \"EventNo/I\");\n");
    else if(varStr[i].find("RunNo") != std::string::npos)      printf("\tinTree_p->Branch(\"RunNo\", &RunNo, \"RunNo/I\");\n");
    else if(varStr[i].find("year") != std::string::npos)       printf("\tinTree_p->Branch(\"year\", &year, \"year/I\");\n");
    else if(varStr[i].find("subDir") != std::string::npos)     printf("\tinTree_p->Branch(\"subDir\", &subDir, \"subDir/I\");\n");
    else if(varStr[i].find("process") != std::string::npos)    printf("\tinTree_p->Branch(\"process\", &process, \"process/I\");\n");
    else if(varStr[i].find("source") != std::string::npos)     printf("\tinTree_p->Branch(\"source\", &source, \"source/I\");\n");
    else if(varStr[i].find("bFlag") != std::string::npos)      printf("\tinTree_p->Branch(\"bFlag\", &bFlag, \"bFlag/I\");\n");
    else if(varStr[i].find("nParticle") != std::string::npos)  printf("\tinTree_p->Branch(\"nParticle\", &nParticle, \"nParticle/I\"); \n");
    else if(varStr[i].find("pid") != std::string::npos)        printf("\tinTree_p->Branch(\"pid\", pid, \"pid[nParticle]/I\");\n");

    // float
    else{
      if(varStr[i].find("Energy") != std::string::npos)                 printf("\tinTree_p->Branch(\"Energy\", &Energy, \"Energy/F\");\n");
      else if(varStr[i].find("particleWeight") != std::string::npos)    printf("\tinTree_p->Branch(\"particleWeight\", &particleWeight, \"particleWeight/F\");\n");
      else if(varStr[i].find("bx") != std::string::npos)                printf("\tinTree_p->Branch(\"bx\", &bx, \"bx/F\");\n");
      else if(varStr[i].find("by") != std::string::npos)                printf("\tinTree_p->Branch(\"by\", &by, \"by/F\");\n");
      else if(varStr[i].find("ebx") != std::string::npos)               printf("\tinTree_p->Branch(\"ebx\", &ebx, \"ebx/F\");\n");
      else if(varStr[i].find("eby") != std::string::npos)               printf("\tinTree_p->Branch(\"eby\", &eby, \"eby/F\");\n");

      else printf("\tinTree_p->Branch(\"%s\", %s, \"%s[nParticle]/F\");\n", varStr[i].c_str(), varStr[i].c_str(), varStr[i].c_str());
    }
  }
}
inline void particleData::_setReducedPrecision()
{
  for (int i = 0; i < nVar; ++i)
  {
    // bool
    if(varStr[i].find("isMC") != std::string::npos)                     continue;
    else if(varStr[i].find("isOnres") != std::string::npos)             continue;
    else if(varStr[i].find("highPurity") != std::string::npos)          continue;
    else if(varStr[i].find("passesArtificAccept") != std::string::npos) continue;

    // short
    else if(varStr[i].find("charge") != std::string::npos)     continue;
    else if(varStr[i].find("pwflag") != std::string::npos)     continue;
    else if(varStr[i].find("ntpc") != std::string::npos)       continue;
    else if(varStr[i].find("nitc") != std::string::npos)       continue;
    else if(varStr[i].find("nvdet") != std::string::npos)      continue;

    // unsigned long long
    else if(varStr[i].find("uniqueID") != std::string::npos)      continue;

    // int    
    else if(varStr[i].find("EventNo") != std::string::npos)    continue;
    else if(varStr[i].find("RunNo") != std::string::npos)      continue;
    else if(varStr[i].find("year") != std::string::npos)       continue;
    else if(varStr[i].find("subDir") != std::string::npos)     continue;
    else if(varStr[i].find("process") != std::string::npos)    continue;
    else if(varStr[i].find("source") != std::string::npos)     continue;
    else if(varStr[i].find("bFlag") != std::string::npos)      continue;
    else if(varStr[i].find("nParticle") != std::string::npos)  continue;
    else if(varStr[i].find("pid") != std::string::npos)        continue;

    // float
    else{
      if(varStr[i].find("Energy") != std::string::npos)                 continue;
      else if(varStr[i].find("particleWeight") != std::string::npos)    continue;
      else if(varStr[i].find("bx") != std::string::npos)                continue;
      else if(varStr[i].find("by") != std::string::npos)                continue;
      else if(varStr[i].find("ebx") != std::string::npos)               continue;
      else if(varStr[i].find("eby") != std::string::npos)               continue;

      else printf("\t%s[i] = reducedPrecision(%s[i]);\n", varStr[i].c_str(), varStr[i].c_str());
    }
  }
}

#endif
