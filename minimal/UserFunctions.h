#include <algorithm>
#include <regex>
#include <TH2F.h>
#include <TRandom2.h>

const float photonsEtaMaxBarrel = 1.4442;
//const float photonsEtaMinEndcap = 1.566;
//const float photonsEtaMaxEndcap = 2.5;
// Changed values to reduce e->g backgronud
const float photonsEtaMinEndcap = 1.6;
const float photonsEtaMaxEndcap = 2.5;

// https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation80XReReco
enum bTaggingEnum { CSVv2L, CSVv2M, CSVv2T };
map<bTaggingEnum,float> bTaggingWorkingPoints = {
  {CSVv2L, 0.5426},
  {CSVv2M, 0.8484},
  {CSVv2T, 0.9535}
};

pair<TVector3,TVector3> megajets( const vector<TVector3>& jets ) {
  // code from https://twiki.cern.ch/twiki/bin/view/CMSPublic/RazorLikelihoodHowTo

  TVector3 j1, j2;
  int N_comb = (int) pow( 2, jets.size() );

  double M_min = numeric_limits<double>::max();
  int j_count;
  for(int i = 1; i < N_comb-1; i++){
    TVector3 j_temp1, j_temp2;
    int itemp = i;
    j_count = N_comb/2;
    int count = 0;
    while(j_count > 0){
      if(itemp/j_count == 1){
        j_temp1 += jets[count];
      } else {
        j_temp2 += jets[count];
      }
      itemp -= j_count*(itemp/j_count);
      j_count /= 2;
      count++;
    }
    double M_temp = j_temp1.Mag2()+j_temp2.Mag2();
    // smallest mass
    if(M_temp < M_min){
      M_min = M_temp;
      j1 = j_temp1;
      j2 = j_temp2;
    }
  }

  if(j2.Pt() > j1.Pt()) return pair<TVector3,TVector3>(j2,j1);
  else                  return pair<TVector3,TVector3>(j1,j2);
}


pair<float,float> razorVariables( const pair<TVector3,TVector3>& megajets, const TVector3& met ) {
  float mr2 = pow( megajets.first.Mag() + megajets.second.Mag(), 2 ) - pow( megajets.first.Z() + megajets.second.Z(), 2 );
  float mrt2 = 0.5 * ( met.Mag()*( megajets.first.Pt() + megajets.second.Pt() ) - met.XYvector() * ( megajets.first.XYvector() + megajets.second.XYvector() ) );
  float r2 = mrt2 / mr2;
  return pair<float,float>(sqrt(mr2),r2);
}

std::ostream& operator << ( std::ostream& os, const TVector3& p ) {
  os << p.Pt() << "\t" << p.Eta() << "\t" << p.Phi();
  return os;
}


string getOutputFilename( const string& inputFileName, const string& appendix="hists" ) {

  // Converts "/path/to/ntuple/QCD_nTuple.root" to "QCD_hists.root"

  auto startPos = inputFileName.rfind("/");
  auto endPos = inputFileName.find("_nTuple.root");
  string outputName = "out.root";
  if( endPos != string::npos ) {
    outputName = inputFileName.substr( startPos+1, endPos-startPos-1 ) + "_"+appendix+".root";
  }

  // for signal scans
  if( inputFileName.find("/T5") != string::npos ) {
    endPos = inputFileName.find(".root");
    if( endPos != string::npos ) {
      outputName = inputFileName.substr( startPos+1, endPos-startPos-1 ) + "_"+appendix+".root";
    }
  }
  return outputName;

}

template <typename VectorClass>
int indexOfMatchedParticle( const tree::Particle& tag, const std::vector<VectorClass>& particles, float deltaR=.1, float relPt=-1 ) {
  for( int i=0; i<(int)particles.size(); ++i ) {
    auto p = particles.at(i)->p;
    if(    ( deltaR<0 || p.DeltaR( tag.p ) < deltaR )
        && ( relPt<0  || fabs(p.Pt()-tag.p.Pt())/tag.p.Pt() < relPt ) ) {
      return i;
    }
  }
  return -1;
}

vector<int> getRunList( const string& filename="/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions15/13TeV/Cert_246908-260627_13TeV_PromptReco_Collisions15_25ns_JSON.txt" ) {
  ifstream t(filename);
  string s((std::istreambuf_iterator<char>(t)),std::istreambuf_iterator<char>());
  smatch m;
  regex e ("\"([\\d]+)\"");
  vector<int> runList;
  while (regex_search (s,m,e)) {
    runList.push_back( stoi( m[1] ) );
    s = m.suffix().str();
  }
  return runList;
}

std::ostream &operator<<(std::ostream &os, TVector3 &p) {
      return os << p.Pt() << "\t" << p.Eta() << "\t" << p.Phi();
}

inline int pseudoRandom( float f ) {
  //return *(int*)&f; // Interprets float as int
  return ((int)(f*1e6)) % 1000;
}

inline bool pseudoRandomSort(const tree::Particle p1, const tree::Particle p2) {
  return pseudoRandom(p1.p.Phi()) > pseudoRandom(p2.p.Phi());
}

class JetSelector {
  public:
    JetSelector( const string& filename, const string& histname ) {
      TFile f( filename.c_str() );
      TH2F* h2 = (TH2F*)f.Get( histname.c_str() );
      if(h2) {
        for(unsigned i=1; i<(unsigned)h2->GetNbinsX()+2; i++) {
          auto h = h2->ProjectionY(("proj"+to_string(i)).c_str(), i, i );
          if(h->Integral()) {
            h->Scale(1./h->Integral());
            hMap[i-1] = *h;
          }
        }
      }
    }

    ~JetSelector(){
    }

    unsigned getJetN( unsigned nJets ){
      if(!hMap.count(nJets)) return 0;
      return (unsigned)std::round(hMap.at(nJets).GetRandom());
    }

  private:
    map<unsigned,TH1D> hMap;
};

float sampleCrossSection(const string& inputFileName) {
  auto startPos = inputFileName.rfind("/");
  auto endPos = inputFileName.find("_nTuple.root");
  auto extPos = inputFileName.find("_ext");
  string sampleName = "";
  if (endPos != string::npos) {
    sampleName = inputFileName.substr(startPos+1, min(endPos,extPos)-startPos-1);
  }

  map<string,float> xs = {
    {"GJets_HT-40To100", 20790},
    {"GJets_HT-100To200", 9238},
    {"GJets_HT-200To400", 2305},
    {"GJets_HT-400To600", 274.4},
    {"GJets_HT-600ToInf", 93.46},
    {"GJets_DR-0p4_HT-40To100", 17420},
    {"GJets_DR-0p4_HT-100To200", 5383},
    {"GJets_DR-0p4_HT-200To400", 1176},
    {"GJets_DR-0p4_HT-400To600", 132.1},
    {"GJets_DR-0p4_HT-600ToInf", 44.32},
    {"QCD_HT50to100", 246400000},
    {"QCD_HT100to200", 27990000},
    {"QCD_HT200to300", 1712000},
    {"QCD_HT300to500", 347700},
    {"QCD_HT500to700", 32100},
    {"QCD_HT700to1000", 6831},
    {"QCD_HT1000to1500", 1207},
    {"QCD_HT1500to2000", 119.9},
    {"QCD_HT2000toInf", 25.24},
    {"WJetsToLNu_HT-100To200", 1345.*1.21},
    {"WJetsToLNu_HT-200To400", 359.7*1.21},
    {"WJetsToLNu_HT-400To600", 48.91*1.21},
    {"WJetsToLNu_HT-600To800", 12.05*1.21},
    {"WJetsToLNu_HT-800To1200", 5.501*1.21},
    {"WJetsToLNu_HT-1200To2500", 1.329*1.21},
    {"WJetsToLNu_HT-2500ToInf", 0.03216*1.21},
    {"WJetsToLNu_HT-600ToInf", 18.77*1.21},
    {"ZJetsToNuNu_HT-100To200", 280.47*1.23},
    {"ZJetsToNuNu_HT-200To400", 77.67*1.23},
    {"ZJetsToNuNu_HT-400To600", 10.73*1.23},
    {"ZJetsToNuNu_HT-600To800", 2.559*1.23},
    {"ZJetsToNuNu_HT-800To1200", 1.1796*1.23},
    {"ZJetsToNuNu_HT-1200To2500", 0.28833*1.23},
    {"ZJetsToNuNu_HT-2500ToInf", 0.006945*1.23},
    {"TTJets-amcatnloFXFX", 831.76},
    {"TTJets-madgraphMLM", 831.76},
    {"TTJets_HT-0to600", 831.76},
    {"TTJets_HT-600to800", 2.6267},
    {"TTJets_HT-800to1200", 1.0817},
    {"TTJets_HT-1200to2500", 0.19579},
    {"TTJets_HT-2500toInf", 0.0023331},
    {"WGJets_MonoPhoton_PtG-40to130", 21.76},
    {"WGJets_MonoPhoton_PtG-130", 1.125},
    {"WGToLNuG-amcatnloFXFX_ext", 512.1},
    {"WGToLNuG_PtG-130-amcatnloFXFX", 1.125},
    {"ZNuNuGJets_MonoPhoton_PtG-40to130", 2.789*1.39},
    {"ZNuNuGJets_MonoPhoton_PtG-130", 0.1832*1.39},
    {"ZGTo2NuG", 27.99},
    {"ZGTo2NuG_PtG-130", 0.2762},
    {"TTGJets", 3.697},
    {"TGJets_amcatnlo_madspin", 2.967},
    {"DYJetsToLL_M-50", 6025.2},
    {"SMS-T5Wg_1600_100", 0.00810078},
    {"SMS-T5Wg_1600_1500", 0.00810078},
    {"SMS-T5Wg_1750_1700", 0.00359842},
    {"SMS-T5Wg_2000_100", 0.000981077},
    {"SMS-T6gg_1750_1650", 0.000646271}
  };
  if (xs.count(sampleName)) {
    return xs.at(sampleName);
  } else {
    cout << "No cross section defined for " << inputFileName << endl;
    return -1;
  }
}

std::pair<float, float> jerScales(float eta) {
  // Copied on 2016-09-14 from https://twiki.cern.ch/twiki/bin/viewauth/CMS/JetResolution
  // for 80X
  eta = std::abs(eta);
  float sf=0, e=0;
  if (eta<.5) {
    sf = 1.122; e = 0.026;
  } else if (eta<.8) {
    sf = 1.167; e = 0.048;
  } else if (eta<1.1) {
    sf = 1.168; e = 0.046;
  } else if (eta<1.3) {
    sf = 1.029; e = 0.066;
  } else if (eta<1.7) {
    sf = 1.115; e = 0.03;
  } else if (eta<1.9) {
    sf = 1.041; e = 0.062;
  } else if (eta<2.1) {
    sf = 1.167; e = 0.086;
  } else if (eta<2.3) {
    sf = 1.094; e = 0.093;
  } else if (eta<2.5) {
    sf = 1.168; e = 0.120;
  } else if (eta<2.8) {
    sf = 1.266; e = 0.132;
  } else if (eta<3.0) {
    sf = 1.595; e = 0.175;
  } else if (eta<3.2) {
    sf = 0.998; e = 0.066;
  } else {
    sf = 1.226; e = 0.145;
  }
  return std::pair<float, float> (sf, e);
}

float smearedPtDataMC(const TVector3& p, const vector<tree::Particle>& genJets, TRandom2& rand) {
  auto sfe = jerScales(p.Eta());
  float oldPt = p.Pt();
  float genPt = 0;
  for (const auto& genJ : genJets) {
    if (p.DeltaR(genJ.p) < .2 && std::abs(oldPt-genJ.p.Pt()) < 3*sfe.second*oldPt) {
      genPt = genJ.p.Pt();
      break;
    }
  }
  float newPt = genPt > 10 ?
      std::max((float)0., genPt + sfe.first*(oldPt-genPt))
    : rand.Gaus(oldPt, sqrt(pow(sfe.first,2)-1) * sfe.second*oldPt);
  return newPt;
}

int genMatchNegativePrompt(const tree::Particle& p, const std::vector<tree::GenParticle>& particles, float maxDr=.15, float maxDpt=.15) {
  for (auto const& genP : particles) {
    auto id = fabs(genP.pdgId);
    auto dr = p.p.DeltaR(genP.p);
    auto dpt = p.p.Pt()/genP.p.Pt();
    if (dr < maxDr && fabs(dpt-1) < maxDpt) {
      if (genP.isPrompt) return id;
      else return -id;
    }
  }
  return 0;
}

int genMatchWZDecay(const tree::Particle& p, const std::vector<tree::IntermediateGenParticle>& particles) {
  // match to daghters of massive particles
  for (auto const& genP : particles) {
    for (auto const & d : genP.daughters) {
      auto id = fabs(d.pdgId);
      auto dr = p.p.DeltaR(d.p);
      auto dpt = p.p.Pt()/d.p.Pt();
      if (dr < 0.4) return id;
    }
  }
  return 0;
}

bool genMatchToId(const tree::Particle& p, const std::vector<tree::GenParticle>& particles, unsigned pdgId, float maxDr=.15, float maxDpt=100) {
  for (auto const& genP : particles) {
    if (pdgId == fabs(genP.pdgId) && p.p.DeltaR(genP.p) < maxDr && fabs(p.p.Pt()/genP.p.Pt()-1)<maxDpt) {
      return true;
    }
  }
  return false;
}

bool genMatchToId(const tree::Particle& p, const std::vector<tree::IntermediateGenParticle>& particles, unsigned pdgId, float maxDr=.15, float maxDpt=100) {
  for (auto const& genP : particles) {
    for (auto const & d : genP.daughters) {
      if (pdgId == fabs(d.pdgId) && p.p.DeltaR(d.p) < maxDr && fabs(p.p.Pt()/d.p.Pt()-1)<maxDpt) {
        return true;
      }
    }
  }
  return false;
}

float topPtReweighting(std::vector<tree::GenParticle>& particles) {
  // scale factors from 8TeV combination https://twiki.cern.ch/twiki/bin/viewauth/CMS/TopPtReweighting
  // Take care that you correct for the additional weight, as this is not fixed to 1
  float a = 0.156, b = -0.00137;
  vector<tree::GenParticle*> selTops;
  for (auto& genP : particles) {
    if (fabs(genP.pdgId)==6) {
      selTops.push_back( &genP );
    }
  }
  float scaleFactor = 1;
  if (selTops.size()==2) {
    float pt1 = selTops.at(0)->p.Pt();
    float pt2 = selTops.at(1)->p.Pt();
    if (pt1>400) pt1=400;
    if (pt2>400) pt2=400;
    scaleFactor = sqrt( exp(a+b*pt1) * exp(a+b*pt2) );
  } else if (selTops.size() == 1 ) {
    float pt1 = selTops.at(0)->p.Pt();
    if (pt1>400) pt1=400;
    scaleFactor = exp(a+b*pt1);
  } else if (!selTops.size()) {
    //cout << "Counted " << selTops.size() << " tops, applying scale factor of 1 " << endl;
  }
  return scaleFactor;
}

string getSignalPointName(unsigned short nBinos, unsigned short m1, unsigned short m2) {
  string out = "";
  switch (nBinos) {
    case 0: out += "WW"; break;
    case 1: out += "Wg"; break;
    case 2: out += "gg"; break;
    default: out += "xx";
  }
  out += "_"+to_string(m1);
  out += "_"+to_string(m2);
  return out;
}

float ZGammaKFactor(std::vector<tree::GenParticle>& particles) {
  // k-factor from http://cms.cern.ch/iCMS/jsp/db_notes/noteInfo.jsp?cmsnoteid=CMS%20AN-2016/078 EXO-16-014
  float kfactor = 1.39; // kfactor for pt<175
  float pt=0;
  for (auto& genP : particles) {
    if (fabs(genP.pdgId)==22) {
      pt = genP.p.Pt();
      break;
    }
  }
  if (pt>400) kfactor = 1.23;
  else if (pt>250) kfactor = 1.30;
  else if (pt>190) kfactor = 1.35;
  return kfactor;
}

float transverseMass(const TVector3& v1, const TVector3& v2) {
    return TMath::Sqrt(2*(v1.Pt()*v2.Pt() - v1.X()*v2.X()-v1.Y()*v2.Y()));
}

float crystalResponse(float energy) {
  auto adc = energy / 0.04;
  // Taken from page 16 of https://indico.cern.ch/event/578798/contributions/2398504/
  float response = 1.;
  if (adc > 16000) {
    response = 0.88;
  } else if (adc > 8500) {
    response = 0.591574 + 0.000127616*adc -1.18194e-08*pow(adc,2) + 3.15593e-13*pow(adc,3)-0.0161;
  } else if (adc > 4000) {
    response = 0.486492 + 0.000304324*adc -5.4362e-08*pow(adc,2) +2.86294e-12*pow(adc,3)-0.01723;
  }
  if (response > 1) response = 1.;
  return response;
}

float getPropagatorPt(const vector<tree::GenParticle>& genParticles) {
  TVector3 propagator;
  int counter = 0;
  for (const auto& p: genParticles) {
    if (p.pdgId == 1000023 or abs(p.pdgId) == 1000024 or p.pdgId == 1000025) {
      propagator += p.p;
      counter++;
    }
  }
  return counter == 2 ? propagator.Pt() : 0;
}

float isrReweightingEWK(float pt, bool err=false) {
  // Taken from https://indico.cern.ch/event/616816/contributions/2489809/attachments/1418579/
  if (err) {
    if (pt>600) return 0.023;
    if (pt>400) return 0.012;
    if (pt>300) return 0.008;
    if (pt>200) return 0.004;
    if (pt>150) return 0.003;
    if (pt>100) return 0.002;
    if (pt>50) return 0.001;
    return 0.0001; // not stated, but should be really small
  } else {
    if (pt>600) return 0.783;
    if (pt>400) return 0.912;
    if (pt>300) return 1.000;
    if (pt>200) return 1.057;
    if (pt>150) return 1.150;
    if (pt>100) return 1.179;
    if (pt>50) return 1.052;
    return 1;
  }
}

float isrReweighting(unsigned nJet, bool err=false) {
  // Taken from https://indico.cern.ch/event/592621/contributions/2398559/
  if (err) {
    switch (nJet) {
      case 0: return 0;
      case 1: return 0.04;
      case 2: return 0.09;
      case 3: return 0.143;
      case 4: return 0.170;
      case 5: return 0.221;
      default: return 0.247;
    }
  } else {
    switch (nJet) {
      case 0: return 1;
      case 1: return 0.92;
      case 2: return 0.821;
      case 3: return 0.715;
      case 4: return 0.662;
      case 5: return 0.561;
      default: return 0.511;
    }
  }
}

bool unMatchedSuspiciousJet(const vector<tree::Jet>& jets, const vector<tree::Particle>& genJets) {
  for (const auto& j: jets) {
    if (fabs(j.p.Eta())<2.5 && j.chf<.1 && !count_if(genJets.begin(), genJets.end(), [&j] (const tree::Particle& p) { return p.p.DeltaR(j.p)<.3;})) return true;
  }
  return false;
}
