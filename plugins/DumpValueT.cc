#include "HLTTrigger/TallinnHLTPFTauAnalyzer/plugins/DumpValueT.h"

typedef DumpValueT<bool> DumpBool;
typedef DumpValueT<double> DumpDouble;
typedef DumpValueT<float> DumpFloat;
typedef DumpValueT<int> DumpInt;

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE(DumpBool);
DEFINE_FWK_MODULE(DumpDouble);
DEFINE_FWK_MODULE(DumpFloat);
DEFINE_FWK_MODULE(DumpInt);





