#ifndef METAMINIMAC_H_INCLUDED
#define METAMINIMAC_H_INCLUDED

#include "StringBasics.h"
#include "HaplotypeSet.h"
#include "simplex.h"
#include <fstream>
#include <sstream>
#include <map>
#include <algorithm>
#include <iterator>
#include <math.h>
#include "MyVariables.h"
#include <iomanip>
#include <bitset>

using namespace std;

template<class Con>
void printcon(const Con& c){
  std::cout.precision(12);
  copy( c.begin(),
	c.end(),
	ostream_iterator<typename Con::value_type>(cout, "  ") );
  cout<<endl;
}

class MetaMinimac
{
	public:
        // Input Prefixes and Other User Variables
        myUserVariables ThisVariable;
        vector<String> InPrefixList;
        int NoInPrefix;
        vector<HaplotypeSet> InputData;


        // Variables for Common Typed Sites
        vector<variant> CommonTypedVariantList;
        int CommonTypedVariantListCounter;
        vector<string> SortedCommonGenoList;
        int NoCommonGenoVariants;


        // Variables for Chunking
        int NoChunks;
        vector<ThisChunk> ChunkList;


        // Variables for output dosage file stream
        IFILE vcfdosepartial, metaWeight;
        char *VcfPrintStringPointer;
        int VcfPrintStringPointerLength;


        // Variables for input dosage file stream and records
        vector<VcfFileReader*> InputDosageStream;
        vector<VcfRecord*> CurrentRecordFromStudy;
        vector<int> StudiesHasVariant;
        int CurrentFirstVariantBp;
        int NoStudiesHasVariant;


        // Variables for Meta imputed dosage
        vector<float> CurrentMetaImputedDosage;
        vector<float> CurrentMetaImputedDosageSum;
        vector<vector<double> > LSQEstimates;
        vector<double> ErrorPerSamplePerChunk;
        vector<double> ErrorPerSample;
        double ErrorSumSq,ErrorSumSqPerChunk;


    // FUNCTIONS


    String      AnalyzeExperiment               (myUserVariables &ThisVariables);
    bool        ParseInputVCFFiles              ();
    bool        CheckSampleNameCompatibility    ();
    bool        ReadEmpVariantAndChunk          ();
    bool        FindCommonGenotypedVariants     ();
    bool        CreateChunks                    ();
    void        GetNumChunks                    ();
    void        PrintChunkInformation           ();
    String      PerformFinalAnalysis            ();
    void        GetMetaImpEstimates             (int Sample, ThisChunk &MyChunk);
    void        OpenStreamInputDosageFiles      ();
    bool        OpenStreamOutputDosageFiles     ();
    void        PrintCurrentVariant             ();
    void        CreateMetaImputedData           ();
    bool        doesExistFile                   (String filename);
    string      GetDosageFileFullName           (String prefix);
    void        PrintMetaImputedData            ();
    int         IsVariantEqual                  (VcfRecord &Rec1, VcfRecord &Rec2);
    void        ReadCurrentDosageData           ();
    void        PrintHaploidDosage              (float &x);
    void        PrintDiploidDosage              (float &x, float &y);
    void        PrintVariant                    (VcfRecord *temp);
    void        AppendtoMainVcfFaster           (int ChunkNo);
    void        AppendtoMainWeightsFile           (int ChunkNo);
    void        UpdateCurrentRecords            ();
    void        Initialize                      ();
    void        FindCurrentMinimumPosition      ();
    void        PrintWeightForHaplotype         (int haploId);
    string      CreateInfo                      ();
                MetaMinimac                     (){}

};





#endif // METAMINIMAC_H_INCLUDED
