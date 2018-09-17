#ifndef HAPLOTYPESET_H_INCLUDED
#define HAPLOTYPESET_H_INCLUDED

#include "StringBasics.h"
#include "VcfFileReader.h"
#include "VcfHeader.h"
#include <fstream>
#include <sstream>
#include <unistd.h>
#include <list>
#include "MyVariables.h"

using namespace std;


class variant
{
public:

    string name;
    int bp;
    string chr;
    string refAlleleString,altAlleleString;

    variant(){};
    variant(string &id,string &CHR,int &BP)
    {
        name=id;
        bp=BP;
        chr=CHR;
    };
};



class ThisChunk
{

public:
    int ChunkNo;
    int NoVariants;
    int StartBp,EndBp;
    int StartWithWindowIndex,EndWithWindowIndex;
    int NoGenoAllStudies;
    vector< vector<int> > ThisChunkInterAllTypedSitesReverseMap;

    ThisChunk(){};
};


class HaplotypeSet
{


public:

    // File Name Variables
    String InfilePrefix;
    string DoseFileName;
    string EmpDoseFileName;

    // Summary Variables
    int         numHaplotypes,numSamples;
    int         numActualHaps;
    int         numMarkers;
    vector<string> individualName;
    vector<int> SampleNoHaplotypes;
    vector<int> CummulativeSampleNoHaplotypes;
    vector<variant> VariantList;
    vector<variant> TypedVariantList;
    int noTypedMarkers;
    string finChromosome;


    // Dosage Data
    vector<float> CurrentHapDosage;
    vector<vector<double> > LooDosage;
    vector<vector<double> > TypedGT;
    vector<vector<int> > FlankLength;
    vector<vector<double> > FlankFrac;


    // FUNCTIONS
    void        LoadEmpVariantList                      ();
    void        LoadVariantList                         (string inFile);
    bool        CheckSampleConsistency                  (int tempNoSamples, vector<string> &tempindividualName, vector<int> tempSampleNoHaplotypes, string File1, string File2);
    void        ReadBasedOnSortCommonGenotypeList       (vector<string> &SortedCommonGenoList);
    void        SortCommonGenotypeList                  (std::unordered_set<string> &CommonGenotypeVariantNameList, vector<string> &SortedCommonGenoList, vector<variant> &CommonTypedVariantList);
    bool        CheckSuffixFile                         (string prefix, const char* suffix, string &FinalName);
    void        LoadHapDoseVariant                      (VcfRecordGenotype &ThisGenotype);


    bool        GetSampleInformation                    (string filename);
    bool        GetSampleInformationfromHDS                    (string filename);
    void        LoadLooVariant                          (VcfRecordGenotype &ThisGenotype,int loonumReadRecords);
    bool        LoadSampleNames                         (string prefix);
    bool        doesExistFile                           (string filename);
};









#endif // HAPLOTYPESET_H_INCLUDED