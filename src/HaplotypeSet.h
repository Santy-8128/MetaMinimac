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
    char refAllele,altAllele;
    string refAlleleString,altAlleleString;
    string MajAlleleString,MinAlleleString;
    bool swapped;




    int NoStudies;
    vector<int> StudyID;
    vector<int> StudyIDIndex;

    vector<vector<int> > StudyMap;

    variant(){};
    variant(string &id,string &CHR,int &BP)
    {
        name=id;
        bp=BP;
        chr=CHR;
    };
    void assignValues(string &id,string &CHR,int &BP)
    {
        name=id;
        bp=BP;
        chr=CHR;
    };
    void assignSwap(bool Swap)
    {
        swapped=Swap;
    };
     void assignRefAlt(string &refe,string &alt)
    {
        refAlleleString=refe;
        altAlleleString=alt;
    };
    void assignMajMin(string &refe,string &alt)
    {
        MajAlleleString=refe;
        MinAlleleString=alt;
    };
    void import(variant &FromVariant,int noStudies)
    {
        name=FromVariant.name;
        bp=FromVariant.bp;
        chr=FromVariant.chr;
        refAlleleString=FromVariant.refAlleleString;
        altAlleleString=FromVariant.altAlleleString;
        NoStudies=noStudies;
        StudyID.resize(NoStudies);
        StudyIDIndex.resize(NoStudies);
    };



};



class ThisChunk
{

	public:
	    int ChunkNo;
	    int NoVariants;
    int StartBp,EndBp;
    int Start,End;
	    int StartIndex,EndIndex;
	    int StartWithWindowIndex,EndWithWindowIndex;
        int NoGenoAllStudies;
        vector< vector<int> > ThisChunkInterAllTypedSitesReverseMap;

        ThisChunk()
        {

        };

};


class HaplotypeSet
{


	public:
		int         numHaplotypes,numSamples;
		int         numMarkers;
//        vector<vector<double> > Dosage;
        String InfilePrefix;
        String InfileName;

        string DoseFileName;
        string EmpDoseFileName;
        vector<string> markerName;
        vector<string> individualName;
    vector<int> SampleNoHaplotypes;
    vector<int> CummulativeSampleNoHaplotypes;



    vector<variant> VariantList;
    vector<variant> TypedVariantList;
    vector<variant> AllVariantList;
    int noTypedMarkers, noAllMarkers;


    vector<vector<double> > LooDosage;
        vector<vector<double> > HapDosage;
        vector<vector<double> > TypedGT;
        vector<vector<double> > MetaDosage;
        vector<vector<double> > RecomSpec;
        vector<vector<int> > FlankLength;
        vector<vector<double> > FlankFrac;
        double SplitLimit;


        vector<double> EstRsq;
        int BufferSize;
        String outFile;
        bool gzip;
        String Method;
        string finChromosome;

//        vector<vector<double> > LooDosage;
        vector<vector<double> > GTDosage;
        vector<int> TypedSites;
        vector<bool> TypedSitesIndicator;
        vector<bool> InterAllTypedSitesIndicator;
        vector<int> InterAllTypedSitesReverseMap;
        int TotalNoReadRecods;
        int TotalNoReadRecodsBegin;

        int NoTypedSites;




    void LoadEmpVariantList();
    void LoadVariantList(string inFile);

    bool CheckSampleConsistency(int tempNoSamples,
                                              vector<string> &tempindividualName,
                                              vector<int> tempSampleNoHaplotypes,
                                string File1, string File2);


    void ReadBasedOnSortCommonGenotypeList(vector<string> &SortedCommonGenoList);
    void SortCommonGenotypeList(
                                              std::unordered_set<string> &CommonGenotypeVariantNameList,
                                              vector<string> &SortedCommonGenoList,
                                              vector<variant> &CommonTypedVariantList);








//        bool DS,GP;
        vector<int>        optEndPoints;
		vector<int>        ScaffoldIndex;
		vector<int>        UnScaffoldIndex;
//		vector<ReducedHaplotypeInfo> ReducedStructureInfo;
		//vector<vector<char> >     haplotypes;
		vector<vector<char> >     haplotypesUnscaffolded;
		vector<vector<double> > alleleFreq;


//		vector<vector<float> > dosage;
		vector<vector<float> > GP1;
		vector<vector<float> > GP2;
		vector<vector<char> > ImputedAlleles;
		int PrintStartIndex,PrintEndIndex;
        bool onlyCompress,filter;
        bool EstimateNcompress;
        bool PseudoAutosomal;
        bool AllMaleTarget;
        String removeSam;
      	vector<char> refAlleleList,major, minor;
		vector<bool> missing, MarkerIndices;

		bool allowMissing, vcfType,m3vcfxType,machType;
        String VcfDose;
        String Info,IdDelimiter;
        String Format;
        String Type;
        vector<string> InPrefixList;
        int NoInPrefix;

        HaplotypeSet()
        {

        }

        HaplotypeSet(String vcfDose,String Outfile,bool Gzip,
                     String format,int bufferSize)
        {
            VcfDose=vcfDose;
            outFile=Outfile;
            gzip=Gzip;
            Format=format;
            BufferSize=bufferSize;

        }

    bool GetSummary(string prefix, myUserVariables &ThisVariable);
    bool CheckSuffixFile(string prefix, const char* suffix, string &FinalName);

        bool        UpdateFlankLength(int index,vector<double> tempRecom);
        bool        ProcessRecombination(int index,vector<double> tempRecom);
        bool        LoadHapDoseVariant(VcfRecordGenotype &ThisGenotype,int &numReadRecords);

        bool GetSampleInformation(string filename);



        void        PrintFlankLength();

        bool        LoadDosageChunkData(ThisChunk &MyChunk,VcfFileReader &ThisVcfFileStream,
                                             VcfRecord &ThisRecord);

        bool        LoadLooVariant(VcfRecordGenotype &ThisGenotype,int loonumReadRecords);

        bool        InitializeDosageChunkData(VcfFileReader &ThisVcfFileStream,
                                             VcfRecord &ThisRecord);
        bool        LoadInfoFile                        (string prefix);
		bool        LoadSampleNames                     (string prefix);
        bool        LoadRecomSpectrum                   (string prefix);
        bool        LoadDosageData                       (string prefix);
        bool        GetSummary                          (string prefix);

        bool        CheckMarkerCount                    (string prefix);

        bool        FastLoadHaplotypes                  ();
        void        PrintOutputVCF                       ();

		string      DetectReferenceFileType             (String filename);
//        bool        WriteMachFile                       ();
        bool        ReadInputFiles                      (string prefix);
        bool        doesExistFile                       (string filename);
        bool        printInfoError                      (int RowNo,string filename);
};









#endif // HAPLOTYPESET_H_INCLUDED
