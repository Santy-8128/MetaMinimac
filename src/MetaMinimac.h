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
        String outFile;
        String inputFile;
        String Method;
        bool gzip;
        vector<string> individualName;
        String Format;
        vector<String> InPrefixList;
        vector<HaplotypeSet> InputData;
        HaplotypeSet OutputData;
        int NoInPrefix;
        int Window,Overlap;
        string finChromosome;
        double SplitLimit;

        int numHaplotypes;
        int NoUnionVariants;
        vector<int> FinalBP;
        vector<int> NoVariantsByCategory;
        vector<vector<int> > ListofVariantsByCategory;
        vector<int> GenoVariantsAllStudies;
        int TotalNoGenoVariantsAllStudies;



        int NoChunks;
        vector<ThisChunk> ChunkList;



        vector<vector<double> > LSQEstimates;
        vector<double> ErrorPerSamplePerChunk;
        vector<double> ErrorPerSample;
        double ErrorSumSq,ErrorSumSqPerChunk;



        vector<int> StudyOneVariants;
        vector<int> StudyTwoVariants;
        vector<int> StudyBothVariants;
//        vector<int> StudyBothVariants;



        MetaMinimac(String InputFile,String Outfile,bool Gzip,
                     String format,int window, int overlap, double Limit, String method)
        {
            inputFile=InputFile;
            outFile=Outfile;
            gzip=Gzip;
            Format=format;
            Window=window;
            Overlap=overlap;
            Method=method;
            SplitLimit=Limit;

        }



void PrintLSQEstimates(int ID);
        bool        CreateChunks                    ();
        bool        AnalyzeChunks                    ();

        bool        CreateMetaRecomDosage                    ();
        bool        ImputeUniqueVariants                    ();
        bool        CreateLooFitDosage                   (ThisChunk &MyChunk);
        bool        CreateInfoR2Dosage                    ();
        bool        MergeVariantList                    ();
        bool        ReadInputFiles                      ();
        bool        MergeDosageData                     ();
        void        CopyParameters                      (HaplotypeSet &HapData,String InPrefix);
        void        CopyParameters                      (HaplotypeSet &FromHapData,HaplotypeSet &ToHapData);

void GetLeastSquareEstimates(int Sample, ThisChunk &MyChunk);
void GetLogOddsLeastSquareEstimates(int Sample, ThisChunk &MyChunk);



void UpdateFlankFractions();


void PrintVCFHeader();

void PrintVCFChunk(ThisChunk &MyChunk);

void AnalyzeLooDosage(ThisChunk &MyChunk);

};





#endif // METAMINIMAC_H_INCLUDED
