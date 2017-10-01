

#include <iostream>
#include <ctime>
#include <stdio.h>
#include "Parameters.h"
#include "StringBasics.h"
#include "MetaMinimac.h"
#include "MyVariables.h"


using namespace std;
void MetaMinimacVersion();
void helpFile();

int main(int argc, char ** argv)
{
	// Parameter Options
	myUserVariables MyVariables;

	ParameterList inputParameters;
	PhoneHome::allThinning = 50;

	BEGIN_LONG_PARAMETERS(longParameterList)
		LONG_PARAMETER_GROUP("Input Files")
		LONG_STRINGPARAMETER("input", &MyVariables.inputFiles)
		LONG_PARAMETER_GROUP("Output Parameters")
		LONG_STRINGPARAMETER("prefix", &MyVariables.outfile)
		LONG_STRINGPARAMETER("format", &MyVariables.format)
		LONG_PARAMETER("nobgzip", &MyVariables.nobgzip)
		LONG_PARAMETER_GROUP("Other Parameters")
		LONG_INTPARAMETER("window", &MyVariables.window)
		LONG_INTPARAMETER("overlap", &MyVariables.overlap)
		LONG_PARAMETER("help", &MyVariables.help)
		LONG_PARAMETER("log", &MyVariables.log)
		LONG_PARAMETER("params", &MyVariables.params)
		LONG_PHONEHOME(VERSION)
		BEGIN_LEGACY_PARAMETERS()
		LONG_STRINGPARAMETER("method", &MyVariables.method)
		LONG_DOUBLEPARAMETER("switchLimit", &MyVariables.switchLimit)
		END_LONG_PARAMETERS();

    int start_time = time(0);
	inputParameters.Add(new LongParameters(" Command Line Options: ",longParameterList));
	FILE *LogFile=NULL;
	dup2(fileno(stdout), fileno(stderr));
	MetaMinimacVersion();
    String compStatus;
	inputParameters.Read(argc, &(argv[0]));

	if(MyVariables.log)
		LogFile=freopen(MyVariables.outfile+".logfile","w",stdout);



	if (MyVariables.help)
	{
		helpFile();
		return(-1);
	}

	inputParameters.Status();


	MetaMinimac inputData;
	String MySuccessStatus="Error";

	MySuccessStatus = inputData.AnalyzeExperiment(MyVariables);


	if(MySuccessStatus!="Success")
	{
		compStatus=MySuccessStatus;
		PhoneHome::completionStatus(compStatus.c_str());
		return(-1);
	}

	int time_tot = time(0) - start_time;

	cout<<"\n ------------------------------------------------------------------------------"<<endl;
	cout<<"                                END OF PROGRAM                                 "<<endl;
	cout<<" ------------------------------------------------------------------------------"<<endl;


	cout << "\n Program Successfully Implemented... \n ";


	printf("\n Total Run completed in %d hours, %d mins, %d seconds.\n",
		   time_tot / 3600, (time_tot % 3600) / 60, time_tot % 60);

	cout<<"\n Thank You for using MetaMinimac !!! "<<endl<<endl;


	compStatus="Success";
	PhoneHome::completionStatus(compStatus.c_str());

	return 0;

//
//    MetaMinimac inputData(inputFiles,outfile,gzip,format,window,overlap,switchLimit,method);


//
//	if (!inputData.MergeDosageData())
//	{
//		cout << "\n Program Exiting ... \n\n";
//		compStatus="Analysis.Error";
//        PhoneHome::completionStatus(compStatus.c_str());
//        return(-1);
//
//	}






}




void MetaMinimacVersion()
{
	printf("\n\n -------------------------------------------------------------------------------- \n");
	printf("          MetaMinimac - Converting Dosage from Minimac3 to other Formats     \n");
	printf(" --------------------------------------------------------------------------------\n");
    printf(" (c) 2016 - Sayantan Das \n");
	cout<< " Version: " << VERSION<< ";\n Built: " << DATE << " by " << USER << std::endl;
}

void helpFile()
{

  printf("\n URL = http://genome.sph.umich.edu/wiki/MetaMinimac\n");
  printf(" GIT = https://github.com/Santy-8128/MetaMinimac\n");


    printf("\n MetaMinimac converts dosage files from minimac3 to other formats.\n");


	printf("\n Usage: ./MetaMinimac  --inDose      TestDataImputedVCF.dose.vcf.gz");
	printf("\n                           --info         TestDataImputedVCF.info");
	printf("\n                           --prefix       OutputFilePrefix");
	printf("\n                           --type         plink OR mach   // depending on output format");
	printf("\n                           --format       DS or GP        // based on if you want to output");
	printf("\n                                                          // dosage (DS) or genotype prob (GP)");
    printf("\n                           --window       1000000           // Number of Markers to import and ");
	printf("\n                                                          // print at a time (valid only for ");
	printf("\n                                                          // MaCH format)");
    printf("\n                           --idDelimiter  _               // Delimiter to Split VCF Sample ID into");
	printf("\n                                                          // FID and IID for PLINK format ");

    printf("\n\n URL = http://genome.sph.umich.edu/wiki/MetaMinimac\n");
    printf(" GIT = https://github.com/Santy-8128/MetaMinimac\n");
  printf("\n Visit website for more details ...\n");

cout<<endl<<endl;
	return;
}


