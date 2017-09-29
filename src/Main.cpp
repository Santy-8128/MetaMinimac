

#include <iostream>
#include <ctime>
#include <stdio.h>
#include "Parameters.h"
#include "StringBasics.h"
#include "MetaMinimac.h"

using namespace std;
void MetaMinimacVersion();
void helpFile();

int main(int argc, char ** argv)
{
	// Parameter Options

	String inputFiles = "";
	String info = "", snps = "",removeSam="";
	String outfile = "MetaMinimac.Output",method="B";
	String format = "DS";
	String type = "mach", errFile = "",chr="",golden="";
	String idDelimiter = "";
//	int max_indiv = 0, max_marker = 0;
    vector<bool> formatVector(3,false);

    int window=5000000,overlap=100000;
    double switchLimit=0.1;

	bool  gzip = true, nobgzip = false;
	bool   help = false, params = false;

	ParameterList inputParameters;
	PhoneHome::allThinning = 50;

	BEGIN_LONG_PARAMETERS(longParameterList)
		LONG_PARAMETER_GROUP("Input Files")
		LONG_STRINGPARAMETER("input", &inputFiles)
		//LONG_PARAMETER("rsid", &rsid)
		LONG_PARAMETER_GROUP("Output Parameters")
		LONG_STRINGPARAMETER("prefix", &outfile)
		LONG_STRINGPARAMETER("format", &format)
		LONG_PARAMETER("nobgzip", &nobgzip)
		LONG_PARAMETER_GROUP("Other Parameters")
		LONG_STRINGPARAMETER("method", &method)
		LONG_INTPARAMETER("window", &window)
		LONG_INTPARAMETER("overlap", &overlap)
		LONG_DOUBLEPARAMETER("switchLimit", &switchLimit)
		LONG_PARAMETER("help", &help)
		LONG_PARAMETER("params", &params)
		LONG_PHONEHOME(VERSION)
		BEGIN_LEGACY_PARAMETERS()
//		LONG_PARAMETER("onlyRefMarkers", &onlyRefMarkers)
//		LONG_STRINGPARAMETER("golden", &golden)
//		LONG_INTPARAMETER("sample", &max_indiv)
//		LONG_INTPARAMETER("marker", &max_marker)
//		LONG_STRINGPARAMETER("remove", &removeSam)
		END_LONG_PARAMETERS();

    int start_time = time(0);

	inputParameters.Add(new LongParameters(" Command Line Options: ",longParameterList));
	MetaMinimacVersion();
    String compStatus;
	inputParameters.Read(argc, &(argv[0]));
	if (help)
	{
		helpFile();
		return(-1);
	}

	inputParameters.Status();
    if(nobgzip)
        gzip=false;


	if (inputFiles == "")
	{
		cout<< " Missing \"--inputFiles\", a required parameter.\n\n";
		cout<< " Try \"--help\" for usage ...\n\n";
		cout<< " Program Exiting ...\n\n";
		compStatus="Command.Line.Error";
		PhoneHome::completionStatus(compStatus.c_str());
		return(-1);
	}

	if(format!="GP" && format!="DS")
    {
		cout<< " Invalid input for \"--format\" parameter : "<<format<<endl;
		cout<< " Parameter must be equal to \"GP\" or \"DS\" ..."<<endl<<endl;
		cout<< " Try \"--help\" for usage ...\n\n";
		cout<< " Program Exiting ...\n\n";
		compStatus="Command.Line.Error";
		PhoneHome::completionStatus(compStatus.c_str());
		return(-1);
	}



    if(window<=0)
    {
        cout<< " Invalid input for parameter \"--window\" : "<<window<<endl;
        cout<< " \"--window\"  can only take Positive Integers ...\n";
        cout<< " Try \"--help\" for usage ...\n\n";
        cout<< " Program Exiting ...\n\n";
        compStatus="Command.Line.Error";
        PhoneHome::completionStatus(compStatus.c_str());
        return(-1);
    }

    if(overlap<=0)
    {
        cout<< " Invalid input for parameter \"--overlap\" : "<<overlap<<endl;
        cout<< " \"--overlap\"  can only take Positive Integers ...\n";
        cout<< " Try \"--help\" for usage ...\n\n";
        cout<< " Program Exiting ...\n\n";
        compStatus="Command.Line.Error";
        PhoneHome::completionStatus(compStatus.c_str());
        return(-1);
    }

    MetaMinimac inputData(inputFiles,outfile,gzip,format,window,overlap,switchLimit,method);

    cout<<" ------------------------------------------------------------------------------"<<endl;
    cout<<"                             INPUT VCF DOSAGE FILE                             "<<endl;
    cout<<" ------------------------------------------------------------------------------"<<endl;

//
	if (!inputData.ReadInputFiles())
	{
		cout << "\n Program Exiting ... \n\n";
		compStatus="Input.VCF.Dose.Load.Error";
        PhoneHome::completionStatus(compStatus.c_str());
        return(-1);



	}

//
//	if (!inputData.MergeDosageData())
//	{
//		cout << "\n Program Exiting ... \n\n";
//		compStatus="Analysis.Error";
//        PhoneHome::completionStatus(compStatus.c_str());
//        return(-1);
//
//	}




    cout<<"\n ------------------------------------------------------------------------------"<<endl;
    cout<<"                                END OF PROGRAM                                 "<<endl;
    cout<<" ------------------------------------------------------------------------------"<<endl;

	int time_tot = time(0) - start_time;

    cout << "\n Program Successfully Implemented... \n ";


	printf("\n Total Run completed in %d hours, %d mins, %d seconds.\n",
		time_tot / 3600, (time_tot % 3600) / 60, time_tot % 60);

    cout<<"\n Thank You for using MetaMinimac !!! "<<endl<<endl;



    compStatus="Success";
    PhoneHome::completionStatus(compStatus.c_str());

	return 0;

}




void MetaMinimacVersion()
{
	printf("\n\n -------------------------------------------------------------------------------- \n");
	printf("          MetaMinimac - Converting Dosage from Minimac3 to other Formats     \n");
	printf(" --------------------------------------------------------------------------------\n");
    printf(" (c) 2015 - Sayantan Das \n");
//	printf(" Version	: Undocumented Release\n");
//	printf(" Built		: sayantan\n\n");
	cout<<"\n Version: " << VERSION<< ";\n Built: " << DATE << " by " << USER << std::endl;
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


