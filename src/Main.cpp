

#include <iostream>
#include <ctime>
#include <unistd.h>
#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <stdio.h>
#include <getopt.h>
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
    bool log = false, help = false, params = false;

	ParameterList inputParameters;
	PhoneHome::allThinning = 50;

    int c;
    static struct option loptions[] =
    {
            {"input",required_argument,NULL,'i'},
            {"output",required_argument,NULL,'p'},
            {"format",required_argument,NULL,'f'},
            {"skipInfo",no_argument,NULL,'s'},
            {"nobgzip",no_argument,NULL,'n'},
            {"chunkLength",required_argument,NULL,'c'},
            {"help",no_argument,NULL,'h'},
            {"log",no_argument,NULL,'l'},
            {"debug",no_argument,NULL,'d'},
            {NULL,0,NULL,0}
    };

    while ((c = getopt_long(argc, argv, "i:p:f:c:snhld",loptions,NULL)) >= 0)
    {
        switch (c) {
            case 'i': MyVariables.inputFiles = optarg; break;
            case 'p': MyVariables.outfile = optarg; break;
            case 'f': MyVariables.formatString = optarg; break;
            case 'c': MyVariables.chunkLength = atoi(optarg); break;
            case 'n': MyVariables.nobgzip=true; break;
            case 's': MyVariables.infoDetails=false; break;
            case 'd': MyVariables.debug=true; break;
            case 'h': help=true; break;
            case 'l': log=true; break;
            case '?': helpFile(); return 1; break;
            default:  printf("[ERROR:] Unknown argument: %s\n", optarg);
        }
    }

    int start_time = time(0);
    MyVariables.CreateCommandLine(argc,argv);

    String compStatus;
    FILE *LogFile=NULL;
    if(log)
        LogFile=freopen(MyVariables.outfile +".logfile","w",stdout);
    dup2(fileno(stdout), fileno(stderr));

    MetaMinimacVersion();
	if (help)
	{
		helpFile();
		return(-1);
	}

	MetaMinimac myAnalysis;
	String MySuccessStatus="Error";

	MySuccessStatus = myAnalysis.AnalyzeExperiment(MyVariables);


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

}




void MetaMinimacVersion()
{
	printf("\n\n -------------------------------------------------- \n");
	printf("                   MetaMinimac     \n");
	printf(" --------------------------------------------------\n");
    printf(" (c) 2018 - Sayantan Das \n");
	cout<< " Version : " << VERSION<< ";\n Built   : " << DATE << " by " << USER << std::endl;
}

void helpFile()
{
    printf( "\n About   : Combine GWAS data imputed against different panels  \n");
    printf( " Usage   : metaMinimac [options] \n");
    printf( "\n");
    printf( " Options :\n");
    printf( "   -i, --input  <prefix1 prefix2 ...>  Prefixes of input data to meta-impute\n");
    printf( "   -o, --output <prefix>               Output prefix [MetaMinimac.Output] \n");
    printf( "   -f, --format <string>               Comma separated FORMAT tags [GT,DS]\n");
    printf( "   -s, --skipInfo                      Skip INFO in output [FALSE] \n");
    printf( "   -n, --nobgzip                       Output unzipped file [FALSE]\n");
    printf( "   -c, --chunkLength                   Length of chunk in mb [10] \n");
    printf( "   -b, --buffer                        Print Buffer [1e8] \n");
    printf( "   -d, --debug                         Debug mode [FALSE] \n");
    printf("\n URL = http://genome.sph.umich.edu/wiki/MetaMinimac\n");
    printf(" GIT = https://github.com/Santy-8128/MetaMinimac\n");
    printf("\n Visit website for more details ...\n");

    cout<<endl<<endl;
	return;
}


