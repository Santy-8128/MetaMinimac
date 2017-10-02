#include "MetaMinimac.h"
#include "HaplotypeSet.h"
#include <math.h>
#include "Minimizers.h"
#include "MyVariables.h"
using BT::Simplex;

String MetaMinimac::AnalyzeExperiment(myUserVariables &ThisVariables)
{
    ThisVariable=ThisVariables;

    if(!ThisVariable.CheckValidity()) return "Command.Line.Error";


    cout<<" ------------------------------------------------------------------------------"<<endl;
    cout<<"                             INPUT VCF DOSAGE FILE                             "<<endl;
    cout<<" ------------------------------------------------------------------------------"<<endl;


    if (!ParseInputVCFFiles())
    {
        cout << "\n Program Exiting ... \n\n";
        return "Command.Line.Error";
    }

    if (!CheckSampleNameCompatibility())
    {
        cout << "\n Program Exiting ... \n\n";
        return "Input.VCF.Dose.Error";
    }

    if (!ReadEmpVariantAndChunk())
    {
        cout << "\n Program Exiting ... \n\n";
        return "Input.VCF.Dose.Error";
    }



    for(int i=0; i<NoChunks; i++)
    {

        for(int j=0;j<InputData[0].numHaplotypes; j++)
            GetMetaImpEstimates(j,ChunkList[i]);
    }

    return "Success";
}


void MetaMinimac::GetMetaImpEstimates(int Sample, ThisChunk &MyChunk)
{

    LogOddsModel ThisSampleAnalysis;
    ThisSampleAnalysis.metaInitialize(Sample,InputData,this,MyChunk);

    int h=0;

//    else if(Method=="B")
//    {
//
//
//        vector<double> init(NoInPrefix, 0.0);
////        LSQEstimates[Sample].resize(NoInPrefix+1);
//
//
//        LSQEstimates[Sample]=Simplex(ThisSampleAnalysis, init);
//        ErrorPerSamplePerChunk[Sample]=ThisSampleAnalysis(LSQEstimates[Sample]);
//        ErrorPerSample[Sample]+=ErrorPerSamplePerChunk[Sample];
//        ErrorSumSqPerChunk+=ErrorPerSamplePerChunk[Sample];
////        logTransform(MiniMizer,LSQEstimates[Sample],NoInPrefix);
//
//    }

//    LeastSquareError ThisSampleAnalysis;
//    ThisSampleAnalysis.initialize(Sample,InputData,this,MyChunk);
//
//    vector<double> init(NoInPrefix, 0.0);
//    //        cout<<"\n WELL \n";
//    LSQEstimates[Sample].resize(NoInPrefix);
//    PositiveTransform(Simplex(ThisSampleAnalysis, init),LSQEstimates[Sample],NoInPrefix);
//


}


bool MetaMinimac::ReadEmpVariantAndChunk()
{
    cout<<"\n Gathering information for chunking ... "<<endl;
    for(int i=0;i<NoInPrefix;i++)
    {
        InputData[i].LoadEmpVariantList();
        cout<<" No Genotyped Sites = "<<InputData[i].noTypedMarkers<<endl;
    }

    cout<<" Successful !!! "<<endl;


    FindCommonGenotypedVariants();

    GetNumChunks();
    CreateChunks();

    return true;
}


void MetaMinimac::GetNumChunks() {

    int NoMarkers = CommonTypedVariantList.size();
    int StartPos = CommonTypedVariantList[0].bp;
    int EndPos = CommonTypedVariantList[NoMarkers - 1].bp;

    NoChunks = ((int) (EndPos - StartPos) / (int) (1000000 * ThisVariable.chunkLength));

    if (NoChunks == 0)
        NoChunks = 1;

    return;
}





bool MetaMinimac::ParseInputVCFFiles()
{

    InPrefixList.clear();
    size_t pos = 0;
    std::string delimiter(ThisVariable.FileDelimiter) ;
    std::string token;
    int Count=0;
    string tempName=ThisVariable.inputFiles.c_str();
    while ((pos = tempName.find(delimiter)) != std::string::npos)
    {
        token = tempName.substr(0, pos);
        cout<<  " Study "<<Count+1<<" Prefix                     : "<<token<<endl;

        InPrefixList.push_back(token.c_str());
        tempName.erase(0, pos + delimiter.length());
        Count++;
    }
    cout<<  " Study "<<Count+1<<" Prefix                     : "<<tempName<<endl;
    InPrefixList.push_back(tempName.c_str());


    NoInPrefix=(int)InPrefixList.size();
    InputData.clear();
    InputData.resize(NoInPrefix);

    cout<<  " Number of Studies to Meta-Impute   : "<<NoInPrefix<<endl;

    if(NoInPrefix<2)
    {
        cout<<"\n ERROR ! Must have at least 2 studies for meta-imputation to work !!! "<<endl;
        cout<<" Program aborting ... "<<endl<<endl;
        return false;
    }
    if(NoInPrefix>4)
    {
        cout<<"\n ERROR ! Must have less than 5 studies for meta-imputation to work !!! "<<endl;
        cout<<" Program aborting ... "<<endl<<endl;
        return false;
    }

    return true;
}


bool MetaMinimac::CheckSampleNameCompatibility()
{
    cout<<"\n Checking Sample Compatibility across files ... "<<endl;

    for(int i=0;i<NoInPrefix;i++)
    {
        if(!InputData[i].LoadSampleNames(InPrefixList[i].c_str()))
            return false;
        if(i>0)
            if(!InputData[i].CheckSampleConsistency(InputData[i-1].numSamples,
                                                    InputData[i-1].individualName,
                                                    InputData[i-1].SampleNoHaplotypes,
                                                    InputData[i-1].DoseFileName,
                                                    InputData[i].DoseFileName))
                return false;

    }

    cout<<" Successful !!! "<<endl;
    return true;
}



bool MetaMinimac::ReadInputFiles()
{
    cout<<"\n ----------------------------------"<<endl;
    cout<<"  READING SAMPLE AND VARIANT NAMES "<<endl;
    cout<<" ----------------------------------"<<endl<<endl;


    for(int i=0;i<NoInPrefix;i++)
    {
        //CopyParameters(InputData[i],InPrefixList[i]);
        if(!InputData[i].LoadSampleNames(InPrefixList[i].c_str()))
            return false;
//        if(i>0)
//            if(!InputData[i].CheckSampleConsistency(InputData[i-1].numSamples,
//                                                       InputData[i-1].individualName,
//                                                       InputData[i-1].SampleNoHaplotypes))
//                return false;

    }

    cout<<"\n --------------------------------"<<endl;
    cout<<"  FINDING COMMON SET OF VARIANTS "<<endl;
    cout<<" --------------------------------"<<endl<<endl;

    CopyParameters(InputData[0],OutputData);
    MergeVariantList();

    numHaplotypes=InputData[0].numHaplotypes;
    finChromosome=InputData[0].VariantList[0].chr;
    individualName=InputData[0].individualName;


    if(Method!="A")
    {

        cout<<"\n --------------------------------------"<<endl;
        cout<<"  READING RECOMBINATION SPECTRUM FILES "<<endl;
        cout<<" --------------------------------------"<<endl<<endl;

        for(int i=0;i<NoInPrefix;i++)
        {
            if(!InputData[i].LoadRecomSpectrum(InPrefixList[i].c_str()))
                return false;

        }
        UpdateFlankFractions();

        for(int i=0;i<NoInPrefix;i++)
        {
            InputData[i].PrintFlankLength();
        }



    }




    LSQEstimates.resize(numHaplotypes);
    ErrorPerSamplePerChunk.resize(numHaplotypes);
    ErrorPerSample.resize(numHaplotypes,0.0);
    ErrorSumSq=0.0;


    cout<<"\n ---------------"<<endl;
    cout<<"  CHUNK SUMMARY "<<endl;
    cout<<" ---------------"<<endl<<endl;

    CreateChunks();
    if(!AnalyzeChunks())
        return false;

    cout<<" Total Sum of Squares for all Chunks = "<<ErrorSumSq<<endl;


    return true;



    for(int i=0;i<NoInPrefix;i++)
    {
        cout<<"\n ----------------------------"<<endl;
        cout<<"       INPUT FILE : ["<<i+1<<"] "<<endl;
        cout<<" ----------------------------"<<endl;

        cout<<" File Prefix        : "<<InPrefixList[i] <<endl;

//        CopyParameters(InputData[i]);

//        if(!InputData[i].ReadInputFiles(InPrefixList[i]))
//            return false;



    }




    return true;
}



bool MetaMinimac::CreateChunks()
{

    int NoMarkers=CommonTypedVariantList.size();

    ThisVariable.chunkLength=(NoMarkers/NoChunks);
    ChunkList.resize(NoChunks);
    ChunkList[0].StartBp=0;
    ChunkList[0].StartWithWindowIndex=0;

    int tempNoChunkMarkers = (NoMarkers/NoChunks);
    int tempNoWindowMarkers = tempNoChunkMarkers/10;


    int counter=0, chunkCounter = 0;

    while(counter<NoCommonGenoVariants)
    {
        ThisChunk &tempChunk = ChunkList[chunkCounter];

        tempChunk.EndWithWindowIndex=counter+tempNoChunkMarkers+tempNoWindowMarkers;
        tempChunk.NoGenoAllStudies=tempChunk.EndWithWindowIndex-tempChunk.StartWithWindowIndex+1;
        tempChunk.ChunkNo=chunkCounter;


        if(chunkCounter>0)
        {
            tempChunk.StartWithWindowIndex=counter - tempNoWindowMarkers;
        }

        if(chunkCounter<NoChunks-1)
        {
            tempChunk.EndBp=CommonTypedVariantList[counter+tempNoChunkMarkers].bp;
            ChunkList[chunkCounter+1].StartBp=tempChunk.EndBp+1;
        }
        else
        {
            tempChunk.EndBp=999999999;
            tempChunk.EndWithWindowIndex=NoCommonGenoVariants-1;
            tempChunk.NoGenoAllStudies=tempChunk.EndWithWindowIndex-tempChunk.StartWithWindowIndex+1;
        }

        counter+=tempNoChunkMarkers;
        chunkCounter++;

    }

    return true;
}




bool MetaMinimac::FindCommonGenotypedVariants()
{
    std::map<string, int > HashUnionVariantMap;
    std::unordered_set<string> CommonGenotypeVariantNameList;

    for(int i=0;i<NoInPrefix;i++)
    {
        for(int j=0;j<InputData[i].noTypedMarkers;j++)
        {
            variant *thisVariant=&InputData[i].TypedVariantList[j];
            HashUnionVariantMap[thisVariant->name]++;
        }
    }

    std::map<string, int >::iterator itStudy;

    for (itStudy=HashUnionVariantMap.begin(); itStudy!=HashUnionVariantMap.end(); ++itStudy)
    {
        if(itStudy->second==NoInPrefix)
        {
            CommonGenotypeVariantNameList.insert(itStudy->first);
        }

    }
    NoCommonGenoVariants = CommonGenotypeVariantNameList.size();



    InputData[0].SortCommonGenotypeList(CommonGenotypeVariantNameList, SortedCommonGenoList, CommonTypedVariantList);

    for(int i=1; i<NoInPrefix; i++)
        InputData[i].ReadBasedOnSortCommonGenotypeList(SortedCommonGenoList);

//
//
//    std::cout << " Total Union Number of Variants                      : " << OutputData.numMarkers << endl;
//
//    int MultipleSitesCount=0;
//
//    for(int i=0;i<NoInPrefix;i++)
//    {
//        std::cout << " Total Number of Variants in "<<i+1<<" Studies               : "
//                  << NoVariantsByCategory[i] << endl;
//        if(i>0)
//            MultipleSitesCount+=NoVariantsByCategory[i];
//    }
//    std::cout << " Total Number of Variants in more than one study     : " << MultipleSitesCount<< endl <<endl;

}




bool MetaMinimac::MergeVariantList()
{
    std::map<int,vector<int> > HashUnionStudyMap;
    std::map<int,vector<int> > HashUnionIndexMap;

    for(int i=0;i<NoInPrefix;i++)
    {
        for(int j=0;j<InputData[i].numMarkers;j++)
        {
            variant *thisVariant=&InputData[i].VariantList[j];
            HashUnionStudyMap[thisVariant->bp].push_back(i);
            HashUnionIndexMap[thisVariant->bp].push_back(j);
        }
    }

    vector<int > StudyClassify;

    OutputData.VariantList.clear();
    int TotalUnionSites=0;

    std::map<int,vector<int> >::iterator itStudy,itIndex=HashUnionIndexMap.begin();
    std::map<string,vector<int> >::iterator itAlt;

    FinalBP.resize(HashUnionStudyMap.size());

    for (itStudy=HashUnionStudyMap.begin(); itStudy!=HashUnionStudyMap.end(); ++itStudy,++itIndex)
    {
        FinalBP[TotalUnionSites]=itStudy->first;
        TotalUnionSites++;

        vector<int> &StudyMapper= itStudy->second;
        vector<int> &IndexMapper= itIndex->second;
        std::map<string,vector<int> > AltAlleleMap;
        string RefString="";


        for(int Counter=0;Counter<(int)StudyMapper.size();Counter++)
        {
            int ThisStudy=StudyMapper[Counter];
            int ThisIndex=IndexMapper[Counter];
            variant *thisVariant=&InputData[ThisStudy].VariantList[ThisIndex];
            if(RefString!="" && thisVariant->refAlleleString!=RefString)
            {
                std::cout << "\n ERROR !!! Mimatching Reference Allele for Position "<<thisVariant->bp << endl;
                return false;
            }

            AltAlleleMap[thisVariant->refAlleleString+"_"+thisVariant->altAlleleString].push_back(ThisStudy);
            AltAlleleMap[thisVariant->refAlleleString+"_"+thisVariant->altAlleleString].push_back(ThisIndex);
        }

        for (itAlt=AltAlleleMap.begin(); itAlt!=AltAlleleMap.end(); ++itAlt)
        {
            vector<int> &AltMapper = itAlt->second;

            variant tempVariant;
            tempVariant.import(InputData[AltMapper[0]].VariantList[AltMapper[1]],(int)(AltMapper.size()/2));

            int Counter=0;
            while(Counter<(int)AltMapper.size())
            {
                tempVariant.StudyID[Counter/2]=AltMapper[Counter];
                tempVariant.StudyIDIndex[Counter/2]=AltMapper[Counter+1];
                Counter+=2;
            }
            OutputData.VariantList.push_back(tempVariant);
        }
    }


    OutputData.numMarkers=(int)OutputData.VariantList.size();
    OutputData.numHaplotypes=InputData[0].numHaplotypes;


    NoVariantsByCategory.resize(NoInPrefix,0);
    ListofVariantsByCategory.resize(NoInPrefix);

    for(int i=0;i<OutputData.numMarkers;i++)
    {
        variant *thisVariant=&OutputData.VariantList[i];
        NoVariantsByCategory[thisVariant->NoStudies-1]++;
        ListofVariantsByCategory[thisVariant->NoStudies-1].push_back(i);


        if(thisVariant->NoStudies==NoInPrefix)
        {
            bool AllGenotyped=true;
            int k=0;
            while(AllGenotyped && k<NoInPrefix)
            {
                int StuId=thisVariant->StudyID[k];
                int MarId=thisVariant->StudyIDIndex[k];
                AllGenotyped=InputData[StuId].TypedSitesIndicator[MarId];
                k++;
            }
            if(AllGenotyped)
            {
                for(k=0;k<NoInPrefix;k++)
                {
                    int StuId=thisVariant->StudyID[k];
                    int MarId=thisVariant->StudyIDIndex[k];
                    InputData[StuId].InterAllTypedSitesIndicator[MarId]=true;
                    InputData[StuId].InterAllTypedSitesReverseMap.push_back(MarId);
                }
                GenoVariantsAllStudies.push_back(i);
            }
        }
    }

    TotalNoGenoVariantsAllStudies=(int)GenoVariantsAllStudies.size();

    std::cout << " Total Union Number of Variants                      : " << OutputData.numMarkers << endl;

    int MultipleSitesCount=0;

    for(int i=0;i<NoInPrefix;i++)
    {
        std::cout << " Total Number of Variants in "<<i+1<<" Studies               : "
                  << NoVariantsByCategory[i] << endl;
        if(i>0)
            MultipleSitesCount+=NoVariantsByCategory[i];
    }
    std::cout << " Total Number of Variants in more than one study     : " << MultipleSitesCount<< endl <<endl;

}




double rosenbrock(vector<double> x)
{
  return -pow(x[0],1/x[0]);
}


void MetaMinimac::PrintVCFHeader()
{

    IFILE vcfdosepartial=NULL;
    vcfdosepartial = ifopen(outFile + ".metaDose.vcf" + (gzip ? ".gz" : ""), "wb", gzip ?InputFile::BGZF:InputFile::UNCOMPRESSED);


    ifprintf(vcfdosepartial,"##fileformat=VCFv4.1\n");
    time_t t = time(0);
    struct tm * now = localtime( & t );
    ifprintf(vcfdosepartial,"##filedate=%d.%d.%d\n",(now->tm_year + 1900),(now->tm_mon + 1) ,now->tm_mday);
    ifprintf(vcfdosepartial,"##source=MetaMinimac\n");
    ifprintf(vcfdosepartial,"##contig=<ID=%s>\n",finChromosome.c_str());


    ifprintf(vcfdosepartial,"##FORMAT=<ID=HDS,Number=1,Type=String,Description=\"Estimated Alternate Allele Dosage : [P(0/1)+2*P(1/1)]\">\n");

    ifprintf(vcfdosepartial,"##FORMAT=<ID=GP,Number=3,Type=Float,Description=\"Estimated Posterior Probabilities for Genotypes 0/0, 0/1 and 1/1 \">\n");


//    ifprintf(vcfdosepartial,"##INFO=<ID=MAF,Number=1,Type=Float,Description=\"Estimated Alternate Allele Frequency\">\n");
//    ifprintf(vcfdosepartial,"##INFO=<ID=R2,Number=1,Type=Float,Description=\"Estimated Imputation Accuracy\">\n");
//    ifprintf(vcfdosepartial,"##INFO=<ID=ER2,Number=1,Type=Float,Description=\"Empirical (Leave-One-Out) R-square (available only for genotyped variants)\">\n");

    ifprintf(vcfdosepartial,"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT");

    for(int hapId=0;hapId<(int)individualName.size();hapId++)
    {
        ifprintf(vcfdosepartial,"\t%s",individualName[hapId].c_str());
    }

    ifprintf(vcfdosepartial,"\n");
    ifclose(vcfdosepartial);

}




void MetaMinimac::PrintVCFChunk(ThisChunk &MyChunk)
{

    IFILE vcfdosepartial=NULL;
    vcfdosepartial = ifopen(outFile + ".metaDose.vcf" + (gzip ? ".gz" : ""), "a", gzip ?InputFile::BGZF:InputFile::UNCOMPRESSED);


    double p1,p2,p3,x,y;
    for(int i=MyChunk.StartIndex;i<=MyChunk.EndIndex;i++)
    {

        int index=i-MyChunk.StartIndex;
        variant &thisVariant=OutputData.VariantList[i];

        ifprintf(vcfdosepartial,"%s\t%d\t%s\t%s\t%s\t.\tPASS\t.\tHDS:GP",
        thisVariant.chr.c_str(),thisVariant.bp,
        thisVariant.name.c_str(),thisVariant.refAlleleString.c_str(),
        thisVariant.altAlleleString.c_str());

        for(int hapId=0;hapId<numHaplotypes;hapId+=2)
        {
            x=OutputData.MetaDosage[index][hapId];
            y=OutputData.MetaDosage[index][hapId+1];
            p3=x*y;
            p2=x*(1-y)+y*(1-x);
            p1=(1-x)*(1-y);
            ifprintf(vcfdosepartial,"\t%.3f|%.3f:%.3f,%.3f,%.3f",x,y,p1,p2,p3);
        }
        ifprintf(vcfdosepartial,"\n");
    }

    ifclose(vcfdosepartial);
}





void PositiveTransform(vector<double> From,vector<double> &To,int NoDimensions)
{

    for(int i=0; i <=NoDimensions; i++)
    {
        To[i]=exp(From[i]);
    }
}

void MetaMinimac::GetLogOddsLeastSquareEstimates(int Sample, ThisChunk &MyChunk)
{

    LogOddsModel ThisSampleAnalysis;
    ThisSampleAnalysis.initialize(Sample,InputData,this,MyChunk);


    if(Method=="A")
    {
        vector<double> init(NoInPrefix-1, 0.0);
        LSQEstimates[Sample].resize(NoInPrefix);

        LSQEstimates[Sample]=Simplex(ThisSampleAnalysis, init);

        ErrorPerSamplePerChunk[Sample]=ThisSampleAnalysis(LSQEstimates[Sample]);
        ErrorPerSample[Sample]+=ErrorPerSamplePerChunk[Sample];
        ErrorSumSqPerChunk+=ErrorPerSamplePerChunk[Sample];

//        logTransform(MiniMizer,LSQEstimates[Sample],NoInPrefix-1);

    }


    else if(Method=="B")
    {


        vector<double> init(NoInPrefix, 0.0);
//        LSQEstimates[Sample].resize(NoInPrefix+1);


        LSQEstimates[Sample]=Simplex(ThisSampleAnalysis, init);
        ErrorPerSamplePerChunk[Sample]=ThisSampleAnalysis(LSQEstimates[Sample]);
        ErrorPerSample[Sample]+=ErrorPerSamplePerChunk[Sample];
        ErrorSumSqPerChunk+=ErrorPerSamplePerChunk[Sample];
//        logTransform(MiniMizer,LSQEstimates[Sample],NoInPrefix);

    }
    else if(Method=="C")
    {

        vector<double> init(2*(NoInPrefix-1),0.0);
//        LSQEstimates[Sample].resize(NoInPrefix+1);
        LSQEstimates[Sample]=Simplex(ThisSampleAnalysis, init);


        ErrorPerSamplePerChunk[Sample]=ThisSampleAnalysis(LSQEstimates[Sample]);
        ErrorPerSample[Sample]+=ErrorPerSamplePerChunk[Sample];
        ErrorSumSqPerChunk+=ErrorPerSamplePerChunk[Sample];

//        logTransform(init,LSQEstimates[Sample],NoInPrefix);

    }
//
//    LeastSquareError ThisSampleAnalysis;
//    ThisSampleAnalysis.initialize(Sample,InputData,this,MyChunk);
//
//    vector<double> init(NoInPrefix, 0.0);
//    //        cout<<"\n WELL \n";
//    LSQEstimates[Sample].resize(NoInPrefix);
//    PositiveTransform(Simplex(ThisSampleAnalysis, init),LSQEstimates[Sample],NoInPrefix);
//


}



void MetaMinimac::GetLeastSquareEstimates(int Sample, ThisChunk &MyChunk)
{



    if(Method=="A")
    {


        LeastSquareError ThisSampleAnalysis;
//        LogOddsModel ThisSampleAnalysis;


        ThisSampleAnalysis.initialize(Sample,InputData,this,MyChunk);

        vector<double> init(NoInPrefix-1, -log((double)NoInPrefix-1));



        LSQEstimates[Sample].resize(NoInPrefix);
        vector<double> MiniMizer=Simplex(ThisSampleAnalysis, init);
        ErrorPerSamplePerChunk[Sample]=ThisSampleAnalysis(MiniMizer);
        ErrorPerSample[Sample]+=ErrorPerSamplePerChunk[Sample];
        ErrorSumSqPerChunk+=ErrorPerSamplePerChunk[Sample];
        logTransform(MiniMizer,LSQEstimates[Sample],NoInPrefix-1);

    }


    else if(Method=="B")
    {

        LeastSquareWithRecombination ThisSampleAnalysis;
//        LogOddsModel ThisSampleAnalysis;


        ThisSampleAnalysis.initialize(Sample,InputData,this,MyChunk);

        vector<double> init(NoInPrefix, -log((double)NoInPrefix));
        LSQEstimates[Sample].resize(NoInPrefix+1);
        vector<double> MiniMizer=Simplex(ThisSampleAnalysis, init);
        ErrorPerSamplePerChunk[Sample]=ThisSampleAnalysis(MiniMizer);
        ErrorPerSample[Sample]+=ErrorPerSamplePerChunk[Sample];
        ErrorSumSqPerChunk+=ErrorPerSamplePerChunk[Sample];
        logTransform(MiniMizer,LSQEstimates[Sample],NoInPrefix);

    }
    else if(Method=="C")
    {

        LeastSquareWithRecombination ThisSampleAnalysis;
        ThisSampleAnalysis.initialize(Sample,InputData,this,MyChunk);

        vector<double> init(NoInPrefix, -log((double)NoInPrefix));
        LSQEstimates[Sample].resize(NoInPrefix+1);
//        vector<double> MiniMizer=Simplex(ThisSampleAnalysis, init);


        ErrorPerSamplePerChunk[Sample]=ThisSampleAnalysis.FixedSSQ();
        ErrorPerSample[Sample]+=ErrorPerSamplePerChunk[Sample];
        ErrorSumSqPerChunk+=ErrorPerSamplePerChunk[Sample];

        logTransform(init,LSQEstimates[Sample],NoInPrefix);

    }
//
//    LeastSquareError ThisSampleAnalysis;
//    ThisSampleAnalysis.initialize(Sample,InputData,this,MyChunk);
//
//    vector<double> init(NoInPrefix, 0.0);
//    //        cout<<"\n WELL \n";
//    LSQEstimates[Sample].resize(NoInPrefix);
//    PositiveTransform(Simplex(ThisSampleAnalysis, init),LSQEstimates[Sample],NoInPrefix);
//


}


void MetaMinimac::AnalyzeLooDosage(ThisChunk &MyChunk)
{

    ErrorSumSqPerChunk=0.0;

    for(int i=0;i<numHaplotypes;i++)
    {



//        GetLeastSquareEstimates(i,MyChunk);
        GetLogOddsLeastSquareEstimates(i,MyChunk);

    }
    cout<<" Total Sum of Squares for this Chunk = "<<ErrorSumSqPerChunk<<endl;

    ErrorSumSq+=ErrorSumSqPerChunk;
}





bool MetaMinimac::CreateLooFitDosage(ThisChunk &MyChunk)
{
    OutputData.MetaDosage.resize(MyChunk.NoVariants);
    for(int i=0;i<MyChunk.NoVariants;i++)
    {
        OutputData.MetaDosage[i].clear();
        OutputData.MetaDosage[i].resize(numHaplotypes,0.0);
    }


    for(int i=MyChunk.StartIndex;i<=MyChunk.EndIndex;i++)
    {

        int index=i-MyChunk.StartIndex;
        variant *thisVariant=&OutputData.VariantList[i];


        if(thisVariant->NoStudies==1)
        {

            int StudyID=thisVariant->StudyID[0];
            int MarkerID=thisVariant->StudyIDIndex[0];

            vector<double> &FromDosage=InputData[StudyID].HapDosage[MarkerID-InputData[StudyID].TotalNoReadRecodsBegin];
            vector<double> &CopyToDosage=(OutputData.MetaDosage[index]);

            CopyToDosage=FromDosage;


        }
        else
        {

            vector<double> &CopyToDosage=(OutputData.MetaDosage[index]);
            vector<double> tempVar(thisVariant->NoStudies);
            vector< vector<double> > Covariates;
            int NoDimensions=thisVariant->NoStudies;


            for(int j=0;j<numHaplotypes;j++)
            {

                double FlankEffect=0.0;


                if(Method=="A")
                {
                    logitTransform(LSQEstimates[j],tempVar,Covariates,thisVariant->NoStudies,Method);
                }

                if(Method=="B" || Method =="C")
                {
                    Covariates.clear();
                    Covariates.resize(1);
                    Covariates[0].resize(NoDimensions-1);


                    for(int StudyIndex=0;StudyIndex<(thisVariant->NoStudies-1);StudyIndex++)
                    {
                        int StudyID=thisVariant->StudyID[StudyIndex];
                        int MarkerID=thisVariant->StudyIDIndex[StudyIndex];
                        Covariates[0][StudyIndex]=(InputData[StudyID].FlankFrac[j][MarkerID]);
                    }

                    logitTransform(LSQEstimates[j],tempVar,Covariates,NoDimensions,Method);

                }



                for(int StudyIndex=0;StudyIndex<thisVariant->NoStudies;StudyIndex++)
                {

                    int StudyID=thisVariant->StudyID[StudyIndex];
                    int MarkerID=thisVariant->StudyIDIndex[StudyIndex];

                    vector<double> &FromDosage=InputData[StudyID].HapDosage[MarkerID-InputData[StudyID].TotalNoReadRecodsBegin];


//
//
//                    if (Method=="B")
//                    {
//                        vector<double> &FromFlankFrac=InputData[StudyID].FlankFrac[j];
//
//
//                        FlankEffect=LSQEstimates[j][NoInPrefix]*(FromFlankFrac[MarkerID]);
//
//
//                    }
//                    else if(Method=="C")
//                    {
//                        vector<double> &FromFlankFrac=InputData[StudyID].FlankFrac[j];
//
//
//                        FlankEffect=(FromFlankFrac[MarkerID]);
//                        LSQEstimates[j][StudyID]=0.0;
//
////                    if(StudyIndex==0)
////                        cout<<"SAMPLE\t"<<j<<"\t"<<MarkerID<<"\t"<<FlankEffect<<endl;
//
//
//                    }
//

                    CopyToDosage[j]+=( FromDosage[j] * (tempVar[StudyIndex]) ) ;

                    if(StudyIndex==0)
                        cout<<"SAMPLE\t"<<j<<"\t"<<MarkerID<<"\t"<<tempVar[StudyIndex]<<endl;
                }

                if(CopyToDosage[j]>1.0000001 || CopyToDosage[j]<0.0)
                {
                    cout<<" LARGE = "<<i<<"\t"<<j<<"\t"<<CopyToDosage[j]<<endl;
//                    abort();


                }

            }








//            cout<<" JESSE ";

//            for(int j=0;j<numHaplotypes;j++)
//            {
//                double sum=0.0;
//                for(int StudyIndex=0;StudyIndex<thisVariant->NoStudies;StudyIndex++)
//                {
//                    int StudyID=thisVariant->StudyID[StudyIndex];
//                    sum+=LSQEstimates[j][StudyID];
////                     if(j==4)
////                        cout<<LSQEstimates[j][StudyID]<<"\t"<<sum<<"\t";
//
//                }
//                CopyToDosage[j]/=sum;
//
////                if(CopyToDosage[j]<0)
////                    cout<<CopyToDosage[j]<<endl;
//////
////                if(CopyToDosage[j]>1)
////                    cout<<" DERP "<<j<<"\t"<<CopyToDosage[j]<<endl;
//
//            }
        }
    }

            return true;

//
//
//    for(int i=0;i<(OutputData.numHaplotypes);i++)
//    {
//
//        vector<double> *FromDosage1=&(InputData[0].LooDosage[i]);
//        vector<double> *FromDosage2=&(InputData[1].LooDosage[i]);
//        vector<double> *FromGT=&(InputData[1].GTDosage[i]);
//        vector<double> *CopyToDosage=&(OutputData.HapDosage[i]);
//        vector<double> *CopyFromDosage1=&(InputData[0].HapDosage[i]);
//        vector<double> *CopyFromDosage2=&(InputData[1].HapDosage[i]);
//
//
//
//        double sumNum=0.0,sumDen=0.0;
//
//        for(int j=0;j<InputData[0].NoTypedSites;j++)
//        {
//            sumDen+=(((*FromDosage1)[j]-(*FromDosage2)[j])*((*FromDosage1)[j]-(*FromDosage2)[j]));
//            sumNum+=(((*FromDosage1)[j]-(*FromDosage2)[j])*((*FromGT)[j]-(*FromDosage2)[j]));
//
//        }
//
//        double weight=sumNum/sumDen;
//
//        cout<<weight<<endl;
//
//        for(int j=0;j<(int)StudyBothVariants.size();j++)
//        {
//            int NewStudyMarkerPos=StudyBothVariants[j];
//
//            variant *thisVariant=&OutputData.VariantList[NewStudyMarkerPos];
//            int OrigStudyMarkerPos0=thisVariant->StudyMap[0][1];
//            int OrigStudyMarkerPos1=thisVariant->StudyMap[1][1];
//
//            (*CopyToDosage)[NewStudyMarkerPos]= (weight*(*CopyFromDosage1)[OrigStudyMarkerPos0])
//                                                    + ((1-weight)*(*CopyFromDosage2)[OrigStudyMarkerPos1]);
//
//        }
//
//    }
//


}





void MetaMinimac::PrintLSQEstimates(int ID)
{
    IFILE temp=NULL;

    if(ID==-1)
    {
        temp = ifopen(outFile + ".LeastSquare" + (gzip ? ".gz" : ""), "wb", gzip ?InputFile::BGZF:InputFile::UNCOMPRESSED);
        ifprintf(temp,"Chunk#\tParam");
        for(int i=0;i<numHaplotypes;i++)
        {
            ifprintf(temp,"\t%s_HAP_%d",individualName[i/2].c_str(),i%2+1);
        }
        ifprintf(temp,"\n");

    }
    else if(ID==-2)
    {
        temp = ifopen(outFile + ".LeastSquare" + (gzip ? ".gz" : ""), "a", gzip ?InputFile::BGZF:InputFile::UNCOMPRESSED);
        ifprintf(temp,"Total\tSumSq",ID);
        for(int i=0;i<numHaplotypes;i++)
        {
              ifprintf(temp,"\t\t%.3f",ErrorPerSample[i]);
        }
        ifprintf(temp,"\n");

    }

    else
    {
        temp=ifopen(outFile + ".LeastSquare" + (gzip ? ".gz" : ""), "a", gzip ?InputFile::BGZF:InputFile::UNCOMPRESSED);
    }


    if(ID>=0)
    {
        ifprintf(temp,"%d\tSumSq",ID);
        for(int i=0;i<numHaplotypes;i++)
        {
            ifprintf(temp,"\t\t%.3f",ErrorPerSamplePerChunk[i]);
        }
        ifprintf(temp,"\n");
        for(int kk=0;kk<(int)LSQEstimates[0].size();kk++)
        {
            ifprintf(temp,"%d\t%d",ID,kk);

            for(int i=0;i<numHaplotypes;i++)
            {
                ifprintf(temp,"\t\t%.3f",LSQEstimates[i][kk]);
            }
            ifprintf(temp,"\n");
        }

    }
    ifclose(temp);

}





bool MetaMinimac::AnalyzeChunks()
{


    list<VcfFileReader> myVcfFileStream;
    list<VcfRecord> MyRecord;

    for(int i=0;i<NoInPrefix;i++)
    {
        myVcfFileStream.emplace_back();
        MyRecord.emplace_back();
        if(!InputData[i].InitializeDosageChunkData(myVcfFileStream.back(), MyRecord.back()))
            return false;
    }


    PrintLSQEstimates(-1);
    PrintVCFHeader();


    for(int i=0;i<NoChunks;i++)
    {
        ThisChunk &tempChunk=ChunkList[i];

        cout<<"\n -------------"<<endl;
        cout<<"  CHUNK NO "<<i+1<<endl;
        cout<<" -------------"<<endl<<endl;
        list<VcfFileReader>::iterator itStream;
        list<VcfRecord>::iterator itRecord=MyRecord.begin();

        int PrefixNo=0;

        for (itStream=myVcfFileStream.begin(); itStream!=myVcfFileStream.end(); ++itStream,++itRecord)
        {
            if(!InputData[PrefixNo++].LoadDosageChunkData(tempChunk,*itStream,*itRecord))
                return false;

        }


        AnalyzeLooDosage(tempChunk);

        PrintLSQEstimates(i);

        CreateLooFitDosage(tempChunk);
        PrintVCFChunk(tempChunk);


    }


    PrintLSQEstimates(-2);
    return true;
}



void TransForm(double &Value)
{


    if(Value>0.7)
        Value=0.98;
    else if (Value <0.3)
        Value = 0.02;
    else
        Value = 0.5;


    Value=log(Value/(1.0-Value));


}

void MetaMinimac::UpdateFlankFractions()
{

    cout<<numHaplotypes<<"\t"<<OutputData.numMarkers<<endl;
    for(int Ind=0;Ind<numHaplotypes;Ind++)
    {
        for(int i=0;i<OutputData.numMarkers;i++)
        {
            variant *thisVariant=&OutputData.VariantList[i];


//            cout<<" I = "<<i<<endl;
            if(thisVariant->NoStudies>1)
            {
                int Sum=0;


                for(int StudyIndex=0;StudyIndex<thisVariant->NoStudies;StudyIndex++)
                {
                    int StudyID=thisVariant->StudyID[StudyIndex];
                    int MarkerID=thisVariant->StudyIDIndex[StudyIndex];


                    Sum += (1+InputData[StudyID].FlankLength[Ind][MarkerID]);
                }

//                cout<<Sum<<endl;
                for(int StudyIndex=0;StudyIndex<thisVariant->NoStudies;StudyIndex++)
                {
                    int StudyID=thisVariant->StudyID[StudyIndex];
                    int MarkerID=thisVariant->StudyIDIndex[StudyIndex];
                    InputData[StudyID].FlankFrac[Ind][MarkerID]=(double)(1+InputData[StudyID].FlankLength[Ind][MarkerID])/(double)Sum;



                    TransForm(InputData[StudyID].FlankFrac[Ind][MarkerID]);







//                    if(InputData[StudyID].FlankFrac[Ind][MarkerID]<=0.0 || InputData[StudyID].FlankFrac[Ind][MarkerID] >=1.0)
//                        cout<<" WELLLLL = "<<InputData[StudyID].FlankFrac[Ind][MarkerID]<<endl;


//                    InputData[StudyID].FlankFrac[Ind][MarkerID]=
//                        (double)(1+InputData[StudyID].FlankLength[Ind][MarkerID])
//                        /(double)(1+InputData[StudyID].FlankLength[Ind][MarkerID]);
                }
            }

        }
    }
}




void MetaMinimac::CopyParameters (HaplotypeSet &HapData,String InPrefix)
{
//    HapData.BufferSize=BufferSize;
    HapData.outFile=outFile;
    HapData.gzip=gzip;
    HapData.finChromosome="";
    HapData.Method=Method;
    HapData.InfilePrefix=InPrefix;
    HapData.SplitLimit=SplitLimit;
}



void MetaMinimac::CopyParameters (HaplotypeSet &FromHapData,HaplotypeSet &ToHapData)
{
//    ToHapData.BufferSize=BufferSize;
    ToHapData.outFile=outFile;
    ToHapData.gzip=gzip;
    ToHapData.finChromosome="";
    ToHapData.individualName=FromHapData.individualName;


}



bool MetaMinimac::ImputeUniqueVariants()
{

    OutputData.HapDosage.resize(OutputData.numHaplotypes);
    for(int i=0;i<(OutputData.numHaplotypes);i++)
        OutputData.HapDosage[i].resize(OutputData.numMarkers);

    for(int i=0;i<(OutputData.numHaplotypes);i++)
    {
        vector<double> *CopyFromDosage1=&(InputData[0].HapDosage[i]);
        vector<double> *CopyFromDosage2=&(InputData[1].HapDosage[i]);
        vector<double> *CopyToDosage=&(OutputData.HapDosage[i]);


        for(int j=0;j<(int)StudyOneVariants.size();j++)
        {
            int NewStudyMarkerPos=StudyOneVariants[j];

            variant *thisVariant=&OutputData.VariantList[NewStudyMarkerPos];

            int OrigStudyMarkerPos=thisVariant->StudyMap[0][1];

            (*CopyToDosage)[NewStudyMarkerPos]=(*CopyFromDosage1)[OrigStudyMarkerPos];
        }


        for(int j=0;j<(int)StudyTwoVariants.size();j++)
        {
            int NewStudyMarkerPos=StudyTwoVariants[j];

            variant *thisVariant=&OutputData.VariantList[NewStudyMarkerPos];

            int OrigStudyMarkerPos=thisVariant->StudyMap[0][1];

            (*CopyToDosage)[NewStudyMarkerPos]=(*CopyFromDosage2)[OrigStudyMarkerPos];
        }
    }
}



bool MetaMinimac::CreateMetaRecomDosage()
{


    int flankLength=100;

    cout<<" Method C : Use Sum of switch probabilities in a Flanking region of "<<flankLength <<" markers ..."<<endl<<endl;



    for(int i=0;i<(OutputData.numHaplotypes);i++)
    {

//        cout<<InputData[0].RecomSpec[2][2000]<<endl;

        vector<double> *FromDosage1=&(InputData[0].LooDosage[i]);
        vector<double> *FromDosage2=&(InputData[1].LooDosage[i]);
        vector<double> *FromGT=&(InputData[1].GTDosage[i]);


        double sumNum=0.0,sumDen=0.0;

        for(int j=0;j<InputData[0].NoTypedSites;j++)
        {
            sumDen+=(((*FromDosage1)[j]-(*FromDosage2)[j])*((*FromDosage1)[j]-(*FromDosage2)[j]));
            sumNum+=(((*FromDosage1)[j]-(*FromDosage2)[j])*((*FromGT)[j]-(*FromDosage2)[j]));

        }

        double weight1=sumNum/sumDen;


        vector<double> *FromRecom1=&(InputData[0].RecomSpec[i]);
        vector<double> *FromRecom2=&(InputData[1].RecomSpec[i]);

        vector<double> Constants1((int)StudyBothVariants.size(),0.0);
        vector<double> Constants2((int)StudyBothVariants.size(),0.0);


        int flankLeft1=0,flankLeft2=0,flankRight1=0,flankRight2=0;


        for(int j=0;j<(int)StudyBothVariants.size();j++)
        {

            int NewStudyMarkerPos=StudyBothVariants[j];
            variant *thisVariant=&OutputData.VariantList[NewStudyMarkerPos];
            int OrigStudyMarkerPos0=thisVariant->StudyMap[0][1];
            int OrigStudyMarkerPos1=thisVariant->StudyMap[1][1];



            if((*FromRecom1)[OrigStudyMarkerPos0]<0.05)
                flankLeft1++;
            else
                flankLeft1=0;


            int kk=OrigStudyMarkerPos0;
            if(flankRight1==0)
            {

                while((*FromRecom1)[kk++]<0.05 && kk<InputData[0].numMarkers)
                    flankRight1++;

            }
            else
                flankRight1--;

            int flank1=(flankLeft1<(flankRight1)? flankLeft1:(flankRight1));





            if((*FromRecom2)[OrigStudyMarkerPos1]<0.05)
                flankLeft2++;
            else
                flankLeft2=0;


            kk=OrigStudyMarkerPos1;
            if(flankRight2==0)
            {

                while((*FromRecom2)[kk++]<0.05 && kk<InputData[1].numMarkers)
                    flankRight2++;

            }
            else
                flankRight2--;

            int flank2=(flankLeft2<(flankRight2)? flankLeft2:(flankRight2));



//            int flank1=0;
//
//            int kk=OrigStudyMarkerPos0;
//            while(kk<InputData[0].numMarkers && (*FromRecom1)[kk]<0.1)
//            {
//                kk++;
//            }
//            flank1=kk-OrigStudyMarkerPos0;
//            kk=OrigStudyMarkerPos0;
//            while(kk>=0&& (*FromRecom1)[kk]<0.05)
//            {
//                kk--;
//            }
//            flank1=(flank1<(OrigStudyMarkerPos0-kk)? flank1:(OrigStudyMarkerPos0-kk));
//
//            int flank2=0;
//
//            kk=OrigStudyMarkerPos1;
//            while(kk<InputData[1].numMarkers && (*FromRecom2)[kk]<0.1)
//            {
//                kk++;
//            }
//            flank2=kk-OrigStudyMarkerPos1;
//            kk=OrigStudyMarkerPos1;
//            while(kk>=0&& (*FromRecom2)[kk]<0.05)
//            {
//                kk--;
//            }
//            flank2=(flank2<(OrigStudyMarkerPos1-kk)? flank2:(OrigStudyMarkerPos1-kk));

            Constants1[j]=flank1;
            Constants2[j]=flank2;

//            cout<<NewStudyMarkerPos<<" "<<flank1<<"\t"<<flank2<<endl;

//            for(int kk=OrigStudyMarkerPos0-flankLength;
//                    kk<=OrigStudyMarkerPos0+flankLength;
//                    kk++)
//            {
//                if(kk>=0 && kk<InputData[0].numMarkers)
//                {
//                    Constants1[j]+=(*FromRecom1)[kk];
//                }
//
//            }
//
//
//
//            for(int kk=OrigStudyMarkerPos1-flankLength;
//
//
//                    kk<=OrigStudyMarkerPos1+flankLength;
//                    kk++)
//            {
//                if(kk>=0 && kk<InputData[1].numMarkers)
//                {
//                    Constants2[j]+=(*FromRecom2)[kk];
//
////                if(i>1)
////                    cout<<(*FromRecom2)[kk]<<endl;
//                }
//            }
//

            vector<double> *CopyToDosage=&(OutputData.HapDosage[i]);
            vector<double> *CopyFromDosage1=&(InputData[0].HapDosage[i]);
            vector<double> *CopyFromDosage2=&(InputData[1].HapDosage[i]);


             double weight;
//            weight1=0;

            if(Constants1[j]==Constants2[j])
                weight=0.5*(weight1+0.5);
            else
                weight=0.5*(weight1+Constants1[j]/(Constants1[j]+Constants2[j]));

//            cout<<weight1<<"\t"<<weight<<endl;

            (*CopyToDosage)[NewStudyMarkerPos]= (weight*(*CopyFromDosage1)[OrigStudyMarkerPos0])
                                                    + ((1-weight)*(*CopyFromDosage2)[OrigStudyMarkerPos1]);

        }

    }



}






bool MetaMinimac::CreateInfoR2Dosage()
{

    cout<<" Method A : Choose Panel with higher Estimated R-Sqaure ..."<<endl<<endl;


    for(int i=0;i<(OutputData.numHaplotypes);i++)
    {

        vector<double> *CopyToDosage=&(OutputData.HapDosage[i]);
        vector<double> *CopyFromDosage2=&(InputData[0].HapDosage[i]);
        vector<double> *CopyFromDosage3=&(InputData[1].HapDosage[i]);

         for(int j=0;j<(int)StudyBothVariants.size();j++)
        {
            int NewStudyMarkerPos=StudyBothVariants[j];

            variant *thisVariant=&OutputData.VariantList[NewStudyMarkerPos];
            int OrigStudyMarkerPos0=thisVariant->StudyMap[0][1];
            int OrigStudyMarkerPos1=thisVariant->StudyMap[1][1];

            if(InputData[0].EstRsq[OrigStudyMarkerPos0] > InputData[1].EstRsq[OrigStudyMarkerPos1])
                (*CopyToDosage)[NewStudyMarkerPos]=(*CopyFromDosage2)[OrigStudyMarkerPos0];
            else
                (*CopyToDosage)[NewStudyMarkerPos]=(*CopyFromDosage3)[OrigStudyMarkerPos1];
        }

    }

}




bool MetaMinimac::MergeDosageData()
{

    cout<<"\n ------------------------------------------------------------------------------"<<endl;
    cout<<"                               MERGE DOSAGE DATA                               "<<endl;
    cout<<" ------------------------------------------------------------------------------"<<endl;



    CopyParameters(InputData[0],OutputData);
    MergeVariantList();

//  cout<<" WOW ="<<InputData[0].RecomSpec[2][2000]<<endl;


    cout<<" Copying Unique Variant from Dosage Data ..."<<endl<<endl;

    ImputeUniqueVariants();

    cout<<" Mering Dosage Data ..."<<endl<<endl;
  if(Method=="A")
    {
        CreateInfoR2Dosage();
    }
    else if(Method=="B")
    {
//        CreateLooFitDosage();
    }
    else if(Method=="C")
    {
        CreateMetaRecomDosage();
    }
    else
    {
        cout<<"\n Improper Value of Parameter Method = "<<Method<<endl;
        cout<<" Please change the value of Method !!!"<<endl;
        return false;

    }

    OutputData.PrintOutputVCF();

    return true;


//
//
//
//    InPrefixList.clear();
//
//    size_t pos = 0;
//    std::string delimiter(":") ;
//    std::string token;
//    int Count=0;
//    string tempName=inputFile.c_str();
//    string temptempName=tempName;
//
//    while ((pos = tempName.find(delimiter)) != std::string::npos)
//    {
//        token = tempName.substr(0, pos);
//        InPrefixList.push_back(token);
////
////        if(Count==0)
////            familyName.push_back(token);
////
////        tempName.erase(0, pos + delimiter.length());
////        if(Count==1)
////        {
////            std::cout << "\n ERROR !!! "<<endl;
////            cout << "\n Program could NOT parse the following sample name with ID delimiter ["<< IdDelimiter<<"] : " << temptempName << endl;
////            cout<<    " More than TWO tokens found : ["<<familyName.back()  <<"] ["
////            <<token << "] ["<<tempName  <<"] " <<endl;
////            cout << "\n Please verify ID Delimitier : [" << IdDelimiter <<"]"<< endl;
////            return false;
////        }
//        tempName.erase(0, pos + delimiter.length());
//        Count++;
//
//
//        }
//
//
//    InPrefixList.push_back(tempName);
//
//
//    NoInPrefix=(int)InPrefixList.size();
//    InputData.clear();
//    InputData.resize(NoInPrefix);
//
//    for(int i=0;i<NoInPrefix;i++)
//    {
//        cout<<"\n ----------------------------"<<endl;
//        cout<<"       INPUT FILE : ["<<i+1<<"] "<<endl;
//        cout<<" ----------------------------"<<endl;
//
//        cout<<" File Prefix        : "<<InPrefixList[i] <<endl;
//
//        InputData[i].BufferSize=BufferSize;
//        InputData[i].outFile=outFile;
//        InputData[i].gzip=gzip;
//
//        InputData[i].ReadInputFiles(InPrefixList[i]);
//
//    }
//
//    InputData[0].PrintOutputVCF();
//

}






//
//    for (std::map<int,vector<int> >::iterator it=HashUnionVarMap.begin(); it!=HashUnionVarMap.end(); ++it)
//    {
//        vector<int> Mapper= it->second;
//        int ThisSize=(int)Mapper.size();
//
//        std::map<string,int> RefAltAlleleMap;
//
//        int Counter=0;
//
//        while(Counter<(int)Mapper.size())
//        {
//            if(Mapper[Counter]==0)
//            {
//
//                variant *thisVariant=&InputData[0].VariantList[Mapper[++Counter]];
//                string tempString=thisVariant->refAlleleString+"_"+thisVariant->altAlleleString;
//               if(RefAltAlleleMap[tempString]>0)
//                {
//                    cout<<"\n Duplicate Variant found in Study 0 "<<thisVariant->bp
//                    <<"_"<<thisVariant->refAlleleString+"_"+thisVariant->altAlleleString<<endl;
//                    cout<<" Please check input file prefix properly ... "<<endl;
//                    return false;
//                }
//                OutputData.VariantList.push_back(*thisVariant);
//                variant *thisVariant2=&OutputData.VariantList.back();
//
//                thisVariant2->NoStudies=1;
//
//                vector<int> tempVector;
//                tempVector.clear();
//                tempVector.push_back(0);
//                tempVector.push_back(Mapper[Counter]);
//                thisVariant2->StudyMap.push_back(tempVector);
//
//                TotalUnionVariants++;
//
//                RefAltAlleleMap[tempString]=TotalUnionVariants;
//                Counter++;
//
//
//            }
//            else
//            {
//                variant *thisVariant=&InputData[1].VariantList[Mapper[++Counter]];
//                string tempString=thisVariant->refAlleleString+"_"+thisVariant->altAlleleString;
//                if(RefAltAlleleMap[tempString]>0)
//                {
//                    variant *thisVariant2=&OutputData.VariantList[RefAltAlleleMap[tempString]-1];
//
//                    thisVariant2->NoStudies++;
//
//                    vector<int> tempVector;
//                    tempVector.clear();
//                    tempVector.push_back(1);
//                    tempVector.push_back(Mapper[Counter]);
//                    thisVariant2->StudyMap.push_back(tempVector);
//
//                }
//                else
//                {
//                    OutputData.VariantList.push_back(*thisVariant);
//
//                    variant *thisVariant2=&OutputData.VariantList.back();
//
//                    thisVariant2->NoStudies=1;
//
//                    vector<int> tempVector;
//                    tempVector.clear();
//                    tempVector.push_back(1);
//                    tempVector.push_back(Mapper[Counter]);
//                    thisVariant2->StudyMap.push_back(tempVector);
//
//                    TotalUnionVariants++;
//                }
//                Counter++;
//            }
//        }
//    }

//
//    OutputData.numMarkers=(int)OutputData.VariantList.size();
//    OutputData.numHaplotypes=InputData[0].numHaplotypes;
//
//    StudyOneVariants.clear();
//    StudyTwoVariants.clear();
//    StudyBothVariants.clear();
//
//
//    for(int i=0;i<OutputData.numMarkers;i++)
//    {
//        variant *thisVariant=&OutputData.VariantList[i];
//
//        if((int)thisVariant->StudyMap.size()==2)
//            StudyBothVariants.push_back(i);
//        else
//        {
//
//            if((int)thisVariant->StudyMap[0][0]==0)
//                StudyOneVariants.push_back(i);
//            else if((int)thisVariant->StudyMap[0][0]==1)
//                StudyTwoVariants.push_back(i);
//            else
//            {
//                 cout<<" ERROR IN MERGING VARIANTS IN METAMINIMAC.cpp. Contact Author"<<endl;
//            }
//
//
//
//        }
//
//    }
//
//
//
//	std::cout << "\n Final Number of Samples   : " << OutputData.numHaplotypes/2 << endl;
//    std::cout << " Final Number of Variants  : " << OutputData.numMarkers << endl;
//    std::cout << " Variants only in Study 1  : " << (int)StudyOneVariants.size() << endl;
//    std::cout << " Variants only in Study 1  : " << (int)StudyTwoVariants.size()  << endl;
//    std::cout << " Variants in both Studies  : " << (int)StudyBothVariants.size()  << endl;
//

