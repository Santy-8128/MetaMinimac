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

    OpenStreamInputDosageFiles();
    if (!OpenStreamOutputDosageFiles())
    {
        cout <<" Please check your write permissions in the output directory\n OR maybe the output directory does NOT exist ...\n";
        cout << "\n Program Exiting ... \n\n";
        return "File.Write.Error";
    }

    PrintChunkInformation();

    return PerformFinalAnalysis();
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
        InPrefixList.push_back(token.c_str());
        tempName.erase(0, pos + delimiter.length());
        Count++;
    }
    InPrefixList.push_back(tempName.c_str());


    NoInPrefix=(int)InPrefixList.size();
    InputData.clear();
    InputData.resize(NoInPrefix);

    cout<<endl<<  " Number of Studies                  : "<<NoInPrefix<<endl;
    for(int i=0;i<NoInPrefix;i++)
    {
        cout<<  " -- Study "<<i+1<<" Prefix                  : "<<InPrefixList[i]<<endl;
    }


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

    cout<<" -- Successful !!! "<<endl;
    return true;
}

bool MetaMinimac::ReadEmpVariantAndChunk()
{
    cout<<"\n Gathering information for chunking ... "<<endl;
    for(int i=0;i<NoInPrefix;i++)
    {
        InputData[i].LoadEmpVariantList();
        cout<<" -- Study "<<i+1<<" #Genotyped Sites = "<<InputData[i].noTypedMarkers<<endl;
    }
    FindCommonGenotypedVariants();
    cout<<" Summary: Found "<< InputData[0].numSamples<<" samples (" << InputData[0].numActualHaps << " haplotypes) across "<<NoInPrefix<<" studies on "<< NoCommonGenoVariants<<" common genotyped sites !" <<endl;

    GetNumChunks();
    CreateChunks();

    return true;
}

void MetaMinimac::FindCommonGenotypedVariants()
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
    CurrentMetaImputedDosage.resize(InputData[0].numHaplotypes);
    CurrentMetaImputedDosageSum.resize(InputData[0].numHaplotypes);


    InputData[0].SortCommonGenotypeList(CommonGenotypeVariantNameList, SortedCommonGenoList, CommonTypedVariantList);

    for(int i=1; i<NoInPrefix; i++)
        InputData[i].ReadBasedOnSortCommonGenotypeList(SortedCommonGenoList);

    cout<<" -- Common  #Genotyped Sites = "<<NoCommonGenoVariants<<endl;
    cout<<" -- Successful !!! "<<endl<<endl;
}

void MetaMinimac::GetNumChunks()
{

    int NoMarkers = CommonTypedVariantList.size();
    int StartPos = CommonTypedVariantList[0].bp;
    int EndPos = CommonTypedVariantList[NoMarkers - 1].bp;

    NoChunks = ((int) (EndPos - StartPos) / (int) (1000000 * ThisVariable.chunkLength));

    if (NoChunks == 0)
        NoChunks = 1;

}

bool MetaMinimac::CreateChunks()
{

    int NoMarkers=CommonTypedVariantList.size();

    ThisVariable.chunkLength=(NoMarkers/NoChunks);
    ChunkList.resize(NoChunks);
    ChunkList[0].StartBp=0;
    ChunkList[0].StartWithWindowIndex=0;

    int tempNoChunkMarkers = (NoMarkers/NoChunks)+1;
    int tempNoWindowMarkers = tempNoChunkMarkers/10;


    int counter=0, chunkCounter = 0;

    while(counter<NoCommonGenoVariants)
    {
        ThisChunk &tempChunk = ChunkList[chunkCounter];

        tempChunk.EndWithWindowIndex=counter+tempNoChunkMarkers+tempNoWindowMarkers;
        tempChunk.ChunkNo=chunkCounter;


        if(chunkCounter>0)
        {
            tempChunk.StartWithWindowIndex=counter - tempNoWindowMarkers;
        }
        tempChunk.NoGenoAllStudies=tempChunk.EndWithWindowIndex-tempChunk.StartWithWindowIndex+1;

        if(chunkCounter<NoChunks-1)
        {
            tempChunk.EndBp=CommonTypedVariantList[counter+tempNoChunkMarkers].bp;
            ChunkList[chunkCounter+1].StartBp=tempChunk.EndBp+1;
        }
        else
        {
            tempChunk.EndBp=999999997;
            tempChunk.EndWithWindowIndex=NoCommonGenoVariants-1;
            tempChunk.NoGenoAllStudies=tempChunk.EndWithWindowIndex-tempChunk.StartWithWindowIndex+1;
        }
        tempChunk.NoVariants=tempChunk.EndWithWindowIndex-tempChunk.StartWithWindowIndex+1;
        counter+=tempNoChunkMarkers;


        chunkCounter++;

    }

    return true;
}

void MetaMinimac::PrintChunkInformation()
{

    cout<<"\n ------------------------------------------------------------------------------"<<endl;
    cout<<"                               CHUKNING SUMMARY                               "<<endl;
    cout<<" ------------------------------------------------------------------------------"<<endl;
    cout<<" No   LeftBuffer      LeftEnd    RightPoint   RightBuffer    #TypedSites"<<endl;
    cout<<" -------------------------------------------------------------------------------"<<endl;

    for(int i=0;i<NoChunks;i++)
    {
        cout<<setw(3)<<i+1<<"  "<<setw(11);
        i==0? cout<<" START":cout<<CommonTypedVariantList[ChunkList[i].StartWithWindowIndex].bp ;
        cout<<"  "<<setw(11);
        i==0? cout<<" START": cout<<ChunkList[i].StartBp;
        cout<<"  "<<setw(12);
        i==NoChunks-1? cout<<" END":cout<<ChunkList[i].EndBp;
        cout<<"  "<<setw(12);
        i==NoChunks-1? cout<<" END":cout<<CommonTypedVariantList[ChunkList[i].EndWithWindowIndex].bp;
        cout<<"  "<<setw(13)<<ChunkList[i].NoVariants<<endl;



    }
    cout<<endl<<endl;
}


void MetaMinimac::Initialize()
{
    LSQEstimates.resize(InputData[0].numHaplotypes);

    for (int j = 0; j < InputData[0].numHaplotypes; j++) {
        LSQEstimates[j].resize(NoInPrefix);
    }

  if(ThisVariable.debug)
  {
      ErrorPerSamplePerChunk.resize(InputData[0].numHaplotypes);
      ErrorPerSample.resize(InputData[0].numHaplotypes, 0.0);
  }

}

String MetaMinimac::PerformFinalAnalysis()
{


    cout << " ------------------------------------------------------------------------------" << endl;
    cout << "                           META-IMPUTATION ANALYSIS                            " << endl;
    cout << " ------------------------------------------------------------------------------" << endl;

    Initialize();

    for (int i = 0; i < NoChunks; i++) {
        cout << "\n -------------------------------------------" << endl;
        cout << " Analyzing Chunk " << i + 1 << "/" << NoChunks << " [" << CommonTypedVariantList[0].chr << ":";
        cout << CommonTypedVariantList[ChunkList[i].StartWithWindowIndex].bp << "-";
        cout << CommonTypedVariantList[ChunkList[i].EndWithWindowIndex].bp << "]" << endl;
        cout << " -------------------------------------------" << endl;

        for (int j = 0; j < InputData[0].numSamples; j++)
        {
            if (InputData[0].SampleNoHaplotypes[j] == 2)
            {
                GetMetaImpEstimates(2 * j, ChunkList[i]);
                GetMetaImpEstimates(2 * j + 1, ChunkList[i]);
            }
            else
                GetMetaImpEstimates(2 * j, ChunkList[i]);
        }

       AppendtoMainVcfFaster(i);

       if(ThisVariable.debug)
           AppendtoMainWeightsFile(i);

    }

    for (int i = 0; i < NoInPrefix; i++) {
        delete InputDosageStream[i];
        delete CurrentRecordFromStudy[i];
    }

    return "Success";

}

void MetaMinimac::GetMetaImpEstimates(int Sample, ThisChunk &MyChunk)
{

    LogOddsModel ThisSampleAnalysis;

    ThisSampleAnalysis.metaInitialize(Sample, InputData, this, MyChunk);

    vector<double> init(NoInPrefix-1, 0.0);
    vector<double> MiniMizer = Simplex(ThisSampleAnalysis, init);
    logitTransform(MiniMizer, LSQEstimates[Sample]);

    if(ThisVariable.debug)
    {
        ErrorPerSamplePerChunk[Sample]=ThisSampleAnalysis(LSQEstimates[Sample]);
        ErrorPerSample[Sample]+=ErrorPerSamplePerChunk[Sample];
        ErrorSumSqPerChunk+=ErrorPerSamplePerChunk[Sample];
    }

    if((Sample+1)%200==0 || Sample==InputData[0].numHaplotypes-1)
        cout<<" Finished Analyzing "<< Sample/2+1 <<" samples ..."<<endl;
}

void MetaMinimac::OpenStreamInputDosageFiles()
{
    InputDosageStream.resize(NoInPrefix);
    CurrentRecordFromStudy.resize(NoInPrefix);
    StudiesHasVariant.resize(NoInPrefix);
    for(int i=0; i<NoInPrefix;i++)
    {
        VcfHeader header;
        InputDosageStream[i] = new VcfFileReader();
        CurrentRecordFromStudy[i]= new VcfRecord();
        InputDosageStream[i]->open( (GetDosageFileFullName(InPrefixList[i])).c_str() , header);

        InputDosageStream[i]->readRecord(*CurrentRecordFromStudy[i]);
        //        cout<<CurrentRecordFromStudy[i]->get1BasedPosition()<<endl;
    }
}

bool MetaMinimac::OpenStreamOutputDosageFiles()
{
    CommonTypedVariantListCounter=0;
    bool gzip=ThisVariable.gzip;
    VcfPrintStringPointer = (char*)malloc(sizeof(char) * (ThisVariable.PrintBuffer));

    vcfdosepartial = ifopen(ThisVariable.outfile + ".metaDose.vcf" + (gzip ? ".gz" : ""), "wb", gzip ?InputFile::BGZF:InputFile::UNCOMPRESSED);
    if(vcfdosepartial==NULL)
    {
        cout <<"\n\n ERROR !!! \n Could NOT create the following file : "<< ThisVariable.outfile + ".dose.vcf" + (gzip ? ".gz" : "") <<endl;
        return false;
    }

    ifprintf(vcfdosepartial,"##fileformat=VCFv4.1\n");
    time_t t = time(0);
    struct tm * now = localtime( & t );
    ifprintf(vcfdosepartial,"##filedate=%d.%d.%d\n",(now->tm_year + 1900),(now->tm_mon + 1) ,now->tm_mday);
    ifprintf(vcfdosepartial,"##source=MetaMinimac.v%s\n",VERSION);
    ifprintf(vcfdosepartial,"##contig=<ID=%s>\n",InputData[0].finChromosome.c_str());
    ifprintf(vcfdosepartial,"##INFO=<ID=AF,Number=1,Type=Float,Description=\"Estimated Alternate Allele Frequency\">\n");
    ifprintf(vcfdosepartial,"##INFO=<ID=MAF,Number=1,Type=Float,Description=\"Estimated Minor Allele Frequency\">\n");
    //ifprintf(vcfdosepartial,"##INFO=<ID=R2,Number=1,Type=Float,Description=\"Estimated Imputation Accuracy (R-square)\">\n");
    //ifprintf(vcfdosepartial,"##INFO=<ID=ER2,Number=1,Type=Float,Description=\"Empirical (Leave-One-Out) R-square (available only for genotyped variants)\">\n");
    //ifprintf(vcfdosepartial,"##INFO=<ID=IMPUTED,Number=0,Type=Flag,Description=\"Marker was imputed but NOT genotyped\">\n");
    //ifprintf(vcfdosepartial,"##INFO=<ID=TYPED,Number=0,Type=Flag,Description=\"Marker was genotyped AND imputed\">\n");
    ifprintf(vcfdosepartial,"##INFO=<ID=TRAINING,Number=0,Type=Flag,Description=\"Marker was used to train meta-imputation weights\">\n");

    if(ThisVariable.infoDetails)
    {
        ifprintf(vcfdosepartial,"##INFO=<ID=NST,Number=1,Type=Integer,Description=\"Number of studies marker was found during meta-imputation\">\n");
        for(int i=0; i<NoInPrefix; i++)
            ifprintf(vcfdosepartial,"##INFO=<ID=S%d,Number=0,Type=Flag,Description=\"Marker was present in Study %d\">\n",i+1,i+1);
    }

    if(ThisVariable.GT)
        ifprintf(vcfdosepartial,"##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n");
    if(ThisVariable.DS)
        ifprintf(vcfdosepartial,"##FORMAT=<ID=DS,Number=1,Type=Float,Description=\"Estimated Alternate Allele Dosage : [P(0/1)+2*P(1/1)]\">\n");
    if(ThisVariable.HDS)
        ifprintf(vcfdosepartial,"##FORMAT=<ID=HDS,Number=2,Type=Float,Description=\"Estimated Haploid Alternate Allele Dosage \">\n");
    if(ThisVariable.GP)
        ifprintf(vcfdosepartial,"##FORMAT=<ID=GP,Number=3,Type=Float,Description=\"Estimated Posterior Probabilities for Genotypes 0/0, 0/1 and 1/1 \">\n");
    if(ThisVariable.SD)
        ifprintf(vcfdosepartial,"##FORMAT=<ID=SD,Number=1,Type=Float,Description=\"Variance of Posterior Genotype Probabilities\">\n");

    ifprintf(vcfdosepartial,"##metaMinimac_Command=%s\n",ThisVariable.CommandLine.c_str());

    ifprintf(vcfdosepartial,"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT");

    for(int Id=0;Id<InputData[0].numSamples;Id++)
    {
        ifprintf(vcfdosepartial,"\t%s",InputData[0].individualName[Id].c_str());
    }
    ifprintf(vcfdosepartial,"\n");
    ifclose(vcfdosepartial);


    if(ThisVariable.debug)
    {
        metaWeight = ifopen(ThisVariable.outfile + ".metaWeights.vcf" + (gzip ? ".gz" : ""), "wb", gzip ?InputFile::BGZF:InputFile::UNCOMPRESSED);
        if(metaWeight==NULL)
        {
            cout <<"\n\n ERROR !!! \n Could NOT create the following file : "<< ThisVariable.outfile + ".metaWeights.vcf" + (gzip ? ".gz" : "") <<endl;
            return false;
        }
        ifprintf(metaWeight,"##fileformat=NA\n");
        time_t t = time(0);
        struct tm * now = localtime( & t );
        ifprintf(metaWeight,"##filedate=%d.%d.%d\n",(now->tm_year + 1900),(now->tm_mon + 1) ,now->tm_mday);
        ifprintf(metaWeight,"##source=MetaMinimac.v%s\n",VERSION);
        ifprintf(metaWeight,"##contig=<ID=%s>\n",InputData[0].finChromosome.c_str());
        ifprintf(metaWeight,"##INFO=<ID=AF,Number=1,Type=Float,Description=\"Estimated Alternate Allele Frequency\">\n");
        ifprintf(metaWeight,"##INFO=<ID=MAF,Number=1,Type=Float,Description=\"Estimated Minor Allele Frequency\">\n");
        ifprintf(metaWeight,"##metaMinimac_Command=%s\n",ThisVariable.CommandLine.c_str());

        ifprintf(metaWeight,"#CHUNK");

        for(int Id=0;Id<InputData[0].numSamples;Id++)
        {
            ifprintf(metaWeight,"\t%s",InputData[0].individualName[Id].c_str());
        }
        ifprintf(metaWeight,"\n");
        ifclose(metaWeight);

    }


    return true;
}

bool MetaMinimac::doesExistFile(String filename)
{
    IFILE ifs = ifopen(filename.c_str(), "r");
    if (ifs)
    {
        ifclose(ifs);
        return true;
    }
    else
    {
        ifclose(ifs);
        return false;
    }
}

string MetaMinimac::GetDosageFileFullName(String prefix)
{
    String filename;
    if(doesExistFile(prefix+".dose.vcf"))
        filename = prefix+".dose.vcf";
    else if(doesExistFile(prefix+".dose.vcf.gz"))
        filename = prefix+".dose.vcf.gz";
    string FileFullName(filename.c_str());
    return FileFullName;
}

//void MetaMinimac::PrintVariant(VcfRecord *temp)
//{
//    cout<<temp->get1BasedPosition()<<":"<<temp->getRefStr()<<"_"<<temp->getAltStr()<<endl;
//}

void MetaMinimac::UpdateCurrentRecords()
{
   // cout<<" UPDATE "<<endl;
    for(int i=0; i<NoStudiesHasVariant;i++)
    {
        int index = StudiesHasVariant[i];
     //   PrintVariant(CurrentRecordFromStudy[index]);
        if(!InputDosageStream[index]->readRecord(*CurrentRecordFromStudy[index]))
            CurrentRecordFromStudy[index]->set1BasedPosition(999999999);
       // PrintVariant(CurrentRecordFromStudy[index]);

    }

    //abort();
}

int MetaMinimac::IsVariantEqual(VcfRecord &Rec1, VcfRecord &Rec2)
{
    if(strcmp(Rec1.getRefStr(),Rec2.getRefStr())!=0)
        return 0;
    if(strcmp(Rec1.getAltStr(),Rec2.getAltStr())!=0)
        return 0;
    return 1;
}

void MetaMinimac::FindCurrentMinimumPosition()
{

    if(NoInPrefix==2)
    {
        int a=CurrentRecordFromStudy[0]->get1BasedPosition();
        int b=CurrentRecordFromStudy[1]->get1BasedPosition();
        CurrentFirstVariantBp=a;
        NoStudiesHasVariant = 1;
        StudiesHasVariant[0]=0;

        if (b == a)
        {
            if(IsVariantEqual(*CurrentRecordFromStudy[0], *CurrentRecordFromStudy[1])==1)
            {
                NoStudiesHasVariant = 2;
                StudiesHasVariant[1] = 1;
            }
        }
        else if (b < a)
        {
            StudiesHasVariant[0]=1;
            CurrentFirstVariantBp=b;
        }

    }
    else if (NoInPrefix>2)
    {

        CurrentFirstVariantBp=CurrentRecordFromStudy[0]->get1BasedPosition();

        for(int i=1;i<NoInPrefix;i++)
            if(CurrentRecordFromStudy[i]->get1BasedPosition() < CurrentFirstVariantBp)
                CurrentFirstVariantBp=CurrentRecordFromStudy[i]->get1BasedPosition();

        NoStudiesHasVariant=0;
        VcfRecord *minRecord;
        minRecord = new VcfRecord();
        int Begin=0;
        for(int i=0;i<NoInPrefix;i++)
        {
            if(CurrentRecordFromStudy[i]->get1BasedPosition() == CurrentFirstVariantBp)
            {
                if(Begin==0)
                {
                    Begin=1;
                    minRecord=CurrentRecordFromStudy[i];
                    StudiesHasVariant[NoStudiesHasVariant] = i;
                    NoStudiesHasVariant++;
                }
                else if (IsVariantEqual(*minRecord, *CurrentRecordFromStudy[i])==1)
                {
                    StudiesHasVariant[NoStudiesHasVariant] = i;
                    NoStudiesHasVariant++;
                }
            }
        }
    }
}

string MetaMinimac::CreateInfo()
{
    stringstream ss;
    double freq=0.0;
    for(int i=0;i<InputData[0].numSamples;i++)
    {
        if(InputData[0].SampleNoHaplotypes[i]==2)
        {
            freq+=(CurrentMetaImputedDosage[2*i]+CurrentMetaImputedDosage[2*i+1]);
        }
    }
    freq/=InputData[0].numHaplotypes;
    double maf=freq > 0.5 ? 1.0 - freq : freq;
    ss<<"AF="<< fixed << setprecision(5) <<freq<<";MAF=";
    ss<<fixed << setprecision(5) << maf;
    ss<<";";

    if(!ThisVariable.infoDetails)
        return ".";
    ss<<"NST=";
    ss<<NoStudiesHasVariant;
    for(int i=0; i<NoStudiesHasVariant; i++)
    {
        ss<<";S";
        ss<<StudiesHasVariant[i]+1;
    }


    VcfRecord* tempRecord = CurrentRecordFromStudy[StudiesHasVariant[0]];
    if(tempRecord->getIDStr()==CommonTypedVariantList[CommonTypedVariantListCounter].name)
    {
        ss<<";TRAINING";
        CommonTypedVariantListCounter++;
    }
    return ss.str();
}

void MetaMinimac::PrintCurrentVariant()
{

    VcfRecord* tempRecord = CurrentRecordFromStudy[StudiesHasVariant[0]];
    VcfPrintStringPointerLength+=sprintf(VcfPrintStringPointer+VcfPrintStringPointerLength, "%s\t%d\t%s\t%s\t%s\t.\tPASS\t%s\t%s",
                                         tempRecord->getChromStr(), tempRecord->get1BasedPosition(),
                                         tempRecord->getIDStr(), tempRecord->getRefStr(),
                                         tempRecord->getAltStr(), CreateInfo().c_str(),
                                         ThisVariable.formatStringForVCF.c_str());


}

void MetaMinimac::PrintMetaImputedData()
{

    for(int Id=0;Id<InputData[0].numSamples;Id++)
    {
        int NoHaps = InputData[0].SampleNoHaplotypes[Id];
        if(NoHaps==2)
            PrintDiploidDosage((CurrentMetaImputedDosage[2*Id]), (CurrentMetaImputedDosage[2*Id+1]));
        else if(NoHaps==1)
            PrintHaploidDosage((CurrentMetaImputedDosage[2*Id]));
        else
            abort();

    }

    VcfPrintStringPointerLength+=sprintf(VcfPrintStringPointer+VcfPrintStringPointerLength,"\n");
    if(VcfPrintStringPointerLength > 0.9 * (float)(ThisVariable.PrintBuffer))
    {
        ifprintf(vcfdosepartial,"%s",VcfPrintStringPointer);
        VcfPrintStringPointerLength=0;
    }

}

void MetaMinimac::PrintDiploidDosage(float &x, float &y)
{

    bool colonIndex=false;
    VcfPrintStringPointerLength+=sprintf(VcfPrintStringPointer+VcfPrintStringPointerLength,"\t");

    if(x<0.0005 && y<0.0005)
    {
        if(ThisVariable.GT)
        {

            VcfPrintStringPointerLength+=sprintf(VcfPrintStringPointer+VcfPrintStringPointerLength,"0|0");
            colonIndex=true;
        }
        if(ThisVariable.DS)
        {

            if(colonIndex)
                VcfPrintStringPointerLength+=sprintf(VcfPrintStringPointer+VcfPrintStringPointerLength,":");
            VcfPrintStringPointerLength+=sprintf(VcfPrintStringPointer+VcfPrintStringPointerLength,"0");
            colonIndex=true;
        }
        if(ThisVariable.HDS)
        {

            if(colonIndex)
                VcfPrintStringPointerLength+=sprintf(VcfPrintStringPointer+VcfPrintStringPointerLength,":");
            VcfPrintStringPointerLength+=sprintf(VcfPrintStringPointer+VcfPrintStringPointerLength,"0,0");
            colonIndex=true;
        }
        if(ThisVariable.GP)
        {

            if(colonIndex)
                VcfPrintStringPointerLength+=sprintf(VcfPrintStringPointer+VcfPrintStringPointerLength,":");
            colonIndex=true;
            VcfPrintStringPointerLength+=sprintf(VcfPrintStringPointer+VcfPrintStringPointerLength,"1,0,0");
        }
        if(ThisVariable.SD)
        {
            if(colonIndex)
                VcfPrintStringPointerLength+=sprintf(VcfPrintStringPointer+VcfPrintStringPointerLength,":");
            colonIndex=true;
            VcfPrintStringPointerLength+=sprintf(VcfPrintStringPointer+VcfPrintStringPointerLength,"0");
        }
        return;
    }


    if(ThisVariable.GT)
    {

        VcfPrintStringPointerLength+=sprintf(VcfPrintStringPointer+VcfPrintStringPointerLength,"%d|%d",(x>0.5),(y>0.5));
        colonIndex=true;
    }
    if(ThisVariable.DS)
    {

        if(colonIndex)
            VcfPrintStringPointerLength+=sprintf(VcfPrintStringPointer+VcfPrintStringPointerLength,":");
        VcfPrintStringPointerLength+=sprintf(VcfPrintStringPointer+VcfPrintStringPointerLength,"%.3f",x+ y);
        colonIndex=true;
    }
    if(ThisVariable.HDS)
    {

        if(colonIndex)
            VcfPrintStringPointerLength+=sprintf(VcfPrintStringPointer+VcfPrintStringPointerLength,":");
        VcfPrintStringPointerLength+=sprintf(VcfPrintStringPointer+VcfPrintStringPointerLength,"%.3f,%.3f",x , y);
        colonIndex=true;
    }
    if(ThisVariable.GP)
    {

        if(colonIndex)
            VcfPrintStringPointerLength+=sprintf(VcfPrintStringPointer+VcfPrintStringPointerLength,":");
        colonIndex=true;
        VcfPrintStringPointerLength+=sprintf(VcfPrintStringPointer+VcfPrintStringPointerLength,"%.3f,%.3f,%.3f",(1-x)*(1-y),x*(1-y)+y*(1-x),x*y);
    }
    if(ThisVariable.SD)
    {
        if(colonIndex)
            VcfPrintStringPointerLength+=sprintf(VcfPrintStringPointer+VcfPrintStringPointerLength,":");
        colonIndex=true;
        VcfPrintStringPointerLength+=sprintf(VcfPrintStringPointer+VcfPrintStringPointerLength,"%.3f", x*(1-x) + y*(1-y));
    }



}

void MetaMinimac::PrintHaploidDosage(float &x)
{
    bool colonIndex=false;
    VcfPrintStringPointerLength+=sprintf(VcfPrintStringPointer+VcfPrintStringPointerLength,"\t");

    if(x<0.0005)
    {
        if(ThisVariable.GT)
        {
            VcfPrintStringPointerLength+=sprintf(VcfPrintStringPointer+VcfPrintStringPointerLength,"0");
            colonIndex=true;
        }
        if(ThisVariable.DS)
        {

            if(colonIndex)
                VcfPrintStringPointerLength+=sprintf(VcfPrintStringPointer+VcfPrintStringPointerLength,":");
            VcfPrintStringPointerLength+=sprintf(VcfPrintStringPointer+VcfPrintStringPointerLength,"0");
            colonIndex=true;
        }
        if(ThisVariable.HDS)
        {

            if(colonIndex)
                VcfPrintStringPointerLength+=sprintf(VcfPrintStringPointer+VcfPrintStringPointerLength,":");
            VcfPrintStringPointerLength+=sprintf(VcfPrintStringPointer+VcfPrintStringPointerLength,"0");
            colonIndex=true;
        }
        if(ThisVariable.GP)
        {

            if(colonIndex)
                VcfPrintStringPointerLength+=sprintf(VcfPrintStringPointer+VcfPrintStringPointerLength,":");
            colonIndex=true;
            VcfPrintStringPointerLength+=sprintf(VcfPrintStringPointer+VcfPrintStringPointerLength,"1,0");
        }
        if(ThisVariable.SD)
        {
            if(colonIndex)
                VcfPrintStringPointerLength+=sprintf(VcfPrintStringPointer+VcfPrintStringPointerLength,":");
            colonIndex=true;
            VcfPrintStringPointerLength+=sprintf(VcfPrintStringPointer+VcfPrintStringPointerLength,"0");
        }
        return;

    }



    if(ThisVariable.GT)
    {
        VcfPrintStringPointerLength+=sprintf(VcfPrintStringPointer+VcfPrintStringPointerLength,"%d",(x>0.5));
        colonIndex=true;
    }
    if(ThisVariable.DS)
    {

        if(colonIndex)
            VcfPrintStringPointerLength+=sprintf(VcfPrintStringPointer+VcfPrintStringPointerLength,":");
        VcfPrintStringPointerLength+=sprintf(VcfPrintStringPointer+VcfPrintStringPointerLength,"%.3f",x);
        colonIndex=true;
    }
    if(ThisVariable.HDS)
    {

        if(colonIndex)
            VcfPrintStringPointerLength+=sprintf(VcfPrintStringPointer+VcfPrintStringPointerLength,":");
        VcfPrintStringPointerLength+=sprintf(VcfPrintStringPointer+VcfPrintStringPointerLength,"%.3f",x );
        colonIndex=true;
    }
    if(ThisVariable.GP)
    {

        if(colonIndex)
            VcfPrintStringPointerLength+=sprintf(VcfPrintStringPointer+VcfPrintStringPointerLength,":");
        colonIndex=true;
        VcfPrintStringPointerLength+=sprintf(VcfPrintStringPointer+VcfPrintStringPointerLength,"%.3f,%.3f",1-x,x);
    }
    if(ThisVariable.SD)
    {
        if(colonIndex)
            VcfPrintStringPointerLength+=sprintf(VcfPrintStringPointer+VcfPrintStringPointerLength,":");
        colonIndex=true;
        VcfPrintStringPointerLength+=sprintf(VcfPrintStringPointer+VcfPrintStringPointerLength,"%.3f", x*(1-x));
    }
}

void MetaMinimac::ReadCurrentDosageData()
{

    for(int i=0; i<NoStudiesHasVariant; i++)
    {
        int index = StudiesHasVariant[i];
        InputData[index].LoadHapDoseVariant(CurrentRecordFromStudy[index]->getGenotypeInfo());
    }
}

void MetaMinimac::CreateMetaImputedData()
{


    if(NoStudiesHasVariant==1)
    {
        CurrentMetaImputedDosage=InputData[StudiesHasVariant[0]].CurrentHapDosage;
        return;
    }


    fill(CurrentMetaImputedDosage.begin(),CurrentMetaImputedDosage.end(),0.0);
    fill(CurrentMetaImputedDosageSum.begin(),CurrentMetaImputedDosageSum.end(),1e-6);
    for(int i=0; i<NoStudiesHasVariant; i++)
    {
        int index=StudiesHasVariant[i];
        for(int j=0; j<InputData[0].numHaplotypes; j++)
        {
            CurrentMetaImputedDosageSum[j]+= LSQEstimates[j][index];
            CurrentMetaImputedDosage[j]+=  LSQEstimates[j][index] * InputData[index].CurrentHapDosage[j];
        }
    }
    for(int j=0; j<InputData[0].numHaplotypes; j++)
    {
        CurrentMetaImputedDosage[j]/=CurrentMetaImputedDosageSum[j];
    }

}

void MetaMinimac::AppendtoMainVcfFaster(int ChunkNo)
{
    VcfPrintStringPointerLength=0;

    vcfdosepartial = ifopen(ThisVariable.outfile + ".metaDose.vcf" + (ThisVariable.gzip ? ".gz" : ""), "a", ThisVariable.gzip ?InputFile::BGZF:InputFile::UNCOMPRESSED);

    do
    {
        FindCurrentMinimumPosition();
        if(CurrentFirstVariantBp>ChunkList[ChunkNo].EndBp)
            break;
        ReadCurrentDosageData();
        CreateMetaImputedData();
        PrintCurrentVariant();
        PrintMetaImputedData();
        UpdateCurrentRecords();
    }while(true);


    ifprintf(vcfdosepartial,"%s",VcfPrintStringPointer);
    ifclose(vcfdosepartial);
}




void MetaMinimac::AppendtoMainWeightsFile(int ChunkNo)
{
    VcfPrintStringPointerLength=0;

    metaWeight = ifopen(ThisVariable.outfile + ".metaWeights.vcf" + (ThisVariable.gzip ? ".gz" : ""), "a", ThisVariable.gzip ?InputFile::BGZF:InputFile::UNCOMPRESSED);
    ifprintf(metaWeight,"Chunk_%d",ChunkNo);

    for (int j = 0; j < InputData[0].numSamples; j++)
    {
        ifprintf(metaWeight,"\t");
        if (InputData[0].SampleNoHaplotypes[j] == 2)
        {
            PrintWeightForHaplotype(2*j);
            ifprintf(metaWeight,"|");
            PrintWeightForHaplotype(2*j+1);
        }
        else
            PrintWeightForHaplotype(2*j);
    }

    ifprintf(metaWeight,"\n");
    ifclose(metaWeight);
}



void MetaMinimac::PrintWeightForHaplotype(int haploId)
{

    ifprintf(metaWeight,"%0.5f",LSQEstimates[haploId][0]);
    for(int i=1;i<NoInPrefix;i++)
        ifprintf(metaWeight,",%0.5f",LSQEstimates[haploId][i]);

}

