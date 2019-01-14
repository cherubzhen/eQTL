#include <unordered_map>
#include <stdlib.h>
#include <iostream>
#include <vector>
#include <sstream>
#include <fstream>
#include <limits.h>
#include <algorithm>
using namespace std;

class lib
{
public:
  static const int CIS_ONLY = 1;
  static const int TRANS_ONLY = 0;
  static const int CIS_TRANS = 2;
  static const int NO_DIRECT = 0;
  static const int MEAN_DIRECT = 1;
  static const int MEDIAN_DIRECT = 2;
  static vector<string> strsplit(string s, char split)
  {
    vector<string> res;
    stringstream sline(s);
    string temp;
    while(getline(sline,temp,split)) res.push_back(temp);
    //if(res.size()==0) res.push_back(s);
    return res;
  }
};

class GP
{
public:
  vector<string> dats;
  ifstream* directFile;
  ifstream* directIdxFile;
  ifstream* genotypingFile; // Imputed Genotyping file
  ifstream geneFile;
  ifstream snpFile;	// SNP list
  ofstream outFile;
  bool Direct;
  double MAF;
  double HWE;
  double FDR;
  int cis_type;
  int direct_type;
  
  GP()
  {
    this->MAF=0;
    this->HWE=0;
    this->FDR=1;
    this->cis_type=lib::CIS_TRANS;
    this->direct_type = lib::NO_DIRECT;
  };
  void Print()
  {
    cout<<"Intersect among: ";
    for(int i = 0; i < dats.size();i++) cout<<dats[i]<<" ";cout<<"\n";
    cout<<"MAF: "<< this->MAF<<"\n";
    cout<<"HWE: "<< this->HWE<<"\n";
    cout<<"FDR: "<< this->FDR<<"\n";
    if(this->cis_type != lib::CIS_TRANS) cout<<"Cis/Trans\t"<<this->cis_type<<"\n";
    if(this->direct_type != lib::NO_DIRECT) cout<<"Direction\t"<<this->cis_type<<"\n";
  }
  void setGeneFile(string Genefile)
  {
    geneFile.open(Genefile.c_str(),ifstream::in);
    if(!geneFile.is_open())
    {
      cout<<"Can not open Gene file\n";
      exit(1);
    }
  }
  void setSNPFile(string SNPfile)
  {
    snpFile.open(SNPfile.c_str(),ifstream::in);
    if(!snpFile.is_open())
    {
      cout<<"Can not open SNP file\n";
      exit(1);
    }
  }
  void setOutFile(string outfile)
  {
    outFile.open(outfile.c_str(),ifstream::out);
    if(!outFile.is_open())
    {
      cout<<"Can not creat output file\n";
      exit(1);
    }
  }
  void setDirectFile(string s)
  {
    this->dats = lib:: strsplit(s,',');
    this->directFile = new ifstream[this->dats.size()];
    this->directIdxFile = new ifstream[this->dats.size()];
    for(int i = 0; i < dats.size(); i++)
    {
      string direct = this->dats[i] + "/"+ this->dats[i]+"_final_eQTL_Direct.txt";
      string directIdx = this->dats[i] + "/"+ this->dats[i]+"_final_eQTL_Direct_index.txt";
      this->directFile[i].open(direct.c_str(),ifstream::in);
      this->directIdxFile[i].open(directIdx.c_str(),ifstream::in);
      if(!directFile[i].is_open())
      {
	cout<<"Can not open direction file\n";
	exit(1);
      }
    }
  }
  void SetGenotypingFile()
  {
    this->genotypingFile = new ifstream[this->dats.size()];
    for(int i = 0; i < dats.size(); i++)
    {
      string stemp = this->dats[i] + "/"+ this->dats[i]+"_final_Genotype.txt";
      genotypingFile[i].open(stemp.c_str(),ifstream::in);
      if(!genotypingFile[i].is_open())
      {
	cout<<this->dats[i]<<":Can not open Genotyping file\n";
	exit(1);
      }
    }
  }
};

class Filter
{
public:
  unordered_map<string, int> Gene;
  unordered_map<string, int> SNP;
  unordered_map<string, string>* eQTL;
  double MAF;
  double HWE;
  double FDR;
  bool Direct;
  int cis_type;
  int direct_type;
  
  Filter(GP* pGP)
  {
    this->MAF=pGP->MAF;
    this->HWE=pGP->HWE;
    this->FDR=pGP->FDR;
    this->cis_type = pGP->cis_type;
    this->Direct = pGP->Direct;
    this->eQTL = new unordered_map<string, string>[pGP->dats.size()+1];
    this->direct_type = pGP->direct_type;
  }
  void setGene(GP* pGP)
  {
    string line;
    int iter = 0;
    try{
      while(!pGP->geneFile.eof())
      {
	getline(pGP->geneFile,line);
	if(line.length()>0){Gene.insert(pair<string, int>(line,1));iter++;}
      }
      pGP->geneFile.close();
    }
    catch(...)
    {
      cout<<"Error in reading gene file\n";
    }
    cout<<iter<<" Genes in Genelist\n";
  }
  void setSNP(GP* pGP)
  {
    string line;
    int iter = 0;
    try
    {
      while(!pGP->snpFile.eof())
      {
	getline(pGP->snpFile,line);
	if(line.length()>0){SNP.insert(pair<string, int>(line,1));iter++;}
      }
      pGP->snpFile.close();
    }
    catch(...)
    {
      cout<<"Error in reading snp file\n";
    }
    cout<<iter<<" SNPs in snplist\n";
  }
  
  bool CreatFilter(GP* pGP)
  {
    bool res = true;
    cout<<"Reading Gene/SNP filtering information\n";
    setGene(pGP);
    setSNP(pGP);
    return res;
  }
  
  void filterSNP(GP* pGP)
  {
    bool allSNP = false;
    if(this->SNP.size()==0)
    { 
      allSNP =true;
      this->SNP.insert(pair<string,int>("rsID",1));
    }

    for(int i = 0; i< pGP->dats.size();i++)
    {
      string line;
      vector<string> scontent;

      while(!pGP->genotypingFile[i].eof())
      {
	getline(pGP->genotypingFile[i],line);
	if(line.length()==0) continue;
	scontent=lib::strsplit(line,'\t');
	if(!allSNP && this->SNP.count(scontent[0])>0)
	{
	  if(atof(scontent[4].c_str()) < this->MAF || atof(scontent[5].c_str()) < this->HWE)
	    this->SNP.erase(scontent[0]);
	}
	if(allSNP && atof(scontent[4].c_str()) >= this->MAF && atof(scontent[5].c_str()) >= this->HWE)
	    this->SNP.insert(pair<string, int>(scontent[0],1));
      }

      pGP->genotypingFile[i].close();
      if(allSNP) allSNP = false;
      cout<<pGP->dats[i]<<":After MAF HWE filtering "<<this->SNP.size()<<" SNPs left\n";
    }
    cout<<"After MAF HWE filtering "<<this->SNP.size()<<" SNPs left\n";
  }
  
  string OrganizeHeader(GP* pGP)
  {
    string header = "SNP\tGene\tis_cis\t";
    for(int i = 0; i < pGP->dats.size(); i++)
    {
      header= header+pGP->dats[i]+"_t-stat\t"+pGP->dats[i]+"_p-value\t"+pGP->dats[i]+"_FDR\t";
      header= header+pGP->dats[i]+"_m0\t"+pGP->dats[i]+"_m1\t"+pGP->dats[i]+"_m2\t";
      if(i < pGP->dats.size()-1)header= header+pGP->dats[i]+"_med0\t"+pGP->dats[i]+"_med1\t"+pGP->dats[i]+"_med2\t";
      else header= header+pGP->dats[i]+"_med0\t"+pGP->dats[i]+"_med1\t"+pGP->dats[i]+"_med2\n";
    }
    return header;
  }
  string OrganizePrefix(vector<string> scontent)
  {
    string res = scontent[0] + "\t"+scontent[1]+"\t"+scontent[5]+"\t";
    return res;
  }
  string OrganizeContent(vector<string> scontent)
  {
    string res = scontent[2] + "\t"+scontent[3]+"\t"+scontent[4]+"\t";
    res += scontent[6] + "\t"+scontent[7]+"\t"+scontent[8]+"\t";
    res += scontent[9] + "\t"+scontent[10]+"\t"+scontent[11]+"\t";
    return res;
  }
  
  vector<int> ReadIdx(GP* pGP, ifstream& in)
  {
    vector<int> res;
    vector<string> scontent;
    string line;
    int iFound = 0;
    while(!in.eof())
    {
      if(iFound == this->Gene.size()) break;
      getline(in,line);
      if(line.length()==0) continue;
      scontent=lib::strsplit(line,'\t');
      if(this->Gene.count(scontent[0]) == 0) continue;
      iFound++;
      for(int i = 1; i < scontent.size();i++)
	res.push_back(atoi(scontent[i].c_str()));
    }
    in.close();
    sort(res.begin(),res.end());
    return res;
  }
  
  void AdjustMemory(GP* pGP, bool* bFinish, int& isFinish)
  {
    cout<<"Adjusting Memory\n";
    unordered_map<string, string>::iterator it;
    isFinish = -1;
    for(int i = 2; i < pGP->dats.size();i++)
    {
      if(bFinish)
      for(it = eQTL[i].begin();it != eQTL[i].end();it++)
      {
	if(eQTL[0].count(it->first)==0)
	  eQTL[i].erase(it);
      }
    }
    cout<<"Finish adjustment\n";
  }
  
  void filtereQTL(GP* pGP)
  {
    string finalcontent;
    unordered_map<string, string>::iterator it,iit;      
    pGP->outFile<<OrganizeHeader(pGP);
    int isFinish = -1, icount = 0;
    bool hasFinished = false;
    bool* bFinish = new bool[pGP->dats.size()];

    #pragma omp parallel for num_threads(4) shared(isFinish, bFinish,pGP,hasFinished)
    for(int i = 0; i < pGP->dats.size(); i++)
    {
      vector<string> scontent;
      string line;
      vector<int> Idx;
      if(pGP->directIdxFile[i].is_open())
      {
	cout<<pGP->dats[i]<<":Reading index...\n";
	Idx = ReadIdx(pGP,pGP->directIdxFile[i]);
	if(Idx.size()==0) Idx.push_back(-1);
      }
      else Idx.push_back(-1);
      cout<<pGP->dats[i]<<":Filtering eQTL...\n";
      int inum = this->SNP.size()*this->Gene.size();
      if(inum == 0) inum=INT_MAX;
      int iter = 0;

      while(!pGP->directFile[i].eof() && inum > 0)
      {
	//if(isFinish>-1)AdjustMemory(pGP,bFinish,isFinish);
	getline(pGP->directFile[i],line);
	if(Idx.size()==0) break;
	if(Idx[0] > 0 && iter++ < Idx[0]) continue;
	else if(Idx[0] > 0) Idx.erase(Idx.begin());
	
	if(line.length()==0) continue;
	scontent=lib::strsplit(line,'\t');
	if(this->SNP.size() > 0 && this->SNP.count(scontent[0])<=0) continue;
	if(this->Gene.size() > 0 && this->Gene.count(scontent[1])<=0) continue;
	inum--;
	if(this->cis_type == lib::TRANS_ONLY && atoi(scontent[5].c_str()) == 1) continue;
	if(this->cis_type == lib::CIS_ONLY && atoi(scontent[5].c_str()) == 0) continue;
	if(i == 0 && this->FDR >= atof(scontent[4].c_str())) 
	{
	  eQTL[i].insert(pair<string,string>(scontent[0]+"_"+scontent[1],OrganizePrefix(scontent)));
	  eQTL[i+1].insert(pair<string,string>(scontent[0]+"_"+scontent[1],OrganizeContent(scontent)));
	}
	else if(i > 0 && this->FDR >= atof(scontent[4].c_str()))
	{
	  if(!hasFinished || hasFinished && eQTL[0].count(scontent[0]+"_"+scontent[1]) > 0)
	    eQTL[i+1].insert(pair<string,string>(scontent[0]+"_"+scontent[1],OrganizeContent(scontent)));
	}
      }
      pGP->directFile[i].close();
      cout<<pGP->dats[i]<<":done\n";
      bFinish[i] = true;
      isFinish = i;
      hasFinished = true;
    }
    
    cout<<"Writing File...\n";
    vector<string> scontent;
    for(it = eQTL[0].begin();it!=eQTL[0].end(); it++)
    {
      int iNum = 0;
      finalcontent=it->second;
      bool isUp = false;
      bool isDown = false;
      for(int i = 0; i < pGP->dats.size();i++)
      {
	iit = eQTL[i+1].find(it->first);
	if(iit != eQTL[i+1].end())
	{
	  finalcontent+=iit->second;
	  iNum++;
	}
      }
      if(iNum == pGP->dats.size() && this->direct_type == lib::NO_DIRECT)
      {
	 icount++;
	finalcontent[finalcontent.length()-1]='\n';
	pGP->outFile<<finalcontent;
      }
      else if(iNum == pGP->dats.size() && this->direct_type != lib::NO_DIRECT)
      {
	scontent=lib::strsplit(finalcontent,'\t');
	if(this->direct_type == lib::MEAN_DIRECT)
	{
	  for(int i = 0; i < pGP->dats.size();i++)
	    if(atof(scontent[6+i*9].c_str())>atof(scontent[8+i*9].c_str()))
	      isDown = true;
	    else
	      isUp = true;
	}
	else
	{
	  for(int i = 0; i < pGP->dats.size();i++)
	    if(atof(scontent[9+i*9].c_str())>atof(scontent[11+i*9].c_str()))
	      isDown = true;
	    else
	      isUp = true;
	}
	if(isDown != isUp) 
	{
	   icount++;
	  finalcontent[finalcontent.length()-1]='\n';
	  pGP->outFile<<finalcontent;
	}
      }
    }
      
    pGP->outFile.close();
    cout<<icount<<" eQTLs left after filtering\n";
  }
};

bool GetOpt(int argc, char** argv, GP* pGP)
{
  string temp;
  bool res = true;
  for(int i = 1; i < argc; i++)
  {
    temp = string(argv[i]);
    if(temp == "-I")
    {
      temp = string(argv[++i]);
      pGP->setDirectFile(temp);
    }
    else if(temp == "-G")
    {
      temp = string(argv[++i]);
      pGP->setGeneFile(temp);
    }
    else if(temp == "-D")
    {
      temp = string(argv[++i]);
      pGP->direct_type =atoi(temp.c_str());
    }
    else if(temp == "-C")
    {
      temp = string(argv[++i]);
      pGP->cis_type = atoi(temp.c_str());
    }
    else if(temp == "-S")
    {
      temp = string(argv[++i]);
      pGP->setSNPFile(temp);
    }
    else if(temp == "-M")
    {
      temp = string(argv[++i]);
      pGP->MAF= atof(temp.c_str());
    }
    else if(temp == "-H")
    {
      temp = string(argv[++i]);
      pGP->HWE= atof(temp.c_str());
    }
    else if(temp == "-O")
    {
      temp = string(argv[++i]);
      pGP->setOutFile(temp);
    }
    else if(temp == "-F")
    {
      temp = string(argv[++i]);
      pGP->FDR= atof(temp.c_str());
    }
    else if(temp == "-h")
    {
      cout<<"Usage Filter: -I dataset(, as separator)\n";
      cout<<"\t -O Output file\n";
      cout<<"\t [option] -G Gene list file\n";
      cout<<"\t [option] -S SNP list file\n";
      cout<<"\t [option] -M MAF cutoff\n";
      cout<<"\t [option] -H HWE cutoff\n";
      cout<<"\t [option] -F FDR cutoff\n";
      cout<<"\t [option] -D direction calculated either by (mean: 1) or (median: 2)\n";
      cout<<"\t [option] -C cis eQTL indicator (0: trans only; 1: cis only; 2[default]: both)\n";
      exit(1);
    }
    else
    {
      cout<<"Unknown parameter \n";
      res = false;
    }
  }
  return res;
}

int main(int argc, char** argv)
{
  GP* pGP = new GP();
  if(!GetOpt(argc,argv, pGP)) exit(1);
  pGP->Print();
  Filter* pF = new Filter(pGP);
  if(!(pF->CreatFilter(pGP))) exit(1);
  if(pGP->HWE>0 || pGP->MAF > 0)
  {
    pGP->SetGenotypingFile();
    pF->filterSNP(pGP);
  }
  pF->filtereQTL(pGP);
  cout<<"done\n";
}