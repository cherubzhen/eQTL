#include <unordered_map>
#include <sstream>
#include <stdlib.h>
#include <fstream>
#include <vector>
#include <iostream>
using namespace std;

vector<string> split(string s, char split)
{
   vector<string> res;
   stringstream sline(s);
   string temp;
   while(getline(sline,temp,split)) res.push_back(temp);
    //if(res.size()==0) res.push_back(s);
   return res;
}
int main(int argc, char** argv)
{
  string temp;
  int index;
  string filename;
  string output;
  for(int i = 0; i < argc; i++)
  {
    temp = string(argv[i]);
    if(temp == "-f")
    {
      temp = string(argv[++i]);
      index = atoi(temp.c_str());
    }
    else if(temp == "-o")
    {
      temp = string(argv[++i]);
      output = temp;
    }
  }
  filename= string(argv[argc-1]);
  ifstream inFile(filename.c_str());
  string line;
  vector<string> scontent;
  unordered_map<string,vector<int> > tree;
  unordered_map<string,vector<int> >::iterator it;
  int iter = 0;
  cout<<"creating index...\n";
  while(!inFile.eof())
  {
    getline(inFile,line);
    if(line.length()>0)
    {
      scontent = split(line,'\t');
      if(tree.count(scontent[1]) == 0)
      {
	vector<int> Idx;
	tree.insert(pair<string,vector<int> >(scontent[index],Idx));
      }
      else
      {
	it = tree.find(scontent[index]);
	it->second.push_back(iter);
      }
    }
    iter++;
  }
  inFile.close();
  ofstream outputFile(output,ios::out);
  cout<<"saving index...\n";
  for(it = tree.begin();it != tree.end();it++)
  {
    outputFile<<it->first;
    for(int i = 0; i < it->second.size();i++)
      outputFile<<"\t"<<it->second[i];
    outputFile<<"\n";
  }
  outputFile.close();
}