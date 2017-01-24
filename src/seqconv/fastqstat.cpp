#include <fstream>
#include <iostream>
#include <cstdlib>
#include <map>

#include <fastq.h>
#include <stddev.h>

using namespace std;
using namespace orpara;

void usage() {
   cerr << "Usage: fastqstat -l <len distribution file> input-fastqfile outputfile\n"
      << "Options:\n"
      << "  -i inputFastqFile input fastq file\n"
      << "  -o outputFile summar statistics result to be written to\n"
      << "  -l lengthDistributionFile length distribution file\n";
   exit(1);
}

void filemap(ostream &ous, const map<int,int>& lens);
string buildOutputName(const string &infile, const string &suffix);

int main(int argc, char* argv[]) {
   string infile, outfile;
   string lenfile; //("length.dis");
   int i=1;
   while (i < argc) {
      if (string(argv[i]) == "-i") infile = argv[++i];
      else if (string(argv[i]) == "-o") outfile = argv[++i];
      else if (string(argv[i]) == "-l") lenfile = argv[++i];
      else if (string(argv[i]) == "--help") usage();
      else {
         infile=argv[i];
         if (i+1 < argc && argv[i+1][0] != '-') {
            ++i;
            outfile = argv[i];
         }
      }
      ++i;
   }
   if (infile.empty()) {
      usage();
   }

   if (lenfile.empty()) {
      lenfile = buildOutputName(infile, ".lendis.tab");
   }
   if (outfile.empty()) {
      outfile = buildOutputName(infile, ".stat");
   }

   ifstream inf(infile.c_str());
   ofstream ouf(outfile.c_str());
   if (inf.fail()) {
      cerr << "Failed to open " << infile << endl;
      return 1;
   }
   if (ouf.fail()) {
      cerr << "Failed to open " << outfile << endl;
   }
   cout << "Calculating simple length stat for fastq sequeces from file: " 
      << infile << endl;
   stddev avgstd;
   Fastq fasq;
   ofstream len(lenfile.c_str());
   if (len.fail()) {
      cerr << "Failed to open lengh distribution file: " << lenfile << endl;
      return 1;
   }
   map<int,int> lencnt;
   while (fasq.read(inf)) {
      avgstd(fasq.length());
      ++lencnt[fasq.length()];
   }
   inf.close();
   filemap(len, lencnt);
   ouf << avgstd << endl;
   cout << avgstd << endl;
   cerr << "Simple length statistics written to " << outfile << endl;
   cerr << "Length distribution written to " << lenfile << endl;

   return 0;
}

string buildOutputName(const string &infile, const string &suffix) {
   string oufile;
   size_t j = infile.rfind('.');
   if (j != string::npos) oufile = infile.substr(0,j);
   else oufile = infile;
   oufile += suffix;
   return oufile;
}

void filemap(ostream &ous, const map<int,int>& lens) {
   map<int,int>::const_iterator it = lens.begin();
   cerr << lens.size() << " map entries to write\n";
   ous << "length\tcount\n";
   while (it != lens.end()) {
      ous << it->first << '\t' << it->second << endl;
      ++it;
   }
}
