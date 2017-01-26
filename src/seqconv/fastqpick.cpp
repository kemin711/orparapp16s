// (C) Kemin Zhou 2012 at orpara.com
#include <fstream>
#include <iostream>
#include <cstdlib>
#include <stdexcept>
#include <set>

#include <fastq.h>
#include <stddev.h>

using namespace std;
using namespace orpara;

void usage() {
      cout << "Usage: fastqpick -l <length cutoff: default 180> input-fastqfile output-fastafile\n";
      cout << "or fastqpick -q <qvalue cutoff: default 20> input-fastqfile output-fastafile\n";
      cout << "or fastqpick -a <accuracy cutoff: default 0.991> input-fastqfile output-fastafile\n";
      cout << "or fastqpick -n or --id-list <idlistFile> infastq outfasq\n";
      exit(1);
}

void filterByQval(const string &infile, const string &outfile, const double Qcut, bool writeFasta);
void filterByLength(const string &infile, const string &outfile, const size_t lencut);
double accuracy2PhredQ(const double accuracy) throw(range_error) {
   if (accuracy > 1 || accuracy < 0.1) {
      throw range_error("accuracy must be between 0 and 1");
   }
   return Fastq::p2q(1-accuracy);
}
set<string> readId(const string &idfile);
// compiler got a bug the function must be defined before main for this one.
// not sure why?
//void pickByIdFile(const string $infile, const string &idlistfile, const string &outfile);
void pickByIdFile(const string &infile, const string &idlistfile, const string &outfile) {
   ifstream inf(infile.c_str());
   ofstream ouf(outfile.c_str());
   if (inf.fail()) {
      cerr << "Failed to open " << infile << endl;
      exit(1);
   }
   if (ouf.fail()) {
      cerr << "Failed to open output file: " << outfile << endl;
      exit(1);
   }
   cerr << "picking fastq sequences by listed ids in file: " << idlistfile 
      << " into : " << outfile << endl;
   set<string> ids = readId(idlistfile);
   set<string>::const_iterator sit;
   Fastq fasq;
   unsigned int j = 0;
   unsigned int total = 0;
   while (fasq.read(inf)) {
      ++total;
      if (ids.find(fasq.getName()) != ids.end()) {
         ouf << fasq;
         ++j;
      }
   }
   cout << j << " of " << total << " fastq sequences picked and written to " 
      << outfile << endl;
}


/**
 * Helper program to copy fastq sequence by
 * given list of id, Phred score, accuracy, or length.
 * Other features of sequence could be added.
 */
int main(int argc, char* argv[]) {
   string infile, outfile, idlistFile;
   size_t lengthcut=180;
   double QscoreCut=20;
   int i=1;
   // action l by length, q by qval, i by idlist
   char action='l';
   bool outputFas = true;
   while (i < argc) {
      if (string(argv[i]) == "-i") infile = argv[++i];
      else if (string(argv[i]) == "-o") outfile = argv[++i];
      else if (string(argv[i]) == "--no-fasta") outputFas = false;
      else if (string(argv[i]) == "-l") {
         lengthcut = atoi(argv[++i]);
      }
      else if (string(argv[i]) == "-q") {
         QscoreCut = atof(argv[++i]);
         action='q';
      }
      else if (string(argv[i]) == "-a") { // accuracy such as 0.999
         // this program will convert it to Qval
         QscoreCut = accuracy2PhredQ(atof(argv[++i]));
         action='q';
      }
      else if (string(argv[i]) == "--id-list" || string(argv[i]) == "-n") {
         idlistFile = argv[++i];
         action='i';
      }
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
   unsigned int j;
   if (action == 'l') {
      if (outfile.empty()) {
         j = infile.rfind('.');
         if (j != string::npos) outfile = infile.substr(0,j);
         else outfile = infile;
         outfile += ".lengthcut.fastq";
      }
      filterByLength(infile, outfile, lengthcut);
   }
   else if (action == 'q') {
      if (outfile.empty()) {
         j = infile.rfind('.');
         if (j != string::npos) outfile = infile.substr(0,j);
         else outfile = infile;
         outfile += ".qcut.fastq";
      }
      filterByQval(infile, outfile, QscoreCut, outputFas);
   }
   else if (action == 'i') {
      if (outfile.empty()) {
         j = infile.rfind('.');
         if (j != string::npos) outfile = infile.substr(0,j);
         else outfile = infile;
         outfile += ".byid.fastq";
      }
      pickByIdFile(infile, idlistFile, outfile);
   }
   else {
      cerr << "action " << action << " not specified\n";
      return 1;
   }
   return 0;
}

set<string> readId(const string &idfile) {
   set<string> result;
   ifstream inf(idfile);
   string id;
   inf >> id;
   while (!inf.eof()) {
      result.insert(id);
      inf >> id;
   }
   return result;
}

void filterByLength(const string &infile, const string &outfile, const size_t lencut) {
   ifstream inf(infile.c_str());
   ofstream ouf(outfile.c_str());
   if (inf.fail()) {
      cerr << "Failed to open " << infile << endl;
      exit(1);
   }
   if (ouf.fail()) {
      cerr << "Failed to open output file: " << outfile << endl;
      exit(1);
   }
   cerr << "picking fastq sequences longer than : " << lencut 
      << " into : " << outfile << endl;

   Fastq fasq;
   unsigned int j = 0;
   unsigned int total = 0;
   while (fasq.read(inf)) {
      ++total;
      if (fasq.length() > lencut) {
         ouf << fasq;
         ++j;
      }
   }
   cout << j << " of " << total << " fastq sequences picked and written to " 
      << outfile << endl;
}


void filterByQval(const string &infile, const string &outfile, const double Qcut,
      bool writeFasta) {
   ifstream inf(infile.c_str());
   ofstream ouf(outfile.c_str());
   ofstream oufas;
   if (inf.fail()) {
      cerr << "Failed to open input file: " << infile << endl;
      exit(1);
   }
   if (ouf.fail()) {
      cerr << "Failed to open output file: " << outfile << endl;
      exit(1);
   }
   string fasfileName;
   if (writeFasta) {
      fasfileName = outfile.substr(0, outfile.rfind('.')) + ".fasta";
      oufas.open(fasfileName.c_str());
      if (oufas.fail()) {
         cerr << "failed to open " << fasfileName << endl;
         exit(1);
      }
   }
   cerr << "picking fastq sequences with average score greater than : " << Qcut 
      << " into : " << outfile << endl;

   Fastq fasq;
   int j = 0;
   int total = 0;
   while (fasq.read(inf)) {
      ++total;
      if (fasq.getAverageQuality() > Qcut) {
         ouf << fasq;
         if (writeFasta) fasq.writeFasta(oufas);
         ++j;
      }
   }
   ouf.close();
   oufas.close();
   cout << j << " out of " << total << " fastq sequences picked and written to " 
      << outfile << endl;
   if (writeFasta)
      cout << "fasta file: " << fasfileName << endl;
}
