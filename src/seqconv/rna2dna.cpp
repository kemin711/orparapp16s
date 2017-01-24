#include <cstring>
#include <iostream>
#include <fstream>

#include <bioseq.h>

using namespace std;
using namespace orpara;

void usage() {
   cerr << "Usage: rna2dna -i input.fasta -o output.fasta\n"
      << " or fastreduceid input.fasta output.fasta\n";
   exit(1);
}
string nameOut(const string &inf) {
   string::size_type i=inf.rfind('.');
   return inf.substr(0,i) + "dna.fas";
}

int main(int argc, char* argv[]) {
   string infile, outfile;
   int i=1;
   while (i<argc) {
      if (!strcmp(argv[i], "-i")) infile=argv[++i];
      else if (!strcmp(argv[i], "-o")) outfile=argv[++i];
      else {
         if (argv[i][0] == '-') {
            cerr << "unsupported option: " << argv[i] << endl;
            return 1;
         }
         infile = argv[i++];
         if (i < argc) {
            outfile = argv[i];
         }
      }
      ++i;
   }
   if (infile.empty()) usage();
   if (outfile.empty()) outfile=nameOut(infile);
   cerr << "input file: " << infile << " output file: " << outfile << endl;

   ifstream inf(infile.c_str());
   if (inf.fail()) {
      cerr << "Failed to open " << infile << endl;
      return 1;
   }
   ofstream ouf(outfile.c_str());
   if (ouf.fail()) {
      cerr << "Failed to open " << outfile << endl;
      return 1;
   }

   int cnt=0;
   DNA seq;
   int every=10;
   while (seq.read(inf)) {
      ++cnt;
      ouf << seq;
      if (cnt % every == 0) {
         cerr << cnt << " sequences processed\n";
         every *= 2;
      }
   }
   cout << cnt << " fasta sequences processed into " << outfile << endl;

   return 0;
}
