// (C) 2011 Kemin Zhou at orpara.com
#include <cstring>
#include <iostream>
#include <fstream>

#include "bioseq.h"

using namespace std;
using namespace orpara;

string nameOut(const string &infname, const string &prefix) {
   size_t i=infname.rfind('.');
   string outname;
   if (i != string::npos) {
      outname = infname.substr(0, i+1) + "shortid" + infname.substr(i);
   }
   else 
      outname = infname + ".shortid.fasta";
   return outname;
}

void usage() {
   cerr << "Usage: fastareduceid --prefix tag -i input.fasta -o output.fasta\n"
      << " or fastreduceid -p newprefix input.fasta output.fasta\n";
   exit(1);
}

/** 
 * Sequencers tends to produce very long sequence names
 * This can make sequence objects occupying large amounts
 * of memory. In making the ids shorter, we can gain
 * performance.
 */
int main(int argc, char* argv[]) {
   string infile, outfile, prefix;
   int i=1;
   while (i<argc) {
      if (!strcmp(argv[i], "-i")) infile=argv[++i];
      else if (!strcmp(argv[i], "-o")) outfile=argv[++i];
      else if (!strcmp(argv[i], "--prefix") || !strcmp(argv[i], "-p")) prefix=argv[++i];
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
   if (prefix.empty()) prefix="TAG";
   if (outfile.empty()) outfile=nameOut(infile, prefix);
   cerr << "input file: " << infile << " output file: " << outfile << endl;

   bioseq seq;
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

   string header;
   int id=1;
   while (seq.read(inf, header)) {
      if (!seq.hasTitle()) {
         seq.setTitle(seq.getName());
      }
      else {
         seq.setTitle(seq.getName() + " " + seq.getTitle());
      }
      seq.setName(prefix + to_string(id++));
      ouf << seq;
   }
   cout << id << " fasta sequences with shorter id written to " 
      << outfile << endl;

   return 0;
}
