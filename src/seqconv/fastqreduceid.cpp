#include <fstream>
#include <iostream>
#include <cstdlib>
#include <map>
#include <string>

#include <fastq.h>
#include <stddev.h>

using namespace std;
using namespace orpara;

void filemap(ostream &ous, const map<int,int>& lens);
string buildOutputName(const string &infile, const string &tag);
void usage() {
   cerr << "Usage: fastqreduceid --prefix ccs input-fastqfile outputfile\n";
   cerr << "Or fastqreduceid -i input-fastqfile -o outputfile\n"
      << "Options:\n"
      << "   -p or --prefix new prefix for the sequence ids\n"
      << "      The prefix can be any string.  Default is ccs\n"
      << "   -i input fastq file\n"
      << "   -o output fastq file. If not given this program will generate one for you.\n"
      << "   --help or ? will print this message\n"
      << "This can also make a file of sequences with redundant ids unique\n";
   exit(1);
}

/**
 * replace the super long ids with a counter after some
 * user give prefix or default prefix.
 * This simulate the behavior of database serial id.
 */
int main(int argc, char* argv[]) {
   string infile, outfile, prefix="ccs";
   int i=1;
   while (i < argc) {
      if (string(argv[i]) == "-i") infile = argv[++i];
      else if (string(argv[i]) == "-o") outfile = argv[++i];
      else if (string(argv[i]) == "--prefix" || 
            string(argv[i]) == "-p") prefix = argv[++i];
      else if (string(argv[i]) == "--help" || argv[i][0] == '?') 
         usage();
      else {
         infile=argv[i];
         if (i+1 < argc && argv[i+1][0] != '-') {
            ++i;
            outfile = argv[i];
         }
      }
      ++i;
   }
   if (infile.empty()) usage();

   if (outfile.empty()) {
      outfile = buildOutputName(infile, "sid");
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
   cerr << "producing fastq file with nice ids from " << infile << endl;

   Fastq fasq;
   string newname;
   int id=1;
   while (fasq.read(inf)) {
      // C++11 has the to_string method in the string library
      //newname = prefix + itos(id++);
      newname = prefix + to_string(id++);
      if (!fasq.hasDescription()) {
         fasq.setDescription(fasq.getName());
         //cout << "fastq has no title\n";
      }
      else {
         //cout << "fastq title: " << fasq.getDescription() << endl;
         fasq.setDescription(fasq.getName() + " " + fasq.getDescription());
      }
      fasq.setName(newname);
      ouf << fasq;
      if (id % 5000 == 0) 
         cerr << "processed " << id << " fastq sequence\n";
   }
   cerr << id << " total sequences processed\n";
   inf.close();
   ouf.close();
   cout << outfile << endl;

   return 0;
}

string buildOutputName(const string &infile, const string &tag) {
   string oufile = infile;
   size_t j = oufile.rfind('.');
   if (j != string::npos) oufile.insert(j, "." + tag);
   else oufile += ("." + tag + ".fastq");
   return oufile;
}

