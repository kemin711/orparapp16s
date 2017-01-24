#include <iostream>
#include <string>
#include <fstream>
#include <cstring>

#include <bioseq.h>
#include <fastq.h>
#include <strformat.h>
#include "dnaqualstore.h"

using namespace std;
using namespace orpara;

// given a fasta file, it will eliminate identical sequences

class Progparam {
   private:
      string kmerRef;
      bool flipread;
      static string metagData;

   public:
      //Progparam() : kmerRef("/remote/DataAnalysis/metag/refseq/silva123Bacteria.kmc"), flipread(true) { }
      Progparam() : kmerRef(metagData + "/refseq/silva123Bacteria.kmc"), flipread(true) { }
      const string& getKmerReference() const { return kmerRef; }
      bool flip() const { return flipread; }
      /**
       * explicitly set up the location of the kmer reference file
       */
      void setKmerReference(const string &file) {
         kmerRef=file; 
      }
      void autoconfig() {
         ifstream inf(kmerRef);
         if (inf.fail()) {
            inf.close();
            kmerRef="/home/kzhou/work/metag/refseq/silva123Bacteria.kmc";
         }
      }
      static void setMetagData(const string &dirname) { metagData = dirname; }
};

void usage() {
   cerr << "fastabund <seqfile.fas or seq.fastq> <output fas file> \n"
      << "  or fastabund -i xxy.fastq -o xxy.fas\n"
      << "  or fastabund xxy.fastq xxy.fas\n"
      << " you can give a sequence in fasta or fastq format, it will\n"
      << " combine all identical sequences from the input file\n"
      << " Then it will write into a new file with extension _cnt.fas\n"
      << "Options:\n"
      << "   --aq FLAG to tell the program to treat fasta as\n"
      << "        fastq using default Q values for every base\n"
      << "        This is for pipelines where fastq format is\n"
      << "        expected but only fasta file is available\n";
   exit(1);
}
int readFasta(const string &fname, map<DNA,int> &seqcnt);
void writeFasta(const map<DNA, int> &seqcnt, const string &outfile, int total);
float processFastq(const string &fname, const string &outf, const Progparam &par);
float processFasta(const string &fname, const string &outf, const Progparam &par);
/*
bool endwith(const string &txt, const string &suffix) {
   if (suffix.length() >= txt.length()) {
      throw out_of_range("suffix longer than text");
   }
   string tmp = txt.substr(txt.length()-suffix.length());
   return tmp == suffix;
}
*/

/**
 * Given a fasta or fastq file, this program will combined
 * duplicated reads into one sequences. It then writes the 
 * result to a *.fasta file (if input is fastq) and a *.dqc file.
 * The fasta file will use seqid_count as the new sequence id.
 * This can be used for swarm or mothur.
 */
int main(int argc, char* argv[]) {
   string infile, outfile;
   //infile="/home/kzhou/work/metag/refseq/addon.uniq.fas"; // for debug
   Progparam param;
   param.autoconfig();
   bool aAsQ = true;
   int i = 1;
   while (i < argc) {
      if (string(argv[i]) == "--help") {
         usage();
      }
      else if (!strcmp(argv[i], "-i")) infile=argv[++i];
      else if (!strcmp(argv[i], "-o")) outfile=argv[++i];
      else if (!strcmp(argv[i], "--aq")) aAsQ=true;
      else {
         infile=string(argv[i]);
         if (i+1 < argc && argv[i+1][0] != '-') {
            outfile=argv[++i];
         }
      }
      ++i;
   }
   if (infile.empty()) usage();
   if (outfile.empty()) {
      outfile = infile.substr(0, infile.rfind('.')) 
         + "_cnt.fas";
   }
   char *md = getenv("METAG_DATA");
   if (md == NULL) {
      cerr << "METAG_DATA environment variable not set\n";
      exit(1);
   }
   Progparam::setMetagData(md);
   param.setKmerReference(string(md) + "/refseq/silva123Bacteria.kmc");

   map<DNA, int> seqcnt;
   int totalNumseq = 0;
   if (endwith(infile, "fasta") || endwith(infile, "fas")) {
      if (aAsQ) {
         processFasta(infile, outfile, param);
      }
      else {
         totalNumseq=readFasta(infile, seqcnt);
         writeFasta(seqcnt, outfile, totalNumseq);
      }
   }
   else if (endwith(infile, "fastq")) {
      processFastq(infile, outfile, param);
   }
   else {
      cerr << "file not ending with proper suffix\n";
      return 1;
   }

   // output the result

   return 0;
}

float processFastq(const string &fname, const string &outf, const Progparam &par) {
   cerr << "merging identical sequences from " << fname << endl;
   DNAQualstore store;
   store.readAverageQuality(fname);
   if (par.flip()) {
      KmerCount refkmer;
      refkmer.open(par.getKmerReference());
      cerr << "flipping direction of the store ...\n";
      store.flipDirection(refkmer);
   }
   store.writeFasta(outf);
   // store the quality and count format
   string dqcFile=fname.substr(0, fname.rfind('.')) + ".dqc";
   store.save(dqcFile);
   cout << "redundancy: " << store.getRedundancy() << endl;
   cout << "dqc file: " << dqcFile << endl;
   cout << "fasta output: " << outf << endl;
   return store.getRedundancy();
}

float processFasta(const string &fname, const string &outf, const Progparam &par) {
   cerr << "treating fasta file as fastq file\n";
   DNAQualstore store; // fasta input file
   store.readFasta(fname);
   if (par.flip()) {
      KmerCount refkmer;
      refkmer.open(par.getKmerReference());
      store.flipDirection(refkmer);
   }
   store.writeFasta(outf);
   // store the quality and count format
   string dqcFile=fname.substr(0, fname.rfind('.')) + ".dqc";
   store.save(dqcFile);
   cout << "redundancy: " << store.getRedundancy() << endl;
   cout << "dqc file: " << dqcFile << endl;
   return store.getRedundancy();
}

/** return the total number of sequences read
 */
int readFasta(const string &fname, map<DNA,int> &seqcnt) {
   DNA bs;
   string header;
   ifstream ifs(fname.c_str());
   if (ifs.fail()) {
      cerr << "Failed to open fasta file: " << fname << endl;
      exit(1);
   }
   int cnt=0; // total number of sequences
   map<DNA,int>::iterator mit;
   while (bs.read(ifs, header)) {
      mit = seqcnt.find(bs);
      if (mit == seqcnt.end()) {
         DNA bsrc=bs.revcompCopy();
         mit = seqcnt.find(bsrc);
         if (mit == seqcnt.end()) {
            seqcnt.insert(make_pair(bs,1));
         }
         else {
            ++(mit->second);
         }
      }
      else
         ++(mit->second);
      ++cnt;
   }
   return cnt;
}

void writeFasta(const map<DNA, int> &seqcnt, const string &outfile,
      int total) {
   map<int,int> countfreq;
   ofstream ofs(outfile);
   if (ofs.fail()) {
      throw runtime_error("Failed to open outfile " + outfile);
   }
   for (auto mit=seqcnt.begin(); mit != seqcnt.end(); ++mit) {
      bioseq newseq=mit->first;
      newseq.setName(newseq.getName() + "_" + to_string(mit->second));
      ofs << newseq;
      DNAQualstore::updateCount(countfreq, mit->second);
   }
   string summaryFile = outfile.substr(0, outfile.rfind('.')) + ".summary";
   DNAQualstore::displayCount(countfreq, summaryFile);
   cout << seqcnt.size() << " unique sequences out of " << total
      << " written to " << outfile << "\n";
   cout << "fasta output: " << outfile << endl;
}

// static data for testing
string Progparam::metagData="/home/kzhou/work/metag";
