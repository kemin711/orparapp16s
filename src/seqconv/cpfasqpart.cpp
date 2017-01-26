// (C) 2012 Kemin Zhou at orapra.com
// for random and fixed sampling of fastq sequences
// from a file
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <string>
#include <random>
#include <vector>
#include <cstring>
#include <chrono>
#include <stdexcept>

#include <strformat.h>
#include <fastq.h>

using namespace std;
using namespace orpara;

void usage() {
   cout << "cpfasqpart -p 0.5 -o outfile inputfile\n"
      << "Will pick the first fraction of sequences from input file\n"
      << "Options\n"
      << "   -p <fraction> a floating number < 1. The copy is done\n"
      << "      by picking the fraction of the total file. No random smple\n"
      << "      copy starts from the begining. This is a lot faster\n"
      << "      than the random sample method\n"
      << "   -s or --sample '123 999 777' random sample without replacement\n"
      << "      The numbers are the number of samples take from the input file\n"
      << "      They must be either single or double quoted\n"
      << "      Each sample is independent of the other. If you want two\n"
      << "      files with 50 and 100 sequences you will do the following\n"
      << "      cpfastqpart -s '50 100' myfastqfile\n"
      << "  -o FILE output file name. This option can be combined with\n"
      << "      the -p option or the -s option with single numbers\n"
      << "      If multiple numbers are specified with the -s option\n"
      << "      then, the output file names are automatically generated\n"
      << "  -a FILE output file name. This option is the same as -o\n"
      << "      except that the result will be appened to the output file\n"
      << "the output file holds the result.\n";
   exit(1);
}

/**
 * Randomly pick a fraction of reads without replacement.
 */
void randomPick(const string &infile, float frac);
double GetUniform();
vector<Fastq> buildStore(const string& infile);
void pickSample(const vector<Fastq> &pop, int numseq, const string &outfile,
      ios_base::openmode mode = ios_base::out);
/**
 * Repeat picking sequence without replacement.
 * @param sampleSize is the samples you need.
 *    The number in the vector specifies the number of
 *    sequences for each drawing.
 * @inputFile is the input file.
 */
void repeatPick(const string &inputFile, const vector<int> &samplesz);
/**
 * pick numseq from inputFile and write into otfile
 *
 * @param append if set true will append to output file.
 */
void singlePick(const string &inputFile, int numseq, const string &outfile, bool append=false);

vector<int> SampleWithoutReplacement(int populationSize, int sampleSize);

/**
 * copy a fraction of the fastq file into a new file
 * None random yet. Needs to implement a random version.
 */
int main(int argc, char* argv[]) {
   if (argc < 2) {
      usage();
   }
   int i=1;
   string infile, outfile;
   float fraction = 0.50; // percent
   vector<int> sampleSize;
   bool append = false;
   while (i < argc) {
      if (string(argv[i]) == "-p") {
         fraction = atof(argv[++i]);
      }
      else if (string(argv[i]) == "-o") {
         outfile = argv[++i];
      }
      else if (string(argv[i]) == "-a") {
         outfile = argv[++i];
         append = true;
      }
      else if (!strcmp(argv[i], "--sample") || !strcmp(argv[i], "-s")) {
         sampleSize=getAllInt(argv[++i]);
      }
      else if (!strcmp(argv[i], "--help") || !strcmp(argv[i], "-h")
            || argv[i][0]=='?') {
         usage();
      }
      else {
         if (argv[i][0] == '-') {
            cerr << "Unsupported option: " << argv[i] << endl;
            usage();
         }
         infile = argv[i];
      }
      ++i;
   }
   if (!sampleSize.empty()) {
      if (sampleSize.size() > 1) {
         repeatPick(infile, sampleSize);
         cout << sampleSize.size() << " randome sampling without replacement done\n";
      }
      else {
         singlePick(infile, sampleSize[0], outfile, append);
      }
      return 0;
   }

   long begin, end;
   ifstream inf(infile.c_str());
   begin = inf.tellg();
   inf.seekg(0, ios::end);
   end = inf.tellg();
   long mid = long(fraction*(begin + end));
   inf.seekg(mid);
   string line;
   getline(inf, line);
   while (!inf.eof() && line != "+") {
      //cout << line << endl;
      getline(inf, line);
   }
   getline(inf, line);
   // line should be the last line now
   string endline=line;
   long cutpos = inf.tellg();
   inf.seekg(0, ios::beg);
   char buffer[2000]; // = new char[1000]; // next genseq cannot read >500 nt
   // still safe
   ofstream ouf(outfile.c_str());
   inf.read(buffer, 2000);
   while (!inf.eof() && inf.tellg() < cutpos-2000) {
      ouf.write(buffer, 2000);
      inf.read(buffer, 2000);
   }
   ouf.write(buffer, 2000);
   getline(inf, line);
   while (inf.tellg() < cutpos && line != endline) {
      ouf << line << endl;
      getline(inf, line);
   }
   ouf << line << endl;
   //cout << "last line " << line << endl;
   //delete [] buffer;
   return 0;
}

// the order of the pick is still deterministic
void repeatPick(const string &inputFile, const vector<int> &samplesz) {
   vector<Fastq> store=buildStore(inputFile);
   for (unsigned int i=0; i<samplesz.size(); ++i) {
      string outfname = inputFile;
      outfname.insert(inputFile.rfind('.'), "_" + to_string(samplesz[i]));
      pickSample(store, samplesz[i], outfname);
   }
}

void singlePick(const string &inputFile, int numseq, const string &outfile, bool append) {
   vector<Fastq> store=buildStore(inputFile);
   if (append) {
      pickSample(store, numseq, outfile, ios_base::app);
   }
   else {
      pickSample(store, numseq, outfile);
   }
}

vector<Fastq> buildStore(const string& infile) {
   // read into a vector
   ifstream inf(infile);
   if (inf.fail()) {
      throw runtime_error("failed to open file: " + infile + " for write");
   }
   vector<Fastq> store;
   Fastq fastq;
   while (fastq.read(inf)) {
      store.push_back(fastq);
   }
   return store;
}

/**
 * sample numseq sequences from fastq input pop
 *
 * @param pop input population
 * @param numseq number of sequences to draw from input
 * @param outfile the output file name to write to.
 */
void pickSample(const vector<Fastq> &pop, int numseq, const string &outfile, 
      ios_base::openmode mode) {
   vector<int> sample = SampleWithoutReplacement(pop.size(), numseq);
   ofstream ouf(outfile, mode);
   if (ouf.fail()) {
      throw runtime_error("failed to open " + outfile);
   }
   for (size_t i = 0; i < sample.size(); i++ ) {
      //cout << sample[i] << "\t";
      ouf << pop[sample[i]];
   }
   //cout << endl;
   ouf.close();
   cout << numseq << " written to " << outfile << endl;
}

double GetUniform()
{
   std::random_device rd;
   std::mt19937 gen(rd());
   //static std::default_random_engine re;
   static std::uniform_real_distribution<double> Dist(0,1);
   //return Dist(re);
   return Dist(gen);
}

// John D. Cook, http://stackoverflow.com/a/311716/15485
/**
 * @param populationSize size of set sampling from
 * @param sampleSize size of each sample
 * @param samples output, zero-offset indicies to selected items
 */
vector<int> SampleWithoutReplacement(int populationSize, int sampleSize)
{ // Use Knuth's variable names
   vector<int> samples(sampleSize);
   int& n = sampleSize;
   int& N = populationSize;
   int t = 0; // total input records dealt with
   int m = 0; // number of items selected so far
   double u;
   while (m < n) {
      u = GetUniform(); // call a uniform(0,1) random number generator
      if ( (N - t)*u >= n - m ) {
         t++;
      } 
      else {
         samples[m] = t;
         t++; m++;
      }
   }
   return samples;
}

