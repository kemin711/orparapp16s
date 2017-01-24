#include <iostream>
#include <string>
#include <fstream>
#include <vector>

#include <fastq.h>
#include <bioseq.h>
#include "dnaqualstore.h"

using namespace orpara;

/** 
 * helper to be used only by this object file
 */
void DNAQualstore::displayCount(const map<int,int> &cf, const string &sf) {
   ofstream sout(sf.c_str());
   if (sout.fail()) {
      cerr << "Failed to open summary file: " << sf << endl;
      exit(1);
   }
   sout << "abundance\tcount\n";
   for (auto it=cf.begin(); it != cf.end(); ++it) {
      sout << it->first << '\t' << it->second << endl;
   }
   sout.close();
   cout << "summary abundance written to " << sf << endl;
}

/**
 * helper to be used by this object file only.
 */
void DNAQualstore::updateCount(map<int,int> &cntf, int c) {
   map<int,int>::iterator i = cntf.find(c);
   if (i == cntf.end()) {
      cntf.insert(make_pair(c,1));
   }
   else {
      ++(i->second);
   }
}

// using the first sequence's quality
//
void DNAQualstore::read() {
   ifstream ifs(inputFile.c_str());
   if (ifs.fail()) {
      cerr << "Failed to open fastq file: " << inputFile << endl;
      exit(1);
   }
   totalseq=0;
   map<DNAQual,int> seqcnt;
   map<DNAQual,int>::iterator mit;
   Fastq fq;
   // we will use the quality score of the first sequence
   // to do the averaging of the quality is too expensive
   // we could do it in the future.
   int shortseqCnt = 0;
   while (fq.read(ifs)) {
      if (fq.length() < minseqlen) {
         cerr << "WARN: sequence length too short\n"
            << fq << endl;
         cerr << __FILE__ ": " << __LINE__ << endl;
         ++shortseqCnt;
      }
      DNAQual bs(fq.getName(), fq.getSequence(), fq.getQuality());
      mit = seqcnt.find(bs);
      if (mit == seqcnt.end()) {
         DNAQual bsrc = bs.revcompCopy();
         mit = seqcnt.find(bsrc);
         if (mit == seqcnt.end()) {
            seqcnt.insert(make_pair(bs,1));
         }
         else { ++(mit->second); }
      }
      else ++(mit->second);
      ++totalseq;
   }
   if (shortseqCnt > 0) {
      cerr << shortseqCnt << " short sequences!\n";
   }
   // you cannot update elements in a set so we have to 
   // go a long way
   seqs.clear();
   for (mit=seqcnt.begin(); mit != seqcnt.end(); ++mit) {
      seqs.insert(make_pair(mit->first.getName(), DNAQualCount(mit->first, mit->second)));
   }
}

void DNAQualstore::readFasta(const string &inf) {
   cerr << "DNAQualstore reading fasta " << inf 
      << " as if it is fastq\n";
   setInputFile(inf);
   ifstream ifs(inputFile);
   if (ifs.fail()) {
      throw runtime_error("Failed to open fastq file: " + inputFile);
   }
   totalseq=0;
   map<DNAQual,int> seqcnt;
   map<DNAQual,int>::iterator mit;
   DNA dna;
   // we will use the quality score of the first sequence
   // to do the averaging of the quality is too expensive
   // we could do it in the future.
   while (dna.read(ifs)) {
      DNAQual bs(dna);
      mit = seqcnt.find(bs);
      if (mit == seqcnt.end()) {
         DNAQual bsrc = bs.revcompCopy();
         mit = seqcnt.find(bsrc);
         if (mit == seqcnt.end()) {
            seqcnt.insert(make_pair(bs,1));
         }
         else { ++(mit->second); }
      }
      else ++(mit->second);
      ++totalseq;
   }
   cerr << totalseq << " in fasta file: " << inf << endl;
   // you cannot update elements in a set so we have to 
   // go a long way
   seqs.clear();
   for (mit=seqcnt.begin(); mit != seqcnt.end(); ++mit) {
      seqs.insert(make_pair(mit->first.getName(), DNAQualCount(mit->first, mit->second)));
   }
   cerr << seqs.size() << " unique names for sequences\n";
}

// add length filter to prevent other components
// from giving up
void DNAQualstore::readAverageQuality() {
   ifstream ifs(inputFile.c_str());
   if (ifs.fail()) {
      cerr << "Failed to open fastq file: " << inputFile << endl;
      exit(1);
   }
   totalseq=0;
   map<DNA, QualArray> seqcnt;
   map<DNA, QualArray>::iterator mit;
   Fastq fq;
   int shortseqCnt=0;
   while (fq.read(ifs)) {
      //cerr << "readAverageQuality() fastq sequence length: " << fq.length() << endl;
      if (fq.length() < minseqlen) {
         cerr << "WARN: sequence length too short\n"
            << fq << endl;
         cerr << __FILE__ ": " << __LINE__ << endl;
         ++shortseqCnt;
         continue;
      }
      DNA bs(fq.getName(), fq.getSequence());
      int *bsq = fq.getQuality();
      mit = seqcnt.find(bs);
      if (mit == seqcnt.end()) {
         DNA bsrc = bs.revcompCopy();
         mit = seqcnt.find(bsrc);
         if (mit == seqcnt.end()) {
            QualArray tmpq(bs.length());
            tmpq.add(bsq);
            seqcnt[bs] = std::move(tmpq);
         }
         else { 
            mit->second.addReverse(bsq);
         }
      }
      else {
         mit->second.add(bsq);
      }
      ++totalseq;
   }
   if (shortseqCnt > 0) {
      cerr << shortseqCnt << " short sequences!\n";
   }
   seqs.clear();
   // the key for the map is not modifiable!
   // you must make a copy
   for (mit=seqcnt.begin(); mit != seqcnt.end(); ++mit) {
      seqs.insert(make_pair(mit->first.getName(), 
            DNAQualCount(std::move(mit->first), 
            mit->second.getCount(), mit->second.getMean())));
   }
}

void DNAQualstore::writeFasta(const string &outfile) const {
   try {
      ofstream ofs(outfile.c_str());
      map<int,int> countfreq;
      for (auto it=seqs.begin(); it != seqs.end(); ++it) {
         ofs << '>' << it->second.getName() << '_' << it->second.getCount() << endl;
         it->second.printFasta(ofs);
         updateCount(countfreq, it->second.getCount());
      }
      string summaryFile = outfile.substr(0, outfile.rfind('.')) + ".abundant.summary";
      displayCount(countfreq, summaryFile);
   }
   catch (system_error e) { // for older version need to include <system_rror>
      cerr << e.what() << " Failed to write Fasta outfile " << outfile << endl;
   }
   cout << seqs.size() << " unique sequences out of " << totalseq 
      << " written to " << outfile << "\n";
}

void DNAQualstore::flipDirection(const KmerCount &refk) {
   map<string, DNAQualCount>::iterator it = seqs.begin();
   map<string, DNAQualCount>::iterator del;
   pair<int,int> loopreg;
   while (it != seqs.end()) {
      KmerCount kc(refk.getK());
      //cerr << __FILE__ << " working on " << it->first << endl;
      // shold be done outside this function
      //Kmert<6> kmer6(it->second.toString());
      //if (kmer6.isPalindrome(loopreg)) {
      //   cout << it->first << " is palindrome\n";
      //}
      //cerr << "working on " << it->first << endl;
      kc(it->second.getSequence());
      if (!refk.sameDirection(kc)) {
         del=it;
         ++it;
         DNAQualCount tmp=del->second.revcompCopy();
         seqs[tmp.getName()]=tmp;
         seqs.erase(del);
      }
      else ++it;
   }
}

vector<const DNAQualCount*> DNAQualstore::getSequences(const vector<string> &ids) const {
   vector<const DNAQualCount*> tmp(ids.size());
   map<string, DNAQualCount>::const_iterator mit;
   for (int i=0; i<ids.size(); ++i) {
      mit = seqs.find(ids[i]);
      if (mit == seqs.end()) {
         cerr << "input file: " << inputFile << endl;
         throw runtime_error("getSequences() could not find sequence: " + ids[i]);
      }
      tmp[i]=&(mit->second);
   }
   return tmp;
}

list<DNAQualCount*> DNAQualstore::getSequencesList(const vector<string> &ids) {
   list<DNAQualCount*> tmp;
   map<string, DNAQualCount>::iterator mit;
   for (int i=0; i<ids.size(); ++i) {
      mit = seqs.find(ids[i]);
      if (mit == seqs.end()) {
         cerr << "input file: " << inputFile << " " << seqs.size()
            << " sequences " << endl;
         throw runtime_error("getSequencesList could not find sequence: " + ids[i]
               + "|");
      }
      tmp.push_back(&(mit->second));
   }
   return tmp;
}

vector<string> DNAQualstore::stripId(const vector<string> &ids) {
   vector<string> tmp(ids.size());
   for (int i=0; i<ids.size(); ++i) {
      tmp[i]=ids[i].substr(0, ids[i].rfind('_'));
   }
   return tmp;
}

DNAQualCount*  DNAQualstore::addSequence(DNAQualCount &&s) { 
   totalseq += s.getCount(); 
   pair<map<string, DNAQualCount>::iterator,bool> rv = seqs.insert(make_pair(s.getName(), std::move(s))); 
   if (!rv.second) {
      cerr << rv.first->first << " already exists in sequence store!\n";
   }
   return &((rv.first)->second);
}

DNAQualCount*  DNAQualstore::addSequence(const DNAQualCount &s) { 
   totalseq += s.getCount(); 
   pair<map<string, DNAQualCount>::iterator,bool> rv = seqs.insert(make_pair(s.getName(), s)); 
   if (!rv.second) {
      cerr << rv.first->first << " already exists in sequence store!\n";
   }
   return &((rv.first)->second);
}

void DNAQualstore::open(const string &storeFile) {
   cerr << "opening dqc file: " << storeFile << endl;
   ifstream inf(storeFile);
   inputFile=storeFile;
   if (inf.fail()) {
      throw runtime_error("Failed to open DNAQualCount file: " + storeFile);
   }
   DNAQualCount tmp;
   inf >> tmp;
   totalseq=0;
   while (!inf.eof())  {
      totalseq += tmp.getCount();
      seqs[tmp.getName()] = tmp;
      //cout << "|" << tmp.getName() << "|" << endl;
      inf >> tmp;
   }
   seqs[tmp.getName()] = tmp;
   totalseq += tmp.getCount();
   //cerr << "last sequence: " << tmp.getName() << endl;
   cerr << totalseq << " total sequences read into memory "
      << getNumberUnique() << " Unique sequences\n";
}

void DNAQualstore::save(const string &storeFile) {
   ofstream ouf(storeFile);
   if (ouf.fail()) {
      throw runtime_error("Failed to open file " + storeFile + " for writing");
   }
   for (auto it=seqs.begin(); it != seqs.end(); ++it) {
      ouf << it->second;
   }
   cerr << seqs.size() << " unique sequences out of " << totalseq 
      << " written to dqc format file: " << storeFile << "\n";
}

void DNAQualstore::discardShort(int lencut) {
   map<string, DNAQualCount>::iterator del;
   map<string, DNAQualCount>::iterator mit=seqs.begin();
   while (mit != seqs.end()) {
      if (mit->second.length() < lencut) {
         totalseq -= mit->second.getCount();
         del = mit;
         ++mit;
         seqs.erase(del);
      }
      else ++mit;
   }
}

