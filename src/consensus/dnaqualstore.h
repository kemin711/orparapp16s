#ifndef DNAQUALSTORE_H
#define DNAQUALSTORE_H

// (C) 2012 Kemin Zhou at orpara.com
// This is a simple soreage calss for
// working with fastq sequences

#include <string>
#include <vector>
#include <system_error>
#include <list>

#include <bioseq.h>
#include <stddev.h>
#include <kmer.h>

using namespace std;
using namespace orpara;

class QualArray {
   public:
      int count;
      vector<stddev> sarr;

      QualArray() : count(0), sarr() { }
      QualArray(int size)
         : count(0), sarr(size) { }
      void add(const int* q) {
         ++count;
         for (size_t i=0; i<sarr.size(); ++i) {
            sarr[i](q[i]);
         }
      }
      void addReverse(const int* q) {
         ++count;
         for (size_t i=0; i<sarr.size(); ++i) {
            sarr[i](q[sarr.size()-i-1]); 
         }
      }
      int getCount() const { return count; }
      vector<int> getMean() const {
         vector<int> tmp(sarr.size());
         for (size_t i=0; i<sarr.size(); ++i) {
            tmp[i]=int(ceil(sarr[i].getMean()));
         }
         return tmp;
      }
};

/**
 * A class to hold fastq files in memory.
 * The data is represented in DNAQualCount format
 * so that they can be used for sequence alignment.
 * Not intended to be copied around.  This is a 
 * heave store.
 */
class DNAQualstore {
   private:
      /**
       * A fastq file with redundant sequences.
       */
      string inputFile;
      /** most frequenct use is index by id
       * using string id, slower. For better performance
       * integer id should be used.
       */
      map<string, DNAQualCount> seqs;
      /**
       * total number of sequences.
       * A faster short cut. You can calculate this
       * from the raw data. It is just takes a little
       * bit of time.
       */
      int totalseq;
      static const int minseqlen=20;

   public:
      /**
       * default constructor.
       * Start an empty store.
       */
      DNAQualstore() : inputFile(), seqs(), totalseq(0) { }
      /** will initialize the object by reading
       * the content in the file into contanner.
       */
      DNAQualstore(const string &infile) 
         : inputFile(infile), seqs(), totalseq(0) { read(); }
      DNAQualstore(const DNAQualstore &store) = delete;
      DNAQualstore& operator=(const DNAQualstore &store)=delete;

      /** 
       * read from the raw inputFile that should be a 
       * fastq formated file.
       * The quality score of multiple identical sequences will 
       * use the first ones. The score is not very important, this
       * method will be faster compared to readAverageQuality()
       * that will average the quality from all identical sequences.
       */ 
      void read();
      /**
       * Read from fasta file that does not have quality info.
       * And treat fasta file as if the sequences have quality
       * score.  This is for stisfaying formating requirement 
       * of fastq files so that fasta files can be used as
       * input.
       */
      void readFasta(const string &inf);
      /**
       * average the quality from all reads.
       * Expensive operation.
       */
      void readAverageQuality();
      /**
       * Read the fastq file and average the quality
       * of identical sequences. 
       */
      void readAverageQuality(const string &fastqf) {
         setInputFile(fastqf);
         readAverageQuality();
      }

      /**
       * Flip sequences in the reverse direction
       * according to a KmerCount object.
       * @param refk should be supplied by the caller.
       */
      void flipDirection(const KmerCount &refk);
      /**
       * remove sequences shorter than lencut.
       */
      void discardShort(int lencut=200);

      /**
       * Save the store in DNAQualCount format.
       * with file extension dqc.
       * Use open to load the dqc file.
       */
      void save(const string &storeFile);
      /**
       * Open from the store of DNAQualCount format.
       * File extension is usually dqc. File is 
       * saved with the save() function.
       */
      void open(const string &storeFile);
      void setInputFile(const string &input) { inputFile=input; }
      /**
       * @return the input file used to generate this store.
       */
      const string& getInputFile() const { return inputFile; }
      /**
       * write a fasta file that can be used as input
       * for swarm program.
       * It will also write a file with *.summary
       * that conatains the summary of different abundances.
       */
      void writeFasta(const string &outfile) const;
      /**
       * Give a list of ids it will return a vector of pointers
       * to objects in this store.
       * Make sure this store should be in the scope.
       */
      vector<const DNAQualCount*> getSequences(const vector<string> &ids) const;
      /**
       * @return a list type for fast editing operations.
       */
      list<DNAQualCount*> getSequencesList(const vector<string> &ids);
      /**
       * move s into this store.
       * @return a pointer to the const DNAQualConst object.
       */
      DNAQualCount* addSequence(DNAQualCount &&s); 
      DNAQualCount* addSequence(const DNAQualCount &s); 
      size_t getTotalSequences() const { return totalseq; }
      float getRedundancy() const { return float(getTotalSequences())/getNumberUnique(); }
      /**
       * @return number of unique sequences.
       */
      size_t getNumberUnique() const { return seqs.size(); }
      
      /**
       * helper function to remove the ID_123 number for count part
       * The number after the under score is used by swarm. This should
       * be removed in the future when swarm is replaced with a good
       * library.
       */
      static vector<string> stripId(const vector<string> &ids);
      /**
       * Helper for playing with abundance => frequency (count)
       */
      static void updateCount(map<int,int> &cntf, int c);
      static void displayCount(const map<int,int> &cf, const string &sf);
};

#endif
