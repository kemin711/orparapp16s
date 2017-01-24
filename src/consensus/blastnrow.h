#ifndef BLASTNROW_H
#define BLASTNROW_H

// This a basic class to work with blastn results
// (C) 2002 Kemin Zhou at orpara.com

#include <map>
#include <iostream>
#include <iterator>

#include <strformat.h>

using namespace orpara;

/**
 * A utility class to hold blastn result that are 
 * specific to this file only.
 * Add one more taxonomy species name.
 * I am using copy and paste for Blastnrow and Blastnrowd;
 * this is very bad programming method! Should
 * have used inheritance.  This one just one off testing.
 * So live with it.
 */
class Blastnrow {
   public:
      Blastnrow() { }
      /**
       * Constructor use the blastn result row as input.
       * The taxon will be looked up from the tax dictionary.
       */
      Blastnrow(const vector<string> &row);
      /**
       * copy constructgor
       */
      Blastnrow(const Blastnrow &other) = default;
      /*
         : qseqid(other.qseqid), sseqid(other.sseqid), pident(other.pident),
            qlen(other.qlen), slen(other.slen), length(other.length), 
            qstart(other.qstart), qend(other.qend), sstart(other.sstart), send(other.send),  
            bitscore(other.bitscore), taxon(other.taxon),
            queryCount(other.queryCount), consDepth(other.consDepth) { }  
            */

      /** helper used by the constructor
       */
      void setTaxon();
      /**
       * The number of sequence that this query represents
       */
      void setCount(int cnt) { queryCount=cnt; }
      /**
       * Set the consensus's depth
       **/
      void setDepth(int dep) { consDepth=dep; }
      int getConsensusDepth() const { return consDepth; }

      /**
       * Assignment operator
       */
      Blastnrow& operator=(const Blastnrow &other) = default;

      friend ostream& operator<<(ostream& ous, const Blastnrow &r);
      /**
       * Compare by cluster name
       */
      bool operator<(const Blastnrow& r) const;
      int getBitscore() const { return bitscore; }
      bool named() const { return isScientificSpeciesName(taxon.first); }
      /**
       * extract the scientific species name without the strain string.
       */
      string getSpecies() const { return taxon.first; }
      string::size_type getSpeciesLength() const { return getSpecies().length(); }
      /**
       * @return species, strain pair. The later can be empty.
       */
      pair<string,string> getTaxon() const { return taxon; }
      string getMappingName() const {
         if (named()) return getSpecies();
         else return taxon.first + " " + taxon.second;
      }  
      /**
       * use strformat.h function getInt()
       */
      int getQidNumber() const { return getInt(qseqid); }
      /**
       * @return the query name
       */
      const string& getQuery() const { return qseqid; }
      /**
       * @return percent identity of the match.
       */
      float getIdentity() const { return pident; }
      float getSubjcov() const { return float(send-sstart+1)/slen; }
      int getAlignLength() const { return length; }
      int getQueryCount() const { return queryCount; }
      void merge(const Blastnrow &row) { queryCount += row.getQueryCount(); 
         consDepth += row.getConsensusDepth(); }
      /**
       *  @return a subset of the fields for plotting purpose
       */
      tuple<string, int, string, float> getReduced() const {
         return make_tuple(getQuery(), getConsensusDepth(), getSpecies(), getIdentity());
      }
      tuple<string, int, string> getSpeciesDepth() const;
      /**
       * Convert the object to a vector of string
       */
      vector<string> toStringVector();

      /** making a copy
       */
      static void setTaxDictionary(const map<string, pair<string,string> > &id2tax) { 
         taxdic=id2tax; }
      static void writeHeader(ostream &ouf);
      static int getDictionarySize() { return taxdic.size(); }
      static const vector<string>& getHeader() { return header; }

   private:
      /**
       * Query sequence id
       */
      string qseqid;
      /**
       * Subject (target) sequence id
       */
      string sseqid;
      /**
       * Percent identity such as 95.85
       */
      float pident;
      int qlen, slen, length, qstart, qend, sstart, send,  bitscore;
      /**
       * Look up result from Silva taxonomy
       */
      pair<string,string> taxon;
      int queryCount;
      int consDepth;
      /**
       * Taxonomy dictionary from Silva. The dictionary is not loaded
       * by default. The caller needs to load the dictionary.
       */
      static map<string, pair<string,string>> taxdic;
      static vector<string> header;
};

#endif
