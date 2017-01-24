#ifndef CONSENSUS_H
#define CONSENSUS_H

// (C) 2012 orpra.com Kemin Zhou orpara.com
// last modified 2016.  This class started with basic 
// consensus, then modified to add more features.
// After sequences are aligned to a representative member,
// This calss simply them to a column-summary of bases
// This class is adopted to use Nucleic acids, in the future
// abstraction should be done to work with protein sequences.
// Also imporved performance with multi-threading

#include <string>
#include <vector>
#include <map>
#include <functional>
#include <cmath>
#include <iostream>
#include <fstream>
#include <thread>
#include <ios>
#include <list>
#include <tuple>
#include <array>
#include <mutex>

#include <bioseq.h>
#include <dynalnt.h>
#include <alninfo.h>
#include <fastq.h>
#include <stddev.h>
#include "dnaqualstore.h"


//#define DEBUG

using namespace std;
using namespace orpara;

float roundPercent(float frac);

class ConsensusSizeException : public invalid_argument {
   public:
      ConsensusSizeException(const string &msg) 
         : invalid_argument(msg) { 
      }
};

template<class T>
class SortPairByInt {
   public:
      bool operator()(const pair<T, int> &p1, const pair<T, int> &p2) {
         return p1.second > p2.second;
      }
};

/**
 * A simple structure to rank alignment quality.
 * This will be used for filtering.
 * This will be used to decide for poor singletons
 * which cluster they belong. 
 */
class AlignQuality {
   private:
      /** at this point we use the repseq to represnt
       * the consensus. In the future we may do 
       * another round of consensus building by using
       * the consensus from the first round as reference.
       * This can achive a tiny bit of improvement.
       */
      string cons;
      /** identical and alnlen decides identity
       */
      int identical;
      int alnlen;
      float repcov;
      float readcov;

   public:
      AlignQuality() : cons(), identical(0), alnlen(0), repcov(0), 
            readcov(0) { }
      AlignQuality(const string &rep, int iden, int len, 
            float covrep, float covrd) 
         : cons(rep), identical(iden), alnlen(len), repcov(covrep), 
            readcov(covrd) { }
      float getIdentity() const { return float(identical)/alnlen; }
      /**
       * order by idental
       */
      bool operator<(const AlignQuality &o) const;
      /**
       * After selecting the best cluster hit, we can lower
       * our cutoff to the following.
       * (1) identity > 0.99 and alnlen > 200 we will pass QC
       * (2) coverage of either one > 0.9 and identity > 0.85
       * The user should adjust this identity cut
       */
      bool passCutoff(float identitycut=0.85) const;
      friend ostream& operator<<(ostream &ous, const AlignQuality &aq) {
         ous << aq.cons << " identical=" << aq.identical 
            << " alnlen=" << aq.alnlen << " repcov=" << aq.repcov 
            << " readcov=" << aq.readcov << " identity="
            << aq.getIdentity();
         return ous;
      }
      const string& getConsensusName() const { return cons; }
};

/**
 * This should be a light weight object using reference to
 * external resources.
 * We will use a DNAQualCount store to work with.
 * This object is movable but not copyable because of its
 * large size.
 */
class Consensus {
   protected:
      ///// global parameters ///////////////
      /** you only need one store. Intialized to 0.
       * Someone has to give it some content before using the
       * consensu tool. 
       * We may need to revcomp some sequences. So can cannot
       * use a const qualifier.
       * For performance, the sequence in the store 
       * should be in the same direction.
       * */
      static DNAQualstore* seqstore;
      /** 
       * In only need a copy but could also have one 
       * per aligner. This will reduce some cost for
       * threaded operations.
       */
      static SimpleScoreMethod scoreMethod; //(10, -9, -20, -3);
      /**
       * This is the machine used in this class.
       * The use of matrix will slow down a little bit.
       * For threading, this should be an array.
       * Externally supplied machines.
       * Must be provided from the outside!
       */
      //Dynaln<SimpleScoreMethod> aligner;
      static vector<Dynaln<SimpleScoreMethod>* > aligner;
      /**
       * Identity cutoff for identity of the alignment.
       * Default is 0.87, the accuracy of the sequencer.
       */
      static double identityCut;
      /**
       * cutoff value for alignment on the reference sequence
       * below this value, algiments will not be counted.
       * Default is 200 for Iontorrent platforms where 400 nt reads
       * are generated. For Illumina
       * this value should be set by the user to about 1/3 of the
       * read length.
       */
      static int alnlengthCut;
      /**
       * For merging clusters, where consensus is used.
       * Consensus usually does not have identifier, so
       * we make up some.
       */
      static int newseqid;
      static int clusterid;

      ///////////// nonstatic section ////////////////////
      /**
       * member of the cluster excluding the repseq.
       * repseq is the seed sequence for building 
       * consensus. The seed sequence is usually
       * the one with the highest counts.
       */
      mutable list<DNAQualCount*> member;
      /** 
       * sum of counts from all member + 1:repseq
       * For efficiency, this number does not have to be
       * recalculated.
       * Both acess and update will be protected
       * by internal mutex
       */
      mutable int totalseq;
      /**
       * representative sequence. The seed sequence
       * for building consensus.
       * This is also the reference sequence for the
       * first round of consensus building.
       */
      mutable DNAQualCount* repseq;

      ////////// computation results ////////////////
      ///// intermediate data structures ///////
      /**
       * For quality control. Not used much right now.
       */
      map<AlignInfo, int> alnsummary;
      /**
       * Single base count. Each position is a A->100, C->10, G->5, T->1
       */
      vector<map<char, int> > bases;
      /**
       * Insertion of bases one after the base relative to the
       * repsequence if the cluster has not been refined.
       * After refinement, the consensus will be used as
       * reference.
       */
      vector<map<string, int> > inserts;
      /**
       * Fastq quality score average, std at each position.
       * This information is not very useful. We keep them
       * for future usage.
       */
      vector<stddev> quality;

      ////// final results ////////////////
      mutable string consensus;
      mutable vector<int> consensusQ;
      mutable bool consensusValid;
      /**
       * Get a uniform unique name for the consensus
       */
      mutable string consensusName;
      /**
       * The base vertical coverage for each consensus position
       * We record 3 numbers:
       * 1 the count of the consensus base
       * 2 the total bases without gap char
       * 3 total base including gap char
       * When there is a insert relative to the reference
       * the total base is taken from the previous position for the
       * total base with gap char.
       */
      mutable vector<tuple<int,int,int> > coverage;
      /**
       * mean,std for each base coverage level
       * over the length of the consensus sequence.
       * This is a global value for this cluster.
       * 1 consensus 2 bases without gap 3 bases including gap
       */
      mutable array<stddev,3> averageCoverage;
      void clearAverage() const {
         for (int i=0; i<3; ++i) averageCoverage[i].clear(); }
      /**
       * for accumulate, updating inserts, bases, quality
       */
      //mutable mutex mutAccumulate;
      /**
       * For synchronizing resources inside the class
       */
      mutable mutex mut;
      /**
       * For synchronizing external resources such as 
       * writing to external files or streams.
       */
      mutable mutex mutExternal;


      ///////// helper functions ////////////
      /**
       * Align the input to reference and see it 
       * has good alignment quality. If forward alignment
       * is not good, it will try the reverse.
       * Only if the match between the seed of the cluster
       * and s passed certain cutoff, then 
       * accumulate is called to integrate s into
       * the cluster.
       * @param s input sequence. It may be flipped 
       *   to the reverse direction.
       */
      bool consume(unsigned alid, DNAQualCount *s);
      /**
       * @param b start position of the block.
       * @param e end position of the block.
       */
      void consumeBlock(unsigned alignerId, list<DNAQualCount*> &chk);
      /**
       * Requires that the alignment has been done between
       * the input and the repseq (seed).
       * This function update the bases, quality,
       * totalseq, and alnsummary.
       * @param alid aligner id this is determined at
       *    the initialize() step.
       * @param read pointer to a const DNAQualCont object.
       *    The read is used as input for getting the count,
       *    quality only.
       * NOTE: this function will NOT update member
       * because this function is shared by callers
       * integrating internal members and external unknown sequence
       * into the member container. For the operator()
       * the members are already in the member containers
       * so you don't have to add.
       */
      void accumulate(unsigned alid, DNAQualCount* read);

      /**
       * Picke the sequences with the highest Count.
       * If there is a tie, then pick the one with the 
       * most frequenct length.
       * After that pick the longer sequence.
       * The repseq pointer will be pointing to this sequence.
       * The member size will be reduced by one.
       */
      void pickRepseq() const;
      /**
       * bring the object from intial stage to action ready state.
       * 1. Pick repsequence
       * 2. Set up size of containner: bases, inserts, quality
       * 3. Put reqseq into bases
       * 4. Set seq1 of aligner to repseq.
       * Condition: Assume that the members has been intialized by 
       *   constructors.
       */
      void initialize();

      //// static method section ////////////////////
      /**
       * Helper only used by this class
       * Most gaps in PacBio should not be closed.
       */
      static void closeStagerGap(string &top, string &bottom);
      static void scanStagerGap(string &top, string &bottom);

   public:
      /**
       * default constructor.
       * Will set the reference sequence from the default container.
       * Number of aligners to use will be deermined
       * later. default two threads.
       */
      Consensus()  
         : member(), totalseq(0), repseq(0),
            alnsummary(), bases(), inserts(),
            quality(), consensus(), consensusQ(),
            consensusValid(false), consensusName(),
            coverage(), mut(), mutExternal() { }
      /**
       * @param seqnames a list of sequence names to start the consensus.
       *   Or sequence ids.
       */
      Consensus(const vector<string> &seqnames) 
         : member(seqstore->getSequencesList(seqnames)), totalseq(0), repseq(0),
           alnsummary(), bases(), inserts(),
           quality(), consensus(), consensusQ(), 
           consensusValid(false), consensusName(), coverage(),
           mut(), mutExternal()
      { initialize(); }
      /** disable copying */
      Consensus(const Consensus &other)=delete;
      Consensus& operator=(const Consensus &other)=delete;
      /**
       * the zero assignment to reseq is not even needed. Just for safety.
       * Movable but not copyable.
       */
      Consensus(Consensus&& o) 
         : member(std::move(o.member)), totalseq(o.totalseq),
         repseq(o.repseq), alnsummary(std::move(o.alnsummary)),
         bases(std::move(o.bases)), inserts(std::move(o.inserts)),
         quality(std::move(o.quality)), consensus(std::move(o.consensus)),
         consensusQ(std::move(o.consensusQ)),
         consensusValid(o.consensusValid), consensusName(move(o.consensusName)),
         coverage(move(o.coverage)), mut(), mutExternal()
      { o.repseq=0; }
      /** 
       * will not build the consensus
       */
      Consensus& operator=(Consensus &&o);
      ~Consensus() { }

      /** 
       * this is where the actual work is done 
       * iterate over all members and align to refsseq
       * It will no update totalseq, this will be calculated
       * when consensus is being built.
       */
      void operator()();

      /**
       * replace the seed sequence with consensus then update
       * bases data structure.
       * And produce variation positions
       * After refine, repseq will nolonger be used.
       */
      void refine();
      /**
       * Alignment is good enought to be used.
       * Very stragent use 0.99 identity cutoff.
       * @param ai aligner index in the aligner pool.
       */
      bool alignPassQC(unsigned ai) const;

      /**
       * @return a copy of the repseq. If the getBestHaplotype has been called
       *   The returned repseq will the best match from the consensus to the 
       *   haplotype.
       */
      DNAQualCount* getRepseq() const {
         return repseq;
      }
      /**
       * @return the length of the reference sequence.
       */
      size_t getRepseqLength() const { return repseq->length(); }
      /** 
       * take the bases variable and generate a consens
       * string and the quality. So far not discarding the 
       * long tails yet.
       * Right now every base is reported. Future version
       * should consider the average vertical coverage, 
       * if head and tail falls below the cutoff then 
       * these segments should be excluded from the final
       * consensus.
       * This method will also produced a name for the 
       * consensus.
       */
      void buildConsensus() const;
      /**
       * When adding new sequences to the consensus,
       * the consensus needs be invalidate with this 
       * function.
       */
      void invalidateConsensus() { consensusValid=false; }
      bool isConsensusValid() const { return consensusValid; }
      /**
       * @return the consensus string.
       */
      const string& getConsensus() const { 
         if (consensus.empty() || !consensusValid) buildConsensus();
         return consensus; 
      }
      const string& getConsensusName() const {
         return consensusName; 
      }
      void setConsensusName(const string &name) { consensusName=name; }
      /**
       * @param name you need to give a name to this new sequence.
       * @return final result with quality from the average quality score.
       *       Count from the totalseq.
       */
      DNAQualCount getResult(const string &name) const { return DNAQualCount(name, getConsensus(), consensusQ, getTotal()); }
      /**
       * version using the internal integer counter to generate names
       * The name of the sequence will be determined by the 
       * getClusterName() class method.
       */
      DNAQualCount getResult() const;
      /**
       * Only cares about the average.
       */
      vector<int> getMeanQuality() const;
      /**
       * Clear all fields for reuse.
       */
      void clear();

      //void run();

      /**
       * Output finctions.
       */
      void printBases(ostream &ous) const;
      /**
       * Print for human comsumption.
       * >ConsensusName [memberId separated by comma]
       * consensus string one line
       * return a sequence object.
       */
      void printConsensus(ostream &ous) const;
      void printQuality(ostream &ous) const;
      void printAlninfo(ostream &ous) const;
      /**
       * print the consensus in fasta format
       * with copy number in the header
       */
      void printFasta(ostream &ous) const {
         getResult().printFastaWithHeader(ous);
      }
      void printCoverage(ostream &ous) const;
      /**
       * return the stddev object for three
       * levels of the coverqage.
       * TODO: add this value to the consensus fasta file
       * and consensus text format.
       */
      array<stddev,3>  getAverageCoverage() const;
      /**
       * @return The base coverage including gaps.
       */
      double getGrossBaseCoverage() const {
         return getAverageCoverage()[2].getMean(); }

      /**
       * print the result of this algorithm
       *  This include: genotype, codons frequency, and consensus sequence.
       *  If not calculating linked mutations frequencies, then
       *  genofile will be empty.
       *  @param codonfiles output file name for codon
       */
      void printResult(ostream &oucons, ostream &oubase, 
            ostream &ouqual, ostream &oualni) const;
      /**
       * overloaded version also prints base vertical coverage
       * This is used for outputing result from one sample
       * with multiple consensus results.
       */
      void printResult(ostream &oucons, ostream &oubase, 
            ostream &ouqual, ostream &oualni, ostream &oucov) const;
      void printResult(const string &consensusfile, const string &basefile, 
            const string &qualfile, ios_base::openmode writeMode=ios::out|ios::app) const;

      /**
      * Global variable use by tabulateCodon for efficiency
      * Start 0-based index of the codon of interest in the reference seq.
      * This could be changed before each analysis
      */

      /**
       * Helper method to convert map to sorted vector by count
       * This method could be used by the public.
       */
      //static vector<pair<string, int> >  map2vectorByCount(const map<string, int>& src); 
#ifdef DEBUG
      void setAlignOutputFile(const string &fname);
#endif
      void printAlninfo(const string &alninfoFile) const;
      //void trimAlignmentTail(string &top, string &bottom);
      /**
       * Number of unique sequences.
       */
      int numseq() const { return member.size()+1; }
      /**
       * @return total number of sequences in this cluster.
       * Also reset totalseq to the sum.
       */
      int computeTotalSequences() const;
      /**
       * Total number of sequences.
       * Each unique sequence may represent multiple members.
       */
      int getTotal() const;
      /**
       * To make sure totalseq is consistent with 
       * the freshly computed value.
       */
      bool debugValidTotal();
      /**
       * Merge this with another Consensus cons.
       * This method makes sure that cons is smaller.
       * If not it will reject with error thrown.
       * It simply use the consensus of the input cluster
       * as input sequence and use the mergeSequence() method
       * to merge the input consensus into the current cluster.
       * @param chk list of input consensi to be merged with this one.
       *    They should have smaller size than this one.
       *    The container is not altered by this function.
       * @param scrap collect the segments from cons that
       *    does not match the current consensus. Input
       *    consensus may contain this consensus or is
       *    a chimera that contains a segment of this consensus.
       * @param toErase record consensus that have been
       *   merged with this consensus. 
       */
      void merge(unsigned alid, const list<Consensus*> &chk, 
            list<DNAQualCount> &scrap, set<Consensus*> &toErase);
      /**
       * Merge this cluster with a single sequence
       * This can be used to merge a singleton.
       * @param sq input sequence that is not known in the seqstore
       *   so this function will add to the sequence store.
       * @param scrap container to collect sequence sgements.
       * @return true if the merge is done.
       *
       * Note: this method is slightly different from mergeKnownSequence
       *   where the input sequence is taken from the seqstore.
       */
      bool mergeSequence(unsigned ai, DNAQualCount &sq, list<DNAQualCount> &scrap);
      /**
       * This may be used for singleton directly.
       * This method differs from mergeSequence in that 
       * the input sequence sq is known to the sequence store.
       * It first try a very strangent cutoff with 0.99 identity cutoff,
       * then it use identity cutoff which should be be around 0.97.
       * @param scrap container to save unused segments.
       * @param sq is the pointer to a sequence in the seqstore.
       */
      bool mergeKnownSequence(unsigned ai, DNAQualCount *sq, list<DNAQualCount> &scrap);
      /**
       * Used by mergeSequence to see seq2 is a chimera.
       * seq2 is uaully a singleton presumably low quality reads.
       * @param ai aligner index in the aligner pool.
       * @param scrap container to accumulate un-absorbed segment.
       * @param seq2 input sequence to be compared to this Consens
       *         If seq2 is chimra, then the segment similar to
       *         the consensus will be absorbed. Unsimilar segments
       *         will be save into scrap.
       */
      bool mergeSegment(unsigned ai, 
            list<DNAQualCount> &scrap, const DNAQualCount& seq2);
      /**
       * Use the sequence as is.  No QC done.
       * Assumer the caller has done the QC.
       * Large cluster may get a better hand. This could lead to bias
       * toward large clusters.
       * @param ai aligner index in the thread pool
       * @param sq input sequence to be merged with the consensus.
       *      sq must be in the store.
       */
      void salvageKnownSequence(unsigned ai, DNAQualCount *sq);
      /**
       * Add a new sequence to the store and to the member. 
       * @return a pointer to this sequence in the store.
       * Note: should not update totalseq.
       * Does not seem to be useful. Tring to eliminate 
       * call to this function in the implementation.
       * Use of this function should be pretected by 
       * the internal mutex.
       */
      DNAQualCount* addMember(const DNAQualCount &c);
      /**
       * This consensus try to swallow all segments
       * provided in a container.
       * Do simple alignment. No further examination of 
       * alignment structure. If match perfectly then
       * merge the segment into the consensus. Otherwise,
       * just ignore.
       * Will reduced the list (segments) size if there segments 
       * absorbed into this cluster.
       * Very strigent selection.
       * @param segments components from chimera reads
       *     The list will shrink if some segment got absorbed into this
       *     consensus.
       * @return seg2cons save the segments that passed
       *  QC. But needs be decided later.
       *  Segment => consensus alignment information
       *  in tabular format.
       */
       void swallow(list<DNAQualCount> &segments,
            map<DNAQualCount*, vector<AlignQuality> > &undecided );  
      /**
       * One thread of work for a given chunk.
       * @param chnk input DNA segments unknown to the store.  It will be shrunk
       *     at the end of the operation. Each thread has its own.
       *     These sequences have not been added to the store.
       * @param seq2cons is the shared segment => consensus table
       *     This will be updated by individual threads.
       *
       * First try 0.99 identity and 0.9 input seq cov.
       * If not pass then will be put to the salvage pathway
       * if identity > cutoff
       */
      void swallowThread(unsigned ai, list<DNAQualCount> &chnk,
                  map<DNAQualCount*, vector<AlignQuality> > &seg2cons);

      /**
       * This consensus will try to swallow singletons in a containner.
       * This method uses swallowSingletonThread to do parallel
       * job.
       */
      void swallowSingleton(list<DNAQualCount*> &sgl, list<DNAQualCount> &seg,
                     map<DNAQualCount*, vector<AlignQuality> > &sgl_consRel);
      // swallow and swallThread should be applied to singleton
      // now I will write a particular version for singleton.
      // TODO: merge the singlton and non-single swallow function
      // The difference is object reference and pointer version
      // right now.
      /**
       * work on a chunk on one CPU
       * @param sglchk singleton chunk. 
       */
      void swallowSingletonThread(unsigned ai, list<DNAQualCount*> &sglchk, 
            list<DNAQualCount> &unmatchedSeg, 
             map<DNAQualCount*, vector<AlignQuality> > &undecided);
      /**
       * Get the quality information from
       * the particular aligner.
       * @param alid aligner index.
       */
      AlignQuality getAlignQuality(unsigned alid) const;

      /**
       * for debug show
       */
      void alignSequencesThread(unsigned ai, 
            const list<DNAQualCount*> &chk, ostream &ous);
      /**
       * To print the alignment of consens to a bunch of 
       * singletons, for debug.
       * @param seqs input sequenes to be aligned to this consensus.
       */
      void alignSequences(const list<DNAQualCount*> &seqs, ostream &ous);
      /**
       * different type for the input container.
       */
      void alignSequences(list<DNAQualCount> &seqs, ostream &ous);
      /**
       * @return the coverage result for the plotting utility.
       *     It is a const reference to protect from modification.
       */
      const vector<tuple<int,int,int> >& getCoverage() { return coverage; }

      ///// public static functions /////
      static void setAlnlengthCut(int alct) { alnlengthCut=alct; }
      /**
       * Length cutoff for alignment.
       * Different versions for convenience.
       */
      static int getAlnlengthCutoff() { return alnlengthCut; }
      static int getAlnlenCutoff() { return alnlengthCut; }
      static int getAlnlengthCut() { return alnlengthCut; }
      /**
       * 1.5x alnlengthcut
       */
      static int getExpandedAlnlengthCutoff() { return ceil(1.5*alnlengthCut); }
      static void setIdentityCut(double icut) { identityCut = icut; }
      static double getIdentityCut() { return identityCut; }
      static double getIdentityCutoff() { return identityCut; }
      /**
       * to get the same score method for other projects to share
       */
      static SimpleScoreMethod& getScoreMethod() { return scoreMethod; }
      static void setStore(DNAQualstore* store) { seqstore=store; }
      /**
       * Add a sequence to the store of the consensus object.
       * @return a pointer to the inserted object.
       */
      static DNAQualCount* add(const DNAQualCount &seq) { 
         return seqstore->addSequence(seq); }
      /**
       * Add a new sequence into the store.
       * Future reference to the pointer will be to the 
       * object in the store.
       */
      static DNAQualCount* add(DNAQualCount &&seq) { 
         return seqstore->addSequence(std::move(seq)); }
      /**
       * Try to use this function. In the future when 
       * multithreading is enable, you need synchronization.
       */
      static string getNewseqName();
      static string getClusterName() { 
         return string("cl") + to_string(clusterid++); 
      }
      static void resetClusterId() { clusterid=1; }
      /**
      * helper method to separate singleton from clusters
      * @param clusters these clusters should be sorted from large to small.
      */
      static list<DNAQualCount*> separateSingleton(vector<vector<string> > &clusters);
      /**
       * This should be called exactly once before using this class.
       */
      static void initAligner(unsigned numthreads);
      static void deallocateAligner();
      static unsigned int numberOfAligners() { return aligner.size(); }
};

/**
 *  from large to samll
 *  For some reason, the lambda version is causing segmentation fault.
 *  There is a bug in the compiler such that I cannot use elegant code.
 */
class sortVectorBySize {
   public:
      bool operator()(const vector<string> &v1, const vector<string> &v2) const {
         return v1.size() > v2.size();
      }
};

class sortConsensusByTotalSequence {
   public :
      bool operator()(const Consensus *c1, const Consensus *c2) const {
         return c1->getTotal() > c2->getTotal();
      }
};

/**
 * Helper function from consensus result to taxonomy mapping
 */
/**
 * Swarm specific
 * @clfile swarm cluster output file.
 * @return the cluster in a tabular format with large clusters
 *   at the front and smaller cluster at the end. Sorted
 *   by cluster size in decending order.
 *
 * Although the cluster is sorted by the number of unique ids,
 * the real number of sequences may differ because one sequence
 * may represent more member than the other!
 */
vector<vector<string> > readCluster(const string &clfile);
/**
 * Read from Mothur clustering result.
 * Mothur cluster level is identified by the first column.
 * cl_dist #cl  actual members separated by comma, cluster separated by space.
 * 0.01
 * 0.02  766   ccs1,ccs4,ccs9 ccs99,ccs678,ccs1 ccs7,ccs2
 * if you specify 0.015, the it will pick 0.01 row
 * if you specify 0.02, it will pick up the 0.02 row
 * @param cutoff the cluster level at which or the largest less than which
 *        to use.
 */ 
vector<vector<string> > readCluster(const string &clfile, float cutoff);
/**
 *  ouf the output stream. Good for multiple calls.
 */
void writeConsensusFasta(const vector<DNAQualCount> &con, ostream &ouf);
void writeConsensusFasta(const vector<DNAQualCount> &con, const string &fasFile);

class ClusterProgparam {
   private:
      string refkmerFile;
      string fastqfile;
      bool salvage;
      static string bindir;
      static string datadir;

   public:
      ClusterProgparam() : refkmerFile(), fastqfile(), salvage(true) { }
      ClusterProgparam(const string &refkfile) : refkmerFile(refkfile),
         fastqfile(), salvage(true) { }
      void setInputFile(const string& file) {
          fastqfile=file; }
      string getKmerReference() const { return datadir + "/" + refkmerFile; }
      string getFstem() const { 
         return fastqfile.substr(0, fastqfile.rfind('.')); }
      void setReference(const string &file) { refkmerFile=file; }
      bool salvageSingle() const { return salvage; }
      void noSalvage() { salvage = false; }
      void doSalvage() { salvage = true; }
      static const string& getBinaryDirectory() { return bindir; }
      static void setBinaryDirectory();
};

/**
 * This is used by the caller in the shell to continue to work with
 * the consesnsus fasta file.
 * First round of consensus building for clusters with
 * more than one member. In second round, singletons are absorbed.
 * Inside this function, writeConsensus() function was used to
 * write the result to several different files.
 * @param clid only cluster ids both clusters and singletons as input.
 * @return the name for the consensus fasta file.
 */
string rebuildConsensus(vector<vector<string> > &clid, const ClusterProgparam &par);
/**
 * Simpler version.
 * Not doing merging or singleton absorbing. This is used for 
 * Generating reference file.
 */
string buildConsensusAll(vector<vector<string> > &clid, const ClusterProgparam &par);
/**
 * Only apply to cluster with 2 or more members.
 * Singletons are excluded from this operation.
 * large cluster will absorb small clusters
 * @param allcons All the cluster consensus with > 1 members.
 * @param salvage flag to do salvage operation or not.
 * @return the scrap fragments that could be in the reverse direction.
 */
list<DNAQualCount> combineConsensus(list<Consensus*> &allcons, bool salvage);
/**
 * write results to file
 * @param scrap the input sequences for this function
 */
void writeScrap(const list<DNAQualCount> &scrap, 
      const string &suffix, const ClusterProgparam &par);

/**
 * compare each consensus to a singleton to see which singleton
 * can be merged into the consensus (cluster). Should implement
 * better single-linkage clustering methods.
 * @param cons are all the consensus clsuter objects.
 * @param sgl list of singleton sequences to be merged into consensus.
 *   For those singleton that did not meet the cutoff
 *   they will be left in the sgl list.
 *   If they are incorporated into the cluster then
 *   they will be removed from the list.
 */
list<DNAQualCount> consensusSwallowSingleton(list<Consensus*> &cons, list<DNAQualCount*> &sgl);
/**
 * This is a helper method used by rebuildConsesnsus().
 * It dump all the contents of the consensus into files of fixed file names.
 * The stalactite graph output should be added to this function.
 * @param cons is the input Sonsensus objects in a list of pointers.
 * @fstem is a main part of the file name, different suffixes will be added
 *     to this main part for different aspects of the consensus file dump.
 */
string writeConsensus(const list<Consensus*>& cons, const string &fstem);
void writeSequencePointer(const list<DNAQualCount*> &sgl, const string &suffix, const ClusterProgparam &par);

/**
 * Align the consensus to the singletons that showed no
 * significant similarity to any cluster.
 * This is for debug and pipeline development.
 * The Consensus of each cluster is used to for alignment
 * not the seed sequence.
 * @param par parameters for the clustering program.
 *   output file name will be decided by par.
 * @param cons is the consensus still standing after the 
 *    combine operation
 * @param single are the singletons still not absorbed.
 * This can be achieved through running external programs
 * by hand.  This should be removed in mature pipelnes.
 * Note: don't use this one for mature pipeline.
 */
void showConsensusSingleAlign(const list<Consensus*> &cons, 
      const list<DNAQualCount*> &single, const ClusterProgparam &par);

/**
 * only difference is the second input is a list of DNAQualCount objects
 * in stead of pointer to these objects.
 */
void showConsensusScrapAlign(const list<Consensus*> &cons, 
      list<DNAQualCount> &single, const ClusterProgparam &par);
/**
 * One segments match to multiple consensus, this algorithm
 * pick the best consensus (cluster) to be merged to.
 * @param sgl2cons segment and consensus relational table.
 */
void salvageLQSequence(list<Consensus*> &cons, 
      map<DNAQualCount*, vector<AlignQuality> > &sgl2cons);



#endif
