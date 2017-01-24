#include <algorithm>
#include <stack>
#include <iomanip>
#include <cassert>
#include <algorithm>
#include <sys/stat.h>
#include <sys/types.h>

#include <bioseq.h>
#include <alnexaminer.h>
#include <codon.h>
#include <strformat.h>
#include "consensus.h"
#include "stalactitegraph.h"

using namespace orpara;

//#define DEBUG

// if this is used for salvage, then there is no need for
// cutoff at the salvage stage only do salvage or not
bool AlignQuality::passCutoff(float identitycut) const {
   if (getIdentity()>0.99 && alnlen>200) return true;
   if ((repcov > 0.9 || readcov > 0.9) && getIdentity()>identitycut) 
      return true;
   return false;
}

bool AlignQuality::operator<(const AlignQuality &o) const 
{ 
   if (getIdentity() < o.getIdentity()) return true; 
   if (getIdentity() > o.getIdentity()) return false; 
   return repcov < o.repcov;
}
/// Consensus satic members ///////
SimpleScoreMethod Consensus::scoreMethod(11, -10, -39, -4);
double Consensus::identityCut=0.87;
int Consensus::alnlengthCut=200;
DNAQualstore* Consensus::seqstore(0); // not intialized
int Consensus::newseqid=1;
int Consensus::clusterid=1;
vector<Dynaln<SimpleScoreMethod>* > Consensus::aligner = {};

//////// helper to be used only in this object file //////
float roundPercent(float frac) {
      return roundf(frac*10000)/100;
}

string Consensus::getNewseqName() {
   return "merged" + to_string(newseqid++);
}

template<class T>
vector<pair<T, int> > map2vectorByCount(const map<T, int>& src) { 
   vector<pair<T,int> > arr;
   for (auto itr=src.begin(); itr != src.end(); ++itr) {
      arr.push_back(pair<T, int>(itr->first, itr->second));
   }
   sort(arr.begin(), arr.end(), SortPairByInt<T>());
   return arr;
}

////////////// member functions /////////////////
// not tested yet, better not to use this one.
Consensus& Consensus::operator=(Consensus &&o) {
   if (this != &o) {
      member=std::move(o.member);
      totalseq=o.totalseq; repseq=o.repseq;
      alnsummary=std::move(o.alnsummary);
      bases=std::move(o.bases); inserts=std::move(o.inserts);
      quality=std::move(o.quality); consensus=std::move(o.consensus);
      consensusQ=std::move(o.consensusQ);
      consensusValid=o.consensusValid;
      consensusName=std::move(o.consensusName);
      coverage=std::move(o.coverage);
   }
   return *this;
}

void Consensus::pickRepseq() const {
   if (member.size() == 1) { //singleton special case
      repseq=member.front();
      member.clear();
      return;
   }
   // member is not sorted
   typedef typename list<DNAQualCount*>::iterator LIT;
   LIT lit=member.begin();
   int maxcnt = (*lit)->getCount();
   vector<LIT> rep(1, lit); // vector of iterators
   ++lit;
   while (lit != member.end()) {
      if ((*lit)->getCount() > maxcnt) {
         maxcnt = (*lit)->getCount();
         rep={lit};
      }
      else if ((*lit)->getCount() == maxcnt) {
         rep.push_back(lit);
      }
      ++lit;
   }
   if (rep.size() == 1) {
      repseq = *(rep[0]);
      member.erase(rep[0]);
      return;
   }
   int maxlen=(*rep[0])->length();
   vector<LIT> newrep={rep[0]};
   for (int i=1; i<rep.size(); ++i) {
      if ((*rep[i])->length() > maxlen) {
         maxlen = (*rep[i])->length();
         newrep={rep[i]};
      }
      else if ((*rep[i])->length() == maxlen) {
         newrep.push_back(rep[i]);
      }
   }
   //if (newrep.size() > 1) {
   //   cerr << ">1 sequences tie on count and length. Picking the first as rep\n";
      // in the futrue you can pick the one with the higher quality score
   //}
   repseq=*(newrep[0]);
   member.erase(newrep[0]);
}

// this should be called in the constructor
void Consensus::initialize() {
   pickRepseq();
   bases.resize(repseq->length());
   inserts.resize(repseq->length());
   quality.resize(repseq->length());
   for (int i=0; i<repseq->length(); ++i) {
      bases[i].insert(make_pair((*repseq)[i], repseq->getCount()));
      quality[i](repseq->Qat(i), repseq->getCount());
   }
   // set up the proper number of thread to do the work
   for (unsigned i=0; i<aligner.size(); ++i) {
      aligner[i]->setSeq1(*repseq);
   }
}

void displayStagerGap(const string &top, const string &bottom, int i, int w) {
   int s=i-w;
   if (s<0) s=0;
   int len=2*w;

   if (s+len >= top.length()) {
      cerr << top.substr(s) << endl
         << bottom.substr(s) << endl;
   }
   else {
      cerr << top.substr(s, 2*w) << endl
         << bottom.substr(s, 2*w) << endl;
   }
}

void Consensus::closeStagerGap(string &top, string &bottom) {
   int width=12;
   int i=width, j, k, b, e;
   while (i< int(top.size())-width) {
      // only single gap
      if (bottom[i] == '-' && bottom[i+1] != '-') {
         b=i-width;
         j=i-1;
         bool stagerLeft=false;
         while (j>b) {
            if (top[j] == '-') { // stager left
               stagerLeft=true; break;
            }
            --j;
         }
         bool stagerRight=false;
         e=i+width;
         k=i+1;
         while (k<e) {
            if (top[k] == '-') {
               stagerRight=true; break;
            }
            ++k;
         }
         if (stagerLeft && stagerRight) {
            cerr << "double stager gap:\n";
            displayStagerGap(top,bottom,i,width);
            bottom.erase(i,1);
            if (i-j < k-i) { 
               top.erase(j,1); 
            }
            else { 
               top.erase(k,1); 
            }
         }
         else if (stagerLeft) {
            cerr << "stager gap left:\n";
            displayStagerGap(top,bottom,i,width);
            top.erase(j,1);
            bottom.erase(i,1);
         }
         else if (stagerRight) {
            cerr << "stager gap right:\n";
            displayStagerGap(top,bottom,i,width);
            top.erase(k,1);
            bottom.erase(i,1);
         }
         i += 2;
      }
      else {
         if (bottom[i] == '-') {
            while (i < int(top.size()) - width && bottom[i] == '-') 
               ++i;
         }
         else ++i;
      }
   }
}

void Consensus::scanStagerGap(string &top, string &bottom) {
   int width=12;
   int i=width, j, k, b, e;
   while (i< int(top.size())-width) {
      // only single gap
      if (bottom[i] == '-' && bottom[i+1] != '-') {
         b=i-width;
         j=i-1;
         bool stagerLeft=false;
         while (j>b) {
            if (top[j] == '-') { // stager left
               stagerLeft=true; break;
            }
            --j;
         }
         bool stagerRight=false;
         e=i+width;
         k=i+1;
         while (k<e) {
            if (top[k] == '-') {
               stagerRight=true; break;
            }
            ++k;
         }
         if (stagerLeft && stagerRight) {
            cerr << "double stager gap:\n";
            displayStagerGap(top,bottom,i,2*width);
         }
         else if (stagerLeft) {
            cerr << "stager gap left:\n";
            displayStagerGap(top,bottom,i,2*width);
         }
         else if (stagerRight) {
            cerr << "stager gap right:\n";
            displayStagerGap(top,bottom,i,2*width);
         }
         i += 2;
      }
      else {
         if (bottom[i] == '-') {
            while (i < int(top.size()) - width && bottom[i] == '-') 
               ++i;
         }
         else ++i;
      }
   }
}

// use the quality info of the aligned part of 
// read, and its count
// only one thread can run this function at one time
// this is synchronized.
void Consensus::accumulate(unsigned alid, DNAQualCount *read) 
{
   string top(aligner[alid]->getTopAln());
   string bottom(aligner[alid]->getBottomAln());
   //scanStagerGap(top, bottom);
   size_t tb = (unsigned int)aligner[alid]->topBeginIndex();
   size_t bb = (unsigned int)aligner[alid]->bottomBeginIndex();
   // i for index in the alignment
   // refi is for index in the repseq
   string::size_type i, refi, readi;
   i=0;
   refi = tb;
   readi= bb;
   unsigned int t = 0;
   bool insideInsert=false;
   string aInsert;
   lock_guard<mutex> lk(mut);
   while (i < top.length()) { // alignment index
      if (top[i] != '-') { // out of insert state
         if (insideInsert) { // save the insert if came from an insert state
            inserts[refi-1][aInsert] += read->getCount();
            insideInsert = false;
         }
         bases[refi][bottom[i]] += read->getCount();
         if (bottom[i] == '-') { // guess the quality of a gap lower
            quality[refi](floor(0.8*read->Qat(readi)), read->getCount()); // del quality
         }
         else {
            quality[refi](read->Qat(readi), read->getCount());
            ++readi; 
         }
         ++refi;
      }
      else { // top is gap char, bottom cannot be gap
         if (!insideInsert) {
            aInsert = string(1, bottom[i]); // clear the old one by new value
            insideInsert = true;
         }
         else {
            aInsert += bottom[i];
         }
         ++readi;
      }
      ++i;
   }
   alnsummary[AlignInfo(aligner[alid]->getScore(), 
         roundPercent(aligner[alid]->getIdentity()),
         roundPercent(aligner[alid]->getCov1()), 
         aligner[alid]->getSeq2Length())] += read->getCount();
   invalidateConsensus();
}

int Consensus::getTotal() const {
   lock_guard<mutex> lk(mut);
   if (totalseq == 0 || !isConsensusValid()) 
      computeTotalSequences(); 
   return totalseq; 
}

int Consensus::computeTotalSequences() const {
   totalseq=repseq->getCount();
   for (auto it=member.begin(); it != member.end(); ++it) {
      totalseq += (*it)->getCount();
   }
   return totalseq;
}

bool Consensus::debugValidTotal() {
   int freshsum=repseq->getCount();
   for (auto it=member.begin(); it != member.end(); ++it) {
      freshsum += (*it)->getCount();
   }
   if (totalseq != freshsum) {
      cerr << getConsensusName() << " totalseq wrong\n"
         << " cluster size: " << numseq()
         << " fresh: " << freshsum << " total: " << totalseq
         << endl;
   }
   return freshsum == totalseq;
}

bool Consensus::alignPassQC(unsigned ai) const {
   if (aligner[ai]->getAlnlen() < alnlengthCut || 
         aligner[ai]->getIdentity() < identityCut)
      return false;
   if (aligner[ai]->getIdentity() > 0.989 && 
         (aligner[ai]->getCov1()>0.8 || aligner[ai]->getCov2()>0.8))
      return true;
   // all gaps single and idcount+gap/alignLen > 99% 
   pair<int,int> numgap = aligner[ai]->numgaps();
   pair<int,int> gaplen = aligner[ai]->gaplen();
   if (gaplen.first-numgap.first < 3 && gaplen.second - numgap.second < 3
         && aligner[ai]->getNogapIdentity() > 0.989
         && (aligner[ai]->getCov1()>0.8 || aligner[ai]->getCov2()>0.8
            || aligner[ai]->getAlnlen() > 1100)) {
      return true;
   }
   return false;
}

// should speed up the computation if the
// direction is known.
// not all members maybe used! still need to do more fine
// counting here.
bool Consensus::consume(unsigned ai, DNAQualCount* s) {
   aligner[ai]->setSeq2(*s);
   aligner[ai]->runlocal();
   if ((aligner[ai]->getIdentity() > identityCut 
            && (aligner[ai]->getCov1() > 0.78 || aligner[ai]->getCov2() > 0.78))
         || (aligner[ai]->getIdentity()> 0.989 && aligner[ai]->getAlnlen() > getAlnlengthCutoff())
         || alignPassQC(ai)) {
      accumulate(ai, s);
      return true;
   }
   if (aligner[ai]->getIdentity() < 0.35 || aligner[ai]->getCov1() < 0.3 
         || aligner[ai]->getCov2()<0.3) {
      s->revcomp();
      aligner[ai]->setSeq2(*s);
      aligner[ai]->runlocal();
      if (aligner[ai]->getIdentity()>identityCut && 
            (aligner[ai]->getCov1()>0.78 || aligner[ai]->getCov2()>0.78)
         || (aligner[ai]->getIdentity() > 0.989 && aligner[ai]->getAlnlen() > getAlnlengthCutoff())) 
      {
         accumulate(ai, s);
         return true;
      }
      if (aligner[ai]->getIdentity() < identityCut || 
            (aligner[ai]->getCov2() < 0.78 && aligner[ai]->getCov1() < 0.78)) {
         cerr << "repseq: " << *repseq << endl;
         cerr << "cannot consume sequence\n" << *s << endl;
         cerr << "inside consume() function of Consensus\n";
         aligner[ai]->printAlign(cerr);
         cerr << "do more investigation for this one\n";
         //throw logic_error("direction of read relative to repseq needs to be sorted out");
      }
   }
   else if (s->length() > getRepseqLength()+200 
         || getRepseqLength() > s->length()+200) {
      // long sequences needs to test reverse complement
      s->revcomp();
      aligner[ai]->setSeq2(*s);
      aligner[ai]->runlocal();
      if (aligner[ai]->getIdentity() > identityCut && 
            (aligner[ai]->getCov1() > 0.75 || aligner[ai]->getCov2() > 0.75)) {
         accumulate(ai, s);
         return true;
      }
   }
   else { // if not used, simply warning, need to do more.
      cerr << "forward align failed QC. This is fine\n";
      aligner[ai]->printAlign(cerr);
   }
   return false;
}

// for threading a block of elements
void Consensus::consumeBlock(unsigned alignerId, list<DNAQualCount*> &chk) {
   for (auto i=chk.begin(); i != chk.end(); ++i) {
      consume(alignerId, *i);
   }
}

// this is essentially the run method
// In this operator, members should not be added
// because they are intialized alrady in member.
void Consensus::operator()() {
   vector<list<DNAQualCount*> > blocks(min(numberOfAligners(), (unsigned int)member.size()));
   int b=0;
   for (auto it = member.begin(); it != member.end(); ++it) {
      blocks[b++ % blocks.size()].push_back(*it);
   }
   vector<thread> threads(blocks.size()); 
   for (unsigned i=0; i < blocks.size(); ++i) {
      threads[i] = thread(&Consensus::consumeBlock, this, i, ref(blocks[i]));
   }
   for_each(threads.begin(), threads.end(), mem_fn(&thread::join));
   if (consensusName.empty()) {
      consensusName=getClusterName();
   }
   else {
      cerr << "consensusName not empty!\n";
      cerr << "name exists for consensus: " << consensusName << endl;
   }
}

// after refine the aligner's first sequence is
// out of scope, if need reuse, need to make
// the consensus a full object, not just a string.
void Consensus::refine() {
   // need to set reference to consensus first
   // then align every member + repseq to the consensus
   // if repseq == consensus, then the alighment can be skipped?
   cerr << "round 2 consensus for " << consensusName << endl;
   buildConsensus();
   DNAQualCount consObj = getResult();
   // assign aligners to use consensus
   // initialize only if refseq != consensus
   // new members have been added
   int b=0;
   if (getGrossBaseCoverage() < 3) {
      cerr << "cluster depth: " <<  getGrossBaseCoverage()
         << " too low to have meaningfull consensus\n";
   }
   else if (repseq->toString() != getConsensus()) {
      cerr << "reqseq and consensus different\n";
      //cerr << consObj << endl << *repseq << endl;
      bases.clear();
      inserts.clear();
      quality.clear();
      bases.resize(getConsensus().size());
      inserts.resize(getConsensus().size());
      quality.resize(getConsensus().size());
      // set up the proper number of thread to do the work
      for (unsigned i=0; i<aligner.size(); ++i) {
         aligner[i]->setSeq1(consObj);
      }
      vector<list<DNAQualCount*> > blocks(min(numberOfAligners(), (unsigned int)member.size()));
      blocks[0].push_back(repseq);
      for (auto it = member.begin(); it != member.end(); ++it) {
         blocks[b++ % blocks.size()].push_back(*it);
      }
      vector<thread> threads(blocks.size()); 
      for (unsigned i=0; i < blocks.size(); ++i) {
         threads[i] = thread(&Consensus::consumeBlock, this, i, ref(blocks[i]));
      }
      for_each(threads.begin(), threads.end(), mem_fn(&thread::join));
      // already got name 
      if (!debugValidTotal()) {
         cerr << "freshly constructed consensus had bad total sequences\n";
      }
      buildConsensus();
   }
   else {
      cerr << "Lucky, the consesus is the same as repseq!\n";
      // nothing needs to be done
   }

}

template<class T>
void clearVectorMap(vector<map<T, int> > &vm) {
   if (!vm.empty()) {
      for (size_t i=0; i<vm.size(); ++i) {
         vm[i].clear();
      }
   }
}

void Consensus::clear() {
   clearVectorMap(bases);
   clearVectorMap(inserts);
}

/*
 * count all nongap characters.
 */
int getTotalValue(const map<char, int> &bc) {
   int sum = 0;
   for (auto itr = bc.begin(); itr != bc.end(); ++itr) {
      sum += itr->second;
   }
   return sum;
}
int getTotalNogapValue(const map<char, int> &bc) {
   int sum = 0;
   for (auto itr = bc.begin(); itr != bc.end(); ++itr) {
      if (itr->first != '-')
         sum += itr->second;
   }
   return sum;
}

/**
 * Should consider if insert is at high frequency.
 * Insert/depth > 0.5 should be used.
 * Reporting a gap char if it is the highest frequency.
 * Tail and end should be trimmed.
 */
void Consensus::buildConsensus() const {
   consensus.clear();
   coverage.clear();
   consensus.reserve(bases.size() + bases.size()/10>10? bases.size()/10 : 10); 
   coverage.reserve(consensus.capacity());
   // at most 10% insertions
   int totalB, totalBnogap;
   for (size_t i=0; i<bases.size(); ++i) {
      vector<pair<char, int> > bcount = map2vectorByCount(bases[i]);
      totalB=getTotalValue(bases[i]);
      totalBnogap = getTotalNogapValue(bases[i]);
      if (bcount[0].first == '-') {
         // need statistical model here, taking into the counts
         // for a probaility, for future work
         if (double(bcount[0].second)/totalB - double(bcount[1].second)/totalB > 0.15) {
            continue; 
            // dont print a gap char
         }
         else {
            consensus.append(1, bcount[1].first);
            consensusQ.push_back(int(ceil(quality[i].getMean())));
         }
      }
      else {
         consensus.append(1, bcount[0].first);
         consensusQ.push_back(int(ceil(quality[i].getMean())));
      }
      coverage.push_back(make_tuple(bcount[0].second, totalBnogap, totalB));
      if (!inserts[i].empty()) {
         vector<pair<string, int> > icount = map2vectorByCount(inserts[i]); 
         if (icount[0].second > 1 && double(icount[0].second)/totalB > 0.48) {
            int totalInsert=0;
            for (size_t i=0; i<icount.size(); ++i) {
               totalInsert += icount[i].second;
            }
            consensus.append(icount[0].first);
            // TODO: need to add more detailed data structure to insert
            int avgQ=(int)ceil(quality[i].getMean()*0.8); // give it a lower score
            for (int k=0; k<icount[0].first.size(); ++k) {
               consensusQ.push_back(avgQ);
            }
            //cerr << "significant insert: " << icount[0].first << " " << icount[0].second 
            //   << " relative to ref: " 
            //   << bcount[0].first << " " << bcount[0].second << endl;
            for (int r=0; r<icount[0].first.size(); ++r) {
               coverage.push_back(make_tuple(icount[0].second, totalInsert, totalB));
            }
         }
      }
   }
   computeTotalSequences();
   consensusValid=true;
   clearAverage();
}

array<stddev,3> Consensus::getAverageCoverage() const {
   if (!consensusValid) {
      buildConsensus();
   }
   if (averageCoverage[0].empty()) {
      for (size_t i=0; i<coverage.size(); ++i) {
         averageCoverage[0](get<0>(coverage[i]));
         averageCoverage[1](get<1>(coverage[i]));
         averageCoverage[2](get<2>(coverage[i]));
      }
   }
   return averageCoverage;
}

vector<int> Consensus::getMeanQuality() const {
   vector<int> result(quality.size());
   for (int i=0; i<quality.size(); ++i) {
      result[i]=int(ceil(quality[i].getMean()));
   }
   return result;
}

void Consensus::printBases(ostream &ous) const {
   ous << getConsensusName() << endl
      << "position\tbase\tcount(frequency)\n";
   for (size_t i=0; i<bases.size(); ++i) {
      int totalBase = getTotalValue(bases[i]);
      ous << i+1;
      vector<pair<char, int> > bcount = map2vectorByCount(bases[i]);
      for (size_t b = 0; b<bcount.size(); ++b) {
         ous << '\t' << bcount[b].first << '\t' << bcount[b].second << '(' 
            << setprecision(6) << double(bcount[b].second)/totalBase << ')';
      }
      // output potentional insertion
      if (!inserts[i].empty()) {
         auto jtr = inserts[i].begin();
         while (jtr != inserts[i].end()) {
            ous << "\tinsert: " << jtr->first << '\t' << jtr->second 
               << '(' << setprecision(6) << double(jtr->second)/totalBase << ')';
            ++jtr;
         }
      }
      ous << endl;
   }
}

void Consensus::printConsensus(ostream &ous) const {
   if (consensus.empty() || !consensusValid) buildConsensus();
   ous << '>' << getConsensusName() << " representative="
      << getRepseq()->getName() << " members=";
   auto it=member.begin();
   ous << (*it)->getName();
   ++it;
   while (it != member.end()) {
      ous << "," << (*it)->getName();
      ++it;
   }
   ous << endl << consensus << endl;
}

void Consensus::printQuality(ostream &ous) const {
   ous << getConsensusName() << endl
      << "position\tqavg\tqstddev\tcount\n";
   for (auto i=0; i<quality.size(); ++i) {
      ous << i+1 << '\t';
      quality[i].print(ous) << endl;
   }
}

void Consensus::printCoverage(ostream &ous) const {
   if (consensus.empty() || !consensusValid) buildConsensus();
   ous << getConsensusName() << endl
      << "position\tconsensusBase\tnogapTotal\ttotal\n";
   for (size_t i=0; i<coverage.size(); ++i) {
      ous << i+1 << '\t' << get<0>(coverage[i])
         << '\t' << get<1>(coverage[i]) 
         << '\t' << get<2>(coverage[i]) <<endl;
   }
}

void Consensus::printResult(ostream &oucons, 
            ostream &oubase, ostream &ouqual, ostream &oualni) const 
{
   printConsensus(oucons);
   printBases(oubase);
   printQuality(ouqual);
   printAlninfo(oualni);
} 

// over loaded version with base vertical coverage
void Consensus::printResult(ostream &oucons, 
            ostream &oubase, ostream &ouqual, ostream &oualni,
            ostream &oucov) const 
{
   printConsensus(oucons);
   printBases(oubase);
   printQuality(ouqual);
   printAlninfo(oualni);
   printCoverage(oucov);
} 

void Consensus::printResult(const string &consensusfile, 
      const string &basefile, const string &qualfile,
      ios_base::openmode writeMode) const 
{
   ofstream ouc(consensusfile, writeMode);
   ofstream oub(basefile, writeMode);
   printConsensus(ouc);
   printBases(oub);
   ouc.close();
   oub.close();
   ofstream ouqual(qualfile, writeMode);
   if (ouqual.fail()) {
      cerr << "Failed to open " << qualfile << endl;
      exit(1);
   }
   printQuality(ouqual);
   ouqual.close();
   printAlninfo("aligninfo.txt");
}

void Consensus::printAlninfo(ostream &ous) const {
   ous << getConsensusName() << endl
      << "score\tidentity\tcoverage\treadlen\tcount\n";
   for (auto itr = alnsummary.cbegin(); itr != alnsummary.cend(); ++itr) {
      ous << itr->first << "\t" << itr->second << endl;
   }
}

void Consensus::printAlninfo(const string &alninfoFile) const {
   ofstream ous(alninfoFile.c_str());
   if (ous.fail()) {
      cerr << "Failed to open alninfo file: " << alninfoFile << endl;
      exit(1);
   }
   ous << "score\tidentity\tcoverage\treadlen\tcount\n";
   for (auto itr = alnsummary.cbegin(); itr != alnsummary.cend(); ++itr) {
      ous << itr->first << "\t" << itr->second << endl;
   }
   ous.close();
}

DNAQualCount Consensus::getResult() const { 
   DNAQualCount tmp(consensusName, getConsensus(), consensusQ, getTotal());
   tmp.appendTitle("depth " + itos((int)round(getAverageCoverage()[2].getMean())), ", ");
   return tmp;
}

// you have to do by work chunks
// working within one thread, no need for synchronization
void Consensus::merge(unsigned alid, const list<Consensus*> &chk, 
      list<DNAQualCount> &scrap, set<Consensus*> &toErase) {
   aligner[alid]->setSeq1(*repseq);
   //cerr << getClusterName() << " " << getTotal() << endl;
   //different threads will compete for cerr and causing trouble
   for (auto it=chk.begin(); it != chk.end(); ++it) {
      if ((*it)->getTotal() > getTotal()) {
         cerr << "this cluster " << getConsensusName() << " size: "
            << getTotal() << " small cluster: " 
            << (*it)->getClusterName() << " size: " << (*it)->getTotal() << endl;
         throw ConsensusSizeException("left cluster size should be larger!");
      }
      DNAQualCount c2=(*it)->getResult(); 
      if (mergeSequence(alid, c2, scrap)) {
         // modify shared data structure
         lock_guard<mutex> lk(mutExternal);
         toErase.insert(*it);
      }
   }
}

// this may be used to merge singleton
// sq is not known in the seqstore.
bool Consensus::mergeSequence(unsigned ai, DNAQualCount &sq, list<DNAQualCount> &scrap) {
   aligner[ai]->setSeq2(sq);
   aligner[ai]->runlocal();
   if (alignPassQC(ai)) {
      unique_lock<mutex> lk(mut);
      DNAQualCount *newseq=addMember(sq);
      lk.unlock();
      accumulate(ai, newseq);
      // consensus  =======
      // sq    -------------------
      // if extra sequence then extract new seq
      if (aligner[ai]->getCov1() > 0.9 && aligner[ai]->getCov2() < 0.7) { 
         // extract extra segment of sq
         if (aligner[ai]->bottomBeginIndex() > alnlengthCut) {
            lock_guard<mutex> guard(mutExternal);
            scrap.push_back(sq.subsequenceWithName(0, aligner[ai]->bottomBeginIndex()));
         }
         if (sq.length() - aligner[ai]->bottomEndIndex() > alnlengthCut) {
            lock_guard<mutex> guard(mutExternal);
            scrap.push_back(sq.subsequenceWithName(aligner[ai]->bottomEndIndex()));
         }
      }
      return true;
   }
   bool seqMerged=false;
   if (aligner[ai]->getIdentity() > identityCut 
         && aligner[ai]->getAlnlen() > 3*alnlengthCut
         && (aligner[ai]->getCov1() > 0.7 || aligner[ai]->getCov2() > 0.7)) 
   {
      seqMerged = mergeSegment(ai, scrap, sq);
   }
   return seqMerged;
}

bool Consensus::mergeKnownSequence(unsigned ai, 
      DNAQualCount *sq, list<DNAQualCount> &scrap)
{
   aligner[ai]->setSeq2(*sq);
   aligner[ai]->runlocal();
   if (alignPassQC(ai)) {
      accumulate(ai, sq);
      lock_guard<mutex> lk(mut);
      member.push_back(sq);
      // debug code can be removed in production version
      cerr << __FILE__ << ":" << __LINE__ << " " << __func__ << " " << getConsensusName() << ": HQ merge\n";
      return true;
   }
   bool seqMerged=false;
   if (aligner[ai]->getIdentity() > identityCut 
         && aligner[ai]->getAlnlen() > 3*alnlengthCut
         && (aligner[ai]->getCov1() > 0.7 || aligner[ai]->getCov2() > 0.7)) 
   {
      cerr << __FILE__ << ":" << __LINE__ << " " << __func__ 
         << " " << getConsensusName() << ": passed identity cut: " << identityCut << " trying mergeSegment\n";
      seqMerged = mergeSegment(ai, scrap, *sq);
   }
   return seqMerged;
} 

void Consensus::salvageKnownSequence(unsigned ai, DNAQualCount *sq) {
   aligner[ai]->setSeq1(*repseq);
   aligner[ai]->setSeq2(*sq);
   aligner[ai]->runlocal();
   //controlled by AlignQuality.passCutoff(identity=0.85)
   accumulate(ai, sq);
   lock_guard<mutex> lk(mut);
   member.push_back(sq);
}

// only used by mergeSequence
bool Consensus::mergeSegment(unsigned ai, list<DNAQualCount> &scrap, const DNAQualCount &seq2) 
{
   Alnexaminer seger;
   vector<Alnseg> zone = seger(aligner[ai]->getMiddleAln());
   if (!seger.isChimera()) return false;
   vector<pair<int,int> > matchIndex=aligner[ai]->getAlnindexVector();
   int numsegSwallowed=0;
   for (size_t i=0; i<zone.size(); ++i) {
      if (zone[i].length() <= getAlnlenCutoff()) {
         continue;
      }
      int zb=zone[i].b;
      int ze=zone[i].e;
      int bb=matchIndex[zb].second;
      while (bb == -1) {
         ++zb;
         bb=matchIndex[zb].second;
      }
      ++bb;
      --ze;
      int ee=matchIndex[ze].second;
      while (ee == -1) {
         --ze;
         ee=matchIndex[ze].second;
      }
      ++ee;
      if (zone[i].identity > 0.99) {
         unique_lock<mutex> lk(mut);
         DNAQualCount *segSeq = add(seq2.subseq(bb,ee));
         lk.unlock();
         if (consume(ai, segSeq)) {
            unique_lock<mutex> lkm(mut);
            member.push_back(segSeq);
            lkm.unlock();
            ++numsegSwallowed;
         }
         else {
            cerr << "warn segment should be consumed but not. Check code\n";
         }
      }
      else {
         lock_guard<mutex> lg(mutExternal);
         scrap.push_back(seq2.subseq(bb,ee));
      }
   }
   if (numsegSwallowed>0) {
      return true;
   }
   return false;
}

// do paralelle job
// moving object around to improve performance
// not using pointers.
void Consensus::swallow(list<DNAQualCount> &segs,
      map<DNAQualCount*, vector<AlignQuality> > &undecided) 
{
   // this is shared among threads
   // divide into chunks
   cerr << "working on " << getConsensusName() 
      << " start with " << undecided.size() << " undecided\n";
   vector<list<DNAQualCount> > chunks(min((unsigned int)segs.size(), numberOfAligners()));
   int b=0;
   for (auto it=segs.begin(); it != segs.end(); ++it) {
      chunks[b++ % chunks.size()].push_back(std::move(*it));
   }
   vector<thread> threads(chunks.size());
   for (unsigned i=0; i<chunks.size(); ++i) {
      threads[i] = thread(&Consensus::swallowThread, this, i, ref(chunks[i]), ref(undecided));
   }
   for_each(threads.begin(), threads.end(), mem_fn(&thread::join));
   // chunks should be shrink now, undecided should be larger
   // rebuild the segments from each chunk
   segs.clear();
   for (unsigned i=0; i<chunks.size(); ++i) {
      if (!chunks[i].empty()) {
         for (auto it=chunks[i].begin(); it != chunks[i].end(); ++it) {
            segs.push_back(move(*it));
         }
      }
   }
   cerr << "after swallow operation has " 
      << segs.size() << " unswallowed segements. "
      << undecided.size() << " undecided sgements\n";
}

// should try reverse complement if does not match
// in the forward direction
// Or should flip all of the segments
void Consensus::swallowThread(unsigned ai, list<DNAQualCount> &chnk,
      map<DNAQualCount*, vector<AlignQuality> > &seg2cons) {
   list<DNAQualCount>::iterator del, it;
   it=chnk.begin();
   aligner[ai]->setSeq1(*repseq); // set up repseq
   while (it != chnk.end()) {
      aligner[ai]->setSeq2(*it);
      aligner[ai]->runlocal();
      if (aligner[ai]->getIdentity()>0.99 && aligner[ai]->getCov2()>0.9) {
         unique_lock<mutex> lk(mut);
         DNAQualCount* newsq = addMember(*it);
         lk.unlock();
         accumulate(ai, newsq);
         del=it;
         ++it;
         chnk.erase(del);
      }
      else {
         AlignQuality Q = getAlignQuality(ai);
         if (Q.passCutoff(getIdentityCut())) {
            unique_lock<mutex> lk2(mut);
            DNAQualCount* newsq = add(*it); // for slavage operation
            lk2.unlock();
            unique_lock<mutex> extl(mutExternal);
            auto mit=seg2cons.find(newsq);
            if (mit == seg2cons.end()) {
               seg2cons.insert(make_pair(newsq, vector<AlignQuality>(1,Q)));
            }
            else {
               mit->second.push_back(Q);
            }
            extl.unlock();
            del=it; ++it;
            chnk.erase(del);
         }
         else ++it;
      }
   }
}

void Consensus::swallowSingletonThread(unsigned ai, list<DNAQualCount*> &sglchk,
      list<DNAQualCount> &unmatchedSeg, map<DNAQualCount*, vector<AlignQuality> > &undecided)
{
   // this consensus will work with a list of singletons
   list<DNAQualCount*>::iterator it, del;
   aligner[ai]->setSeq1(*repseq);
   it = sglchk.begin();
   int mergeCnt=0;
   int salvageCnt=0;
   while (it != sglchk.end()) {
      if (mergeKnownSequence(ai, *it, unmatchedSeg)) {
         del = it; ++it;
         sglchk.erase(del);
         ++mergeCnt;
      }
      else { 
         AlignQuality Q = getAlignQuality(ai);
         if (Q.passCutoff(getIdentityCut())) { // go to salvage pathway
            unique_lock<mutex> lk(mutExternal);
            auto mit = undecided.find(*it);
            if (mit == undecided.end()) {
               undecided.insert(make_pair(*it, vector<AlignQuality>(1,Q)));
            }
            else {
               mit->second.push_back(Q);
            }
            lk.unlock();
            del = it; ++it;
            sglchk.erase(del);
            ++salvageCnt;
         }
         else ++it; 
      }
   }
   // the following code is for debugging
   if (mergeCnt > 0) {
      unique_lock<mutex> lk(mutExternal);
      cerr << __FILE__ << ":" << __func__ << ":" << __LINE__ << " " 
         << getConsensusName() << " absorbed " << mergeCnt << " singletons\n";
      lk.unlock();
   }
   if (salvageCnt > 0) {
      unique_lock<mutex> lk(mutExternal);
      cerr << __FILE__ << ":" << __func__ << ":" << __LINE__ << " " 
         << getConsensusName() << " salvaged " << salvageCnt << " singletons\n";
      lk.unlock();
   }
}

void Consensus::swallowSingleton(list<DNAQualCount*> &sgl, 
      list<DNAQualCount> &seg, 
      map<DNAQualCount*, vector<AlignQuality> > &sgl_consRel) 
{
   list<DNAQualCount*>::iterator it, del;
   map<DNAQualCount*, vector<AlignQuality> >::iterator mit;
      // break up the task into smaller jobs
   vector<list<DNAQualCount*> > chunks(min(numberOfAligners(), (unsigned int)sgl.size()));
   unsigned int i=0;
   for (it = sgl.begin(); it != sgl.end(); ++it) {
      chunks[i++ % chunks.size()].push_back(*it);
   }
   vector<thread> threads(chunks.size());
   for (i=0; i<chunks.size(); ++i) {
      threads[i] = thread(&Consensus::swallowSingletonThread, this, i, 
            ref(chunks[i]), ref(seg), ref(sgl_consRel));
   }
   for_each(threads.begin(), threads.end(), mem_fn(&thread::join));
   // reconstruct un-used singletons
   sgl.clear();
   for (i=0; i<chunks.size(); ++i) {
      if (!chunks[i].empty()) {
         sgl.insert(sgl.end(), chunks[i].begin(), chunks[i].end());
      }
   }
}

// update total seq
DNAQualCount* Consensus::addMember(const DNAQualCount &c) {
   DNAQualCount* ns = add(c);
   member.push_back(ns);
   return ns;
}

// static functions 
list<DNAQualCount*> Consensus::separateSingleton(vector<vector<string> > &clusters) {
   cerr << clusters.size() << " input clusters\n";
   int i=0;
   while (clusters[i].size() > 1) ++i;
   int E=i;
   vector<string> single;
   while (i<clusters.size()) {
      single.push_back(clusters[i][0]);
      ++i;
   }
   clusters.resize(E);
   cerr << clusters.size() << " nonsingles "
      << single.size() << " singles\n";
   return seqstore->getSequencesList(single);
}

AlignQuality Consensus::getAlignQuality(unsigned alid) const {
   return AlignQuality(getConsensusName(),
         aligner[alid]->getIdentical(), aligner[alid]->getAlnlen(),
         aligner[alid]->getCov1(), aligner[alid]->getCov2());
}

void Consensus::initAligner(unsigned numthreads) {
   if (numthreads > thread::hardware_concurrency()) 
      numthreads = thread::hardware_concurrency();
   if (!aligner.empty()) {
      deallocateAligner();
      aligner.clear();
   }
   for (unsigned i=0; i<numthreads; ++i) {
      aligner.push_back(new Dynaln<SimpleScoreMethod>(scoreMethod));
   }
}

void Consensus::deallocateAligner() {
   for (unsigned i=0; i<aligner.size(); ++i) 
      delete aligner[i];
}

void Consensus::alignSequences(list<DNAQualCount> &seqs, ostream &ous) {
   list<DNAQualCount*> sptr;
   for (auto it=seqs.begin(); it != seqs.end(); ++it) {
      sptr.push_back(&(*it));
   }
   alignSequences(sptr, ous);
}

void Consensus::alignSequences(const list<DNAQualCount*> &seqs, ostream &ous)
{
   int b=0;
   vector<list<DNAQualCount*> > chunk(min(numberOfAligners(), (unsigned)seqs.size()));
   for (auto it=seqs.begin(); it != seqs.end(); ++it) {
      chunk[b++ % chunk.size()].push_back(*it);
   }
   vector<thread> threads(chunk.size());
   for (int i=0; i<chunk.size(); ++i) {
      threads[i] = thread(&Consensus::alignSequencesThread, this, i, ref(chunk[i]), ref(ous));
   }
   for_each(threads.begin(), threads.end(), mem_fn(&thread::join));
}
      
void Consensus::alignSequencesThread(unsigned ai, const list<DNAQualCount*> &chk, ostream &ous) 
{
   aligner[ai]->setSeq1(getResult());
   for (auto it=chk.begin(); it != chk.end(); ++it) {
      aligner[ai]->setSeq2(**it);
      aligner[ai]->runlocal();
      lock_guard<mutex> lk(mutExternal);
      aligner[ai]->printAlign(ous);
   }
}

//////////////////// helper function for post processing of consensus ////////
void writeConsensusFasta(const vector<DNAQualCount> &con, ostream &ouf) {
   for (int i=0; i<con.size(); ++i) {
      con[i].printFastaWithHeader(ouf);
   }
}

void writeConsensusFasta(const vector<DNAQualCount> &con, const string &fasFile) {
   ofstream ouf(fasFile);
   for (int i=0; i<con.size(); ++i) {
      con[i].printFastaWithHeader(ouf);
   }
}

// for swarm
vector<vector<string> > readCluster(const string &clfile) {
   ifstream inf(clfile.c_str());
   if (inf.fail()) {
      cerr << "Failed to open input file " << clfile << endl;
      exit(1);
   }
   vector<vector<string> > result;
   string line;
   getline(inf, line);
   while (!inf.eof())  {
      //vector<string> row=split(line, ' ');
      vector<string> row=DNAQualstore::stripId(split(line, ' '));
      result.push_back(row);
      getline(inf, line);
   }
   sort(result.begin(), result.end(), sortVectorBySize());
   return result;
}

// mothur specific
vector<vector<string> > readCluster(const string &clfile, float cutoff) {
   ifstream inf(clfile.c_str());
   if (inf.fail()) {
      throw runtime_error("Failed to open input file " + clfile);
   }
   vector<vector<string> > result;
   string line;
   getline(inf, line); // header
   getline(inf, line); // unique line
   getline(inf, line); // first data line
   vector<string> lastrow, row;
   while (!inf.eof()) {
      row = split(line, '\t');
      float dist = stof(row[0]);
      if (dist == cutoff) break;
      if (dist > cutoff) {
         row=lastrow; break;
      }
      lastrow=row;
      getline(inf, line);
   }
   // erase the first two elements
   if (!row.empty()) {
      cerr << row.size()-2 << " clusters\n";
      for (size_t i=2; i<row.size(); ++i) {
         result.push_back(DNAQualstore::stripId(split(row[i], ',')));
      }
      if (result.empty()) {
         cerr << "WARN: There is no cluster found at " << cutoff << " cutoff!\n";
      }
      else {
         sort(result.begin(), result.end(), sortVectorBySize());
      }
   }
   return result;
}

/// for cluster program parameter class
//string ClusterProgparam::bindir = "/remote/RSU/sw-cache/metag/bin";
string ClusterProgparam::bindir = "/isilon/Apps/inf_dis/metag/bin";
//string ClusterProgparam::datadir = "/home/kzhou/work/metag/refseq";
string ClusterProgparam::datadir = "/isilon/Data/inf_dis/metag/refseq";
void ClusterProgparam::setBinaryDirectory() {
   char *metagHome = getenv("METAG_HOME");
   if (metagHome == NULL) {
      metagHome=getenv("HOME");
   }
   bindir = string(metagHome) + "/bin";
   char* metagData = getenv("METAG_DATA");
   if (metagData != NULL) {
      datadir=string(metagData) + "/refseq";
   }
}
///////////////////////////////////////

string rebuildConsensus(vector<vector<string> > &clid, const ClusterProgparam &par) {
   list<DNAQualCount*> singleton = Consensus::separateSingleton(clid);
   cerr << clid.size() << " cluster with more than one member\n";
   cerr << "doing consensus building ...\n";
   list<Consensus*> conslist;
   for (int i=0; i<clid.size(); ++i) {
      Consensus *cons = new Consensus(clid[i]);
      (*cons)();
      conslist.push_back(cons);
   }
   cerr << conslist.size() << " consensi build from clusters\n";
   list<DNAQualCount> scrapSeg = combineConsensus(conslist, par.salvageSingle());
   cerr << scrapSeg.size() << " scraps after consensus merge\n";
   writeScrap(scrapSeg, "consensus", par);
   list<DNAQualCount> scrapOfSingle = consensusSwallowSingleton(conslist, singleton);
   cerr << scrapOfSingle.size() << " scraps from singleton after merge\n";
   writeScrap(scrapOfSingle, "single", par);
   string consensusFasFile = writeConsensus(conslist, par.getFstem());
   writeSequencePointer(singleton, "badsingle", par); // singletons not absorbed
   showConsensusSingleAlign(conslist, singleton, par);
   showConsensusScrapAlign(conslist, scrapOfSingle, par);
   // deallocate memory
   // testing the second round of conseus
   // put graphics results in a separate directory
   StalactiteGraph stalac(PNG);
   const int dir_err = mkdir("caveplots", S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
   if (dir_err == -1) {
      cerr << "Failed to create caveplots directory\n";
      exit(1);
   }
   stalac.setOutputDir("caveplots");
   for (auto it=conslist.begin(); it != conslist.end(); ++it) {
      stalac.plot((*it)->getCoverage(), (*it)->getConsensusName());
      //(*it)->refine(); //TODO: in the process of developing
      delete (*it);
   }
   return consensusFasFile;
}

string buildConsensusAll(vector<vector<string> > &clid, const ClusterProgparam &par) {
   list<DNAQualCount*> singleton = Consensus::separateSingleton(clid);
   cerr << clid.size() << " cluster with more than one member\n";
   cerr << "Generating consensus from clusters with >1 members ...\n";
   list<Consensus*> conslist;
   for (int i=0; i<clid.size(); ++i) {
      Consensus *cons = new Consensus(clid[i]);
      (*cons)();
      conslist.push_back(cons);
   }
   cerr << conslist.size() << " consensi build from clusters\n";
   string consensusFasFile = writeConsensus(conslist, par.getFstem());
   writeSequencePointer(singleton, "singleton", par); 
   for (auto it=conslist.begin(); it != conslist.end(); ++it) {
      delete (*it);
   }
   return consensusFasFile;
}

list<DNAQualCount> combineConsensus(list<Consensus*> &allcons, bool salvage) {
   cerr << "combining similar clusters ...\n";
   allcons.sort(sortConsensusByTotalSequence());
   list<Consensus*>::iterator delit, it, it2;
   it = allcons.begin();
   list<DNAQualCount> unmatchedSeg;
   unsigned numthread=Consensus::numberOfAligners();
   int round=1;
   while (it != allcons.end()) {
      //cerr << "round " << round++ << " working on consensus " 
      //   << (*it)->getConsensusName() << " clsize: " << (*it)->getTotal()
      //   << " "
      //   << distance(it, allcons.end()) << " smaller clusters" << endl;
      it2=it; ++it2;
      if (it2 == allcons.end()) break;
      // break job into smaller chunks
      vector<list<Consensus*> > chunk(min(numthread, (unsigned)allcons.size()-1));
      int b = 0;
      while (it2 != allcons.end()) { // only if size difference > 3x
         //cerr << "right cluster " << (*it2)->getConsensusName()
         //   << " size: " << " " << (*it2)->getTotal() << endl;
         if ((*it)->getTotal() > (*it2)->getTotal()*3 || 
               ((*it)->getTotal() > (*it2)->getTotal() && (*it2)->getTotal() < 6)) {
            chunk[b++ % chunk.size()].push_back(*it2);
            //cerr << "candidate for merging\n";
         }
         ++it2;
      }
      //cerr << chunk.size() << " chunks\n";
      // make sure all chunks are not empty
      for (unsigned i=0; i<chunk.size(); ++i) {
         if (chunk[i].empty()) {
            //cerr << "inside combineConsensus: got empty chunk, need to shinkg vector\n";
            chunk.resize(i);
            break;
         }
      }
      // each thread reduces the chunk only, not touching the allcons!
      vector<thread> threads(chunk.size());
      set<Consensus*> mergedCons;
      for (unsigned i=0; i < chunk.size(); ++i) {
         threads[i] = thread(&Consensus::merge, *it, i, ref(chunk[i]), 
               ref(unmatchedSeg), ref(mergedCons));
      }
      for_each(threads.begin(), threads.end(), mem_fn(&thread::join));
      // remove those that merged
      if (!mergedCons.empty()) {
         //cerr << mergedCons.size() << " smaller clusters merged\n";
         it2=it; ++it2;
         while (it2 != allcons.end()) {
            if (mergedCons.find(*it2) != mergedCons.end()) {
               delit = it2;
               ++it2;
               delete (*delit);
               allcons.erase(delit);
            }
            else {
               ++it2;
            }
         }
      }
      ++it;
   }
   cerr << allcons.size() << " clusters at end of merge\n";
   cerr << unmatchedSeg.size() << " scrap fragments\n";
   // swallow the segments from previous operation
   // in rare cases segment may be in the reverse orientation
   map<DNAQualCount*, vector<AlignQuality> > seg_consRel;
   for (it=allcons.begin(); it != allcons.end(); ++it) {
      if (unmatchedSeg.empty()) {
         cerr << (*it)->getConsensusName() << " has no more segment to swallow\n";
         break;
      }
      (*it)->swallow(unmatchedSeg, seg_consRel);
   }
   if (!unmatchedSeg.empty()) {
      cerr << unmatchedSeg.size() << " scraps left after swallow operation\n";
   }
   if (salvage && !seg_consRel.empty()) {
      salvageLQSequence(allcons, seg_consRel);
      cerr << seg_consRel.size() << " consensus segments salvaged\n";
   }
   else {
      cerr << "no salvage to do\n";
   }
   /*
   cout << "Validate totalseq after combineConsensus operation\n";
   for (it=allcons.begin(); it != allcons.end(); ++it) {
      (*it)->debugValidTotal();
   }
   */
   return unmatchedSeg;
}

void writeScrap(const list<DNAQualCount> &scrap, const string &suffix, const ClusterProgparam &par) {
   if (scrap.empty()) {
      cerr << "there is no scrap to work with\n";
      return;
   }
   string scrapFile=par.getFstem() + suffix + ".scrap.fas";
   ofstream ouf(scrapFile);
   if (ouf.fail()) {
      throw runtime_error("failed to open file: " + scrapFile);
   }
   for (auto it=scrap.begin(); it != scrap.end(); ++it) {
      it->printFastaWithHeader(ouf);
   }
   cerr << scrap.size() << " scrap written to file: " << scrapFile << endl;
}

// used for the first round of selection. Used higher cutoff
// the left over sgl need to go through another round of 
// selection.
list<DNAQualCount> consensusSwallowSingleton(list<Consensus*> &cons, 
      list<DNAQualCount*> &sgl) 
{
   cerr << __FILE__ << ":" << __LINE__
      << cons.size() << " consensus swallowing " << sgl.size() << " singletons\n";
   list<DNAQualCount> unmatchedSeg;
   map<DNAQualCount*, vector<AlignQuality> > sgl_consRelation;
   for (auto it=cons.begin(); it != cons.end(); ++it) {
      (*it)->swallowSingleton(sgl, unmatchedSeg, sgl_consRelation);
   }
   cerr << sgl.size() << " singltons after merging with consensus\n";
   cerr << unmatchedSeg.size() << " scrap fragments\n";
   // swallow the segments from above operation
   map<DNAQualCount*, vector<AlignQuality> > seg_consRelation;
   for (auto it=cons.begin(); it != cons.end(); ++it) {
      (*it)->swallow(unmatchedSeg, seg_consRelation);
   }
   cerr << unmatchedSeg.size() << " scraps after Cluster swallow single operation\n";
   cerr << "trying to salvate singleton from " << sgl_consRelation.size() << " relationships\n";
   salvageLQSequence(cons, sgl_consRelation);
   cerr << "trying to salvage singleton segments from " << seg_consRelation.size() << " relationships\n";
   salvageLQSequence(cons, seg_consRelation);
   /*
    * Because the totalseq is computed when it is requested
    * The validation is nolonger needed.
   cout << "VVVVV Validate totalseq after consensus swallow singleton\n";
   for (auto it=cons.begin(); it != cons.end(); ++it) {
      (*it)->debugValidTotal();
   }
   */
   return unmatchedSeg;
}

string writeConsensus(const list<Consensus*>& cons, const string &fstem) {
   string conFile, baseFile, qualFile, alninfoFile, covFile;
   conFile = fstem + ".consensus.txt";
   baseFile = fstem + ".base.txt";
   qualFile = fstem + ".qual.txt";
   alninfoFile = fstem + ".alninfo.txt";
   covFile = fstem + ".vertcov.txt";
   ofstream ocon(conFile);
   ofstream obas(baseFile);
   ofstream oqua(qualFile);
   ofstream oalni(alninfoFile);
   ofstream ocov(covFile);
   string conFasFile = fstem + "consensus.fas";
   ofstream oconfas(conFasFile);
   cerr << "writing cluster consensus to file ...\n";
   for (auto it=cons.begin(); it != cons.end(); ++it) {
      (*it)->printResult(ocon, obas, oqua, oalni, ocov);
      (*it)->printFasta(oconfas);
   }
   cerr << cons.size() << " consensus written to multiple files\n";
   // send information to the caller in the shell
   cout << "Consensus fasta file: " << conFasFile << endl;
   return conFasFile;
}

void writeSequencePointer(const list<DNAQualCount*> &sgl, 
      const string &suffix, const ClusterProgparam &par) 
{
   string fasfile= par.getFstem() + suffix + ".fas";
   ofstream ops(fasfile);
   if (ops.fail()) {
      throw runtime_error("cannot open file: " + fasfile);
   }
   for (auto it=sgl.begin(); it != sgl.end(); ++it) {
      (*it)->printFastaWithHeader(ops);
   }
   cerr << sgl.size() << " " << suffix << " written to file: " << fasfile << endl;
}

// for debug
void showConsensusSingleAlign(const list<Consensus*> &cons, 
      const list<DNAQualCount*> &single, const ClusterProgparam &par) {
   string file= par.getFstem() + "consensus_badsingleton.aln";
   ofstream ouf(file);
   if (ouf.fail()) {
      throw runtime_error("cannot open " + file);
   }
   ouf << "This file is for debug and improvement of the pipeline\n";
   for (auto it=cons.begin(); it != cons.end(); ++it) {
      ouf << "Consensus: " << (*it)->getResult() << endl;
      (*it)->alignSequences(single, ouf);
   }
   cerr << "Consensus to unabsorbed singleton alignments written to "
      << file << endl;
}

// now only showing the scrap derived from singleton
void showConsensusScrapAlign(const list<Consensus*> &cons, 
      list<DNAQualCount> &seq, const ClusterProgparam &par) {
   string file= par.getFstem() + "consensus_scrap.aln";
   ofstream ouf(file);
   if (ouf.fail()) {
      throw runtime_error("cannot open " + file);
   }
   ouf << "This file is for debug and improvement of the pipeline\n";
   for (auto it=cons.begin(); it != cons.end(); ++it) {
      ouf << "Consensus: " << (*it)->getResult() << endl;
      (*it)->alignSequences(seq, ouf);
   }
   cerr << "Consensus to scrap alignments written to " << file << endl;
}

void salvageLQSequence(list<Consensus*> &cons, 
      map<DNAQualCount*, vector<AlignQuality> > &sgl2cons) 
{
   cerr << __FILE__ << ":" << __func__ << " " << __LINE__
      << " salvaging " << sgl2cons.size() << " sequences\n";
   // build a look up map from consensusName to consensus
   map<string, Consensus*> consusdic;
   for (auto i=cons.begin(); i != cons.end(); i++) {
      consusdic.insert(make_pair((*i)->getConsensusName(), (*i)));
   }
   // seq already aligned to consensus and filtered before
   // seq => consensus list
   unsigned numthreads=Consensus::numberOfAligners();
   map<DNAQualCount*, vector<AlignQuality> > sgl_consRelation;
   auto it=sgl2cons.begin();
   while (it != sgl2cons.end()) {
      int b=0;
      vector<thread> threads;
      // spawn new threads
      while (it != sgl2cons.end() && b < numthreads) {
         threads.push_back(thread(&Consensus::salvageKnownSequence, consusdic.find(it->second.back().getConsensusName())->second, b, it->first));
         ++it; ++b;
      }
      for_each(threads.begin(), threads.end(), mem_fn(&thread::join));
      // wait for all threads to finish
   }
}

