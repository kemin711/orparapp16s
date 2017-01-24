#include "blastnrow.h"

////// Blastnrow calss ////////
// start with a dictionary
// This needs to be updated in the future
map<string, pair<string,string> > Blastnrow::taxdic={};
vector<string> Blastnrow::header={"clusterid", "clsize", "cldepth", "hitid", "species", "strain", "percent_ident", "consensesu_len", "hit_len", "align_len", "qstart", "qend", "sstart", "send", "bitscore"};

void Blastnrow::writeHeader(ostream &ouf) {
   auto it = header.begin();
   ouf << *it;
   ++it;
   while (it != header.end()) {
      ouf << "\t" << *it;
      ++it;
   }  
}

Blastnrow::Blastnrow(const vector<string> &row)
   : qseqid(row[0]), sseqid(row[1]),
     pident(stof(row[2])), 
     qlen(stoi(row[3])), slen(stoi(row[4])), length(stoi(row[5])),
     qstart(stoi(row[6])), qend(stoi(row[7])),
     sstart(stoi(row[8])), send(stoi(row[9])),
     bitscore(stoi(row[10])),
     taxon(), queryCount(0), consDepth(0)
{
   setTaxon();
}

void Blastnrow::setTaxon() {
   map<string, pair<string,string> >::const_iterator mit=taxdic.find(sseqid);
   if (mit == taxdic.end()) {
      cerr << "failed to find " << sseqid << " in tax dictioanry\n";
      throw runtime_error(sseqid + " not in taxdictionary");
   }
   taxon=mit->second;
}

ostream& operator<<(ostream& ous, const Blastnrow& r) {
   static char sep='\t';
   ous << r.qseqid << sep << r.queryCount << sep << r.getConsensusDepth() << sep
      << r.sseqid << sep << r.taxon.first << sep;
   if (r.taxon.second.empty()) {
      ous << "\\N"; // for mysql loading
   }
   else ous << r.taxon.second;
   ous << sep
      << r.pident << sep << r.qlen << sep
      << r.slen << sep << r.length << sep
      << r.qstart << sep << r.qend << sep
      << r.sstart << sep << r.send << sep << r.bitscore;
   return ous;
}
bool Blastnrow::operator<(const Blastnrow& r) const {
   if (getQidNumber() < r.getQidNumber()) {
      return true;
   }
   if (getQidNumber() > r.getQidNumber()) {
      return false;
   }
   // same id
   if (getBitscore() > r.getBitscore()) {
      return true;
   }
   return false;
}

tuple<string, int, string> Blastnrow::getSpeciesDepth() const {
   if (taxon.second.empty())
      return make_tuple(getQuery(), getConsensusDepth(), taxon.first);
   return make_tuple(getQuery(), getConsensusDepth(), taxon.first + " " + taxon.second);
}

vector<string> Blastnrow::toStringVector() {
   vector<string> res;
   res.push_back(qseqid); res.push_back(to_string(queryCount)); res.push_back(to_string(consDepth));
   res.push_back(sseqid); res.push_back(taxon.first); res.push_back(taxon.second);
   res.push_back(ftos(pident, 3));
   res.push_back(to_string(qlen));
   res.push_back(to_string(slen));
   res.push_back(to_string(length));
   res.push_back(to_string(qstart));
   res.push_back(to_string(qend));
   res.push_back(to_string(sstart));
   res.push_back(to_string(send));
   res.push_back(to_string(bitscore));
   return res;
}

