#ifndef _SPLICING_GRAPH_HPP_
#define _SPLICING_GRAPH_HPP_

#include <iostream>
#include <fstream>
#include <list>
#include <utility>
#include <string>
#include <map>
#include <unordered_map>
#include <algorithm>
#include <regex>
#include <set>

#include "sdsl/bit_vectors.hpp"
#include "utils.hpp"

// struttura per la definizione di una feature
struct Feature {
  // indice per muoversi sui campi di una linea del file gtf 
  int i;
  // bho
  int index;

  // campi del file gtf
  std::string seqid;
  std::string type;
  std::string featString;
  int start;
  int end;
  bool strand; //+: true; -:false
  std::string id;
  std::string transcript_id;
  // costruttore default 
  Feature() {
    i = 1;
  }

  // costruttore che prende in input una stringa
  Feature(std::string line) {
    featString = line;
    i = 1;
    std::string token;

    // nel file gtf ho i campi divisi dai tab
    std::string delimiter = "\t";
    std::size_t pos;
    
    while ((pos = line.find(delimiter)) != std::string::npos) {
      // estrapolo i singoli token della riga, elimino quella parte dalla 
      // stringa e effettuto l'add alla struttura
      token = line.substr(0, pos);
      line.erase(0, pos + delimiter.length());
      add(token);
    }
    // aggiungo anche l'ultima parte contenente vari ID
    add(line);
  }

  // aggiunta dei singoli token
  void add(std::string elem) {
    // uso un contatior i che viene aumentato per sapere che campo sto
    // analizzando 
    switch(i) {
      // primo caso: id della sequenza
    case 1:
      seqid = elem;
      break;
      
      // secondo caso ignorato (provenienza dato ?)
      // terzo caso: tipo (exon, cds, etc...)
    case 3:
      type = elem;
      break;
      
      // quarto caso: indice di partenza
    case 4:
      start = std::stoi(elem);
      break;
      
      // quinto caso: indice di fine
    case 5:
      end = std::stoi(elem);
      break;
      
      // sesto caso saltato 
      // settimo caso: direzione strand in bool: "+" -> true, "-" -> false
    case 7:
      if(elem.compare("+") == 0) {
	strand = true;
      } else {
	strand = false;
      }
      break;

      // ottavo caso saltato
      // nono caso: vari id 
    case 9:
      bool flag = true;
      //GFF
      //std::string delimiter = ";";
      //std::string string_to_search = "ID=";

      //GTF
      std::string delimiter = "; ";
      std::string string_to_search;


      // a seconda del tipo si cerca l'id giusto 
      if(type.compare("gene") == 0) {
	string_to_search = "gene_id \"";
      } else if(type.compare("transcript") == 0) {
	string_to_search = "transcript_id \"";
      } else if(type.compare("exon") == 0) {
	string_to_search = "exon_id \"";
      }

      // si procede cercando nella stringa, grazie al delimitatore ";",
      // e rimuovendo quanto non necessario
      std::size_t pos;
      std::string token;
      std::string subtoken;
      while((pos = elem.find(delimiter)) != std::string::npos) {
	token = elem.substr(0, pos);
	if(token.substr(0, string_to_search.size()).compare(string_to_search)
	   == 0) {
	  //GTF
	  id = token.substr(string_to_search.size(),
			    token.size() - string_to_search.size() - 1);
	  //GFF
	  //id = token.substr(string_to_search.size(),token.size()
	  //     - string_to_search.size());
	  flag = false;
	  break;
	}
	elem.erase(0, pos + delimiter.length());
      }
      if(flag &&
	 elem.substr(0,string_to_search.size()).compare(string_to_search)
	 == 0) {
	
	//GTF
	id = elem.substr(string_to_search.size(),
			 elem.size() - string_to_search.size() - 2);
	//GFF
	//id = elem.substr(string_to_search.size(),elem.size()
	//     - string_to_search.size()-1);
	flag = false;
      }
      if(flag) {
	id = ".";
      }
    }
    i++;
  }

  std::string toStr() {
    return seqid + " " +
      type + " " +
      std::to_string(start) + " " +
      std::to_string(end) + " " +
      std::to_string(strand) + " " +
      id;
  }

  std::string getTId(){
    std::regex transcriptRegex("(transcript_id \"(.*?)\";)");
    std::smatch match;
    std::regex_search(featString, match, transcriptRegex);
    return match[2];
  }
};

// definizione della classe per lo splicing graph
class SplicingGraph {
private:
  std::string T;
  std::string TIntrons;
  std::string reference;
  int refLen;
  bool strand;
  std::vector<std::pair<int, int> > ExonsPos;
  std::vector<std::list<int> > parents;
  std::vector<std::list<int> > sons;
  int exsN;
  std::vector<std::list<int> > parentsIntron;
  std::vector<std::list<int> > sonsIntron;
  
  std::map<std::string, std::vector<std::pair<unsigned long long int, unsigned long long int> > > transcriptExons;
  std::map<std::string, int> id2index;

  std::set<std::pair<unsigned long long int, unsigned long long int> > IntronPos;
  std::map<std::pair<unsigned long long int, unsigned long long int>, int > intronsId;
  std::map<std::pair<unsigned long long int, unsigned long long int>,
   	   std::set<std::pair<unsigned long long int, unsigned long long int> > > idIntronsList;
  std::map<std::pair<unsigned long long int, unsigned long long int>,
   	   std::vector<unsigned long long int> > idIntronsMap;
  
  std::vector<std::vector<int> > edges;
  sdsl::rrr_vector<> bitVector;
  sdsl::rrr_vector<>::select_1_type selectBV;
  sdsl::rrr_vector<>::rank_1_type rankBV;

  sdsl::rrr_vector<> bitVectorIntrons;
  sdsl::rrr_vector<>::select_1_type selectBVIntrons;
  sdsl::rrr_vector<>::rank_1_type rankBVIntrons;

  void setupBitVector();
  void setupBitVectorIntrons();

  void save(const std::string);
  void saveIntron(const std::string);
  bool verbose;
  //void load(const std::string);
public:
  //SplicingGraph(const std::string&);
  SplicingGraph(const char*, const std::string&, const bool&, const bool&);
  std::string getText() const;
  std::string getTextIntrons() const;
  
  std::list<int> getParents(const int&) const;
  std::list<int> getSons(const int&) const;
  std::list<int> getParentsIntron(const int&) const;
  std::list<int> getSonsIntron(const int&) const;
  std::string getExon(const int&) const;
  std::string getIntron(const int&) const;

  int rank(const int&) const;
  int select(const int&) const;
  int rankIntrons(const int&) const;
  int selectIntrons(const int&) const;
  
  bool contain(const int&, const int&) const;
  bool isNew(const int&, const int&) const;
  void print() const;
  int getExonsNumber() const;
  
  // std::map<std::pair<int, int>, std::pair<int,int> > getIntronsPos() const;
  // std::map<std::pair<int, int>, int > getIntronsId() const;
  std::map<std::pair<unsigned long long int, unsigned long long int>, std::vector< unsigned long long int> > getIntronMap() const;
  // std::vector<std::pair<std::string, std::string> >
  // getIntrons(std::pair<int, int>, std::map<std::pair<std::string, std::string>,
  // std::vector<std::pair<int,int> > >) const;
};

#endif
