#include "SplicingGraph.hpp"

// da coppia di indici restuisce la stringa <id1>:<id2>
std::string getExonID(unsigned long long int s, unsigned long long int e) {
  return std::to_string(s) + ":" + std::to_string(e);
}

bool sortinrev(const std::pair<unsigned long long int, unsigned long long int> &a,  
               const std::pair<unsigned long long int, unsigned long long int> &b) { 
       return (a.first > b.first); 
} 
// costruttore del grafo, prendo in input la stringa completa del genoma e il
// nome del file gtf con le annotazioni 

SplicingGraph::SplicingGraph(const char *genomic,
                             const std::string &gtf,
			     const bool& verbose,
			     const bool& intron_calc) {

  // lunghezza del gene
  refLen = strlen(genomic);

  // stream in ingresso del file gtf
  std::ifstream gtfFile;

  // inizializzione matrice adiacenza con 0
  edges.push_back({0});
  edges.push_back({0});

  // inizizzione liste dei parents e dei sons di un nodo
  parents.push_back({});
  sons.push_back({});

  // inizializzo linearizzazione che inizia col delimitatore 
  T = "|";

  // inizializzazione id esone di partenza
  int exID = 1;
  

  // inizializzazione indici
  int curr_i = 1;
  int last_i = -1;

  // inizializzazione stringa che riceverà in input ogni singola riga
  // del file gtf
  std::string line;

  // apertura del file 
  gtfFile.open(gtf);
  
  if(gtfFile.is_open()) {
    // per ogni riga del file gtf si studia il necessario
    while(getline(gtfFile,line)) {
      // costruisco la feature feat
      Feature feat (line);

      // se si ha una feature gene aggiorno reference e strand di conseguenza
      if(feat.type.compare("gene") == 0) {
	reference = feat.seqid;
	strand = feat.strand;
      }
      // se si ha una feature transcript o mRNA si "azzera" l'indice di fine
      else if(feat.type.compare("transcript") == 0 ||
	      feat.type.compare("mRNA") == 0) {
	last_i = -1;
      }
      // se si ha una feature exon
      else if(feat.type.compare("exon") == 0) {
	// credo la stringa con "[a:b]" con a indice iniziale esone e b finale
	std::string posID = getExonID(feat.start, feat.end);
	transcriptExons[feat.getTId()].push_back(std::make_pair(feat.start,
								feat.end));
	
	//std::cout << feat.getTranscript() << " \n";
	// se l'elemento nuovo è l'ultimo della mappa (quindi è nuovo ???)
	if (id2index.find(posID) == id2index.end()) {
	  // assegno l'indice dell'esone in considerazione
	  id2index[posID] = exID;
	  
	  //genero la stringa corrispondente alle posizioni dell'esone
	  std::string currExString(genomic + feat.start - 1,
				   feat.end - feat.start + 1);
	  // std::cout <<  feat.start<< " " << feat.end << ": "<<
	  // currExString << std::endl;
	  // aggiorno la sequenza esonica
	  T += currExString + "|";

	  // aggiungo la coppia inizio-fine al vettore delle coppie esoniche
	  ExonsPos.push_back(std::make_pair(feat.start, feat.end));

	  // è il primo elemento, parents e sons vengono solo inizializzati
	  parents.push_back({});
	  sons.push_back({});

	  // si incrementa il contatore degli id
	  ++exID;
	}

	// si aggiorna l'indice di partenza
	curr_i = id2index[posID];
	
	// procedo sse si è analizzato almeno un elemento e sse non si sta
	// analizzando due volte di fila lo stesso elemento (?)
	if(last_i != -1) {
	  if(last_i != curr_i) {

	    // se ho un last_i maggiore del numero di archi 
	    // si ridimensiona (?)
	    if(last_i >= (int)edges.size()) {
	      int i = edges.size();
	      while(i <= last_i + 1) {
		edges.push_back({0});
		++i;
	      }
	    }

	    // se ho un curr_i maggiore del numero di vertici di adiacenza
	    // dell'ultimo si ridimensiona tale lista (?)
	    if(curr_i >= (int)edges[last_i].size()) {
	      int i = edges[last_i].size();
	      while(i<=curr_i+1) {
		edges[last_i].push_back(0);
		++i;
	      }
	    }
	    //edges[last_i] vector<vector<int>> rightValidVariants(1, {i});

	    // si mette un arco tra l'ultimo vertice e quello appena aggiunto
	    edges[last_i][curr_i] = 1;

	    // si ridimensiona la lista dei parents
	    if(curr_i >= (int)parents.size()) {
	      int i = 0;
	      while(i <= curr_i + 1) {
		parents.push_back({});
		++i;
	      }
	    }

	    // alla lista dei parents del vertice appena aggiunto si mette il
	    // vecchio ultimo vertice
	    parents[curr_i].push_back(last_i);

	    // si ridimensiona la lista dei sons
	    if(last_i >= (int)sons.size()) {
	      int i = 0;
	      while(i <= last_i + 1) {
		sons.push_back({});
		++i;
	      }
	    }
	    // alla lista dei sons dell'ultimo vertice si aggiunge quello
	    // corrente
	    sons[last_i].push_back(curr_i);
	  }
	}
	// l'ultimo vertice diventa quello corrente
	last_i = curr_i;
      }
    }
  }

  // Transitive closure on the graph
  // si attacca ogni vertice a quello seguente
  int i = 1;
  // si itera su ogni coppia nel vector di std::pair
  for(const std::pair<int,int>& p1 : ExonsPos) {

    // se col contatore si supera il numero di vertici (grandezza colonna
    // matrice di adiacenza) ridimensiono
    if(i >= (int)edges.size()) {
      int i_ = edges.size();
      while(i_ <= i+1) {
	edges.push_back({0});
	++i_;
      }
    }
    
    int j = 1;
    // si cicla nuovamente sui vertici (ciclando sulle varie pair esoniche)
    for(const std::pair<int,int>& p2 : ExonsPos) {
      
      // se questa nuova coppia è successiva a quella del primo for
      if(p1.second <= p2.first) {
	
	// se necessario si ridimensiona
	if(j >= (int)edges[i].size()) {
	  int j_ = edges[i].size();
	  while(j_ <= j+1) {
	    edges[i].push_back(0);
	    ++j_;
	  }
	}

	// se già non si ha un nodo lo si crea mettendo però 2 nella matrice
	// di adiacenza per segnalare il tipo diverso di arco
	// e si aggiornano parents e sons
	if(edges[i][j] == 0) {
	  edges[i][j] = 2;
	  parents[j].push_back(i);
	  sons[i].push_back(j);
	}
      }
       
      ++j;
    }
    ++i;
  }
  if(verbose){
    int count = 0;
    std::cout << "parents:\n";
    for(auto i: parents){
      std::cout << count << ": ";
      for(auto o: i){
	std::cout << o << " ";
      }
      count++;
      std::cout << "\n";
    }
  }
   if(verbose){
    int count = 0;
    std::cout << "sons:\n";
    for(auto i: sons){
      std::cout << count << ": ";
      for(auto o: i){
	std::cout << o << " ";
      }
      count++;
      std::cout << "\n";
    }
  }
  // set contenente, in modo unico, tutti possibili introni su ogni trascritto
  // un esone presente in un trascritto ma non nell'altro è considerato
  // possibile introne
  std::set<std::pair<unsigned long long int, unsigned long long int> > possible_introns;
  // mappa che associa a due esoni gli introni tra essi contenuti
  std::map<std::pair<std::string, std::string>,
   	   std::vector<std::pair<unsigned long long int,
				 unsigned long long int> > > introns;
  std::map<std::string, std::vector<std::pair<unsigned long long int,
					      unsigned long long int>
				    > > transcriptExonsNew;

  if(intron_calc){
    // inizializzo la linerizzazione degli introni
    TIntrons = "|";
    // creo la lista dei possibili introni considerando anche gli esoni di altri
    // trascritti, inserendo in possible_introns
    for(auto t: transcriptExons){
      std::sort(t.second.begin(), t.second.end(), sortinrev);
      transcriptExonsNew[t.first] = t.second;
    }
    // sistemo con elementi ordinati
    transcriptExons = transcriptExonsNew;

    for(auto t: transcriptExons){
      for(int i = t.second.size() - 1; i >= 0; i--){
	possible_introns.insert(std::make_pair(t.second[i].first,
					       t.second[i].second));
	if(i != 0 )
	  possible_introns.insert(std::make_pair(t.second[i].second,
						 t.second[i - 1].first));
      }
    }
    // carico la mappa che associa id esoni introni
    for(auto t: transcriptExons){
      for(unsigned long i =  t.second.size() - 1 ; i > 0; i--){

	// genero la stringa formata dagli indici dei due esoni consecutivi
	// in un trascritto
	std::string exons = getExonID(t.second[i].first, t.second[i].second)
	  + " " + getExonID(t.second[i - 1].first, t.second[i - 1].second);
	// carico mappa introni usando indici per transcritto
	for(auto intr: possible_introns){
	  if(intr.first >= t.second[i].second &&
	     intr.second <= t.second[i - 1].first){

	    introns[std::make_pair(t.first, exons)].push_back(intr);
	  }
	}
	// rimuovo la distrinzione dei trascritti creando una nuova mappa
	// salvo in un set gli intronin (come coppia di indici)
      	auto exon1 = getExonID(t.second[i].first, t.second[i].second);
	auto exon2 = getExonID(t.second[i - 1].first, t.second[i - 1].second);
	for(auto o: introns[std::make_pair(t.first, exons)]){
	  
	  idIntronsList[std::make_pair(id2index[exon1],
				       id2index[exon2])].insert(o);
	  IntronPos.insert(o);
	}
            
      }
    }

    // cosruisco la linearizzazione e associo un id agli introni
    int idI = 1;
    for(auto i: IntronPos){
      //std::cout << i.first << " "<<i.second <<  "\n";
      std::string intronCurr(genomic + i.first - 1,
			     i.second - i.first + 1);
      //std::cout << intronCurr << "\n";
      if(verbose){
	std::cout << i.first << " " << i.second << ": "
		  << intronCurr << std::endl;
      }
      TIntrons += intronCurr + "|";
      intronsId[std::make_pair(i.first, i.second)] = idI;
      idI++;
    }

    // costruisco la mappa finale con gli id
    for(auto i: idIntronsList){
      for(auto o: i.second)
	idIntronsMap[i.first].push_back(intronsId[o]);
    }
  }
  // stampa informazioni
  if(verbose){
    std::cout << "exons data:\n";
    int j = 1; 
    for(auto i: ExonsPos){
      std::cout << j <<" "<< i.first<< " "<< i.second << ":\n";
      std::string intronCurr(genomic + i.first - 1,
			     i.second - i.first + 1);

      std::cout << intronCurr << std::endl;
      j++;
    }
    std::cout << std::endl;
  }
  if(intron_calc){
    //std::cout << IntronPos.size();
    for(int i = 0; i < (int)IntronPos.size() + 1; i++){
      sonsIntron.push_back({});
      parentsIntron.push_back({});
    }
    for(int i = 1; i < (int)IntronPos.size() + 1; i++){
      std::list<int> tmpS;
      std::list<int> tmpP;
      for(int o = 1; o < i; o++){
	 //for(int o = 1; o < (int)IntronPos.size(); o++){
	tmpP.push_back(o);
      }
      for(int o = i ; o < (int)IntronPos.size() ; o++){
      //for(int o = 1; o < (int)IntronPos.size(); o++){
	tmpS.push_back(o+1);
      }
      sonsIntron[i] = tmpS;
      parentsIntron[i] = tmpP;
    }
  }
  if(verbose && intron_calc){
    std::cout << "parentsIntrons:\n";
    int count = 0;
    for(auto i: parentsIntron){
      std::cout << count << ": ";
      for(auto o: i){
	std::cout << o << " ";
      }
      count++;
      std::cout << std::endl;
    }
    count = 0;
    std::cout << "sonsIntrons:\n";
    for(auto i: sonsIntron){
      std::cout << count << ": ";
      for(auto o: i){
	std::cout << o << " ";
      }
      count++;
      std::cout << std::endl;
    }
  }
  if(verbose && intron_calc){
    std::cout << "introns map construction, step 0:\n";
    int count_intron = 1;
    for(auto i: IntronPos){
      std::cout << i.first << " " << i.second << ": "<< count_intron << std::endl;
      count_intron++;
    }
    std::cout << "introns map construction, step 1:\n";
    for(auto i: introns){
      std::cout << i.first.first << " " << i.first.second<< ": ";
      for(auto o: i.second){
	std::cout << o.first << " "  << o.second << " | ";
      }
      std::cout << std::endl;
    }
    std::cout << std::endl;
    std::cout << "introns map construction, step 2:\n";
    for(auto i: idIntronsList){
      std::cout << i.first.first << " " << i.first.second << " :";
      for(auto o: i.second){
	std::cout << o.first << ", " << o.second << " | ";
      }
      std::cout<<std::endl;
    }
    std::cout << std::endl;
    std::cout << "final id introns map:\n";
    for(auto i: idIntronsMap){
      std::cout << i.first.first << " " << i.first.second<< ": ";
      for(auto o: i.second){
	std::cout << o << " | ";
      }
      std::cout << std::endl;
    }
    std::cout << std::endl;
    std::cout << "Introns linearization:\n";
    std::cout << TIntrons << "\n";
  }

  gtfFile.close();
  setupBitVector();
  setupBitVectorIntrons();
  save(gtf);
  saveIntron(gtf);
}

// metodo per il bitvector
void SplicingGraph::setupBitVector() {
  // si inizializza a 0 il bitvector lungo come la stringa della linearizzazione
  sdsl::bit_vector BV (T.size(), 0);
  
  unsigned int i = 0;
  while(i < T.size()) {
    // si mette 1 in corrispondenza del carattere speciale 
    if(T[i] == '|') {
      BV[i] = 1;
    }
    i++;
  }
  // rrr_vector è un h_0 compressed vector
  bitVector = sdsl::rrr_vector<>(BV);
  
  // setup di rank e select 
  selectBV = sdsl::rrr_vector<>::select_1_type(&bitVector);
  rankBV = sdsl::rrr_vector<>::rank_1_type(&bitVector);
}

// getter della linearizzazione
std::string SplicingGraph::getText() const {
  return T;
}


// metodo per il bitvector introni
void SplicingGraph::setupBitVectorIntrons() {
  // si inizializza a 0 il bitvector lungo come la stringa della linearizzazione
  sdsl::bit_vector BV (TIntrons.size(), 0);
  
  unsigned int i = 0;
  while(i < TIntrons.size()) {
    // si mette 1 in corrispondenza del carattere speciale 
    if(TIntrons[i] == '|') {
      BV[i] = 1;
    }
    i++;
  }
  // rrr_vector è un h_0 compressed vector
  bitVectorIntrons = sdsl::rrr_vector<>(BV);
  
  // setup di rank e select introni
  selectBVIntrons = sdsl::rrr_vector<>::select_1_type(&bitVectorIntrons);
  rankBVIntrons = sdsl::rrr_vector<>::rank_1_type(&bitVectorIntrons);
}

// getter della linearizzazione introni
std::string SplicingGraph::getTextIntrons() const {
  return TIntrons;
}



// getter della lista parens di un nodo
std::list<int> SplicingGraph::getParents(const int& i) const {
  return parents[i];
}

// getter della lista sons di un nodo
std::list<int> SplicingGraph::getSons(const int& i) const {
  return sons[i];
}

// getter dell'esone sfruttando la select (che ritorna l'indice di inzio
// dell'i-simo 1) per poi ritornare la sottostringa
std::string SplicingGraph::getExon(const int& i) const {
  int s = selectBV(i) + 1;
  int e = selectBV(i + 1);
  std::string exonText = T.substr(s, e - s);
  return exonText;
}

std::string SplicingGraph::getIntron(const int& i) const {
  int s = selectBVIntrons(i) + 1;
  int e = selectBVIntrons(i + 1);
  std::string intronText = TIntrons.substr(s, e - s);
  return intronText;
}
// metodo per la rank
int SplicingGraph::rank(const int& i) const {
  return rankBV(i);
}

// metodo per la select
int SplicingGraph::select(const int& i) const {
  return selectBV(i);
}

// rank introns
int SplicingGraph::rankIntrons(const int& i) const {
  return rankBVIntrons(i);
}

// select introns
int SplicingGraph::selectIntrons(const int& i) const {
  return selectBVIntrons(i);
}

// metodo che ritorna true se ho un arco tra due nodi
bool SplicingGraph::contain(const int& x, const int& y) const {
  if(edges[x][y] >= 1) {
    return true;
  } else {
    return false;
  }
}

// metodo che ritorna true se ho un arco novel tra due nodi 
bool SplicingGraph::isNew(const int& x, const int& y) const {
  if(edges[x][y] > 1) {
    return true;
  } else {
    return false;
  }
}

// metodo per la stampa
void SplicingGraph::print() const {
  std::cout << T << std::endl;
  std::cout << std::endl;
  unsigned int i = 0;
  while(i < bitVector.size()) {
    std::cout << bitVector[i];
    i++;
  }
  std::cout << std::endl;
  std::cout << std::endl;
  for(std::vector<int> v : edges) {
    for(const int& e : v) {
      std::cout << e << " ";
    }
    std::cout << std::endl;
  }
  std::cout << std::endl;
  for(auto e : parents) {
    for(auto p : e) {
      std::cout << p << " ";
    }
    std::cout << std::endl;
  }
  std::cout << std::endl;
  for(auto e : sons) {
    for(auto s : e) {
      std::cout << s << " ";
    }
    std::cout << std::endl;
  }
}

// metodo per il slvataggio su file
void SplicingGraph::save(const std::string path) {
  std::ofstream ofile;
  ofile.open(path + ".sg");
  ofile << reference << " " << refLen << " ";
  if(strand)
    ofile << "+";
  else
    ofile << "-";
  ofile << "\n";
  ofile << T << "\n";
  ofile << exsN << "\n";
  for(const std::vector<int>& v : edges) {
    for(const int& e : v) {
      ofile << e << " ";
    }
    ofile << "; ";
  }
  ofile << "\n";
  for(const std::pair<int,int>& p : ExonsPos) {
    ofile << p.first << "," << p.second << " ";
  }
  ofile << "\n";
  // for(const std::string& name : ExonsName) {
  //     ofile << name << " ";
  // }
  ofile << "\n";
  ofile.close();
}

int SplicingGraph::getExonsNumber() const {
  return exsN;
}

bool sortbysec(const std::pair<int,int> &a, 
	       const std::pair<int,int> &b) { 
  return (a.second < b.second); 
}

// std::map<std::pair<int, int>, int > SplicingGraph::getIntronsId() const{
//   return intronsId;
// }

std::map<std::pair<unsigned long long int, unsigned long long int>,
	 std::vector<unsigned long long int> > SplicingGraph::getIntronMap()
  const{
  return idIntronsMap;
}

  
// metodo per il slvataggio su file
void SplicingGraph::saveIntron(const std::string path) {
  std::ofstream ofile;
  ofile.open(path + ".intron.sg");
  ofile << reference << " " << refLen << " ";
  if(strand)
    ofile << "+";
  else
    ofile << "-";
  ofile << "\n";
  ofile << TIntrons << "\n";
  for(const std::pair<unsigned long long int, unsigned long long int>& p : IntronPos) {
    ofile << p.first << "," << p.second << " ";
  }
  ofile << "\n";
  ofile.close();
}

std::list<int> SplicingGraph::getParentsIntron(const int& i) const {
  // std::list<int> res;
  // unsigned int size = parentsIntron.size();
  // if(parentsIntron[i].size() < size/5){
  //   res = parentsIntron[i];
  // }else{
  //   for(unsigned int o = (parentsIntron[i].size() - size/5) + 1; o <= parentsIntron[i].size(); o++){
  //     res.push_back(o);
  //   }
  // }
  // return res;
  return parentsIntron[i];
}

std::list<int> SplicingGraph::getSonsIntron(const int& i) const {
  // std::list<int> res;
  // unsigned int size = sonsIntron.size();
  // unsigned int count = 0;
  // for(auto o: sonsIntron[i]){
  //   if(count == size/5)
  //     break;
  //   res.push_back(o);
  //   count++;
  // }
  // return res;
  return sonsIntron[i];
}
