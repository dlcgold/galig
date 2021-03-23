#include "MEMsGraph.hpp"


bool equalMems(const Mem& a, const Mem& b){
  return ((a.t == b.t )&& (a.p == b.p) && (a.l == b.l));
}


bool operator<(const Mem &a, const Mem &b){
  return a.t + a.l < b.t + b.l;
}



MemsGraph::MemsGraph(const std::string& read_,
                     const int& L_,
                     const int& eps_,
                     const int& exsN_,
                     const bool& verbose_,
		     const bool& intron_calc_,
		     const bool& longIntron) : AnnNodesMap(AnnGraph),
						 AnnEdgesMap(AnnGraph),
						 NovNodesMap(NovGraph),
						 NovEdgesMap(NovGraph){
  // read e lunghezza della stessa
  read = read_;
  m = read.size();

  // parametri
  L = L_;
  LI = (longIntron) ? m / 2 : L;
  eps = eps_;
  float t1 = eps * m * L;
  float t2 = eps * m;

  // ?
  K0 = ceil(t1 / 100);

  // lunghezza gap
  K1 = ceil(t1 / 100) - L + 1;

  // errore gap
  K2 = ceil(t2 / 100);
  //std::cout << m << " " << t2 << " " << K2 << std::endl;
  if(intron_calc)
    K2 *= 4;
  // numero esoni
  exsN = exsN_;

  // opzione verbose
  verbose = verbose_;

  // opzione per introni
  intron_calc = intron_calc_;
  
  // nodi iniziali e finali del grafo degli annotati e dei novel
  AnnStart = AnnGraph.addNode();
  AnnEnd = AnnGraph.addNode();
  NovStart = NovGraph.addNode();
  NovEnd = NovGraph.addNode();
  AnnNodesMap[AnnStart] = Mem(0, 0, 0, false);
  AnnNodesMap[AnnEnd] = Mem(-1, -1, -1, false);
  NovNodesMap[NovStart] = Mem(0, 0, 0, false);
  NovNodesMap[NovEnd] = Mem(-1, -1, -1, false);
}

std::pair<bool, int> MemsGraph::checkMEMs(const SplicingGraph& sg,
                                          const Mem& m1,
                                          const Mem& m2) {
  
  // si usa la rank per capire da dove arrivano le mem, di cui uso l'indice di
  // partenza sulla linearizzazione
  int id1 = sg.rank(m1.t - 1);
  int id2 = sg.rank(m2.t - 1);

  // si estrapola le seq dei due esoni (posto che siano due)
  std::string exon1_text = sg.getExon(id1);
  std::string exon2_text = sg.getExon(id2);
  
  int err = -1;
  bool type = true;
  
  
  if(verbose) {
    std::cout << "Extending " << m1.toStrCat() << " with "
	      << m2.toStrCat() << std::endl;
  }
  
  // primo caso: si è sullo stesso esone (casi 1, 2, 3)
  if(id1 == id2) { //m1 and m2 in the same exon
    //if(m2.p+m2.l>m1.p+m1.l && m1.t<m2.t && m2.t<m1.t+m1.l+K1
    // && m1.t+m1.l<m2.t+m2.l) {
    // seconda mem completamente a destra
    if(verbose)
      std::cout << "check mem single\n";
    if(m2.p + m2.l > m1.p +m1.l &&
       m1.t < m2.t && m1.t + m1.l < m2.t + m2.l) {
      if(verbose) {
	std::cout << "same exon" << std::endl;
      }
      // if(validStart(sg, m2).first){
      // 	startIntrons(m2, sg);
      // }else if (validEnd(sg, m1).first){
      // 	endIntrons(m1, sg);
      // }
      // gap sulla read
      
      int gapP = m2.p - m1.p - m1.l;

      // gap sull'esone 
      int gap_E = m2.t - m1.t - m1.l;

      // casi 1  
      if(gapP >= 0 && gap_E >= 0) {
	
	// prima parte caso 1
	if(gapP == 0) {
	  // se il gap sull'esone non è trascurabile probabilmente si ha intorne
	  if(gap_E > K2) {
	    //Possible intron
	    if(verbose) {
	      std::cout << "Possible intron without overlap" << std::endl;
	    }
	    // peso a 0
	    err = 0;
	    type = false;
	  } else {
	    //Errors
	    if(verbose) {
	      std::cout << "Nothing" << std::endl;
	    }
	    // errore pari al gap (non avendo gap sulla read la distanza di
	    // edit è pari alla lunghezza del gap sull'esone)
	    err = gap_E;
	    type = true;
	  }
	} else if(abs(gapP - gap_E) <= K2) {
	  // seconda parte caso 1 
	  // Possible SNV (single nucleotide variant)
	  if(verbose) {
	    std::cout << "Nothing" << std::endl;
	  }

	  // si estrapola sequenza del gap su read ed esone
	  // std::string::substr(inizio, lunghezza)
	  std::string sub_P = read.substr(m1.p + m1.l - 1,
					  m2.p - m1.p - m1.l);
	  std::string sub_E = exon1_text.substr(m1.t + m1.l -
						sg.select(id1) - 1 - 1,
						m2.t - m1.t - m1.l);

	  // errore (peso arco) come la distanza di edit tra i due gap
	  err = editDistance(sub_P, sub_E);
	  type = true;
	}
      } else if(gapP <= 0 && gap_E <= 0) {
	// caso 2
	if(verbose) {
	  std::cout << "Nothing" << std::endl;
	}
	// errore (peso arco) come modulo tra le lunghezze dei due gap
	err = abs(gapP - gap_E);
	type = true;
	
      } else if(gapP <= 0 && gap_E > K2) {
	// caso 3
	// Possible intron
	if(verbose) {
	  std::cout << "Possible intron with overlap" << std::endl;
	}
	// peso a 0
	err = 0;
	type = false;
      } else {
	if(verbose) {
	  std::cout << "Nothing" << std::endl;
	}

	// altrimenti peso come somma dei gap
	err = abs(gapP) + abs(gap_E);
	type = true;
      }
    }
    // caso in cui i due mems sono da esoni diversi
    // casi 4, 5, 6
  } else { //m1 and m2 in different exons
    // check che gli esoni vengano da due esoni collegati nello splicing graph
    if(verbose)
      std::cout << "check different\n";
    if(sg.contain(id1, id2)) {
      if(verbose) {
	std::cout << "different exons" << std::endl;
      }
      // mem2 deve essere completamente a destra del mem1 sulla read
      if(m2.p + m2.l > m1.p + m1.l) {

	// if(validStart(sg, m2).first){
	//   startIntrons(m2, sg);
	// }else if (validEnd(sg, m1).first){
	//   endIntrons(m1, sg);
	// }
	
	// lunghezza gap tra i mem sulla read
	int gapP = m2.p - m1.p - m1.l;
	
	// gapE1: lunghezza suff_Z(m1)
	int gapE1 = sg.select(id1 + 1) + 1 - m1.t - m1.l;

	// gapE2: lunghezza pref_Z(m2)
	int gapE2 = m2.t - sg.select(id2) - 1 - 1;
	// caso 4
	if(gapP <= 0) {
	  err = 0; //abs(gapP);
	  if(verbose) {
	    std::cout <<"gap negativo: "<<id1 << " " << id2 << " " << sg.isNew(id1, id2)
		      << " " << gapE1 << " " << gapE2 << std::endl;
	  }
	  
	  // se inoltre non ho un arco novel ed entrambi i gap esonici nulli
	  if(!sg.isNew(id1, id2) && gapE1 == 0 && gapE2 == 0) {
	    type = true;
	  }
	  // altrimenti guardo la qualità dell'errore
	  else if(err <= K2) {
	    //Possible Competing
	    type = false;
	  } else
	    err = -1;
	  
	} else {
	  // caso 5 
	  if(gapE1 == 0 && gapE2 == 0) {
	    //Possible insertion (only if annotated edge)
	    // se ho un arco non novel potrei avere un introne
	    if(!sg.isNew(id1, id2)) {
	      err = 0;
	      type = false;
	    }
	  } else {
	    // caso 6
	    if(abs(gapP - (gapE1 + gapE2)) <= K2) {
	      //Possible SNV (single nucleotide variant)
	      if(verbose) {
		std::cout << "SNV" << std::endl;
	      }
	      
	      // si estrapolano le sequenze utili
	      // sequenza relativa la gap sulla read
	      std::string subP = read.substr(m1.p + m1.l - 1, gapP);

	      // sequenza relativa a suff_Z(m1)
	      std::string subE1 = exon1_text.substr(m1.t + m1.l -
						    sg.select(id1) - 1 - 1,
						    gapE1);

	      // sequenza relativa a pref_Z(m2)
	      std::string subE2 = exon2_text.substr(0, gapE2);

	      // concatenazione di suff_Z(m1) e pref_Z(m2) 
	      std::string subE = subE1 + subE2;

	      // calcolo distanza di edit 
	      err = editDistance(subP, subE);
	      
	      if(!sg.isNew(id1, id2)) {
		type = true;
	      } else {
		type = false;
	      }
	    }
	  }
	}
      }
    } else {
      if(verbose) {
	std::cout << "no edge" << std::endl;
      }
    }
  }
  if(err > K2) {
    err = -1;
  }
  if(verbose) {
    std::cout << type << " " << err << std::endl;
  }

  // si restituisce la coppia del boolenao e dell'errore
  return std::make_pair(type, err);
}

// metodo per stabilire se un mem è un mem di partenza
std::pair<bool, int> MemsGraph::validStart(const SplicingGraph& sg,
					   const Mem& Mem) {
  
  // si controlla di essere nel range dello start proximity (K0)
  if(Mem.p <= K0) {
    if(verbose)
      std::cout << "minore di " << K0 << "\n";
    // so inizializza l'errore a partire da K2 (errore sul gap permesso)
    int err = K2 + 1;

    // se si è all'inizio della read non si ha però errore
    if(Mem.p == 1) {
      err = 0;
    } else {
      // si usa la rank per identificare l'esone
      int id = sg.rank(Mem.t - 1);

      // si recupera la sequenza dell'esone 
      std::string exon_text = sg.getExon(id);

      // si identifica la lunghezza del prefisso sulla read e la sequenza
      int l = Mem.p - 1;
      std::string sub_P = read.substr(0, l);

      // si inizializza per sottostringa esonica
      std::string sub_E;

      // calcolo della lunghezza del prefisso dell'esone
      int exon_pref_len = Mem.t - sg.select(id) - 1 - 1;

      // caso in cui il prefisso sull'esone è più corto di quello della read
      if(exon_pref_len < l) {
	// calcolo della lunghezza del prefisso comune a esone e read
	int shared_pref_len = l - exon_pref_len;

	// si estrae la sottosequenza del prefisso dell'esone
	std::string exon_pref = exon_text.substr(0, exon_pref_len);

	// si setta l'errore alla lunghezza del prefisso sulla read
	err = l;

	// si prende la lista coi parents dell'esone dallo splicing graph
	std::list<int> parents = sg.getParents(id);
	if(verbose){
	  for(auto i: parents){
	    std::cout << i << " ";
	  }
	}
	// si itera su ogni nodo nella lista parents per stabilire se ho un mem
	// iniziale
	for(std::list<int>::iterator it=parents.begin();
	    it!=parents.end(); ++it) {
	  /**
	   * We look only at the father of the node,
	   * IF he is long enough, we get its suffix;
	   * ELSE we get all its label (without going further checking all 
	   * its parents)
	   **/
	  // si estrae la sequenza del nodo preso in considerazione
	  int par = *it;
	  std::string par_text = sg.getExon(par);
	  // si guarda se il nodo padre copre la differenza tra il prefisso
	  // di read ed esone
	  if(sg.select(par + 1) - shared_pref_len - sg.select(par) - 1 >= 0) {

	    // se la copre prendiamo il suffisso del padre che permette di
	    // coprire il gap e lo sommiamo al prefisso dll'esone 
	    sub_E = par_text.substr(sg.select(par+1) - shared_pref_len -
				    sg.select(par) - 1, shared_pref_len)
	      + exon_pref;
	  } else {
	    // altrimenti si prende tutta la sequenza del nodo parent in analisi
	    // e la si somma al prefisso dell'esone
	    sub_E = par_text + exon_pref;
	  }
	  // if(verbose){
	  //   std::cout << "subP: " << sub_P << "\n";
	  //   std::cout << "subE: " << sub_E << "\n";
	  // }
	  // si aggiorna l'errore tenendo il minimo
	  int curr_err = editDistance(sub_P, sub_E);
	  if(curr_err < err) {
	    // if(verbose){
	    //   std::cout << "min\n";
	    //   std::cout << "subP: " << sub_P << "\n";
	    //   std::cout << "subE: " << sub_E << "\n";
	    // }
	    err = curr_err;
	  }
	}
      } else {
	// se invece il prefisso dell'esone è più lungo di quello della read
	// sub_E è diventa il prefisso dell'esone che eccede oltre quello della
	// read
	sub_E = exon_text.substr(Mem.t - l - sg.select(id) - 1 - 1, l);
	err = editDistance(sub_P, sub_E);
      }
    }
    if(verbose){
      std::cout << "err: " << err << "\n";
      std::cout << "k2: " << K2 << "\n";
    }
    // se alla fine l'errore è accettabile ritorno che è un valid start
    if(err <= K2) {
      return std::make_pair(true, err);
    }
  }
  // si ritorna la coppia indicante che non è un valid start (con un errore
  // sentinella)
  return std::make_pair(false, K2 + 1);
}

// metodo per verificare se un mem è un valido mem finale
std::pair<bool, int> MemsGraph::validEnd(const SplicingGraph& sg,
					 const Mem& Mem) {

  // se la mem sulla read termina in prossimità (prossimità K0) della fine della
  // read 
  if(Mem.p + Mem.l >= m - K0) {
    // errore impostato con K2
    int err = K2 + 1;

    // se la mem sulla read supera la fine errore a 0
    if(Mem.p + Mem.l == m + 1) {
      err = 0;
    } else {
      // si recupera l'id dell'esone
      int id = sg.rank(Mem.t - 1);

      // si recupera la sequenza dell'esone
      std::string exon_text = sg.getExon(id);

      // lunghezza del suffisso della read oltre l'ultimo mem
      int l = m - (Mem.p + Mem.l ) + 1;

      // estrazione della sequenza di tale suffisso
      std::string sub_P = read.substr(Mem.p+Mem.l-1, l);

      // inizializzazione per la sottostringa sull'esone
      std::string sub_E;

      // calcolo della lunghezza del suffisso sull'esone a partire dalla fine
      // del mem (sfrutto l'inizio dell'esone successivo)
      int exon_suff_len = sg.select(id + 1) - (Mem.t + Mem.l) + 1;

      // se il suffisso sull'esone è minore di quello sulla read
      if(exon_suff_len < l) {

	// si estrapola la lista dei sons
	std::list<int> sons = sg.getSons(id);

	// si calcola la lunghezza del suffisso condiviso tra read esone
	int shared_suff_len = l - exon_suff_len;

	// si inizalizza la stringa per il suffisso dell'esone
	std::string exon_suff;
	if(exon_suff_len == 0) {
	  // ovviamente se la lunghezza di tale suffisso è nulla si ha la
	  // stringa vuota
	  exon_suff = "";
	} else {
	  // altrimenti prendo la sequenza di tale suffisso 
	  exon_suff = exon_text.substr(Mem.t + Mem.l -sg.select(id) - 1 - 1,
				       exon_suff_len);
	}
	// si imposta l'errore alla lunghezza del suffisso sulla read
	err = l;

	// si itera sui sons per vedere se posso stabilire che ho un mem finale
	for(std::list<int>::iterator it=sons.begin(); it!=sons.end(); ++it) {

	  // si estrapola la sequenza del son in considerazione
	  int son = *it;
	  std::string son_text = sg.getExon(son);
	  // si calcola il nuovo suffisso partendo da quello dell'esone
	  sub_E = exon_suff + son_text.substr(0, shared_suff_len);
	  int curr_err = editDistance(sub_P, sub_E);

	  // si calcola l'errore con la distanza di edit
	  if(curr_err < err) {
	    err = curr_err;
	  }
	}
      } else {
	// se la lunghezza del suffisso sulla read è minore
	// sub_E diventa il suffisso dell'esone a partire dall'ultima mem di
	// lunghezza pari al suffisso della read
	sub_E = exon_text.substr(Mem.t + Mem.l - sg.select(id) - 1 - 1, l);

	// si usa la distanza di edit per calclare l'errore
	err = editDistance(sub_P, sub_E);
      }
    }
    //std::cout << "K2: " << K2 << "\n";

    if(err <= K2) {
      //std::cout << "K2: " << K2 << "\n";
      // se l'errore alla fine è accettabile so che è un mem finale con
      // quell'errore
      return std::make_pair(true, err);
    }
  }
  // se non è un mem finale lo segnalo
  return std::make_pair(false,K2+1);
}


// metodo per la costruzione del memGraph 
void MemsGraph::build(const SplicingGraph& sg, const std::list<Mem>& MEMs_,
		      BackwardMEM& bmi) {
  

  if(!MEMs_.empty()){
    std::vector<std::forward_list<Mem> > MEMs (m+1);
    
    // ciclo che itera sui mems e carica nella posizone p del vector di lists
    // mems il mem che inizia all'indice p sulla read
    // nella lista ho i mem di varia lunghezza inizianti in quell'indice
    for(const Mem& m : MEMs_) {
      int p = m.p;
      //std::cout << m.t<< " "<<p << "\n";
      MEMs[p].push_front(m);
    }
    //std::cout << "test\n";
    // std::cout << "size introns: " << MEMsI_.size() << "\n";
    // si itera sul vector di lists
    bool addStart = false;
    for(std::vector<std::forward_list<Mem> >::iterator it=MEMs.begin();
	it!=MEMs.end(); ++it) {
      
      // itero sulle mems nella list
      for(Mem& m1 : *it) {
	// estrapolazione della posizione di inizio della mem sulla read
	int p1 = m1.p;
	// inizializzazione dei due nodi, uno per il grafo annotato e uno per
	// quello novel
	Node AnnNode1;
	Node NovNode1;
	// estrapolazione dell'id del primo esone e della sua sequenza
	int id1 = sg.rank(m1.t - 1);
	std::string exon1_text = sg.getExon(id1);
	//std::cout << "id1 esone: " << id1 << "\n";
	/*********
	 * Start *
	 *********/
	// valuto se il mem in questione è uno valido per essere definito inizale
	std::pair<bool, int> startInfo = validStart(sg, m1);
	bool startFlag = startInfo.first;
	int err = startInfo.second;
	if(verbose)
	  std::cout << "start for " << m1.toStrCat() << " "
		    << startFlag << " " << err << "\n";
	
	if(intron_calc && ((!startFlag && !addStart && m1.p > L) ||
			   (!addStart && startFlag && m1.p > L))){
	  auto readI = read.substr(0, m1.p);
	  // std::cout << "inizio calcolo mem iniziali\n";
	  auto MEMsI_ = bmi.getMEMs(readI, L);
	  // std::cout << "fine calcolo mem iniziali\n";

	  if(verbose){
	    std::cout << "start intron size: " << MEMsI_.size() << "\n";
	    for(auto i: MEMsI_){
	      std::cout << i.toStrCat() << "\n";
	    }
	    std::cout << "\n";
	  }
	  if(!MEMsI_.empty() && !MEMsI_.empty()){
	    // std::cout << "inizio startIntrons\n";
	    auto starti = startIntrons(m1, sg, MEMsI_, readI.size());
	    // std::cout << "fine startIntrons\n";
	    // std::cout << "inizio grafo start\n";
	    if(!starti.empty()){
	      for(auto vec: starti){
		if(!vec.empty()){
		  if(m1.isNew){
		    AnnNode1 = AnnGraph.addNode();
		    m1.setAnnNode(AnnNode1);
		    AnnNodesMap[AnnNode1] = m1;
		  }else{
		    AnnNode1 = m1.AnnNode;
		  }
		  std::vector<Node> startAnnIntronNodes;
		  for(auto e: vec){
		    Node tempAnn = AnnGraph.addNode();
		    e.setAnnNode(tempAnn);
		    AnnNodesMap[tempAnn] = e;
		    startAnnIntronNodes.push_back(tempAnn);
		  }
	      
		  for(unsigned int i = 0; i < startAnnIntronNodes.size() - 1; i++){
		    Arc arc = AnnGraph.addArc(startAnnIntronNodes[i],
					      startAnnIntronNodes[i + 1]);
		    // AnnEdgesMap[arc] = std::abs(AnnNodesMap[startAnnIntronNodes[i + 1]].p -
		    // 			    (AnnNodesMap[startAnnIntronNodes[i]].p +
		    // 			     AnnNodesMap[startAnnIntronNodes[i]].l));
		    AnnEdgesMap[arc] = 0;
		  }
		  addStart = true;
		  Arc arcIF = AnnGraph.addArc(AnnStart,
					      startAnnIntronNodes.front());
		  AnnEdgesMap[arcIF] = (AnnNodesMap[startAnnIntronNodes.front()].p - 1);
		  //AnnEdgesMap[arcIF] = 1;
		  Arc arcEI = AnnGraph.addArc(startAnnIntronNodes.back(),
					      AnnNode1);
		  // AnnEdgesMap[arcEI] = std::abs(AnnNodesMap[AnnNode1].p -
		  // 				    (AnnNodesMap[startAnnIntronNodes.back()].p +
		  // 				     AnnNodesMap[startAnnIntronNodes.back()].l));
		  AnnEdgesMap[arcEI] = 1;
		  // Arc arcEE = AnnGraph.addArc(AnnNode1, AnnEnd);
		  // // AnnEdgesMap[arcEE] = std::abs((int)read.size() -
		  // // 				    (AnnNodesMap[AnnNode1].p +
		  // // 				     AnnNodesMap[AnnNode1].l - 1));
		  // AnnEdgesMap[arcEE] = 1;
		  startAnnIntronNodes.clear();
		}
	      }
	    }
	    // std::cout << "fine grafo start\n";
	  }
	}
	//std::cout <<"test start: " << startFlag << " " << err <<"\n";
	// se effettivamente il mem è valido come iniziale
	if(startFlag) {
	  //std::cout << "valid start\n";
	  if(verbose)
	    std::cout << "start "<<err<<std::endl;
	  // studio prefisso read
	  //startIntrons(m1, sg);
	  if(m1.isNew) {
	    // si creano e si aggiungono i due nodi ai due grafi
	    AnnNode1 = AnnGraph.addNode();
	    NovNode1 = NovGraph.addNode();

	    // setter per memorizzare i nodi nella struttura stessa del mem 
	    m1.setAnnNode(AnnNode1);
	    m1.setNovNode(NovNode1);

	    // si settano anche le due mappe
	    AnnNodesMap[AnnNode1] = m1;
	    NovNodesMap[NovNode1] = m1;
	  
	    // if the MEM has a father, we don't link it to the START
	    // (even if it is a valid candidate as starting mem)
	    // si aggiunge un arco tra il nodo start e il nodo in considerazione
	    // in entrambi i grafi con peso pari a err
	    addStart = true;
	    Arc arc = AnnGraph.addArc(AnnStart, AnnNode1);
	    AnnEdgesMap[arc] = err;
	    arc = NovGraph.addArc(NovStart, NovNode1);
	    NovEdgesMap[arc] = err;
	  
	  } else {
	    // altrimenti recupero i nodi
	    AnnNode1 = m1.AnnNode;
	    NovNode1 = m1.NovNode;
	  }
	  
	} else {
	  
	  //std::cout << "start new " << m1.isNew << "\n";
	  // se non è iniziale e non è nuovo recupero i nodi
	  // if(intron_calc && m1.p > L){
	  //   auto readI = read.substr(0, m1.p);
	  //   auto MEMsI_ = bmi.getMEMs(readI, L);
	  //   if(verbose){
	  //     std::cout << "start intron size: " << MEMsI_.size() << "\n";
	  //     for(auto i: MEMsI_){
	  // 	std::cout << i.toStrCat() << "\n";
	  //     }
	  //     std::cout << "\n";
	  //   }
	  //   if(!MEMsI_.empty() && !MEMsI_.empty()){
	  //     auto starti = startIntrons(m1, sg, MEMsI_, readI.size());
	  //     if(!starti.empty()){
	  // 	for(auto vec: starti){
	  // 	  if(!vec.empty()){
	  // 	    if(m1.isNew){
	  // 	      AnnNode1 = AnnGraph.addNode();
	  // 	      m1.setAnnNode(AnnNode1);
	  // 	      AnnNodesMap[AnnNode1] = m1;
	  // 	    }else{
	  // 	      AnnNode1 = m1.AnnNode;
	  // 	    }
	  // 	    std::vector<Node> startAnnIntronNodes;
	  // 	    for(auto e: vec){
	  // 	      Node tempAnn = AnnGraph.addNode();
	  // 	      e.setAnnNode(tempAnn);
	  // 	      AnnNodesMap[tempAnn] = e;
	  // 	      startAnnIntronNodes.push_back(tempAnn);
	  // 	    }
	      
	  // 	    for(unsigned int i = 0; i < startAnnIntronNodes.size() - 1; i++){
	  // 	      Arc arc = AnnGraph.addArc(startAnnIntronNodes[i],
	  // 					startAnnIntronNodes[i + 1]);
	  // 	      // AnnEdgesMap[arc] = std::abs(AnnNodesMap[startAnnIntronNodes[i + 1]].p -
	  // 	      // 			    (AnnNodesMap[startAnnIntronNodes[i]].p +
	  // 	      // 			     AnnNodesMap[startAnnIntronNodes[i]].l));
	  // 	      AnnEdgesMap[arc] = 0;
	  // 	    }
	  // 	    Arc arcIF = AnnGraph.addArc(AnnStart,
	  // 					startAnnIntronNodes.front());
	  // 	    // AnnEdgesMap[arcIF] = (AnnNodesMap[startAnnIntronNodes.front()].p - 1) + 1;
	  // 	    AnnEdgesMap[arcIF] = 1;
	  // 	    Arc arcEI = AnnGraph.addArc(startAnnIntronNodes.back(),
	  // 					AnnNode1);
	  // 	    // AnnEdgesMap[arcEI] = std::abs(AnnNodesMap[AnnNode1].p -
	  // 	    // 				    (AnnNodesMap[startAnnIntronNodes.back()].p +
	  // 	    // 				     AnnNodesMap[startAnnIntronNodes.back()].l));
	  // 	    AnnEdgesMap[arcEI] = 1;
	  // 	    Arc arcEE = AnnGraph.addArc(AnnNode1, AnnEnd);
	  // 	    // AnnEdgesMap[arcEE] = std::abs((int)read.size() -
	  // 	    // 				    (AnnNodesMap[AnnNode1].p +
	  // 	    // 				     AnnNodesMap[AnnNode1].l - 1));
	  // 	    AnnEdgesMap[arcEE] = 0;
	  // 	    startAnnIntronNodes.clear();
	  // 	  }
	  // 	}
	  //     }
	  //   }
	  // }
	  
	  if(m1.isNew){
	    continue;
	  
	  } else {
	    AnnNode1 = m1.AnnNode;
	    NovNode1 = m1.NovNode;
	  }
	}    
	
	
	/*************
	 * Extending *
	 *************/
	bool AnnExt = false;
	bool NovExt = false;
	int p2 = p1+1;
      
	//std::cout <<"extend p2: "<< p2 << "\n";
	// si itera fino ad uscire dalla lunghezza della read (al più di L
	// lunghezza minima della mem)
	while(p2 <= m - L + 1) {
	  //std::cout << "extend" << p2 <<"\n";
	  // si itera su tutti i mems che iniziano dopo m1 (quelli almeno
	  // parzialmente a destra) in posizione p2 secondo le varie lunghezze
	  for(Mem& m2 : MEMs[p2]) {
	    // inizializzo nodi e faccio il check mem valutando il pair di ritorno
	    Node AnnNode2;
	    Node NovNode2;
	    std::pair<bool, int> linkageInfo = checkMEMs(sg, m1, m2);
	    bool flag = linkageInfo.first;
	    int err = linkageInfo.second;
	    if(err>=0) {
	    
	      // se m2 è nuovo
	      if(m2.isNew) {
		// si creano i nodi e si aggiungono
		AnnNode2 = AnnGraph.addNode();
		NovNode2 = NovGraph.addNode();
		m2.setAnnNode(AnnNode2);
		m2.setNovNode(NovNode2);
		AnnNodesMap[AnnNode2] = m2;
		NovNodesMap[NovNode2] = m2;
	      } else {
		// altrimenti si recuperano i nodi
		AnnNode2 = m2.AnnNode;
		NovNode2 = m2.NovNode;
	      }
	      // inoltre se il flag è a true si aggiunge anche l'arco tra i due
	      // mem in entrambi i grafi con err come peso e si segnala che
	      // non sono nodi iniziali ma fanno parte di quelli che estendono (?)
	      if(flag) {
	        
		Arc arc = AnnGraph.addArc(AnnNode1, AnnNode2);
		AnnEdgesMap[arc] = err;
		arc = NovGraph.addArc(NovNode1, NovNode2);
		NovEdgesMap[arc] = err;
		AnnExt = true;
		NovExt = true;
	      } else {
		// se il flag è false si aggiunge l'arco solo nel grafo novel
		Arc arc = NovGraph.addArc(NovNode1, NovNode2);
		NovEdgesMap[arc] = err;
		NovExt = true;
	
	      }
	    }
	    if(intron_calc && (!flag || err <0)){
	      if(m2.isNew) {
		// si creano i nodi e si aggiungono
		AnnNode2 = AnnGraph.addNode();
		NovNode2 = NovGraph.addNode();
		m2.setAnnNode(AnnNode2);
		m2.setNovNode(NovNode2);
		AnnNodesMap[AnnNode2] = m2;
		NovNodesMap[NovNode2] = m2;
	      } else {
		// altrimenti si recuperano i nodi
		AnnNode2 = m2.AnnNode;
		NovNode2 = m2.NovNode;
	      }
	      int gapM = m2.p - (m1.p + m1.l);
	     
	      if(gapM > L){
		auto readM = read.substr(m1.p + m1.l, gapM);
		// std::cout << "inizio calcolo mem middle\n";

		auto MEMsI_ = bmi.getMEMs(read, L);
		// std::cout << "fine calcolo mem middle\n";
			  
		if(verbose && !MEMsI_.empty()){
		  std::cout << "middle intron size: " << MEMsI_.size() << "\n";
		}
		// std::cout << "inizio grafo middle\n";
		if(!MEMsI_.empty()){
		  // std::cout << "inizio middleIntrons\n";
		  auto middlei = middleIntrons(m1, m2, sg, MEMsI_, read.size());
		  // std::cout << "fine middleIntrons\n";
		  if(!middlei.empty()){
		    for(auto vec: middlei){
		      if(!vec.empty()){
			std::vector<Node> middleAnnIntronNodes;
			for(auto e: vec){
			  Node tempAnn = AnnGraph.addNode();
			  e.setAnnNode(tempAnn);
			  AnnNodesMap[tempAnn] = e;
			  middleAnnIntronNodes.push_back(tempAnn);
			}
			for(unsigned int i = 0; i < middleAnnIntronNodes.size() - 1;
			    i++){
			  Arc arc = AnnGraph.addArc(middleAnnIntronNodes[i],
						    middleAnnIntronNodes[i + 1]);
			  // AnnEdgesMap[arc] = std::abs(AnnNodesMap[middleAnnIntronNodes[i + 1]].p -
			  // 				  (AnnNodesMap[middleAnnIntronNodes[i]].p +
			  // 				   AnnNodesMap[middleAnnIntronNodes[i]].l));
			  AnnEdgesMap[arc] = 0;
			}
			Arc arcEI = AnnGraph.addArc(AnnNode1,
						    middleAnnIntronNodes.front());
			// AnnEdgesMap[arcEI] = std::abs(AnnNodesMap[middleAnnIntronNodes.front()].p -
			// 				  (AnnNodesMap[AnnNode1].p +
			// 				   AnnNodesMap[AnnNode1].l));
			AnnEdgesMap[arcEI] = 1;
			Arc arcIB = AnnGraph.addArc(middleAnnIntronNodes.back(),
						    AnnNode2);
			// AnnEdgesMap[arcIB] = std::abs(AnnNodesMap[AnnNode2].p -
			// 				  (AnnNodesMap[middleAnnIntronNodes.back()].p +
			// 				   AnnNodesMap[middleAnnIntronNodes.back()].l)) + 1;
			AnnEdgesMap[arcIB] = 1;
			middleAnnIntronNodes.clear();
		      }
		      AnnExt = true;
		    }
		  }
		}
		// std::cout << "fine grafo start\n";
	      }
	    }
	  }
	  // vado nella lista di mem successiva
	  ++p2;
	}
      
	/*******
	 * End *
	 *******/
	// if the MEM doesn't have a son, we link it to the END
	// if it is a valid candidate as ending mem
	// si valuta se si sta studiando un mem finale
	// i campi AnnExt e NovExt restano false
	if(!AnnExt && !NovExt) {
	  // si recuperano i dati di validEnd
	
	  std::pair<bool, int> endInfo = validEnd(sg, m1);
	  bool endFlag = endInfo.first;
	  err = endInfo.second;
	
	  // se si ha a che fare con un end valido
	  if(endFlag) {
	    //std::cout << " end \n";
	    //std::cout << "valid end\n";
	    if(verbose)
	      std::cout << "end " << err << std::endl;
	  
	    // si crea un arco al nodo end di peso err in entrambi i grafi
	    Arc arc = AnnGraph.addArc(AnnNode1, AnnEnd);
	    AnnEdgesMap[arc] = err;
	    arc = NovGraph.addArc(NovNode1,NovEnd);
	    NovEdgesMap[arc] = err;
	  }
	  if(intron_calc){
	    if((!endFlag && ((int)read.size() - (m1.p + m1.l -1) >= L))
	       || (endFlag && (m1.p + m1.l -1) <= (int)read.size() &&
		   (int)read.size() - (m1.p + m1.l-1) >= L)){
	      //if((!endFlag && ((int)read.size() - (m1.p + m1.l) > L))){
	      // std::cout << "\nend (" << m1.toStr() << ": ";
	      auto readE = read.substr(m1.p + m1.l);
	      // std::cout << readE.size() << "\n";
	      // std::cout << "inizio calcolo mem end\n";

	      auto MEMsI_ = bmi.getMEMs(read, L);
	      // std::cout << "fine calcolo mem end\n";

	      if(verbose && !MEMsI_.empty()){
		std::cout << "ending intron size: " << MEMsI_.size() << "\n";
	      }
	      
	      if(!MEMsI_.empty()){
		// std::cout << "inizio endIntrons\n";
		auto endi = endIntrons(m1, sg, MEMsI_, read.size());
		// std::cout << "fine endIntrons\n";
		// std::cout << "inizio grafo end\n";
		if(!endi.empty()){
		  for(auto vec: endi){
		    if(!vec.empty()){
		      std::vector<Node> endAnnIntronNodes;
		      for(auto e: vec){
			Node tempAnn = AnnGraph.addNode();
			e.setAnnNode(tempAnn);
			AnnNodesMap[tempAnn] = e;
			endAnnIntronNodes.push_back(tempAnn);
		      }
		      for(unsigned int i = 0; i < endAnnIntronNodes.size() - 1; i++){
			Arc arc = AnnGraph.addArc(endAnnIntronNodes[i],
						  endAnnIntronNodes[i + 1]);
			// AnnEdgesMap[arc] = std::abs(AnnNodesMap[endAnnIntronNodes[i + 1]].p -
			// 			    (AnnNodesMap[endAnnIntronNodes[i]].p +
			// 			     AnnNodesMap[endAnnIntronNodes[i]].l));
			AnnEdgesMap[arc] = 0;
		      }

		      // Arc arcEII = AnnGraph.addArc(AnnStart, AnnNode1);
		      // AnnEdgesMap[arcEII] = AnnNodesMap[AnnNode1].p - 1;
		      Arc arcEI = AnnGraph.addArc(AnnNode1, endAnnIntronNodes.front());
		      // AnnEdgesMap[arcEI] = std::abs(AnnNodesMap[endAnnIntronNodes.front()].p -
		      // 				    (AnnNodesMap[AnnNode1].p +
		      // 				     AnnNodesMap[AnnNode1].l - 1) - 1);
		      AnnEdgesMap[arcEI] = 1;
		      // std::cout << std::abs(AnnNodesMap[endAnnIntronNodes.front()].p -
		      // 				  (AnnNodesMap[AnnNode1].p +
		      // 				   AnnNodesMap[AnnNode1].l)) << std::endl;
		      Arc arcIB = AnnGraph.addArc(endAnnIntronNodes.back(), AnnEnd);
		      AnnEdgesMap[arcIB] = std::abs((int)read.size() -
						    (AnnNodesMap[endAnnIntronNodes.back()].p +
						     AnnNodesMap[endAnnIntronNodes.back()].l - 1)) + 1;
		      // AnnEdgesMap[arcIB] = 1;
		      // std::cout <<std::abs((int)read.size() -
		      // 				  (AnnNodesMap[endAnnIntronNodes.back()].p +
		      // 				   AnnNodesMap[endAnnIntronNodes.back()].l - 1)) << std::endl;

		      endAnnIntronNodes.clear();
		      // std::cout << std::endl;
		    
		    }
		  }
		}
		// std::cout << "fine grafo end\n";
	      }
	    } 
	  }   
	}
      }
    }
  }
  
  
  // if(verbose) {
  //   save("PartGraph.dot");
  // }
  // std::cout << "inizo arcToDEleteAnn\n";
  // Transitive reduction
  // rimuovo gli archi transitivi superflui (pesano più della somma degli altri)
  // 2 archi da grafo annotato
  std::list<InArc> arcsToDeleteAnn;
  // per ogni nodo nel grafo annotato
  for(NodeIt x (AnnGraph); x!=lemon::INVALID; ++x) {
    // per ogni arco uscente da tale nodo
    for(OutArc XY (AnnGraph, x); XY!=lemon::INVALID; ++XY) {
      // y è il nodo verso cui va l'arco
      Node y = AnnGraph.target(XY);
      // per tutti gli archi uscenti da y
      for(OutArc YZ (AnnGraph, y); YZ!=lemon::INVALID; ++YZ) {
	// z è il nodo verso cui punta y
	Node z = AnnGraph.target(YZ);
	// per tutti gli archi entranti in z
	for(InArc XZ (AnnGraph, z); XZ!=lemon::INVALID; ++XZ) {
	  // x_ è il nodo da cui arriva tale nodo
	  Node x_ = AnnGraph.source(XZ);
	  // se x è uguale a x_
	  if(AnnGraph.id(x) == AnnGraph.id(x_)) {
	    // se il peso tra x e z è minore a quello di x e y + y e z
	    // aggiungo l'arco tra x e z a quelli da eliminare
	    if(AnnEdgesMap[XY] + AnnEdgesMap[YZ] <= AnnEdgesMap[XZ]) {
	      arcsToDeleteAnn.push_back(XZ);
	    }
	  }
	}
      }
    }
  }
  // std::cout << "fine arcToDEleteAnn\n";
  // std::cout << "inizo arcToDEleteNov\n";
  // come sopra ma sul grafo novel
  std::list<InArc> arcsToDeleteNov;
  for(NodeIt x (NovGraph); x!=lemon::INVALID; ++x) {
    for(OutArc XY (NovGraph, x); XY!=lemon::INVALID; ++XY) {
      Node y = NovGraph.target(XY);
      for(OutArc YZ (NovGraph, y); YZ!=lemon::INVALID; ++YZ) {
	Node z = NovGraph.target(YZ);
	for(InArc XZ (NovGraph, z); XZ!=lemon::INVALID; ++XZ) {
	  Node x_ = NovGraph.source(XZ);
	  if(NovGraph.id(x) == NovGraph.id(x_)) {
	    if(NovEdgesMap[XY]+NovEdgesMap[YZ]<=NovEdgesMap[XZ]) {
	      arcsToDeleteNov.push_back(XZ);
	    }
	  }
	}
      }
    }
  }
    // std::cout << "fine arcToDEleteNov\n";

  // cancello gli archi delle due liste ToDelete
  for(const InArc& a : arcsToDeleteAnn) {
    if(AnnGraph.valid(a)) {
      AnnGraph.erase(a);
    }
  }
  for(const InArc& a : arcsToDeleteNov) {
    if(NovGraph.valid(a)) {
      NovGraph.erase(a);
    }
  }
  int countAnnNode = 0;
  for(lemon::ListDigraph::NodeIt n(AnnGraph); n != lemon::INVALID; ++n){
    countAnnNode++;
  }
  int countNovNode = 0;
  for(lemon::ListDigraph::NodeIt n(NovGraph); n != lemon::INVALID; ++n){
    countNovNode++;
  }
  int finalNov = lemon::countInArcs(NovGraph, NovEnd);
  int finalAnn = lemon::countInArcs(AnnGraph, AnnEnd);
  if(intron_calc && (MEMs_.empty() || (countNovNode == 2 && countAnnNode == 2)
		     || (finalNov == 0 && finalAnn == 0))){
    auto MEMs_I =  bmi.getMEMs(read, LI);
    if(!MEMs_I.empty()){
      if(verbose){
	std::cout << "intron mems: ";
	for(auto i: MEMs_I){
	  std::cout << i.toStrCat() << " ";
	}
      }
      std::vector<std::forward_list<Mem> > MEMs (m+1);
      for(const Mem& m : MEMs_I) {
	int p = m.p;
	MEMs[p].push_front(m);
      }
      for(std::vector<std::forward_list<Mem> >::iterator it=MEMs.begin();
	  it!=MEMs.end(); ++it) {
	for(Mem& m1 : *it) {
	  if(verbose)
	    std::cout << "studing: " << m1.toStrCat() << "\n";
	  int p1 = m1.p;
	  Node AnnNode1;
	  //Node NovNode1;
	  int id1 = sg.rankIntrons(m1.t-1);
	  std::string intron1_text = sg.getIntron(id1);
	  /*********
	   * Start *
	   *********/
	  // std::cout << "inizo validStartIntrons\n";
	  std::pair<bool, int> startInfo = validStartIntrons(sg, m1);
	  // std::cout << "fine validStartIntrons\n";
	  if(verbose){
	    std::cout << "start " <<  startInfo.first << " "
		      <<  startInfo.second << "\n";
	  }
	  bool startFlag = startInfo.first;
	  int err = startInfo.second;
	  if(startFlag) {
	    if(verbose)
	      std::cout << "start "<<err<<std::endl;
	    if(m1.isNew) {
	      AnnNode1 = AnnGraph.addNode();
	      //NovNode1 = NovGraph.addNode();
	      m1.setAnnNode(AnnNode1);
	      //m1.setNovNode(NovNode1);
	      AnnNodesMap[AnnNode1] = m1;
	      //NovNodesMap[NovNode1] = m1;
	      //if the MEM has a father, we don't link it to the START (even if it is a valid candidate as starting mem)
	      Arc arc = AnnGraph.addArc(AnnStart, AnnNode1);
	      AnnEdgesMap[arc] = err;
	      // arc = NovGraph.addArc(NovStart,NovNode1);
	      // NovEdgesMap[arc] = err;
	      //
	    }else {
	      AnnNode1 = m1.AnnNode;
	      //NovNode1 = m1.NovNode;
	    }
	    // Arc arc = AnnGraph.addArc(AnnStart, AnnNode1);
	    // AnnEdgesMap[arc] = err;
	    // arc = NovGraph.addArc(NovStart,NovNode1);
	    // NovEdgesMap[arc] = err;
	  } else {
	    if(m1.isNew) {
	      continue;
	    } else {
	      AnnNode1 = m1.AnnNode;
	      //NovNode1 = m1.NovNode;
	    }
	  }
	  /*************
	   * Extending *
	   *************/
	  bool AnnExt = false;
	  //bool NovExt = false;
	  int p2 = p1+1;
	  //	std::cout << "inizial p2 " << p2 << " and K2" << K2 <<"\n";
	  while(p2 <= m-L+1) {
	    for(Mem& m2 : MEMs[p2]) {
	      Node AnnNode2;
	      //Node NovNode2;
	      // std::cout << "inizo checkMEMsIntrons\n";
	      std::pair<bool, int> linkageInfo = checkMEMsIntrons(sg, m1, m2);
	      // std::cout << "fine checkMEMsIntrons\n";
	      bool flag = linkageInfo.first;
	      int err = linkageInfo.second;
	      if(verbose){
		std::cout << flag << "|" << err << "\n";
	      }
	      //std::cout << "qui\n";
	      if(err>=0) {
		if(m2.isNew) {
		  AnnNode2 = AnnGraph.addNode();
		  //NovNode2 = NovGraph.addNode();
		  m2.setAnnNode(AnnNode2);
		  //m2.setNovNode(NovNode2);
		  AnnNodesMap[AnnNode2] = m2;
		  //NovNodesMap[NovNode2] = m2;
		
		} else {
		  AnnNode2 = m2.AnnNode;
		  //NovNode2 = m2.NovNode;
		}
		if(flag) {
		  Arc arc = AnnGraph.addArc(AnnNode1,AnnNode2);
		  AnnEdgesMap[arc] = err;
		  //arc = NovGraph.addArc(NovNode1,NovNode2);
		  //NovEdgesMap[arc] = err;
		  AnnExt = true;
		  //NovExt = true;
		} else {
		  //std::cout << "qui\n";
		  // Arc arc = NovGraph.addArc(NovNode1,NovNode2);
		  // NovEdgesMap[arc] = err;
		  // NovExt = true;
		  Arc arc = AnnGraph.addArc(AnnNode1,AnnNode2);
		  AnnEdgesMap[arc] = err;
		  //arc = NovGraph.addArc(NovNode1,NovNode2);
		  //NovEdgesMap[arc] = err;
		  AnnExt = true;
		  //NovExt = true;
		  //std::cout << "qui\n";
		  
		}
	      }
	    }
	    ++p2;
	  }
	  /*******
	   * End *
	   *******/
	  //if the MEM doesn't have a son, we link it to the END if it is a valid candidate as ending mem
	  //if(!AnnExt && !NovExt) {
	  if(!AnnExt) {
	    // std::cout << "inizo validEndIntrons\n";
	    std::pair<bool, int> endInfo = validEndIntrons(sg, m1);
	    // std::cout << "fine validEndIntrons\n";
	    bool endFlag = endInfo.first;
	    err = endInfo.second;
	    if(verbose)
	      std::cout << "end for" << m1.toStrCat()<< " -> "
			<< endFlag << " " << err <<"\n";
	    if(endFlag) {
	      Arc arc = AnnGraph.addArc(AnnNode1, AnnEnd);
	      AnnEdgesMap[arc] = err;
	      // arc = NovGraph.addArc(NovNode1,NovEnd);
	      // NovEdgesMap[arc] = err;
	    }
	  }
	}
      }
    }
  }
  arcsToDeleteAnn.clear();
  // per ogni nodo nel grafo annotato
  for(NodeIt x (AnnGraph); x!=lemon::INVALID; ++x) {
    // per ogni arco uscente da tale nodo
    for(OutArc XY (AnnGraph, x); XY!=lemon::INVALID; ++XY) {
      // y è il nodo verso cui va l'arco
      Node y = AnnGraph.target(XY);
      // per tutti gli archi uscenti da y
      for(OutArc YZ (AnnGraph, y); YZ!=lemon::INVALID; ++YZ) {
	// z è il nodo verso cui punta y
	Node z = AnnGraph.target(YZ);
	// per tutti gli archi entranti in z
	for(InArc XZ (AnnGraph, z); XZ!=lemon::INVALID; ++XZ) {
	  // x_ è il nodo da cui arriva tale nodo
	  Node x_ = AnnGraph.source(XZ);
	  // se x è uguale a x_
	  if(AnnGraph.id(x) == AnnGraph.id(x_)) {
	    // se il peso tra x e z è minore a quello di x e y + y e z
	    // aggiungo l'arco tra x e z a quelli da eliminare
	    if(AnnEdgesMap[XY] + AnnEdgesMap[YZ] <= AnnEdgesMap[XZ]) {
	      arcsToDeleteAnn.push_back(XZ);
	    }
	  }
	}
      }
    }
  }
  for(const InArc& a : arcsToDeleteAnn) {
    if(AnnGraph.valid(a)) {
      AnnGraph.erase(a);
    }
  }
  if(verbose) {
    save("Graph.dot");
  }
  
}

// metodo per la visita del grafo
std::list<std::pair<int, std::list<Mem> > > MemsGraph::visit(const SplicingGraph& sg) {

  if(verbose){
    
    int numberN = 0;
    for (lemon::ListDigraph::NodeIt n(AnnGraph); n != lemon::INVALID; ++n){
      numberN++;	
    }
    std::cout << "Number of annotated nodes: " << numberN << std::endl;  
    int numberA = 0;
    for (lemon::ListDigraph::ArcIt arc(AnnGraph);
	 arc != lemon::INVALID; ++arc){
      numberA++;
    }
    std::cout << "Number of arcs: " << numberA << std::endl;
    for (lemon::ListDigraph::ArcIt arc(AnnGraph);
	 arc != lemon::INVALID; ++arc){
      auto idNodeS = AnnGraph.id(AnnGraph.source(arc));
      auto idNodeT = AnnGraph.id(AnnGraph.target(arc));
      auto idArc = AnnGraph.id(arc);
      std::cout << "Arc " << idArc << " goes from node "
		<< AnnGraph.id(AnnGraph.source(arc)) << " "
		<< AnnNodesMap[AnnGraph.nodeFromId(idNodeS)].toStrCat()
		<< " to node "
		<< AnnGraph.id(AnnGraph.target(arc)) << " "
		<< AnnNodesMap[AnnGraph.nodeFromId(idNodeT)].toStrCat()
		<< ", weight = " << AnnEdgesMap[arc] << std::endl;
    }
  }
  // lista dei path da ritornare
  std::list<std::pair<int, std::list<Mem> > > paths;

  // liste di mems di appoggio
  std::list<Mem> AnnPath1;
  std::list<Mem> AnnPath2;
  std::list<Mem> NovPath;

  // check se annotato
  bool FoundAnnotated = false;

  // pesi
  int AnnW1 = K2 + 1;
  int AnnW2 = K2 + 1;
  int NovW = K2 + 1;

  //Visiting Annotated Graph
  // setup per cammino minimo con Dijkstra e heap di fibonacci
  lemon::Dijkstra<Graph, lemon::ListDigraph::ArcMap<int> >
    ::SetStandardHeap<FibH>
    ::SetHeap<FibH, FibM>
    ::Create AnnDijkstra (AnnGraph, AnnEdgesMap);

  // mappa nodi
  FibM AnnHCR (AnnGraph);
  // heap di fibonacci
  FibH AnnHeap (AnnHCR);
  // esecuzione di Dijkstra
  AnnDijkstra.heap(AnnHeap, AnnHCR);
  AnnDijkstra.run(AnnStart, AnnEnd);

  // se dijkstra arriva al nodo end
  if(AnnDijkstra.reached(AnnEnd)) {
    
    // calcolo della distanza tra la root e il nodo end
    AnnW1 = AnnDijkstra.dist(AnnEnd);
    // std::cout << "AnnW1 "<< AnnW1 << std::endl;
    // se tale distanza è accettabile
    if(AnnW1 <= K2) {

      // si marca il booleano segnalando di avere un path
      FoundAnnotated = true;

      // il path viene memorizzato in p
      Path p = AnnDijkstra.path(AnnEnd);

      // check per segnalare che è in analisi il nodo start col suo arco uscente 
      bool first = true;
      
      // si itera sugli archi del path
      for(Path::ArcIt it(p); it != lemon::INVALID; ++it) {

	// eliminazione arco dal nodo start
	if(first) {
	  AnnGraph.erase(it);
	  first = false;
	}
	
	// salvataggio dell'arco in it
	Arc e = it;

	// salvataggio della destinazione dell'arco
	Node target = AnnGraph.target(e);

	// estrazione della mem del nodo destinazione
	Mem m = AnnNodesMap[target];

	// finché non si arriva al nodo end si salvano le mem in una lista
	// di mem
	if(NovGraph.id(target) != NovGraph.id(AnnEnd)) {
	  AnnPath1.push_back(m);
	}
      }
      // si aggiunge alla lista dei path la coppia con la distanza complessiva
      // e il cammino
      paths.push_back(std::make_pair(AnnW1, AnnPath1));
    }
  }

  // seconda esecuzione senza foundannotated (???)
  AnnDijkstra.run(AnnStart, AnnEnd);
  if(AnnDijkstra.reached(AnnEnd)) {
    AnnW2 = AnnDijkstra.dist(AnnEnd);
    if(AnnW2 <= K2) {
      Path p = AnnDijkstra.path(AnnEnd);
      for(Path::ArcIt it(p); it != lemon::INVALID; ++it) {
	Arc e = it;
	Node target = AnnGraph.target(e);
	Mem m = AnnNodesMap[target];
	if(NovGraph.id(target) != NovGraph.id(AnnEnd)) {
	  AnnPath2.push_back(m);
	}
      }
      paths.push_back(std::make_pair(AnnW2, AnnPath2));
    }
  }
  if(verbose){
    for(auto p: paths){
      for(auto e: p.second){
	std::cout << e.toStrCat() << " ";
      }
      std::cout << "\n";
    }
  }
  //Visiting Novel Graph
  lemon::Dijkstra<Graph, lemon::ListDigraph::ArcMap<int> >
    ::SetStandardHeap<FibH>
    ::SetHeap<FibH,FibM>
    ::Create NovDijkstra (NovGraph, NovEdgesMap);
  FibM NovHCR (NovGraph);
  FibH NovHeap (NovHCR);
  NovDijkstra.heap(NovHeap, NovHCR);

  // studio con dijkstra del novel graph 
  NovDijkstra.run(NovStart,NovEnd);
  if(NovDijkstra.reached(NovEnd)) {
    NovW = NovDijkstra.dist(NovEnd);

    // sse il booleano non è marcato
    if(!FoundAnnotated) {
      if(NovW <= K2) {
	Path p = NovDijkstra.path(NovEnd);
	for(Path::ArcIt it(p); it != lemon::INVALID; ++it) {
	  Arc e = it;
	  Node target = NovGraph.target(e);
	  Mem m = NovNodesMap[target];
	  if(NovGraph.id(target) != NovGraph.id(NovEnd)) {
	    NovPath.push_back(m);
	  }
	}
      }
      // (????)
      if(paths.size() > 1)
	paths.pop_back();
      paths.push_front(std::make_pair(NovW, NovPath));
    }
  }
  return paths;
}

// salvataggio su file 
void MemsGraph::save(const std::string& s) {
  std::ofstream myfile;

  myfile.open("Ann"+s);

  std::string dot = "digraph G {\n graph [splines=true overlap=false]\n node [shape=ellipse, width=0.3, height=0.3]\n";
  for(NodeIt n (AnnGraph); n != lemon::INVALID; ++n) {
    dot += " " + std::to_string(AnnGraph.id(n)) + " [label=\"" +
      AnnNodesMap[n].toStrCat() + "\"];\n";
  }
  for(ArcIt a (AnnGraph); a != lemon::INVALID; ++a) {
    dot += " " + std::to_string(AnnGraph.id(AnnGraph.source(a))) + " -> " +
      std::to_string(AnnGraph.id(AnnGraph.target(a))) + "[label=\"" +
      std::to_string(AnnEdgesMap[a]) + "\"];\n";
  }
  dot += "}";

  myfile << dot;
  myfile.close();
    
  myfile.open("Nov"+s);

  dot = "digraph G {\n graph [splines=true overlap=false]\n node  [shape=ellipse, width=0.3, height=0.3]\n";
  for (NodeIt n (NovGraph); n != lemon::INVALID; ++n) {
    dot += " " + std::to_string(NovGraph.id(n)) + " [label=\"" +
      NovNodesMap[n].toStrCat() + "\"];\n";
  }
  for(ArcIt a (NovGraph); a != lemon::INVALID; ++a) {
    dot += " " + std::to_string(NovGraph.id(NovGraph.source(a))) + " -> " +
      std::to_string(NovGraph.id(NovGraph.target(a))) + "[label=\""
      + std::to_string(NovEdgesMap[a]) + "\"];\n";
  }
  dot += "}";

  myfile << dot;
  myfile.close();
}

std::vector<std::vector<Mem> > MemsGraph::startIntrons(const Mem& m, const
						       SplicingGraph& sg,
						       std::list<Mem>& listMem,
						       const int& start) {
 
  int id = sg.rank(m.t - 1);
  if(verbose)
    std::cout << "study starting introns for exon: "<< id << "\n";
  auto idIList = sg.getIntronMap();
  std::string intronText = sg.getTextIntrons();
  std::list<Mem> goodMemList;
  std::set<Mem> setMem;
  for(auto mem: listMem){
    if(mem.p + mem.l > m.p){
      mem.l = m.p - mem.p;
      setMem.insert(mem);
    }else{
      setMem.insert(mem);
    }
  }
  
  //std::cout << "s: " << start << "\n";
  std::vector<std::vector<Mem> > tot(listMem.size());
  std::vector<std::vector<Mem> > res;
  std::vector<int> counts(listMem.size(), 0);
  int count = 1;
  for(auto i: idIList){
    if((int)i.first.second == id){
      for(auto mem: setMem){
	for(auto intr: i.second){
	  if(sg.rankIntrons(mem.t - 1) == (int)intr){
	    if(count == 1){
	      if(mem.l > start){
		tot[0].push_back(Mem(mem.t, mem.p, start - 2, false));
		counts[0] = mem.l;
		count++;
	      }else{
		tot[0].push_back(mem);
		counts[0] = mem.l;
		count++;
	      }
	    }else{
	      for(unsigned int i = 0; i < counts.size(); i++){
		if(tot[i].size() == 0){
		  counts[i] += mem.l;
		  tot[i].push_back(mem);
		  count ++;
		  break;
		}else{
		  if(counts[i] + mem.l <= start &&
		     tot[i].back().p + tot[i].back().l <= mem.p ){
		    tot[i].push_back(mem);
		    counts[i] += mem.l;
		    count++;
		  }
		}
	      }
	    }
	  }
	}
      }
    }
  }
  // std::cout << "set\n";
 
  // for(auto vec:tot){
  //   for(auto m: vec){
  //     std::cout << m.toStrCat() << " ";
  //   }
  //   std::cout << "\n";
  // }
  
  // calcolo il massimo di copertura sulla read
  int part = 0;
  int max = 0;
  for(int i = 0; i < (int)tot.size(); i++){
    part = 0;
    for(int o = 0; o < (int)tot[i].size(); o++){
      part += tot[i][o].l;
    }
    if(part >= max)
      max = part;
  }
  // tengo solo le seguenze di mem migliori (in res)
  for(int i = 0; i < (int)tot.size(); i++){
    part = 0;
    for(int o = 0; o < (int)tot[i].size(); o++){
      part += tot[i][o].l;
    }
    if(part >= max)
      res.push_back(tot[i]);
  }

  
  return res;
}

std::vector<std::vector<Mem> > MemsGraph::middleIntrons(const Mem& m1,
							const Mem& m2,
							const SplicingGraph& sg,
							std::list<Mem>& listMem,
							const int& mid) {
  int id1 = sg.rank(m1.t - 1);
  int id2 = sg.rank(m2.t - 1);
  if(verbose)
    std::cout << "study middle introns for exons: "<< id1 << " " << id2 << "\n";
  std::set<Mem> setMem;

  // calcolo id con rank della posizione della mem sull'esone
  for(auto mem: listMem){
    if(mem.p < m1.p + m1.l && mem.p + mem.l > m1.p + m1.l){
      int iL = mem.l;
      int iP = mem.p;
      int shift = (m1.p + m1.l) - iP ;
      mem.p = m1.p + m1.l;
      mem.t += shift; 
      if(iL - shift < m2.p)
	mem.l = iL - shift;
      else
	mem.l = m2.p - 1;
      setMem.insert(mem);
    }else{
      setMem.insert(mem);
    }
  }
  auto idIList = sg.getIntronMap();
  std::string intronSequence = "";
  std::string intronText = sg.getTextIntrons();
  std::list<Mem> goodMemList;
  std::vector<std::vector<Mem> > tot(listMem.size());
  std::vector<std::vector<Mem> > res;
  std::vector<int> counts(listMem.size(), 0);
  int count = 1;

  if(id1 != id2){
    for(auto i: idIList){
      if((int)i.first.first == id1 && (int)i.first.second == id2){
	for(auto intr: i.second){
	  for(auto mem: setMem){
	    if(sg.rankIntrons(mem.t - 1) == (int)intr){
	      if(count == 1){
		if(mem.l > mid){
		  tot[0].push_back(Mem(mem.t, mem.p, mid - 2, false));
		  counts[0] = mem.l;
		  count++;
		}else{
		  tot[0].push_back(mem);
		  counts[0] = mem.l;
		  count++;
		}
	      }else{
		for(unsigned int i = 0; i < counts.size(); i++){
		  if(tot[i].size() == 0){
		    counts[i] += mem.l;
		    tot[i].push_back(mem);
		    count++;
		    break;
		  }else{
		    if(counts[i] + mem.l <= mid &&
		       tot[i].back().p + tot[i].back().l <= mem.p ){
		      tot[i].push_back(mem);
		      counts[i] += mem.l;
		      count++;
		    }
		  }
		}
	      }
	    }
	  }
	}
      }
    }
  }
  int part = 0;
  int max = 0;
  for(int i = 0; i < (int)tot.size(); i++){
    part = 0;
    for(int o = 0; o < (int)tot[i].size(); o++){
      part += tot[i][o].l;
    }
    if(part >= max)
      max = part;
  }
  // tengo solo le seguenze di mem migliori (in res)
  for(int i = 0; i < (int)tot.size(); i++){
    part = 0;
    for(int o = 0; o < (int)tot[i].size(); o++){
      part += tot[i][o].l;
    }
    if(part >= max)
      res.push_back(tot[i]);
  }
  return res;
}


std::vector<std::vector<Mem> > MemsGraph::endIntrons(const Mem& m,
						     const SplicingGraph& sg,
						     std::list<Mem>& listMem,
						     const int& end) {

  int id = sg.rank(m.t - 1);
  if(verbose)
    std::cout << "study ending introns for exon: "<< id << "\n";
  // richiamo la lista degli id degli introni e inizializzo
  auto idIList = sg.getIntronMap();
  std::string intronSequence = "";
  std::string intronText = sg.getTextIntrons();
  //std::cout << "end: " << intronText << std::endl;
  std::list<Mem> goodMemList;
  std::set<Mem> setMem;
  for(auto mem: listMem){
    if(mem.p < m.p + m.l && mem.p + mem.l > m.p + m.l){
      int iL = mem.l;
      int iP = mem.p - 1;
      int shift = (m.p + m.l) - iP - 1;
      mem.p = m.p + m.l;
      mem.t += shift; 
      mem.l = iL - shift;
      setMem.insert(mem);
    }else{
      setMem.insert(mem);
    }
  }
  // std::cout << "set:\n";
  // for(auto m: setMem){
  //   std::cout << m.toStrCat() << " ";
  // }
  // std::cout << "\n";
  std::vector<std::vector<Mem> > tot(setMem.size());
  std::vector<std::vector<Mem> > res;
  std::vector<int> counts(setMem.size(), 0);
  int count = 1;
  for(auto i: idIList){
    if((int)i.first.first == id){
      for(auto intr: i.second){
	for(auto mem: setMem){
	  if(sg.rankIntrons(mem.t - 1) == (int)intr && mem.p >= m.p + m.l){
	    if(count == 1){
	      if(mem.l > end){
		tot[0].push_back(Mem(mem.t, mem.p, end - 2, false));
		counts[0] += mem.l;
		count++;
	      }else{
		tot[0].push_back(mem);
		counts[0] += mem.l;
		count++;
	      }
	    }else{
	      for(unsigned int i = 0; i < counts.size(); i++){
		if(tot[i].size() == 0){
		  counts[i] += mem.l;
		  tot[i].push_back(mem);
		  count ++;
		  break;
		}else{
		  if(counts[i] + mem.l <= end &&
		     tot[i].back().p + tot[i].back().l <= mem.p){
		    tot[i].push_back(mem);
		    counts[i] += mem.l;
		    count++;
		  }
		}
	      }
	    }
	  }
	}
      }
    }
  }
  int part = 0;
  int max = 0;
  for(int i = 0; i < (int)tot.size(); i++){
    part = 0;
    for(int o = 0; o < (int)tot[i].size(); o++){
      part += tot[i][o].l;
    }
    if(part >= max)
      max = part;
  }
  // tengo solo le seguenze di mem migliori (in res)
  for(int i = 0; i < (int)tot.size(); i++){
    part = 0;
    for(int o = 0; o < (int)tot[i].size(); o++){
      part += tot[i][o].l;
    }
    if(part >= max)
      res.push_back(tot[i]);
  }

  return res;
}


std::vector<std::vector<Mem> > MemsGraph::totIntrons(std::list<Mem>& listMem) {
  std::vector<Mem> firstMemList;
  //std::vector<std::vector<Mem> > tot(listMem.size());
  std::vector<std::vector<Mem> > res;
  // ordino per indice su read
  if(listMem.size() > 1){
    listMem.sort([](const Mem& a, const Mem& b) -> bool{
      return a.p < b.p;
    });
  }
  int countfirst = 0;
  int firstp = listMem.front().p;
  for(auto m: listMem){
    if(m.p == firstp){
      countfirst++;
      firstMemList.push_back(m);
    }
  }
 
  //std::cout << countfirst << "\n";
  std::vector<std::vector<Mem> > tot(listMem.size() * countfirst);
  int count = 1;
  // salvo la copertura parziale (totale dei l dei mem) in counts
  std::vector<int> counts(listMem.size() * countfirst, 0);
  // set per lavorare su mem unici
  std::set<Mem> setMem(listMem.begin(), listMem.end());
  if(verbose){
    for(auto m: setMem){
      std::cout << m.toStr() << " ";
    }
    std::cout << std::endl;
  }
  // ciclo su ogni singola mem
  for(auto mem: listMem){
    if(verbose)
      std::cout << "analyze " << mem.toStr() << std::endl;
    // se è la prima la metto all'inizio di ogni possibile sequenza
    if(count == 1){
      int index = 0;
      unsigned int checkf = 0;
      for(unsigned int o = 0; o < listMem.size() * countfirst; o++){
	auto m = firstMemList[index];
	tot[o].push_back(m);
	counts[o] = m.l + m.p;
	if(countfirst != 1)
	  checkf++;
	if(checkf == listMem.size()){
	  index++;
	  checkf = 0;
	}
      }
      count++;
    }else{
      // ciclo su ogni sequenza
      for(unsigned int i = 0; i < counts.size(); i++){
	auto tmpmem = mem;
	// se la sequenza non è satura e se il mem in analisi è a destra
	// dell'ultimo della sequenza lo inserisco
	// e aggiorno in counts la lunghezza completa della sequenza
	if(counts[i] + mem.l < (int)read.size() &&
	   tot[i].back().p + tot[i].back().l <= mem.p){
	  // std::cout << "qui\n";
	  int l = tot[i].size();
	  if(i != 0 && !equalMems(tot[i-1][l-2], mem)
	     && !equalMems(tot[i-1].back(), mem)
	     && mem.p >= tot[i].back().p + tot[i].back().l -1){
	    
	    //tot[i].push_back(mem);
	    int p = 0;
	    if(mem.p - (tot[i].back().p + tot[i].back().l) > 0)
	      p = mem.p -( tot[i].back().p + tot[i].back().l);
	    tot[i].push_back(mem);
	    counts[i] += (mem.l + p);
	  }
	}else{
	  // se invece sforo la lunghezza della read aggiungendo il mem
	  // aggiungo il mem di lunghezza fino a saturare
	  if(counts[i] + mem.l >= (int)read.size() &&
	     tot[i].back().p + tot[i].back().l <= mem.p){
	    int gap = (int)read.size() - counts[i];
	    if(gap >= L){
	      auto nm = Mem(mem.t, mem.p, gap, true);
	      int l = tot[i].size();
	      if(i != 0 && !equalMems(tot[i-1][l-2], tot[i].back())
		 && !equalMems(tot[i-1].back(), nm)
		 && nm.p >= tot[i].back().p + tot[i].back().l -1){
		
	        int p = 0;
		if(nm.p - (tot[i].back().p + tot[i].back().l) > 0)
		  p = nm.p -( tot[i].back().p + tot[i].back().l);
		if(nm.l + tot[i].back().l == (int)read.size() - 1)
		  nm.l += 1;
		else if(nm.l + tot[i].back().l == (int)read.size())
		  nm.l -= 1;
		tot[i].push_back(nm);
		counts[i] += (nm.l + p);
		//mem = tmpmem;
	      }
	    }
	    // se invece non è completamente a sinistra
	  }else if(tot[i].back().p + tot[i].back().l > mem.p){
	    // calcolo l'overlap per capire di quanto shiftare
	    int shift = tot[i].back().p + tot[i].back().l - mem.p;
	    if(verbose)
	      std::cout << "shift: " << shift << "\n";
	    // scelgo chi ridimensionare in base alla lunghezza del mem
	    // scegliendo il più corto
	    // se è più lungo l'ultimo della sequenza
	    if(tot[i].back().l > mem.l && counts[i] < (int)read.size()){
	      if(verbose)
		std::cout << "start mem: " << mem.toStr() << " -> ";
	      // shifto il mem, lo accorcio e lo aggiungo
	      if(mem.l - shift >= L){
		mem.p += (shift);
		mem.t += (shift);
		mem.l -= shift;
	      }else{
		int shiftLessL = mem.l - L;
		
		if(counts[i] + mem.l <= (int)read.size()){
		  //std::cout << "" << "\n";
		  mem.l = L;
		  mem.p += (shiftLessL + 1);
		  mem.t += (shiftLessL + 1);

		  tot[i].back().l -= (shiftLessL);
		}else{
		  // std::cout << "check\n";
		  // std::cout <<  tot[i].back().toStr() << "\n";
		  int l = mem.l;
		  mem.l = l - (counts[i] + l - (int)read.size() - 1);
		  int readshift = tot[i].back().p + tot[i].back().l - mem.p;
		  //std::cout <<"red shi " << readshift<< "\n";
		  //int s = tot[i].back().l - (tot[i].back().l - mem.p -1);
		  tot[i].back().l -= readshift;
		  counts[i] -= readshift;
		}
	      }
	      if(verbose){
		std::cout << "mem: " << mem.toStr() << "\n ";
		std::cout << "index " << i << "\n";
	      }
	      if(i != 0 && verbose){
		std::cout << "i-1: " << tot[i-1].back().toStr() << "\n";
	      }
	      int l = tot[i].size();
	      if(((i != 0 && (!equalMems(tot[i-1][l-2], tot[i].back())
			      && !equalMems(tot[i-1].back(), mem))) || i == 0)
		 && mem.p >= tot[i].back().p + tot[i].back().l -1){
		
		int p = 0;
	        if(mem.p - (tot[i].back().p + tot[i].back().l) > 0)
		  p = mem.p -( tot[i].back().p + tot[i].back().l);
		tot[i].push_back(mem);
		counts[i] += (mem.l + p);
		mem = tmpmem;
	      }
	      // segnalo di riprendere il ciclo
	      
	      // altrimenti accorcio solo
	    }else if(tot[i].back().l < mem.l && counts[i] < (int)read.size()){
	      if(verbose)
		std::cout << "end mem: " << tot[i].back().toStr() << " -> ";
	      if(tot[i].back().l - shift >= L){
		tot[i].back().l -= shift;
		counts[i] -= shift;
	      }else{
		tot[i].back().l = L;
		int t = L + 1 - mem.p;
		int p = mem.p; 
		mem.p = (L+1);
		mem.t += t;
		mem.l -= (mem.p - p);
	      }
	      int l = tot[i].size();
	      if(i != 0 && (!equalMems(tot[i-1][l - 2], tot[i].back())
			    && !equalMems(tot[i-1].back(), mem))
		 && mem.p >= tot[i].back().p + tot[i].back().l -1) {	 
		
		int p = 0;
		if(mem.p - (tot[i].back().p + tot[i].back().l) > 0)
		  p = mem.p -( tot[i].back().p + tot[i].back().l);
		tot[i].push_back(mem);
		counts[i] += (mem.l + p);
		mem = tmpmem;
	      }
	      if(verbose)
		std::cout << "mem: " << tot[i].back().toStr() << " \n ";
	     
	    }
	    // restore della mem
	    mem = tmpmem;
	    if(verbose)
	      std::cout << "end analyze:\n";    
	  }
	}
      }
    }
    
    if(verbose){
      std::cout << "part: \n";
      int i = 0;
      for(auto o: tot){
	for(auto m: o){
	  std::cout << m.toStr() << " ";
	}
	std::cout << counts[i];
	std::cout << std::endl;
	i++;
      }
    }
    if(verbose)
      std::cout << "endf analyze:\n";    
  }
  // calcolo il massimo di copertura sulla read per ogni possile sequenza
  int part = 0;
  int max = 0;
  for(int i = 0; i < (int)tot.size(); i++){
    part = 0;
    for(int o = 0; o < (int)tot[i].size(); o++){
      part += tot[i][o].l;
    }
    if(part >= max)
      max = part;
  }
  if(verbose){
    std::cout << "tot: \n";
    for(auto o: tot){
      for(auto m: o){
	std::cout << m.toStr() << " ";
      }
      std::cout << std::endl;
    }
  }
  // tengo solo le seguenze di mem migliori (in res) (copertura pari a max)
  std::set<std::vector<Mem> > resSet(res.begin(), res.end());
  for(int i = 0; i < (int)tot.size(); i++){
    part = 0;
    for(int o = 0; o < (int)tot[i].size(); o++){
      part += tot[i][o].l;
    }
    if(part >= max)
      resSet.insert(tot[i]);
    //res.push_back(tot[i]);
  }
  if(verbose){
    std::cout << "res: \n";
    for(auto o: res){
      for(auto m: o){
	std::cout << m.toStr() << " ";
      }
      std::cout << std::endl;
    }
  }
  
  for(auto vec: resSet)
    res.push_back(vec);
  return res;
}


std::pair<bool, int> MemsGraph::checkMEMsIntrons(const SplicingGraph& sg,
						 const Mem& m1,
						 const Mem& m2) {
  
  // si usa la rank per capire da dove arrivano le mem, di cui uso l'indice di
  // partenza sulla linearizzazione
  int id1 = sg.rankIntrons(m1.t - 1);
  int id2 = sg.rankIntrons(m2.t - 1);

  // si estrapola le seq dei due esoni (posto che siano due)
  std::string intron1_text = sg.getIntron(id1);
  std::string intron2_text = sg.getIntron(id2);
  
  int err = -1;
  bool type = true;
  
  
  if(verbose) {
    std::cout << "Extending " << m1.toStrCat() << " with "
	      << m2.toStrCat() << std::endl;
  }
  
  // primo caso: si è sullo stesso esone (casi 1, 2, 3)
  if(id1 == id2) { //m1 and m2 in the same exon
    if(m2.p + m2.l > m1.p +m1.l &&
       m1.t < m2.t && m1.t + m1.l < m2.t + m2.l) {
      if(verbose) {
	std::cout << "same intron" << std::endl;
      }

      // gap sulla read
      
      long long int gapP = m2.p - m1.p - m1.l;

      // gap sull'esone 
      long long int gap_E = m2.t - m1.t - m1.l;
      if(verbose)
	std::cout << gapP << " " << gap_E << "\n";
      // casi 1  
      if(gapP >= 0 && gap_E >= 0) {
	
	// prima parte caso 1
	if(gapP == 0) {
	  // se il gap sull'esone non è trascurabile probabilmente si ha intorne
	  if(gap_E > K2) {
	    //Possible intron
	    if(verbose) {
	      std::cout << "Possible intron without overlap" << std::endl;
	    }
	    // peso a 0
	    err = 0;
	    type = false;
	  } else {
	    //Errors
	    if(verbose) {
	      std::cout << "Nothing" << std::endl;
	    }
	    // errore pari al gap (non avendo gap sulla read la distanza di
	    // edit è pari alla lunghezza del gap sull'esone)
	    err = gap_E;
	    type = true;
	  }
	} else if(abs(gapP - gap_E) <= K2) {
	  // seconda parte caso 1 
	  // Possible SNV (single nucleotide variant)
	  if(verbose) {
	    std::cout << "Nothing" << std::endl;
	  }

	  // si estrapola sequenza del gap su read ed esone
	  // std::string::substr(inizio, lunghezza)
	  std::string sub_P = read.substr(m1.p + m1.l - 1,
					  m2.p - m1.p - m1.l);
	  std::string sub_E = intron1_text.substr(m1.t + m1.l -
						  sg.selectIntrons(id1) - 1 - 1,
						  m2.t - m1.t - m1.l);

	
	  // errore (peso arco) come la distanza di edit tra i due gap
	  err = editDistance(sub_P, sub_E);
	  if(verbose){
	    std::cout << sub_P << "\n";
	    std::cout << sub_E << "\n";
	    std::cout << err << "\n";
	  }
	  type = true;
	}
      } else if(gapP <= 0 && gap_E <= 0) {
	// caso 2
	if(verbose) {
	  std::cout << "Nothing" << std::endl;
	}
	// errore (peso arco) come modulo tra le lunghezze dei due gap
	err = abs(gapP - gap_E);
	type = true;
	
      } else if(gapP <= 0 && gap_E > K2) {
	// caso 3
	// Possible intron
	if(verbose) {
	  std::cout << "Possible intron with overlap" << std::endl;
	}
	// peso a 0
	err = 0;
	type = false;
      } else {
	if(verbose) {
	  std::cout << "Nothing" << std::endl;
	}

	// altrimenti peso come somma dei gap
	err = abs(gapP) + abs(gap_E);
	type = true;
      }
    }
    // caso in cui i due mems sono da esoni diversi
    // casi 4, 5, 6
  } else { //m1 and m2 in different exons
    if(id2 > id1){
      // check che gli esoni vengano da due esoni collegati nello splicing graph
      if(verbose) {
	std::cout << "different introns" << std::endl;
      }
      // mem2 deve essere completamente a destra del mem1 sulla read
      if(m2.p + m2.l > m1.p + m1.l) {
	// lunghezza gap tra i mem sulla read
	long long int gapP = m2.p - m1.p - m1.l;
	
	// gapE1: lunghezza suff_Z(m1)
	long long int gapE1 = sg.selectIntrons(id1 + 1) + 1 - m1.t - m1.l;

	// gapE2: lunghezza pref_Z(m2)
	long long int gapE2 = m2.t - sg.selectIntrons(id2) - 1 - 1;
	// caso 4
	if(gapP <= 0) {
	  err = 0; //abs(gapP);
	  if(verbose) {
	    std::cout <<"gap negativo: "<<id1 << " " << id2 << " " 
		      << " " << gapE1 << " " << gapE2 << std::endl;
	  }
	  
	  // se inoltre non ho un arco novel ed entrambi i gap esonici nulli
	  if(gapE1 == 0 && gapE2 == 0) {
	    type = true;
	  }else if(err <= K2) {
	    //Possible Competing
	    type = false;
	  } else{
	    err = -1;
	  }
	} else {
	  // caso 5 
	  if(gapE1 == 0 && gapE2 == 0) {
	    //Possible insertion (only if annotated edge)
	    // se ho un arco non novel potrei avere un introne
	  
	    err = 0;
	    type = false;
	  
	  } else {
	    // caso 6
	    if(abs(gapP - (gapE1 + gapE2)) <= K2) {
	      //Possible SNV (single nucleotide variant)
	      if(verbose) {
		std::cout << "SNV" << std::endl;
	      }
	      
	      // si estrapolano le sequenze utili
	      // sequenza relativa la gap sulla read
	      std::string subP = read.substr(m1.p + m1.l - 1, gapP);

	      // sequenza relativa a suff_Z(m1)
	      std::string subE1 = intron1_text.substr(m1.t + m1.l -
						      sg.selectIntrons(id1) - 1 - 1,
						      gapE1);

	      // sequenza relativa a pref_Z(m2)
	      std::string subE2 = intron2_text.substr(0, gapE2);

	      // concatenazione di suff_Z(m1) e pref_Z(m2) 
	      std::string subE = subE1 + subE2;

	      // calcolo distanza di edit 
	      err = editDistance(subP, subE);
	      
	    
	      type = true;

	    } 
	  }
	}
      }
    }
  }
  if(verbose) {
    std::cout << type << " " << err << std::endl;
  }
  if(err > K2) {
    err = -1;
  }
  if(verbose) {
    std::cout << type << " " << err << std::endl;
  }

  // si restituisce la coppia del boolenao e dell'errore
  return std::make_pair(type, err);
}

std::pair<bool, int> MemsGraph::validStartIntrons(const SplicingGraph& sg,
						  const Mem& Mem) {

  // si controlla di essere nel range dello start proximity (K0)
  if(Mem.p <= K0) {
    // so inizializza l'errore a partire da K2 (errore sul gap permesso)
    int err = K2 + 1;

    // se si è all'inizio della read non si ha però errore
    if(Mem.p == 1) {
      err = 0;
    } else {
      // si usa la rank per identificare l'esone
      int id = sg.rankIntrons(Mem.t - 1);
      
      // si recupera la sequenza dell'esone 
      std::string intron_text = sg.getIntron(id);

      // si identifica la lunghezza del prefisso sulla read e la sequenza
      int l = Mem.p - 1;
      std::string sub_P = read.substr(0, l);

      // si inizializza per sottostringa esonica
      std::string sub_E;

      // calcolo della lunghezza del prefisso dell'esone
      long long int intron_pref_len = Mem.t - sg.selectIntrons(id) - 1 - 1;
      // caso in cui il prefisso sull'esone è più corto di quello della read
      if(intron_pref_len < l) {
	
	// calcolo della lunghezza del prefisso comune a esone e read
	long long int shared_pref_len = l - intron_pref_len;

	// si estrae la sottosequenza del prefisso dell'esone
	std::string intron_pref = intron_text.substr(0, intron_pref_len);

	// si setta l'errore alla lunghezza del prefisso sulla read
	err = l;
	
	// si prende la lista coi parents dell'esone dallo splicing graph
	std::list<int> parents = sg.getParentsIntron(id);

	// si itera su ogni nodo nella lista parents per stabilire se ho un mem
	// iniziale
	for(std::list<int>::iterator it=parents.begin();
	    it!=parents.end(); ++it) {
	  /**
	   * We look only at the father of the node,
	   * IF he is long enough, we get its suffix;
	   * ELSE we get all its label (without going further checking all 
	   * its parents)
	   **/
	  // si estrae la sequenza del nodo preso in considerazione
	  int par = *it;
	  std::string par_text = sg.getIntron(par);
	  // si guarda se il nodo padre copre la differenza tra il prefisso
	  // di read ed esone
	  if(sg.selectIntrons(par + 1) - shared_pref_len - sg.selectIntrons(par) - 1 >= 0) {

	    // se la copre prendiamo il suffisso del padre che permette di
	    // coprire il gap e lo sommiamo al prefisso dll'esone 
	    sub_E = par_text.substr(sg.selectIntrons(par+1) - shared_pref_len -
				    sg.selectIntrons(par) - 1, shared_pref_len)
	      + intron_pref;
	  } else {
	    // altrimenti si prende tutta la sequenza del nodo parent in analisi
	    // e la si somma al prefisso dell'esone
	    sub_E = par_text + intron_pref;
	  }

	  // si aggiorna l'errore tenendo il minimo
	  int curr_err = editDistance(sub_P, sub_E);
	  if(curr_err < err) {
	    err = curr_err;
	  }
	}
      } else {
	// se invece il prefisso dell'esone è più lungo di quello della read
	// sub_E è diventa il prefisso dell'esone che eccede oltre quello della
	// read
	sub_E = intron_text.substr(Mem.t - l - sg.selectIntrons(id) - 1 - 1, l);
	err = editDistance(sub_P, sub_E);
      }
    }
    // se alla fine l'errore è accettabile ritorno che è un valid start
    if(err <= K2) {
      return std::make_pair(true, err);
    }
  }
  // si ritorna la coppia indicante che non è un valid start (con un errore
  // sentinella)
  return std::make_pair(false, K2 + 1);
}

// metodo per verificare se un mem è un valido mem finale
std::pair<bool, int> MemsGraph::validEndIntrons(const SplicingGraph& sg,
						const Mem& Mem) {

  // se la mem sulla read termina in prossimità (prossimità K0) della fine della
  // read 
  if(Mem.p + Mem.l >= m - K0) {
    // errore impostato con K2
    int err = K2 + 1;

    // se la mem sulla read supera la fine errore a 0
    if(Mem.p + Mem.l == m + 1) {
      err = 0;
    } else {
      // si recupera l'id dell'esone
      int id = sg.rankIntrons(Mem.t - 1);

      // si recupera la sequenza dell'esone
      std::string intron_text = sg.getIntron(id);

      // lunghezza del suffisso della read oltre l'ultimo mem
      long long int l = m - (Mem.p + Mem.l ) + 1;

      // estrazione della sequenza di tale suffisso
      std::string sub_P = read.substr(Mem.p+Mem.l-1, l);

      // inizializzazione per la sottostringa sull'esone
      std::string sub_E;

      // calcolo della lunghezza del suffisso sull'esone a partire dalla fine
      // del mem (sfrutto l'inizio dell'esone successivo)
      long long int intron_suff_len = sg.selectIntrons(id + 1) - (Mem.t + Mem.l) + 1;

      // se il suffisso sull'esone è minore di quello sulla read
      if(intron_suff_len < l) {

	// si estrapola la lista dei sons
	std::list<int> sons = sg.getSonsIntron(id);

	// si calcola la lunghezza del suffisso condiviso tra read esone
	long long int shared_suff_len = l - intron_suff_len;

	// si inizalizza la stringa per il suffisso dell'esone
	std::string intron_suff;
	if(intron_suff_len == 0) {
	  // ovviamente se la lunghezza di tale suffisso è nulla si ha la
	  // stringa vuota
	  intron_suff = "";
	} else {
	  // altrimenti prendo la sequenza di tale suffisso 
	  intron_suff = intron_text.substr(Mem.t + Mem.l -sg.selectIntrons(id) - 1 - 1,
					 intron_suff_len);
	}
	// si imposta l'errore alla lunghezza del suffisso sulla read
	err = l;

	// si itera sui sons per vedere se posso stabilire che ho un mem finale
	for(std::list<int>::iterator it=sons.begin(); it!=sons.end(); ++it) {

	  // si estrapola la sequenza del son in considerazione
	  int son = *it;
	  std::string son_text = sg.getIntron(son);
	  // si calcola il nuovo suffisso partendo da quello dell'esone
	  sub_E = intron_suff + son_text.substr(0, shared_suff_len);
	  int curr_err = editDistance(sub_P, sub_E);

	  // si calcola l'errore con la distanza di edit
	  if(curr_err < err) {
	    err = curr_err;
	  }
	}
      } else {
	// se la lunghezza del suffisso sulla read è minore
	// sub_E diventa il suffisso dell'esone a partire dall'ultima mem di
	// lunghezza pari al suffisso della read
	sub_E = intron_text.substr(Mem.t + Mem.l - 
				   sg.selectIntrons(id) - 1 - 1, l);

	// si usa la distanza di edit per calclare l'errore
	err = editDistance(sub_P, sub_E);
      }
    }
    //std::cout << "K2: " << K2 << "\n";

    if(err <= K2) {
      //std::cout << "K2: " << K2 << "\n";
      // se l'errore alla fine è accettabile so che è un mem finale con
      // quell'errore
      return std::make_pair(true, err);
    }
  }
  // se non è un mem finale lo segnalo
  return std::make_pair(false,K2+1);
}
