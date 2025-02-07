#include <iostream>
#include <string>
#include <fstream>
#include <algorithm>
#include <utility>
#include <list>
#include <filesystem>

#include <zlib.h>
#include <stdio.h>

#include "kseq.h"
#include "SplicingGraph.hpp"
#include "bMEM.hpp"
#include "utils.hpp"
#include "MEMsGraph.hpp"

KSEQ_INIT(gzFile, gzread)

void printHelp() {
  std::cout << "Usage: ASGAL [options] (required: -g -a -s -o)\n" << std::endl;
  std::cout << "Options:" << std::endl;
  std::cout << "  -g, --genome <path>" << std::endl;
  std::cout << "  -a, --annotation <path>" << std::endl;
  std::cout << "  -s, --sample <path>" << std::endl;
  std::cout << "  -o, --output <path>: output file" << std::endl;
  std::cout << "  -l, --L <int>: minimum lenght of MEMs used to build the alignments (default: 15)" << std::endl;
  std::cout << "  -e, --eps <int>: error rate, a value from 0 to 100 (default: 3)" << std::endl;
  std::cout << "  -h, --help: show this help message and exit" << std::endl;
  //std::cout << "  -v, --verbose: explain what is being done and save .dot" << std::endl;
  std::cout << "  -i, --intron: search new exons in introns (experimental)" << std::endl;
  std::cout << "  -v, --verbose: explain what is being done" << std::endl;
}

std::pair<char, std::list<std::pair<int, std::list<Mem> > > > analyzeRead(BackwardMEM& bm,
									  BackwardMEM& bmi,
									  const SplicingGraph& sg,
									  const std::string& read,
									  const int& L,
									  const int& eps,
									  const int& exsN,
									  const std::string head,
									  const std::string out,
									  const bool& intron_calc,
									  const bool& verbose) {
  std::list<Mem> mems = bm.getMEMs(read, L);
  // std::cout << "calcolo mem esoni\n";
  // if(read == test){
  //     std::cout << "check\n";
  //     for(auto m: memsI)
  // 	std::cout << m.toStr()  << " ";
  //   }
  // std::list<Mem> memsIntron;
  // FBtr0300326_e_3878_X_294812
  //bool test = false;
  // if(head == "FBtr0342963_e_1467_X_287023"){
  //   std::cout << "test\n";
  //   test = true;
  // }
  if(verbose && mems.size() != 0){
    std::cout << "size exon mem list: "<< mems.size() << std::endl;
  }
  // if(verbose && memsIntron.size() != 0){
  //   std::cout << "size intron mem list: "<< memsIntron.size() << std::endl;
  // }
  // if(verbose && memsIntron.size() != 0){
  //   std::cout << "normal introns mem rank: ";
  //   for(auto mem: memsIntron){
  //     std::cout << sg.rankIntrons(mem.t - 1) << " ";
  //   }
  //   std::cout << std::endl;
  // }
  bool longIntron = false;
  if(sg.getTextIntrons().size() > 500000){
    longIntron = true;
  }
  if(verbose) {
    for(const Mem& m : mems) {
      std::cout << m.toStrCat() << " ";
    }
    std::cout << std::endl;
  }
  if(verbose){
    for(const Mem& m : mems) {
      std::cout << sg.rank(m.t -1) << " ";
    }
    std::cout << std::endl;
  }
  std::list<std::pair<int, std::list<Mem> > > paths; // Path: [(weight, [mems])]
  MemsGraph mg (read, L, eps, exsN, verbose, intron_calc, longIntron);
  mg.build(sg, mems, bmi);
  paths = mg.visit(sg);
  
 
  // Reversed-and-complemented read
 
  std::string readRC = reverseAndComplement(read);
  std::list<Mem> memsRC = bm.getMEMs(readRC, L);
  // std::cout << "mem esoni RC\n";
  if(verbose && memsRC.size() != 0){
    std::cout << "size exon mem listRC: "<< memsRC.size();
    std::cout << std::endl;
  }
 
  if(verbose) {
    for(const Mem& m : memsRC) {
      std::cout << m.toStrCat() << " ";
    }
    std::cout << std::endl;
  }
  if(verbose){
    for(const Mem& m : memsRC) {
      std::cout << sg.rank(m.t -1) << " ";
    }
    std::cout << std::endl;
  }
  std::list<std::pair<int, std::list<Mem> > > pathsRC; // Path: [(weight, [mems])]
  
  MemsGraph mgRC (readRC, L, eps, exsN, verbose, intron_calc, longIntron);
  mgRC.build(sg, memsRC, bmi);
  pathsRC = mgRC.visit(sg);
  
  bool empty = paths.empty();
  bool emptyRC = pathsRC.empty();
  char strand = '/';
  if(!empty || !emptyRC) {
    if(!empty && emptyRC) {
      strand = '+';
    } else if(empty && !emptyRC) {
      paths = pathsRC;
      strand = '-';
    } else {
      if(paths.front().first <= pathsRC.front().first) {
	strand = '+';
      } else {
	paths = pathsRC;
	strand = '-';
      }
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
  return std::make_pair(strand, paths);
}

int main(int argc, char* argv[]) {
  std::string genomic;
  std::string annotation;
  std::string rna_seqs;
  int L = 0;
  int eps = -1;
  std::string out;
  bool verbose = false;
  bool intron_calc = false;
  // - Collecting command line parameters
  // -------------------------------------
  int c;
  while (1) {
    static struct option long_options[] =
      {
	{"genomic", required_argument, 0, 'g'},
	{"annotation", required_argument, 0, 'a'},
	{"sample", required_argument, 0, 's'},
	{"L",  required_argument, 0, 'l'},
	{"erate", required_argument, 0, 'e'},
	{"output", required_argument, 0, 'o'},
	{"help", no_argument, 0, 'h'},
	{"verbose", no_argument, 0, 'v'},
	{"intron", no_argument, 0, 'i'},
	{0, 0, 0, 0}
      };

    int option_index = 0;
    c = getopt_long(argc, argv, "g:a:s:l:e:o:hvi", long_options, &option_index);

    if (c == -1) {
      break;
    }

    switch(c) {
    case 'g':
      genomic = optarg;
      if(!std::filesystem::exists(genomic)){
	std::cerr << optarg << " not exist" << std::endl;
	exit(-1);
      }	
      break;
    case 'a':
      annotation = optarg;
      if(!std::filesystem::exists(annotation)){
	std::cerr << optarg << " not exist" << std::endl;
	exit(-1);
      }
      break;
    case 's':
      rna_seqs = optarg;
      if(!std::filesystem::exists(rna_seqs)){
	std::cerr << optarg << " not exist" << std::endl;
	exit(-1);
      }
      break;
    case 'l':
      L = std::stoi(optarg);
      break;
    case 'e':
      eps = std::stoi(optarg);
      break;
    case 'o':
      out = std::string(optarg);
      break;
    case 'v':
      verbose = true;
      break;
    case 'i':
      intron_calc = true;
      break;
    case 'h':
      printHelp();
      exit(EXIT_SUCCESS);
    default:
      printHelp();
      exit(EXIT_FAILURE);
    }
  }
  if(L == 0) {
    L = 15;
  }
  if(eps < 0 || eps > 100) {
    eps = 3;
  }

  // - Building splicing graph
  // --------------------------
  gzFile fastain = gzopen(genomic.c_str(), "r");
  kseq_t *reference = kseq_init(fastain);
  kseq_read(reference);

  SplicingGraph sg (reference->seq.s, annotation, verbose, intron_calc);
 
  if(verbose && intron_calc){
    auto test = sg.getIntronMap();
    for(auto i: test){
      std::cout << i.first.first << " " << i.first.second << ": ";
      for(auto o: i.second){
	std::cout << o << " ";
      }
      std::cout<<std::endl;
    }
  }
  kseq_destroy(reference);
  gzclose(fastain);
  std::string intronText;
  int exsN = sg.getExonsNumber();
  intronText = (intron_calc) ? sg.getTextIntrons() : "";
  
  // // - Setting up MEMs Index and out file
  // // ---------------------------------------------------
  BackwardMEM bm (sg.getText(), genomic, false);
  BackwardMEM bmi (intronText, genomic, true);
  
  std::ofstream outFile;
  outFile.open(out);
  fastain = gzopen(rna_seqs.c_str(), "r");
  std::pair<char, std::list<std::pair<int, std::list<Mem> > > > paths;

  kseq_t *seqs = kseq_init(fastain);
  int l;
  std::string head;
  std::string read;
  int i = 1;
  //std::cout << intronText.size() << "\n";
  // - Main loop: one iteration, one read
  // ---------------------------------------
  while ((l = kseq_read(seqs)) >= 0) {
    head = seqs->name.s;
    read = seqs->seq.s;
    
    if(verbose)
      std::cout << "\nread " << i - 1 << " ("<< read.size()<<") " << head << " " << read <<":\n";
    paths = analyzeRead(bm, bmi, sg, read, L, eps, exsN, head, out, intron_calc,
			verbose);
    if(paths.first != '/') {
      for(std::pair<int, std::list<Mem> > path : paths.second) {
	if(!path.second.empty()) {
	  int err = path.first;
	  outFile << paths.first << " " << head << " " << err << " ";
	  for(std::list<Mem>::iterator m=path.second.begin();
	      m!=path.second.end(); ++m) {
	    outFile << m->toStrCat() << " ";
	  }
	  if(paths.first == '+')
	    outFile << read;
	  else
	    outFile << reverseAndComplement(read);
	  outFile << "\n";
	}
      }
    }
    if(i%10000 == 0)
      std::cout << "Processed " << i << " genes." << std::endl;    
    ++i;
  }
  outFile.close();
  kseq_destroy(seqs);
  gzclose(fastain);
  
  return 0;
}
