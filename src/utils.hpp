#ifndef _UTILS_HPP_
#define _UTILS_HPP_

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>

#include <lemon/list_graph.h>

// struttura per definire un MEM
struct Mem {
  // tre campi della tripla che descrive la mem
  int t;
  int p;
  int l;

  
  int k;
  bool isNew;

  // bool per identificare la caregoria (esone/introne)
  bool intron;
  lemon::ListDigraph::Node AnnNode;
  lemon::ListDigraph::Node NovNode;
  // lemon::ListDigraph::Node IntronNode;
  // lemon::ListDigraph::Node FinalNode;

  // costruttori
  Mem() {
    t = 0;
    p = 0;
    l = 0;
    k = 0;
    isNew = true;
    intron = false;
  }

  Mem(int t_, int p_, int l_, bool i_) {
    t = t_;
    p = p_;
    l = l_;
    k = 0;
    isNew = true;
    intron = i_;
  }

  bool operator==(const Mem& m) const {
    return (m.t == t && m.p == p && m.l == l);
  }

  // setter di un nodo annotato
  void setAnnNode(lemon::ListDigraph::Node node_) {
    AnnNode = node_;
    isNew = false;
  }

  // setter di un nodo novel
  void setNovNode(lemon::ListDigraph::Node node_) {
    NovNode = node_;
    isNew = false;
  }

  // setter di un nodo introne
  // void setIntronNode(lemon::ListDigraph::Node node_) {
  //   IntronNode = node_;
  //   isNew = false;
  // }

  // void setFinalNode(lemon::ListDigraph::Node node_) {
  //   FinalNode = node_;
  //   isNew = false;
  // }
  // metodo to string
  std::string toStrCat() const {
    std::string category;
    if(intron)
      category = "1";
    else
      category = "0";
    return "(" + std::to_string(t) + "," +
      std::to_string(p) + "," +
      std::to_string(l) + "," +
      category + ")";
  }
  std::string toStr() const {
    return "(" + std::to_string(t) + "," +
      std::to_string(p) + "," +
      std::to_string(l) + ")";
  }
  // bool operator==(Mem const& b) const {
  //   return ((t == b.t) && (p == b.p) && (l == b.l)
  // 	    && (k == b.k) && (isNew == b.isNew) && (intron == b.intron));
  // }
  
};

int editDistance(const std::string&, const std::string&);
std::string reverseAndComplement(const std::string&);
bool compareMEMs(const Mem&, const Mem&);
bool compareMEMsLength(const Mem&, const Mem&);
#endif
