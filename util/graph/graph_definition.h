/*****************************************************************************

    This file is part of Modelica C Compiler.

    Modelica C Compiler is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Modelica C Compiler is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Modelica C Compiler.  If not, see <http://www.gnu.org/licenses/>.

******************************************************************************/

/*! \file graph_definition.h
*   General library to work with set based graphs.
*   The main goal of this module was to abstract
*   classes that describe mathematical entities,
*   such as sets, graphs, etc.
*   If implementation has to change, a new concrete
*   class can be defined, but the client modules
*   of graph_definition will call abstract classes,
*   so no change is needed in the rest of the compiler.
*   Through typename template the abstract class uses
*   the new concrete implementation, therefore, requiring
*   little change appart from concrete class implementation.
*/

#ifndef GRAPH_DEFINITION_
#define GRAPH_DEFINITION_

#include <iostream>
#include <list>
#include <map>
#include <utility>

#include <boost/config.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/optional/optional.hpp>
#include <boost/unordered_set.hpp>

#include <ast/ast_types.h>
#include <ast/equation.h>
#include <util/table.h>

using namespace std;

#define Inf numeric_limits<int>::max()


/*-----------------------------------------------------------------------------------------------*/
/*-----------------------------------------------------------------------------------------------*/
// Abstract classes
/*-----------------------------------------------------------------------------------------------*/
/*-----------------------------------------------------------------------------------------------*/

int gcd(int a, int b){
  int c;

  do{
    c = a % b;
    if(c > 0){
      a = b;
      b = c;
    }
  } while (c != 0);

  return b;
}

int lcm(int a, int b){
  if(a < 0 || b < 0)
    return -1;

  return (a * b) / gcd(a, b);
}

template<template<typename Value, typename Hash = boost::hash<Value>, 
                  typename Pred = std::equal_to<Value>, 
                  typename Alloc = std::allocator<Value>> class CT>
struct IntervalImp1{
  int lo;
  int step;
  int hi;
  bool empty;

  IntervalImp1(){};
  IntervalImp1(bool isEmpty){ 
    lo = -1;
    step = -1;
    hi = -1;
    empty = isEmpty;
  };
  IntervalImp1(int vlo, int vstep, int vhi, bool vempty){ 
    if(vlo >= 0 && vstep > 0 && vhi >= 0){
      empty = vempty;
      lo = vlo;
      step = vstep;
     

      if(vlo <= vhi && !vempty && vhi < Inf){
        int rem = std::fmod(vhi - vlo, vstep); 
        hi = vhi - rem; 
      }

      else if(vlo <= vhi && !vempty && vhi == Inf){
        hi = Inf;
      }

      else{
        cerr << "Wrong values for subscript (check low <= hi)" << endl;
        empty = true;
      }
    }

    else{ 
      cerr << "Subscripts should be positive" << endl;
      lo = -1;
      step = -1;
      hi = -1;
      empty = true;
    }
  }

  int lo_(){
    return lo;
  }

  int step_(){
    return step;
  }

  int hi_(){
    return hi;
  }

  bool empty_() const{
    return empty;
  }

  bool isIn(int x){
    if(x < lo || x > hi || empty)
      return false;

    float aux = fmod(x - lo, step);
    if(aux == 0)
      return true;

    return false;
  } 

  IntervalImp1 cap(IntervalImp1 &inter2){
    int maxLo = max(lo, inter2.lo), newLo = -1;
    int newStep = lcm(step, inter2.step);
    int newEnd = min(hi, inter2.hi);

    if(!empty && !inter2.empty)
      for(int i = 0; i < newStep; i++){
        int res1 = maxLo + i;

        if(isIn(res1) && inter2.isIn(res1)){
          newLo = res1;
          break;
        }
      }

    else
      return IntervalImp1(true);


    if(newLo < 0)
      return IntervalImp1(true);

    return IntervalImp1(newLo, newStep, newEnd, false);
  }

  CT<IntervalImp1> diff(IntervalImp1 &i2){
    CT<IntervalImp1> res;
    IntervalImp1 capres = cap(i2);

    if(capres.empty){
      res.insert(*this);
      return res;
    }

    if(capres == *this){
      IntervalImp1 aux(true);
      res.insert(aux);
      return res;
    }

    // "Before" intersection
    if(lo < capres.lo){
      IntervalImp1 aux = IntervalImp1(lo, 1, capres.lo - 1, false);
      IntervalImp1 left = cap(aux);
      res.insert(left);
    }

    // "During" intersection
    if(capres.step <= (capres.hi - capres.lo)){
      int nInters = capres.step / step;
      for(int i = 1; i < nInters; i++){
        IntervalImp1 aux = IntervalImp1(capres.lo + i * step, capres.step, capres.hi, false);
        res.insert(aux);
      }  
    }

    // "After" intersection
    if(hi > capres.hi){
      IntervalImp1 aux = IntervalImp1(capres.hi + 1, 1, hi, false);
      IntervalImp1 right = cap(aux);
      res.insert(right);
    }
  
    return res;
  }

  int minElem(){
    return lo;
  }

  bool operator==(const IntervalImp1 &other) const{
    return (lo == other.lo) && (step == other.step) && (hi == other.hi) &&
           (empty == other.empty);
  }

  bool operator!=(const IntervalImp1 &other) const{
    return (lo != other.lo) || (step != other.step) || (hi != other.hi) ||
           (empty != other.empty);
  }

  size_t hash(){
    return lo;
  }
};

template<template<typename Value, typename Hash = boost::hash<Value>, 
                  typename Pred = std::equal_to<Value>, 
                  typename Alloc = std::allocator<Value>> class CT>
size_t hash_value(IntervalImp1<CT> inter){
  return inter.hash();
}

template<template<typename Value, typename Hash = boost::hash<Value>, 
                  typename Pred = std::equal_to<Value>, 
                  typename Alloc = std::allocator<Value>> class CT,
         typename IntervalImp, typename NumImp>
struct IntervalAbs{
  IntervalAbs(){}; 
  IntervalAbs(IntervalImp inter){
    i = inter;
  }
  IntervalAbs(NumImp lo, NumImp step, NumImp hi, bool emp){
    i = IntervalImp(lo, step, hi, emp);
  }
  IntervalAbs(bool emp){
    i = IntervalImp(emp);
  }
  
  NumImp lo_(){
    return i.lo_();
  }

  NumImp step_(){
    return i.step_();
  }
 
  NumImp hi_(){
    return i.hi_();
  }

  bool empty_(){
    return i.empty_();
  }

  bool isIn(NumImp x){
    return i.isIn(x);
  }

  IntervalAbs cap(IntervalAbs i2){
    return IntervalAbs(i.cap(i2.i));
  }

  CT<IntervalAbs> diff(IntervalAbs i2){
    CT<IntervalAbs> res;
    CT<IntervalImp> diffres = i.diff(i2.i);

    for(IntervalImp inter : diffres){
      IntervalAbs aux(inter);
      res.insert(aux);
    }

    return res;
  }

  NumImp minElem(){
    return i.minElem();
  }

  bool operator==(const IntervalAbs &other) const{
    return i == other.i;
  }
  
  bool operator!=(const IntervalAbs &other) const{
    return i != other.i;
  }

  size_t hash(){
    return i.lo_(); 
  }
 
  private:
  IntervalImp i; 
};

template<template<typename Value, typename Hash = boost::hash<Value>, 
                  typename Pred = std::equal_to<Value>, 
                  typename Alloc = std::allocator<Value>> class CT,
         typename IntervalImp, typename NumImp>
inline size_t hash_value(IntervalAbs<CT, IntervalImp, NumImp> inter){
  return inter.hash();
}

// MultiIntervals ---------------------------------------------------------------------------------

template<template<typename T, typename = allocator<T>> class CT1,
         template<typename Value, typename Hash = boost::hash<Value>, 
                  typename Pred = std::equal_to<Value>, 
                  typename Alloc = std::allocator<Value>> class CT2,
         typename IntervalImp, typename NumImp>
struct MultiInterImp1{
  CT1<IntervalImp> inters;
  typedef typename CT1<IntervalImp>::iterator intImpIt;

  MultiInterImp1(){
    CT1<IntervalImp> emptyRes;
    inters = emptyRes;  
  }
  MultiInterImp1(CT1<IntervalImp> is){
    inters = is;
  }

  CT1<IntervalImp> inters_(){
    return inters;
  }

  void addInter(IntervalImp i){
    inters.insert(inters.end(), i);
  }

  bool empty(){
    intImpIt it = inters.begin();

    if(inters.empty())
      return true;

    while(it != inters.end()){
      if(!(*it).empty_())
        return false;

      ++it;
    }
  
    return true;
  }

  MultiInterImp1 cap(MultiInterImp1 &mi2){
    CT1<IntervalImp> res;
    intImpIt itres = res.begin();

    intImpIt it1 = inters.begin();
    intImpIt it2 = mi2.inters.begin();
    int minLength = min(inters.size(), mi2.inters.size());
    for(int i = 0; i < minLength; ++i){
      IntervalImp capres = (*it1).cap(*it2);
     
      if(capres.empty_() && !(*it1).empty_() && !(*it2).empty_()){
        CT1<IntervalImp> aux;
        return MultiInterImp1(aux);
      }
 
      itres = res.insert(itres, capres);
      ++itres;

      ++it1;
      ++it2;    
    }

    return MultiInterImp1(res);
  }

  CT2<MultiInterImp1> diff(MultiInterImp1 &mi2){
    MultiInterImp1 capres = cap(mi2);

    CT2<MultiInterImp1> resmi;
    typename CT2<MultiInterImp1>::iterator itresmi = resmi.begin();

    if(capres.empty()){
      resmi.insert(*this);
      return resmi;
    }

    if(inters == capres.inters)
      return resmi;

    intImpIt it1 = inters.begin();
    intImpIt itcap = capres.inters.begin();
    CT1<CT2<IntervalImp>> diffs;
    typename CT1<CT2<IntervalImp>>::iterator itdiffs = diffs.begin();

    while(it1 != inters.end()){
      itdiffs = diffs.insert(itdiffs, (*it1).diff(*itcap));
 
      ++it1;
      ++itcap;
      ++itdiffs;
    }

    it1 = inters.begin();
    ++it1;
    itdiffs = diffs.begin();
    
    while(itdiffs != diffs.end()){
      CT2<IntervalImp> aux = *itdiffs;
      typename CT2<IntervalImp>::iterator itaux = aux.begin();

      int i = distance(diffs.begin(), itdiffs);

      while(itaux != aux.end()){
        CT1<IntervalImp> resi;
        intImpIt itresi = resi.begin();

        itcap = capres.inters.begin();

        if(i > 0){
          for(int j = 0; j < i; j++){
            itresi = resi.insert(itresi, *itcap);         
            ++itresi; 

            ++itcap;
          }
        }

        itresi = resi.insert(itresi, *itaux);
        ++itresi; 

        intImpIt auxit1 = it1;
        while(auxit1 != inters.end()){
          itresi = resi.insert(itresi, *auxit1);
          ++itresi;
        
          ++auxit1;
        }

        itresmi = resmi.insert(itresmi, MultiInterImp1(resi));
        ++itresmi;

        ++itaux;
      }

      ++it1;
      ++itdiffs;
    }

    return resmi;
  }

  CT1<NumImp> minElem(){
    intImpIt it = inters.begin();

    CT1<NumImp> res;
    typename CT1<NumImp>::iterator itres = res.begin();

    while(it != inters.end()){
      itres = res.insert(itres, (*it).minElem());
      ++itres;

      ++it;
    }    

    return res;
  }

  bool operator==(const MultiInterImp1 &other) const{
    return inters == other.inters;
  }

  bool operator!=(const MultiInterImp1 &other) const{
    return inters != other.inters;
  }

  size_t hash(){
    return inters.size(); 
  }
};

template<template<typename T, typename = allocator<T>> class CT1,
         template<typename Value, typename Hash = boost::hash<Value>, 
                  typename Pred = std::equal_to<Value>, 
                  typename Alloc = std::allocator<Value>> class CT2,
         typename IntervalImp, typename NumImp>
size_t hash_value(MultiInterImp1<CT1, CT2, IntervalImp, NumImp> mi){
  return mi.hash();
}

template<template<typename T, typename = allocator<T>> class CT1,
         template<typename Value, typename Hash = boost::hash<Value>, 
                  typename Pred = std::equal_to<Value>, 
                  typename Alloc = std::allocator<Value>> class CT2,
          typename MultiInterImp, typename IntervalImp, typename NumImp>
struct MultiInterAbs{
  MultiInterAbs(){
    CT1<IntervalImp> aux;
    multiInterImp = aux;
  };
  MultiInterAbs(MultiInterImp mi){
    multiInterImp = mi;
  }
  MultiInterAbs(CT1<IntervalImp> ints){
    multiInterImp = MultiInterImp(ints);
  };

  bool empty(){
    return multiInterImp.empty();
  };

  CT1<IntervalImp> inters_(){
    return multiInterImp.inters_();
  }

  void addInter(IntervalImp i){
    multiInterImp.addInter(i);
  }

  MultiInterAbs cap(MultiInterAbs &mi2){
    return MultiInterAbs(multiInterImp.cap(mi2.multiInterImp));
  }

  CT2<MultiInterAbs> diff(MultiInterAbs &mi2){
    CT2<MultiInterImp> diffRes = multiInterImp.diff(mi2.multiInterImp);
    typename CT2<MultiInterImp>::iterator it = diffRes.begin();
    CT2<MultiInterAbs> res;

    while(it != diffRes.end()){
      MultiInterAbs aux(*it);
      res.insert(aux);

      ++it;
    }

    return res;
  };

  CT1<NumImp> minElem(){
    return multiInterImp.minElem();
  }

  bool operator==(const MultiInterAbs &other) const{
    return multiInterImp == other.multiInterImp;
  }

  bool operator!=(const MultiInterAbs &other) const{
    return multiInterImp != other.multiInterImp;
  }

  size_t hash(){
    return multiInterImp.hash();
  }

  private:
  MultiInterImp multiInterImp;
};

template<template<typename T, typename = allocator<T>> class CT1,
         template<typename Value, typename Hash = boost::hash<Value>, 
                  typename Pred = std::equal_to<Value>, 
                  typename Alloc = std::allocator<Value>> class CT2,
          typename MultiInterImp, typename IntervalImp, typename NumImp>
size_t hash_value(MultiInterAbs<CT1, CT2, MultiInterImp, IntervalImp, NumImp> mi){
  return mi.hash();
}

// Atomic sets ------------------------------------------------------------------------------------

/*
template <template<typename T1, typename = allocator<T1>> class CT,
          typename MultiInterImp, typename NumImp>
struct AtomSetImp1{
  MultiInterImp aset;

  AtomSetImp1(){
    MultiInterImp emptyRes;
    aset = emptyRes;
  }
  AtomSetImp1(MultiInterImp as){
    aset = as;
  }

  MultiInterImp aset_(){
    return aset;
  }

  bool empty(){
    return aset.empty();  
  }

  AtomSetImp1 cap(AtomSetImp1 &aset2){
    AtomSetImp1 aux(aset.cap(aset2.aset));
    return aux;
  }

  CT<AtomSetImp1<CT, MultiInterImp, NumImp>> diff(AtomSetImp1 &aset2){
    CT<AtomSetImp1<CT, MultiInterImp, NumImp>> res;
    typename CT<AtomSetImp1<CT, MultiInterImp, NumImp>>::iterator itres = res.begin();

    CT<MultiInterImp> atomicDiff = aset.diff(aset2.aset);

    if(atomicDiff.empty()){
      CT<AtomSetImp1<CT, MultiInterImp, NumImp>> emptyRes;
      return emptyRes;
    } 

    else{
      typename CT<MultiInterImp>::iterator it = atomicDiff.begin();
      while(it != atomicDiff.end()){
        itres = res.insert(itres, AtomSetImp1(*it)); 
        ++itres;

        ++it;
      }
    }

    return res;
  }

  CT<NumImp> minElem(){
    return aset.minElem();
  }

  bool operator==(const AtomSetImp1 &other) const{
    return aset == other.aset;
  }
};

template <template<typename T, typename = allocator<T>> class CT,
          typename ASetImp, typename MultiInterImp, typename NumImp>
struct AtomSetAbs{
  AtomSetAbs(){
    ASetImp aux;
    as = aux;
  };
  AtomSetAbs(ASetImp ass){
    as = ass;
  }
  AtomSetAbs(MultiInterImp mi){
    ASetImp aux(mi);
    as = aux;
  }

  MultiInterImp aset_(){
    return as.aset_();
  }

  bool empty(){
    return as.empty();  
  }

  AtomSetAbs cap(AtomSetAbs &aset2){
    AtomSetAbs aux(as.cap(aset2.as));
    return aux;
  }

  CT<AtomSetAbs<CT, ASetImp, MultiInterImp, NumImp>> diff(AtomSetAbs &aset2){
    CT<ASetImp> diffRes = as.diff(aset2.as);
    typename CT<ASetImp>::iterator it = diffRes.begin();
    CT<AtomSetAbs> res;
    typename CT<AtomSetAbs>::iterator itRes = res.begin();

    while(it != diffRes.end()){
      itRes = res.insert(itRes, AtomSetAbs(*it));
      ++itRes;

      ++it;
    }

    return res;
  }
 
  CT<NumImp> minElem(){
    return as.minElem();
  }

  bool operator==(const AtomSetAbs &other) const{
    return as == other.as;
  }

  private:
  ASetImp as;
};

template<template<class T, class Alloc = allocator<T>> class CT, 
         typename ASetImp>
struct SetImp1{
  typedef CT<ASetImp> setType;
  setType asets;
 
  SetImp1(){
    setType aux;
    asets = aux;
  }
  SetImp1(setType ss){
    asets = ss;
  }

  setType asets_(){
    return asets;
  }

  bool empty(){
    typename setType::iterator it = asets.begin();

    if(asets.empty())
      return true;

    while(it != asets.end()){
      if(!(*it).empty())
        return false;

      ++it;
    }
  
    return true;
  }

  void addAtomSet(ASetImp &aset2){
    if(!aset2.empty())
      asets.insert(asets.end(), aset2);
  }

  void addAtomSets(setType &sets2){
    typename setType::iterator it2 = sets2.begin();

    while(it2 != sets2.end()){
      addAtomSet(*it2);
      ++it2;
    }
  }

  SetImp1 cap(SetImp1 &set2){
    if(empty() || set2.empty()){
      setType emptyRes;
      return emptyRes; 
    }
    
    setType res;
    typename setType::iterator itres = res.begin();
    typename setType::iterator it1 = asets.begin();

    while(it1 != asets.end()){
      typename setType::iterator it2 = set2.asets.begin();

      while(it2 != set2.asets.end()){
        itres = res.insert(itres, (*it1).cap(*it2));
        ++itres;

        ++it2;
      }
     
      ++it1;
    }

    return res;
  }

  SetImp1 diff(SetImp1 &set2){
    SetImp1 res;
    setType capres = cap(set2).asets; 

    if(!capres.empty()){
      typename setType::iterator it1 = asets.begin();

      while(it1 != asets.end()){
        setType aux(1, *it1);
 
        typename setType::iterator it2 = capres.begin();
        while(it2 != capres.end()){
          SetImp1 newSets;

          typename setType::iterator itaux = aux.begin();
          while(itaux != aux.end()){
            setType diffres = (*itaux).diff(*it2);
            newSets.addAtomSets(diffres);

            ++itaux;
          }

          aux = newSets.asets;

          ++it2;
        }

        res.addAtomSets(aux);

        ++it1;
      }
    }

    else
      res.addAtomSets(asets);

    return res; 
  }

  SetImp1 cup(SetImp1 &set2){
    SetImp1 resDiff = diff(set2);

    if(!resDiff.empty())
      addAtomSets(resDiff);

    return SetImp1(asets);
  }

  void simplify(){
    typename setType::iterator it1 = asets.begin();

    while(it1 != asets.end()){
      typename setType::iterator it2 = asets.begin();

      while(it2 != asets.end()){
      }
    }
  }

  bool operator==(const SetImp1 &other) const{
    SetImp1 aux1 = *this;
    SetImp1 aux2 = other;
    auto it1 = asets.begin();

    if(aux1.empty() && aux2.empty())
      return true;

    while(it1 != asets.end()){
      bool aux = false;

      auto it2 = other.asets.begin(); 
      while(it2 != other.asets.end()){
        if(*it1 == *it2)
          aux = true; 
      
        ++it2;
      }

      if(!aux)
        return false;

      ++it1;
    }

    return true;
  }
};

template <template<class T, class Alloc = allocator<T>> class CT, 
          typename SetImp, typename ASetImp>
struct SetAbs{
  SetAbs(){
    SetImp aux;
    set = aux;
  }
  SetAbs(SetImp ss){
    set = ss;
  }
  SetAbs(CT<ASetImp> ass){
    SetImp aux(ass);
    set = aux;
  }

  bool empty(){
    return set.empty(); 
  }

  CT<ASetImp> asets_(){
    return set.asets_();
  }

  void addAtomSet(ASetImp &aset2){
    set.addAtomSet(aset2); 
  }

  void addAtomSets(CT<ASetImp> &sets2){
    set.addAtomSets(sets2);
  }

  SetAbs cap(SetAbs &set2){
    return SetAbs(set.cap(set2.set));
  }

  SetAbs diff(SetAbs &set2){
    return SetAbs(set.diff(set2.set)); 
  }

  SetAbs cup(SetAbs &set2){
    return SetAbs(set.cup(set2.set));
  }

  bool operator==(const SetAbs &other) const{
    return set == other.set;
  }

  private:
  SetImp set;
};

template <template<typename T, typename = std::allocator<T>> class CT>
struct LExprImp1{
  float m;
  float h;

  LExprImp1(){};
  LExprImp1(float mm, float hh){
    m = mm;
    h = hh;
  }

  float m_(){
    return m;
  } 

  float h_(){
    return h;
  }

  LExprImp1 compose(LExprImp1 &le2){
    float newM = le2.m * m;
    float newH = m * le2.h + h;

    return LExprImp1(newM, newH);
  }

  LExprImp1 invExpr(){
    float newM, newH;

    if(m != 0){
      newM = m;
      newH = -h / m;
      return LExprImp1(newM, newH);
    }

    cerr << "The expression is not injective";
    return LExprImp1();
  }
};

template<typename LExprImp, typename NumImp>
struct LExprAbs{
  LExprAbs(){};
  LExprAbs(LExprImp lexpr){
    le = lexpr;
  }

  NumImp m_(){
    return le.m_();
  }

  NumImp h_(){
    return le.h_();
  }

  LExprAbs compose(LExprAbs &le2){
    return LExprAbs(le.compose(le2.le));
  } 

  LExprAbs invExpr(){
    return LExprAbs(le.invExpr());
  }

  private:
  LExprImp le;
};

template <template<typename T, typename = std::allocator<T>> class CT,
          typename IntervalImp, typename LExprImp>
struct LMIntImp1{
  IntervalImp dom;
  LExprImp expr;

  LMIntImp1(){};
  LMIntImp1(IntervalImp d, LExprImp e){
    dom = d;
    expr = e;
  }

  IntervalImp dom_(){
    return dom;
  }

  LExprImp expr_(){
    return expr;
  }

  IntervalImp image(IntervalImp &set){
    IntervalImp capRes = dom.cap(set);

    if(capRes.empty())
      return IntervalImp(true);

    int newLo = expr.m_() * capRes.lo_() + expr.h_();
    int newStep = expr.m_() * capRes.step_();
    int newHi = expr.m_() * capRes.hi_() + expr.h_();

    IntervalImp res(newLo, newStep, newHi, false);
    return res;    
  }

  IntervalImp preImage(IntervalImp &set){
    IntervalImp capRes = image(dom).cap(set);
    IntervalImp res;

    if(expr.m_() == 0){
      if(capRes.isIn(expr.h_()))
        res = dom; 
    }

    else{
      LMIntImp1 aux(capRes, expr.invExpr());
      res = aux.image();
    }

    return res;
  }

  LMIntImp1 compose(LMIntImp1 &lm2){
    IntervalImp newDom;
    LExprImp newExpr = expr.compose(lm2.expr);

    if(lm2.m_() == 0){
      if(dom.isIn(lm2.h_()))
        newDom = lm2.dom;
    }

    else{
      IntervalImp capRes = dom.cap(lm2.image());
      LMIntImp1 invLMInt(capRes, lm2.invExpr()); 
      newDom = invLMInt.dom();
    }

    return LMIntImp1(newDom, newExpr);
  }

  LMIntImp1 miniInv(){
    IntervalImp newDom;
    LExprImp newExpr;

    if(expr.m_() == 0){
      newDom = image(dom);
      newExpr = LExprImp(0, dom.lo_()); 
    }

    else{
      newDom = image(dom);
      newExpr = expr.invExpr();
    }

    return LMIntImp1(newDom, newExpr);
  }
};

template<typename LMIntImp, typename IntervalImp, typename LExprImp>
struct LMIntAbs{
  LMIntAbs(){};
  LMIntAbs(LMIntImp lmap){
    lm = lmap;
  }

  IntervalImp dom_(){
    return lm.dom_();
  }

  LExprImp expr_(){
    return lm.expr_();
  }

  IntervalImp image(IntervalImp &set){
    return lm.image(set);
  }  

  IntervalImp preImage(IntervalImp &set){
    return lm.preImage(set);
  } 

  LMIntAbs compose(LMIntAbs &lm2){
    return LMIntAbs(lm.compose(lm2.lm));
  }

  LMIntAbs miniInv(){
    return LMIntAbs(lm.miniInv());
  } 

  private:
  LMIntImp lm;
};

template<template<typename T, typename = allocator<T>> class CT,
         typename MultiInterImp, typename IntervalImp, typename LExprImp, typename LMIntImp>
struct LMMultiIntImp1{
  MultiInterImp dom;
  CT<LExprImp> expr; 

  LMMultiIntImp1(){};
  LMMultiIntImp1(MultiInterImp d, CT<LExprImp> e){
    if(e.size() >= d.ints_().size()){
      dom = d;
      expr = e;
    }

    else
      cerr << "Expression dimension should be larger than domain dimension";
  }

  MultiInterImp dom_(){
    return dom;
  }

  CT<LExprImp> expr_(){
    return expr;
  }

  MultiInterImp image(MultiInterImp &set){
    typename CT<IntervalImp>::iterator itdom = dom.ints_().begin();
    typename CT<LExprImp>::iterator itexpr = expr.begin();

    CT<IntervalImp> res;
    typename CT<IntervalImp>::iterator itres = res.begin();

    while(itdom != dom.ints().end()){
      LMIntImp aux(*itdom, *itexpr);
      itres = res.insert(itres, aux.image(set));
      ++itres;

      ++itdom;
      ++itexpr;
    }

    return MultiInterImp(res);
  } 

  MultiInterImp preImage(MultiInterImp &set){
    typename CT<IntervalImp>::iterator itdom = dom.ints_().begin();
    typename CT<LExprImp>::iterator itexpr = expr.begin();

    CT<IntervalImp> res;
    typename CT<IntervalImp>::iterator itres = res.begin();

    while(itdom != dom.ints().end()){
      LMIntImp aux(*itdom, *itexpr);
      itres = res.insert(itres, aux.preImage(set));
      ++itres;

      ++itdom;
      ++itexpr;
    }

    return MultiInterImp(res);
  }

  LMMultiIntImp1 compose(LMMultiIntImp1 &lm2){
    typename CT<IntervalImp>::iterator itdom1 = dom.ints_().begin();
    typename CT<IntervalImp>::iterator itdom2 = lm2.dom.ints_().begin();
    typename CT<LExprImp>::iterator itexpr1 = expr.begin();
    typename CT<LExprImp>::iterator itexpr2 = lm2.expr.begin();

    CT<IntervalImp> resDom;
    typename CT<IntervalImp>::iterator itresdom = resDom.begin();
    CT<LExprImp> resExpr;
    typename CT<IntervalImp>::iterator itresexpr = resExpr.begin();

    while(itdom2 != lm2.dom.ints().end()){
      LMIntImp aux1(*itdom1, *itexpr1);
      LMIntImp aux2(*itdom2, *itexpr2);

      LMIntImp res = aux1.compose(aux2);
      itresdom = resDom.insert(itresdom, res.dom_());
      ++itresdom;
      itresexpr = resExpr.insert(itresexpr, res.expr_());
      ++itresexpr;

      ++itdom1;
      ++itdom2;
      ++itexpr1;
      ++itexpr2;
    }

    MultiInterImp aux(resDom);
    return LMMultiIntImp1(resDom, resExpr);
  } 

  LMMultiIntImp1 miniInv(){
    typename CT<IntervalImp>::iterator itdom = dom.ints_().begin();
    typename CT<LExprImp>::iterator itexpr = expr.begin();

    CT<IntervalImp> resDom;
    typename CT<IntervalImp>::iterator itresdom = resDom.begin();
    CT<LExprImp> resExpr;
    typename CT<LExprImp>::iterator itresexpr = resExpr.begin();

    while(itdom != dom.ints_().end()){
      LMIntImp res(*itdom, *itexpr);
      LMIntImp inv = res.miniInv();

      itresdom = resDom.insert(itresdom, inv.dom_());
      ++itresdom;
      itresexpr = resExpr.insert(itresexpr, inv.expr_());
      ++itresexpr;

      ++itdom;
      ++itexpr;
    }

    MultiInterImp aux(resDom);
    return LMMultiIntImp1(aux, resExpr);
  }
};

template<template<typename T, typename = allocator<T>> class CT,
         typename LMMultiIntImp, typename MultiInterImp, typename LExprImp>
struct LMMultiIntAbs{
  LMMultiIntAbs(){};
  LMMultiIntAbs(LMMultiIntImp lmap){
    lm = lmap;
  }

  MultiInterImp dom_(){
    return lm.dom_();
  }

  CT<LExprImp> expr_(){
    return lm.expr_();
  }

  MultiInterImp image(MultiInterImp &set){
    return lm.image();
  } 

  MultiInterImp preImage(MultiInterImp &set){
    return lm.preImage();
  }

  LMMultiIntAbs compose(LMMultiIntAbs &lm2){
    return LMMultiIntAbs(lm.compose(lm2.lm));
  } 

  LMMultiIntAbs miniInv(){
    return LMMultiIntAbs(lm.miniInv());
  }

  private:
  LMMultiIntImp lm;
};

template<template<typename T, typename = allocator<T>> class CT,
         typename SetImp, typename ASetImp, typename MultiInterImp, typename LExprImp>
struct LMSetImp1{
  SetImp dom;
  CT<CT<LExprImp>> expr;

  LMSetImp1(){};
  LMSetImp1(SetImp d, CT<LExprImp> e){
    if(e.size() >= d.ints_().size()){
      dom = d;
      expr = e;
    }

    else
      cerr << "Expression dimension should be larger than domain expression";
  }

  SetImp image(SetImp &set){
    CT<ASetImp> adom = dom.asets_(); 
    typename CT<ASetImp>::iterator itadom = adom.begin();
    typename CT<CT<LExprImp>>::iterator itexpr = expr.begin();

    while(itadom != adom.end()){
      MultiInterImp mi = (*itadom).aset_();
    }
  }
};
*/

/*-----------------------------------------------------------------------------------------------*/
/*-----------------------------------------------------------------------------------------------*/
// Concrete classes
/*-----------------------------------------------------------------------------------------------*/
/*-----------------------------------------------------------------------------------------------*/

typedef int NI;

typedef IntervalImp1<boost::unordered_set> IntervalImp;
typedef IntervalAbs<boost::unordered_set, IntervalImp, NI> Interval;

ostream &operator<<(ostream &out, Interval &i){
  out << "[" << i.lo_() << ":" << i.step_() << ":" << i.hi_() << "]";
  return out;
}

typedef MultiInterImp1<list, boost::unordered_set, Interval, NI> MultiInterImp;
typedef MultiInterAbs<list, boost::unordered_set, MultiInterImp, Interval, NI> MultiInterval;

ostream &operator<<(ostream &out, MultiInterval &mi){
  list<Interval> is = mi.inters_();
  list<Interval>::iterator it = is.begin();

  if(is.size() == 0)
    return out;

  if(is.size() == 1){
    out << *it;
    return out;
  }

  out << *it;
  ++it;
  while(next(it, 1) != is.end()){
    out << "x" << *it;
    ++it;
  }
  out << "x" << *it;

  return out;
}

/*
typedef AtomSetImp1<list, MultiInterval, NI> AtomSetImp;
typedef AtomSetAbs<list, AtomSetImp, MultiInterval, NI> AtomSet;

ostream &operator<<(ostream &out, AtomSet &as){
  MultiInterval mi = as.aset_();

  out << "{" << mi << "}";

  return out;
}

typedef SetImp1<list, AtomSet> SetImp;
typedef SetAbs<list, SetImp, AtomSet> Set;

ostream &operator<<(ostream &out, Set &s){
  list<AtomSet> as = s.asets_();
  list<AtomSet>::iterator it = as.begin();

  if(as.size() == 0){
    out << "{}";
    return out;
  }

  if(as.size() == 1){
    out << "{" << *it << "}";
    return out;
  }

  out << *it;
  ++it;
  while(next(it, 1) != as.end()){
    out << "U" << *it;
    ++it;
  }
  out << "U" << *it;

  return out;
}
typedef LExprImp1<list> LExprImp;
typedef LExprAbs<LExprImp, float> LExpr;

ostream &operator<<(ostream &out, LExpr &le){
  float m = le.m_();
  float h = le.h_();

  out << m << " * x + " << h;

  return out;
}

typedef LMIntImp1<list, Interval, LExpr> LMIntImp;
typedef LMIntAbs<LMIntImp, Interval, LExpr> LMInt;

ostream &operator<<(ostream &out, LMInt &lm){
  Interval d = lm.dom_();
  LExpr e = lm.expr_();

  out << "(" << d << ", " << e << ")";  

  return out;
}

typedef LMMultiIntImp1<list, MultiInterval, Interval, LExpr, LMInt> LMMultIntImp;
typedef LMMultiIntAbs<list, LMMultiIntImp1<list, MultiInterval, 
                                           Interval, LExpr, LMInt>, MultiInterval, LExpr> LMMultiInt;


ostream &operator<<(ostream &out, LMMultiInt &lm){
  MultiInterval d = lm.dom_();
  list<LExpr> e = lm.expr_();
  list<LExpr>::iterator it = e.begin();

  out << "(" << d << ", (";
  while(next(it, 1) != e.end()){
    out << *it << ",";
    ++it;
  }
  out << "))";

  return out;
}

typedef LMMultiInt LMap;

// ------ Graph definition ------ //

struct SetVertex{
  SetVertex(){};
  SetVertex(string nm, Set v); 
  SetVertex(string nm, int i, Set v, int ind); 

  member_(string, name);
  member_(int, id);
  member_(Set, vs);
  member_(int, index);

  printable(SetVertex);
  comparable(SetVertex);
};

struct SetEdge{
  SetEdge(){};
  SetEdge(string nm, LMap e1, LMap e2);
  SetEdge(string nm, int i, LMap e1, LMap e2, int ind);

  member_(string, name);
  member_(int, id);
  member_(LMap, es1);
  member_(LMap, es2);
  member_(int, index);

  printable(SetEdge);
  comparable(SetEdge);
};

member_imp(SetVertex, string, name);
member_imp(SetVertex, int, id);
member_imp(SetVertex, Set, vs);
member_imp(SetVertex, int, index);

SetVertex::SetVertex(string nm, Set v) : name_(nm), vs_(v){
  set_id(-1);
  set_index(-1);
}

SetVertex::SetVertex(string nm, int i, Set v, int ind) : name_(nm), id_(i), vs_(v), index_(ind){}

ostream &operator<<(ostream &out, SetVertex &sv){
  out << sv.name();

  return out;
}

bool SetVertex::operator==(const SetVertex &other) const {
  return (id() == other.id()); 
}

member_imp(SetEdge, string, name);
member_imp(SetEdge, int, id);
member_imp(SetEdge, LMap, es1);
member_imp(SetEdge, LMap, es2);
member_imp(SetEdge, int, index);

SetEdge::SetEdge(string nm, LMap e1, LMap e2) : name_(nm), es1_(e1), es2_(e2){
  set_id(-1);
  set_index(-1);
}
SetEdge::SetEdge(string nm, int i, LMap e1, LMap e2, int ind) 
  : name_(nm), id_(i), es1_(e1), es2_(e2), index_(ind){}

ostream &operator<<(ostream &out, SetEdge &se){
  out << se.name();

  return out;
}

bool SetEdge::operator==(const SetEdge &other) const {
  return (id() == other.id()); 
}

typedef boost::adjacency_list<boost::listS, boost::listS, boost::undirectedS, SetVertex, SetEdge>
 SetBasedGraph;
typedef SetBasedGraph::vertex_descriptor SetVertexDesc;
typedef SetBasedGraph::edge_descriptor SetEdgeDesc;
*/
#endif

// TODO: operaciones de igualdad sobre resultado de diferencias
