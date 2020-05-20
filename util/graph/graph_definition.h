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
#include <util/debug.h>
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
        WARNING("Wrong values for subscript (check low <= hi)");
        empty = true;
      }
    }

    else{ 
      WARNING("Subscripts should be positive");
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
  typedef typename CT1<IntervalImp>::iterator IntImpIt;

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
    IntImpIt it = inters.begin();

    if(inters.empty())
      return true;

    while(it != inters.end()){
      if(!(*it).empty_())
        return false;

      ++it;
    }
  
    return true;
  }

  bool isIn(CT1<NumImp> elem){
    typename CT1<NumImp>::iterator it1 = elem.begin();
    IntImpIt it2 = inters.begin();

    if(elem.size() != inters.size())
      return false;

    while(it1 != elem.end()){
      if(!(*it2).isIn(*it1))
        return false;      

      ++it1;
      ++it2;
    }

    return true;
  }

  MultiInterImp1 cap(MultiInterImp1 &mi2){
    CT1<IntervalImp> res;
    IntImpIt itres = res.begin();

    IntImpIt it1 = inters.begin();
    IntImpIt it2 = mi2.inters.begin();
    int sz = inters.size();
    if(inters.size() == mi2.inters.size()){
      for(int i = 0; i < sz; ++i){
        IntervalImp capres = (*it1).cap(*it2);
     
        if(!(*it1).empty_() && !(*it2).empty_()){
          if(capres.empty_()){
            CT1<IntervalImp> aux;
            return MultiInterImp1(aux);
          }

          itres = res.insert(itres, capres);
          ++itres;
        }

        else{
          CT1<IntervalImp> aux;
          return MultiInterImp1(aux);
        }

        ++it1;
        ++it2;    
      }
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

    IntImpIt it1 = inters.begin();
    IntImpIt itcap = capres.inters.begin();
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
        IntervalImp auxaux = *itaux;

        if(!auxaux.empty_()){
          CT1<IntervalImp> resi;
          IntImpIt itresi = resi.begin();
  
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

          IntImpIt auxit1 = it1;
          while(auxit1 != inters.end()){
            itresi = resi.insert(itresi, *auxit1);
            ++itresi;
        
            ++auxit1;
          }

          itresmi = resmi.insert(itresmi, MultiInterImp1(resi));
          ++itresmi;
        }

        ++itaux;
      }

      ++it1;
      ++itdiffs;
    }

    return resmi;
  }

  MultiInterImp1 crossProd(MultiInterImp1 mi2){
    IntImpIt it2 = mi2.inters.begin();

    while(it2 != mi2.inters.end()){
      inters.insert(inters.end(), *it2);
 
      ++it2;
    }

    return *this;
  }

  CT1<NumImp> minElem(){
    IntImpIt it = inters.begin();

    CT1<NumImp> res;
    typename CT1<NumImp>::iterator itres = res.begin();

    while(it != inters.end()){
      if((*it).empty_()){
        CT1<NumImp> aux;
        return aux;
      }

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
  }
  MultiInterAbs(MultiInterImp mi){
    multiInterImp = mi;
  }
  MultiInterAbs(CT1<IntervalImp> ints){
    multiInterImp = MultiInterImp(ints);
  }

  bool empty(){
    return multiInterImp.empty();
  }

  CT1<IntervalImp> inters_(){
    return multiInterImp.inters_();
  }

  bool isIn(CT1<NumImp> elem){
    return multiInterImp.isIn(elem);
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
  }

  MultiInterAbs crossProd(MultiInterAbs &mi2){
    return MultiInterAbs(multiInterImp.crossProd(mi2.multiInterImp));
  }

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

template<template<typename T, typename = allocator<T>> class CT1, 
         template<typename Value, typename Hash = boost::hash<Value>, 
                  typename Pred = std::equal_to<Value>, 
                  typename Alloc = std::allocator<Value>> class CT2,
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

  bool isIn(CT1<NumImp> elem){
    return aset.isIn(elem);
  }

  AtomSetImp1 cap(AtomSetImp1 &aset2){
    AtomSetImp1 aux(aset.cap(aset2.aset));
    return aux;
  }

  CT2<AtomSetImp1> diff(AtomSetImp1 &aset2){
    CT2<AtomSetImp1> res;
    typename CT2<AtomSetImp1>::iterator itres = res.begin();

    CT2<MultiInterImp> atomicDiff = aset.diff(aset2.aset);

    if(atomicDiff.empty()){
      CT2<AtomSetImp1> emptyRes;
      return emptyRes;
    } 

    else{
      typename CT2<MultiInterImp>::iterator it = atomicDiff.begin();
      while(it != atomicDiff.end()){
        itres = res.insert(itres, AtomSetImp1(*it)); 
        ++itres;

        ++it;
      }
    }

    return res;
  }

  AtomSetImp1 crossProd(AtomSetImp1 &aset2){
    return AtomSetImp1(aset.crossProd(aset2.aset));
  }

  CT1<NumImp> minElem(){
    return aset.minElem();
  }

  bool operator==(const AtomSetImp1 &other) const{
    return aset == other.aset;
  }

  bool operator!=(const AtomSetImp1 &other) const{
    return aset != other.aset;
  }

  size_t hash(){
    return aset.hash();
  }
};

template<template<typename T, typename = allocator<T>> class CT1, 
         template<typename Value, typename Hash = boost::hash<Value>, 
                  typename Pred = std::equal_to<Value>, 
                  typename Alloc = std::allocator<Value>> class CT2,
         typename MultiInterImp, typename NumImp>
size_t hash_value(AtomSetImp1<CT1, CT2, MultiInterImp, NumImp> as){
  return as.hash();
}

template<template<typename T, typename = allocator<T>> class CT1,
         template<typename Value, typename Hash = boost::hash<Value>, 
                  typename Pred = std::equal_to<Value>, 
                  typename Alloc = std::allocator<Value>> class CT2,
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

  bool isIn(CT1<NumImp> elem){
    return as.isIn(elem);
  }

  AtomSetAbs cap(AtomSetAbs &aset2){
    AtomSetAbs aux(as.cap(aset2.as));
    return aux;
  }

  CT2<AtomSetAbs> diff(AtomSetAbs &aset2){
    CT2<ASetImp> diffRes = as.diff(aset2.as);
    typename CT2<ASetImp>::iterator it = diffRes.begin();
    CT2<AtomSetAbs> res;
    typename CT2<AtomSetAbs>::iterator itRes = res.begin();

    while(it != diffRes.end()){
      itRes = res.insert(itRes, AtomSetAbs(*it));
      ++itRes;

      ++it;
    }

    return res;
  }

  AtomSetAbs crossProd(AtomSetAbs &aset2){
    return AtomSetAbs(as.crossProd(aset2.as));
  }
 
  CT1<NumImp> minElem(){
    return as.minElem();
  }

  bool operator==(const AtomSetAbs &other) const{
    return as == other.as;
  }

  bool operator!=(const AtomSetAbs &other) const{
    return as != other.as;
  }

  size_t hash(){
    return as.hash();
  }

  private:
  ASetImp as;
};

template<template<typename T, typename = allocator<T>> class CT1,
         template<typename Value, typename Hash = boost::hash<Value>, 
                  typename Pred = std::equal_to<Value>, 
                  typename Alloc = std::allocator<Value>> class CT2,
         typename ASetImp, typename MultiInterImp, typename NumImp>
size_t hash_value(AtomSetAbs<CT1, CT2, ASetImp, MultiInterImp, NumImp> as){
  return as.hash();
}

// Sets --------------------------------------------------------------------------------------------

template<template<typename T, typename = allocator<T>> class CT1,
         template<typename Value, typename Hash = boost::hash<Value>, 
                  typename Pred = std::equal_to<Value>, 
                  typename Alloc = std::allocator<Value>> class CT2,
         typename ASetImp, typename NumImp>
struct SetImp1{
  typedef CT2<ASetImp> SetType;
  typedef typename SetType::iterator SetIt;

  SetType asets;
 
  SetImp1(){
    SetType aux;
    asets = aux;
  }
  SetImp1(SetType ss){
    asets = ss;
  }

  SetType asets_(){
    return asets;
  }

  bool empty(){
    ASetImp aux;
    SetIt it = asets.begin();

    if(asets.empty())
      return true;

    while(it != asets.end()){
      aux = *it;

      if(!aux.empty())
        return false;

      ++it;
    }
  
    return true;
  }

  bool isIn(CT1<NumImp> elem){
    ASetImp aux;
    SetIt it = asets.begin();

    while(it != asets.end()){
      aux = *it;
      if(aux.isIn(elem))
        return true;
    }

    return false;
  }

  void addAtomSet(ASetImp &aset2){
    if(!aset2.empty())
      asets.insert(aset2);
  }

  void addAtomSets(SetType &sets2){
    ASetImp aux;

    SetIt it = sets2.begin();

    while(it != sets2.end()){
      aux = *it;
      addAtomSet(aux);

      ++it;
    }
  }

  SetImp1 cap(SetImp1 &set2){
    ASetImp aux1, aux2;

    if(empty() || set2.empty()){
      SetType emptyRes;
      return emptyRes; 
    }
    
    SetType res;
    SetIt it1 = asets.begin();

    while(it1 != asets.end()){
      aux1 = *it1;

      SetIt it2 = set2.asets.begin();

      while(it2 != set2.asets.end()){
        aux2 = *it2;
        res.insert(aux1.cap(aux2));

        ++it2;
      }
     
      ++it1;
    }

    return res;
  }

  SetImp1 diff(SetImp1 &set2){
    ASetImp aux1, aux2, aux3;

    SetImp1 res;
    SetType capres = cap(set2).asets; 

    if(!capres.empty()){
      SetIt it1 = asets.begin();

      while(it1 != asets.end()){
        aux1 = *it1;
        SetType aux;
        aux.insert(aux1);
 
        SetIt it2 = capres.begin();
        while(it2 != capres.end()){
          aux2 = *it2;

          SetImp1 newSets;

          SetIt itaux = aux.begin();
          while(itaux != aux.end()){
            aux3 = *itaux;

            SetType diffres = (aux3).diff(aux2);
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

  SetImp1 crossProd(SetImp1 &set2){
    ASetImp aux1;
    ASetImp aux2;
    SetIt it1 = asets.begin();
    SetIt it2 = set2.asets.begin();

    SetType res;

    while(it1 != asets.end()){
      aux1 = *it1;

      while(it2 != set2.asets.end()){
        aux2 = *it2;
 
        ASetImp auxres = aux1.crossProd(aux2);
        res.addAtomSet(auxres);

        ++it2;
      }
 
      ++it1;
    }

    return res;
  }

  CT1<NumImp> minElem(){
    SetIt it = asets.begin();
    ASetImp aux1;

    CT2<CT1<NumImp>> mins;
    typename CT2<CT1<NumImp>>::iterator itmins = mins.begin();

    while(it != asets.end()){
      aux1 = *it;

      itmins = mins.insert(itmins, aux1.minElem());
      ++itmins;

      ++it;
    }

    bool hasValue = false;
    CT1<NumImp> res;
    CT1<NumImp> aux2;
  
    itmins = mins.begin();
    while(itmins != mins.end()){
      aux2 = *itmins;

      if(!hasValue && !aux2.empty()){
        res = aux2;
        hasValue = true;
      }

      if(hasValue && !aux2.empty()){
        typename CT1<NumImp>::iterator it1 = res.begin();
        typename CT1<NumImp>::iterator it2 = aux2.begin();

        while(it1 != res.end()){
          if(*it2 < *it1)
            res = aux2;

          else if(*it1 < *it2)
            break;

          ++it1;
          ++it2;
        }
      }      

      ++itmins;
    }

    return res;
  }

  bool operator==(const SetImp1 &other) const{
    return asets == other.asets;
  }

  bool operator!=(const SetImp1 &other) const{
    return asets != other.asets;
  }

  size_t hash(){
    return asets.size();
  }
};

template<template<typename T, typename = allocator<T>> class CT1, 
         template<typename Value, typename Hash = boost::hash<Value>, 
                  typename Pred = std::equal_to<Value>, 
                  typename Alloc = std::allocator<Value>> class CT2,
         typename ASetImp, typename NumImp>
size_t hash_value(SetImp1<CT1, CT2, ASetImp, NumImp> s){
  return s.hash();
}

template<template<typename T, typename = allocator<T>> class CT1,
         template<typename Value, typename Hash = boost::hash<Value>, 
                  typename Pred = std::equal_to<Value>, 
                  typename Alloc = std::allocator<Value>> class CT2,
          typename SetImp, typename ASetImp, typename NumImp>
struct SetAbs{
  SetAbs(){
    SetImp aux;
    set = aux;
  }
  SetAbs(SetImp ss){
    set = ss;
  }
  SetAbs(CT2<ASetImp> ass){
    SetImp aux(ass);
    set = aux;
  }

  bool empty(){
    return set.empty(); 
  }
 
  bool isIn(CT1<NumImp> elem){
    return set.isIn(elem);
  }

  CT2<ASetImp> asets_(){
    return set.asets_();
  }

  void addAtomSet(ASetImp &aset2){
    set.addAtomSet(aset2); 
  }

  void addAtomSets(CT2<ASetImp> &sets2){
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

  SetAbs crossProd(SetAbs &set2){
    return SetAbs(set.crossProd(set2.set));
  }

  CT1<NumImp> minElem(){
    return set.minElem();
  }

  bool operator==(const SetAbs &other) const{
    return set == other.set;
  }

  bool operator!=(const SetAbs &other) const{
    return set != other.set;
  }

  size_t hash(){
    return set.hash();
  }

  private:
  SetImp set;
};

template<template<typename T, typename = allocator<T>> class CT1,
         template<typename Value, typename Hash = boost::hash<Value>, 
                  typename Pred = std::equal_to<Value>, 
                  typename Alloc = std::allocator<Value>> class CT2,
          typename SetImp, typename ASetImp, typename NumImp>
size_t hash_value(SetAbs<CT1, CT2, SetImp, ASetImp, NumImp> s){
  return s.hash();
}

// LinearMaps ---------------------------------------------------------------------------------------

template <template<typename T, typename = std::allocator<T>> class CT,
          typename NumImp>
struct LMapImp1{
  typedef CT<NumImp> CTNum;
  typedef typename CTNum::iterator CTNumIt;

  CTNum gain;
  CTNum offset;

  LMapImp1(){
    CTNum aux;
    gain = aux;
    offset = aux;
  }
  LMapImp1(CTNum g, CTNum o){
    if(g.size() == o.size()){
      gain = g;
      offset = o;
    }

    else{
      WARNING("Offset and gain should be of the same size");

      CTNum aux;
      gain = aux;
      offset = aux;
    }  
  }

  CTNum gain_(){
    return gain;
  }

  CTNum off_(){
    return offset;
  }

  bool empty(){
    if(gain.empty() && offset.empty())
      return true;

    return false;
  }

  LMapImp1 compose(LMapImp1 &lm2){
    CTNum resg;
    CTNumIt itresg = resg.begin();
    CTNum reso;
    CTNumIt itreso = reso.begin();

    CTNumIt itg1 = gain.begin();
    CTNumIt ito1 = offset.begin();
    CTNumIt itg2 = lm2.gain.begin();
    CTNumIt ito2 = lm2.offset.begin();

    if(gain.size() == lm2.gain.size()){
      while(itg1 != gain.end()){
        itresg = resg.insert(itresg, (*itg1) * (*itg2));
        ++itresg;
        itreso = reso.insert(itreso, (*ito2) * (*itg1) + (*ito1));
        ++itreso;

        ++itg1;
        ++ito1;
        ++itg2;
        ++ito2;
      } 
    }

    else{
      WARNING("Linear maps should be of the same size");
      LMapImp1 aux;
      return aux;
    }

    return LMapImp1(resg, reso);
  }

  LMapImp1 invLMap(){
    CTNum resg;
    CTNumIt itresg = resg.begin();
    CTNum reso;
    CTNumIt itreso = reso.begin();

    CTNumIt itg1 = gain.begin();
    CTNumIt ito1 = offset.begin();

    while(itg1 != gain.end()){
      if((*itg1) != 0){
        itresg = resg.insert(itresg, 1 / (*itg1));
        ++itresg;

        itreso = reso.insert(itreso, -(*ito1) / (*itg1));
        ++itreso;
      }

      else{
        itresg = resg.insert(itresg, Inf);
        itresg++;
      
        itreso = reso.insert(itreso, -Inf);
        ++itreso;
      }

      ++itg1;
      ++ito1;
    }

    return LMapImp1(resg, reso);
  }

  bool operator==(const LMapImp1 &other) const{
    return gain == other.gain && offset == other.offset;
  }
};

template<template<typename T, typename = std::allocator<T>> class CT,
         typename LMapImp, typename NumImp>
struct LMapAbs{
  typedef CT<NumImp> CTNum;
  typedef typename CTNum::iterator CTNumIt;

  LMapAbs(){
    LMapImp aux;
    lm = aux;
  }
  LMapAbs(LMapImp lmap){
    lm = lmap;
  }
  LMapAbs(CTNum g, CTNum o){
    lm = LMapImp(g, o);
  }
  
  CTNum gain_(){
    return lm.gain_();
  }

  CTNum off_(){
    return lm.off_();
  }

  bool empty(){
    return lm.empty();
  }

  LMapAbs compose(LMapAbs &lm2){
    return LMapAbs(lm.compose(lm2.lm));
  }

  LMapAbs invLMap(){
    return LMapAbs(lm.invLMap());
  } 

  bool operator==(const LMapAbs &other) const{
    return lm == other.lm;
  }

  private:
  LMapImp lm;
};

/*
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

typedef int NI1;
typedef float NI2;

template<typename T, class = allocator<T>>
using OrdCT = list<T>;

template<typename Value, typename Hash = boost::hash<Value>, 
         typename Pred = std::equal_to<Value>, 
         typename Alloc = std::allocator<Value>>
using UnordCT = boost::unordered_set<Value>;


typedef IntervalImp1<UnordCT> IntervalImp;
typedef IntervalAbs<UnordCT, IntervalImp, NI1> Interval;

ostream &operator<<(ostream &out, Interval &i){
  out << "[" << i.lo_() << ":" << i.step_() << ":" << i.hi_() << "]";
  return out;
}

typedef MultiInterImp1<OrdCT, UnordCT, Interval, NI1> MultiInterImp;
typedef MultiInterAbs<OrdCT, UnordCT, MultiInterImp, Interval, NI1> MultiInterval;

ostream &operator<<(ostream &out, MultiInterval &mi){
  OrdCT<Interval> is = mi.inters_();
  OrdCT<Interval>::iterator it = is.begin();

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

typedef AtomSetImp1<OrdCT, UnordCT, MultiInterval, NI1> AtomSetImp;
typedef AtomSetAbs<OrdCT, UnordCT, AtomSetImp, MultiInterval, NI1> AtomSet;

ostream &operator<<(ostream &out, AtomSet &as){
  MultiInterval mi = as.aset_();

  out << "{" << mi << "}";

  return out;
}

typedef SetImp1<OrdCT, UnordCT, AtomSet, NI1> SetImp;
typedef SetAbs<OrdCT, UnordCT, SetImp, AtomSet, NI1> Set;

ostream &operator<<(ostream &out, Set &s){
  UnordCT<AtomSet> as = s.asets_();
  UnordCT<AtomSet>::iterator it = as.begin();
  AtomSet aux;

  if(as.size() == 0){
    out << "{}";
    return out;
  }

  if(as.size() == 1){
    aux = *it;
    out << "{" << aux << "}";
    return out;
  }

  aux = *it;
  out << aux;
  ++it;
  while(next(it, 1) != as.end()){
    aux = *it;
    out << "U" << aux;
    ++it;
  }
  aux = *it;
  out << "U" << aux;

  return out;
}

typedef LMapImp1<OrdCT, NI2> LMapImp;
typedef LMapAbs<OrdCT, LMapImp, NI2> LMap;

ostream &operator<<(ostream &out, LMap &lm){
  OrdCT<NI2> g = lm.gain_();
  OrdCT<NI2>::iterator itg = g.begin();
  OrdCT<NI2> o = lm.off_();  
  OrdCT<NI2>::iterator ito = o.begin();
  
  while(itg != g.end()){
    out << *itg << " * x + " << *ito << "\n";

    ++itg;
    ++ito;
  }

  return out;
}

/*
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

