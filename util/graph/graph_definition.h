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
#include <set>
#include <utility>

#include <boost/config.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/icl/discrete_interval.hpp>
#include <boost/optional/optional.hpp>
#include <boost/variant.hpp>
#include <boost/variant/get.hpp>

#include <ast/ast_types.h>
#include <ast/equation.h>
#include <util/table.h>

namespace ICL = boost::icl;
using boost::variant;
using namespace std;

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
  return (a * b) / gcd(a, b);
}

template <template<typename T, typename = allocator<T>> class CT>
struct IntervalImp1{
  bool empty;
  int lo;
  int step;
  int hi;

  IntervalImp1(){};
  IntervalImp1(bool isEmpty){ 
    empty = isEmpty;
  };
  IntervalImp1(int vlo, int vstep, int vhi, bool vempty){ 
    empty = vempty;
    lo = vlo;
    step = vstep;

    if(vlo <= vhi && !empty){
      int rem = std::fmod(vlo - vhi, step); 
      hi = vhi - rem; 
    }

    else{ 
      cerr << "Wrong values for subscript" << endl;
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
    if(x > lo || x < hi)
      return false;

    float aux = (x - lo) / step;
    if(aux == (int) aux)
      return true;

    return false;
  } 

  IntervalImp1 cap(const IntervalImp1 &inter2){
    int newStep = lcm(step, inter2.step), newLo = -1;
    int minEnd = std::min(hi, inter2.hi);

    if(!empty && !inter2.empty)
      for(int i = 0; i < newStep; i += step){
        int res1 = lo + i;
        float res2 = (res1 - inter2.lo) / inter2.step;

        if(res2 == (int) res2){
          newLo = res1;
          break;
        }
      }

    else
      return IntervalImp1(true);

    if(newLo < 0)
      return IntervalImp1(true);

    while(inter2.lo > newLo)
      newLo += newStep;

    int newEnd = newLo;
    while(newEnd < minEnd)
      newEnd += newStep;

    return IntervalImp1(newLo, newStep, newEnd, false);
  }

  CT<IntervalImp1> diff(const IntervalImp1 &i2){
    CT<IntervalImp1> res;
    typename CT<IntervalImp1>::iterator itRes = res.begin();

    IntervalImp1 i1 = *this;
    IntervalImp1 capRes = i1.cap(i2);

    if(capRes.empty()){
      CT<IntervalImp1> emptyRes;
      return emptyRes;
    }

    // "Before" intersection
    if(i1.lo() < capRes.lo){
      IntervalImp1 aux = IntervalImp1(i1.lo, 1, capRes.lo - 1, false);
      IntervalImp1 left = i1.cap(aux); 
      itRes = res.insert(itRes, left);
    }

    // "During" intersection
    if(capRes.step > capRes.hi - capRes.lo){
      CT<IntervalImp1> emptyRes;
      return emptyRes;
    }

    else{
      int nInters = capRes.step / i1.step - 1;
      for(int i = 1; i < nInters; i++){
        IntervalImp1 aux = IntervalImp1(capRes.lo + i, capRes.step, capRes.hi, false);
        itRes = res.insert(itRes, aux);
      }  
    }

    // "After" intersection
    if(i1.hi() > capRes.hi()){
      IntervalImp1 aux = IntervalImp1(capRes.hi + 1, 1, i1.hi, false);
      IntervalImp1 right = i1.cap(aux);
      itRes = res.insert(itRes, right);
    }

    return res; 
  }
};

template<template<typename T, typename = allocator<T>> class CT,
         typename IntervalImp, typename NumImp>
struct IntervalAbs{
  IntervalAbs(){}; 
  IntervalAbs(NumImp lo, NumImp step, NumImp hi, bool emp){
    i = IntervalImp(lo, step, hi, emp);
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

  IntervalAbs cap(IntervalAbs &i2){
    return IntervalAbs(i.cap(i2.i));
  }

  CT<IntervalAbs> diff(IntervalAbs &i2){
    CT<IntervalImp> resDiff = i.diff(i2.i);
    typename CT<IntervalImp>::iterator it = resDiff.begin();
    CT<IntervalAbs> res;
    typename CT<IntervalAbs>::iterator itres = res.begin();
    
    while(it != res.end()){
      itres = res.insert(itres, IntervalAbs(*it));
      ++itres;

      ++it;
    }

    return res;
  }
 
  private:
  IntervalImp i; 
};

// Represents the cartesian product of its elements
template <template<typename T, typename = allocator<T>> class CT,
          typename IntervalImp>
struct MultiInterImp1{
  CT<IntervalImp> inters;
  typedef typename CT<IntervalImp>::iterator intImpIt;

  MultiInterImp1(){
    CT<IntervalImp> emptyRes;
    inters = emptyRes;  
  }
  MultiInterImp1(CT<IntervalImp> is){
    inters = is;
  }

  CT<IntervalImp> inters_(){
    return inters;
  }

  void addInter(IntervalImp i){
    inters.insert(inters.end(), i);
  }

  bool empty(){
    return inters.empty();
  }

  MultiInterImp1 cap(MultiInterImp1 &mi2){
    CT<IntervalImp> res;
    intImpIt itres = res.begin();

    intImpIt it1 = inters.begin();
    intImpIt it2 = mi2.inters.begin();
    int minLength = std::min(inters.size(), mi2.inters.size());
    for(int i = 0; i < minLength; i++){
      IntervalImp capres = (*it1).cap(*it2);

      if(capres.empty()){
        CT<IntervalImp> emptyRes;
        return MultiInterImp1<CT, IntervalImp>(emptyRes);
      }
      
      itres = res.insert(itres, capres);
      ++itres;

      ++it1;
      ++it2;    
    }

    return MultiInterImp1<CT, IntervalImp>(res);
  }

  CT<MultiInterImp1<CT, IntervalImp>> diff(MultiInterImp1 &mi2){
    MultiInterImp1<CT, IntervalImp> capres = cap(mi2);

    if(capres.empty()){
      CT<MultiInterImp1<CT, IntervalImp>> emptyRes;
      return emptyRes;
    }

    CT<MultiInterImp1<CT, IntervalImp>> res1;
    typename CT<MultiInterImp1<CT, IntervalImp>>::iterator itres1 = res1.begin();

    intImpIt it1 = inters.begin();
    intImpIt itcap = capres.inters.begin();

    for(int i = 0; i < inters.size(); i++){
      CT<IntervalImp> res2;
      intImpIt itres2 = res2.begin();

      intImpIt it11 = inters.begin();

      for(int j = 0; j < inters.size(); j++){
        int i1 = std::distance(inters.begin(), it1), i11 = distance(inters.begin(), it11);

        if(i1 < i11)
          it11 = res2.insert(itres2, *itcap);

        else if(i1 == i11){
          CT<IntervalImp> diffres = diff(*it1, *itcap);
          typename CT<IntervalImp>::iterator taux = diffres.begin();
          while(taux != diffres.end())
            it11 = res2.insert(itres2, *taux);
        }

        else if(i1 > i11)
          it11 = res2.insert(itres2, *it11);

        ++it11;
      }

      itres1 = res1.insert(itres1, res2);
      ++itres1;

      ++itcap;
    }

    return res1;
  }
};

template <template <typename T, typename = allocator<T>> class CT,
          typename MultiInterImp, typename IntervalImp>
struct MultiInterAbs{
  MultiInterAbs(){
    CT<IntervalImp> aux;
    multiInterImp = aux;
  };
  MultiInterAbs(CT<IntervalImp> ints){
    multiInterImp = MultiInterImp(ints);
  };

  bool empty(){
    return multiInterImp.empty();
  };

  CT<IntervalImp> inters_(){
    return multiInterImp.inters_();
  }

  void addInter(IntervalImp i){
    multiInterImp.addInter(i);
  }

  MultiInterAbs cap(MultiInterAbs &mi2){
    return MultiInterAbs(multiInterImp.cap(mi2.multiInterImp));
  };

  CT<MultiInterAbs<CT, MultiInterImp, IntervalImp>> diff(MultiInterAbs &mi2){
    CT<MultiInterImp> diffRes = multiInterImp.diff(mi2.multiInterImp);
    typename CT<MultiInterImp>::iterator it = diffRes.begin();
    CT<MultiInterImp> res;
    typename CT<MultiInterImp>::iterator itres = res.begin();

    while(it != diffRes.end()){
      itres = res.insert(itres, MultiInterAbs(*it));
      ++itres;

      ++it;
    }

    return res;
  };

  private:
  MultiInterImp multiInterImp;
};

template <template<typename T1, typename = allocator<T1>> class CT,
          typename MultiInterImp>
struct AtomSetImp1{
  MultiInterImp aset;

  AtomSetImp1(){
    MultiInterImp emptyRes();
    aset = emptyRes;
  };
  AtomSetImp1(MultiInterImp as){
    aset = as;
  }

  MultiInterImp aset_(){
    return aset;
  }

  bool empty(){
    return aset.empty();  
  }

  AtomSetImp1 cap(const AtomSetImp1 &aset2){
    return AtomSetImp1(aset.cap(aset2.aset));
  }

  CT<AtomSetImp1<CT, MultiInterImp>> diff(const AtomSetImp1 &aset2){
    CT<AtomSetImp1<CT, MultiInterImp>> res;
    typename CT<AtomSetImp1<CT, MultiInterImp>>::iterator itres = res.begin();

    CT<MultiInterImp> atomicDiff = aset.diff(aset2);

    if(atomicDiff.empty()){
      CT<AtomSetImp1<CT, MultiInterImp>> emptyRes;
      return emptyRes;
    } 

    else{
      typename CT<MultiInterImp>::iterator it = atomicDiff.begin();
      for(int i = 0; i < atomicDiff.size(); i++)
        itres = res.insert(AtomSetImp1(*it)); 
    }

    return res;
  }
};

template <template<typename T, typename = allocator<T>> class CT,
          typename ASetImp, typename MultiInterImp>
struct AtomSetAbs{
  AtomSetAbs(){};
  AtomSetAbs(ASetImp as){
    aset = as;
  }

  MultiInterImp aset_(){
    return aset.aset_();
  }

  bool empty(){
    return aset.empty();  
  }

  AtomSetAbs cap(const AtomSetAbs &aset2){
    return AtomSetAbs(aset.cap(aset2.aset));
  }

  CT<AtomSetAbs<CT, ASetImp, MultiInterImp>> diff(const AtomSetAbs &aset2){
    CT<ASetImp> diffRes = aset.diff(aset2.aset);
    typename CT<ASetImp>::iterator it = diffRes.begin();
    CT<ASetImp> res;
    typename CT<ASetImp>::iterator itRes = res.begin();

    while(it != diffRes.end()){
      itRes = res.insert(itRes, AtomSetAbs(*it));

      ++it;
    }

    return res;
  }

  private:
  ASetImp aset;
};

template <template<typename T1, typename = std::allocator<T1>> class CT,
          typename ASetImp>
struct SetImp1{
  typedef CT<ASetImp> setType;
  setType asets;
 
  SetImp1(){};
  SetImp1(setType ss){
    asets = ss;
  };

  setType asets_(){
    return asets;
  }

  bool empty(){
    typename setType::iterator it = asets.begin();

    if(asets.empty())
      return true;

    while(it != asets.end())
      if(!(*it).empty())
        return false;
  
    return true;
  }

  SetImp1 addAtomSet(const ASetImp &aset2){
    typename setType::iterator itasets = asets.begin();

    asets.insert(itasets, aset2);
    return SetImp1(asets);
  }

  SetImp1 addAtomSets(const setType &sets2){
    setType res;
    typename setType::iterator it2 = sets2.begin();

    while(it2 != sets2.end()){
      res = res.addAtomSet(*it2);
      ++it2;
    }

    return res;
  }

  SetImp1 cap(const SetImp1 &set2){
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
        itres = res.insert(itres, *it1.cap(*it2));
        ++it2;
      }
     
      ++it1;
    }

    return res;
  }

  SetImp1 diff(const SetImp1 &set2){
    CT<ASetImp> emptyCT;
    SetImp1 res(emptyCT);
    setType setCap = cap(set2); 

    if(!setCap.empty()){
      typename setType::iterator it1 = asets.begin();

      while(it1 != asets.end()){
        setType asets2;
        it1 = insert(it1, *it1);
        typename setType::iterator itasets2 = asets2.begin();
 
        typename setType::iterator it2 = setCap.begin();
        while(it2 != setCap.end()){
          setType newSet;
          typename setType::iterator itnew = newSet.begin();

          while(itasets2 != asets2.end()){
            itnew = newSet.insert(itnew, *itasets2.diff(*it2));

            ++itasets2;
          }

          asets2 = newSet;

          ++it2;
        }

        res.addAtomSets(asets2);

        ++it1;
      }
    }

    return res; 
  }

  SetImp1 cup(const SetImp1 &set2){
    SetImp1 res = *this;
    SetImp1 resDiff = diff(set2);

    if(!resDiff.empty())
      res.addAtomSets(resDiff);

    return res;
  }
};

template <template<typename T1, typename = std::allocator<T1>> class CT,
          typename SetImp, typename ASetImp>
struct SetAbs{
  SetAbs(){};
  SetAbs(SetImp ss){
    set = ss;
  };

  bool empty(){
    return set.empty(); 
  }

  CT<ASetImp> asets_(){
    return set.asets_();
  }

  SetAbs addAtomSet(const ASetImp &aset2){
    return SetAbs(set.addAtomSet(aset2)); 
  }

  SetAbs addAtomSets(const CT<ASetImp> &sets2){
    return SetAbs(set.addAtomSets(sets2));
  }

  SetAbs cap(const SetImp &set2){
    return SetAbs(set.cap(set2));
  }

  SetAbs diff(const SetImp &set2){
    return SetAbs(set.diff(set2)); 
  }

  SetAbs cup(const SetImp &set2){
    return SetAbs(set.cup(set2));
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

/*
  SetImp image(SetImp &set){
    CT<ASetImp> adom = dom.asets_(); 
    typename CT<ASetImp>::iterator itadom = adom.begin();
    typename CT<CT<LExprImp>>::iterator itexpr = expr.begin();

    while(itadom != adom.end()){
      MultiInterImp mi = (*itadom).aset_();
    }
  }
*/
};

/*-----------------------------------------------------------------------------------------------*/
/*-----------------------------------------------------------------------------------------------*/
// Concrete classes
/*-----------------------------------------------------------------------------------------------*/
/*-----------------------------------------------------------------------------------------------*/

typedef IntervalImp1<list> IntervalImp;
typedef IntervalAbs<list, IntervalImp, int> Interval;

ostream &operator<<(ostream &out, Interval &i){
  out << "[" << i.lo_() << ":" << i.step_() << ":" << i.hi_() << "]";
  return out;
}

typedef MultiInterImp1<list, Interval> MultiInterImp;
typedef MultiInterAbs<list, MultiInterImp, Interval> MultiInterval;

ostream &operator<<(ostream &out, MultiInterval &mi){
  list<Interval> is = mi.inters_();
  list<Interval>::iterator it = is.begin();

  while(next(it, 1) != is.end()){
    out << *it << "x";
    ++it;
  }

  return out;
}

typedef AtomSetImp1<list, MultiInterval> AtomSetImp;
typedef AtomSetAbs<list, AtomSetImp, MultiInterval> AtomSet;

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

  while(next(it, 1) != as.end()){
    out << *it << "U";
    ++it;
  }

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

#endif
