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

// TODO:
// Printable instances for all classes!

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

// Interval implementation is directly given because the methods (cap, diff)
// rely too much on its numbered-based implementation. 
// It is needed here so MultiInterAbs knows the struct.
template <template<typename T, typename = std::allocator<T>> class ContainerType>
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


  bool isEmpty() const{
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
    int newStep = lcm(step(), inter2.step()), newLo = -1;
    int minEnd = std::min(hi(), inter2.hi());

    if(!empty() && !inter2.empty())
      for(int i = 0; i < newStep; i += step()){
        int res1 = lo() + i;
        float res2 = (res1 - inter2.lo()) / inter2.step();

        if(res2 == (int) res2){
          newLo = res1;
          break;
        }
      }

    else
      return IntervalImp1(true);

    if(newLo < 0)
      return IntervalImp1(true);

    while(inter2.lo() > newLo)
      newLo += newStep;

    int newEnd = newLo;
    while(newEnd < minEnd)
      newEnd += newStep;

    return IntervalImp1(newLo, newStep, newEnd, false);
  }

  ContainerType<IntervalImp1> diff(const IntervalImp1 &i2){
    ContainerType<IntervalImp1> res;
    typename ContainerType<IntervalImp1>::iterator itRes = res.begin();

    IntervalImp1 i1 = *this;
    IntervalImp1 capres = i1.cap(i2);

    if(capres.empty()){
      ContainerType<IntervalImp1> emptyRes;
      return emptyRes;
    }

    // "Before" intersection
    if(i1.lo() < capres.lo()){
      IntervalImp1 aux = IntervalImp1(i1.lo(), 1, capres.lo() - 1, false);
      IntervalImp1 left = i1.cap(aux); 
      itRes = res.insert(itRes, left);
    }

    // "During" intersection
    if(capres.step() > capres.hi() - capres.lo()){
      ContainerType<IntervalImp1> emptyRes;
      return emptyRes;
    }

    else{
      int nInters = capres.step() / i1.step() - 1;
      for(int i = 1; i < nInters; i++){
        IntervalImp1 aux = IntervalImp1(capres.lo() + i, capres.step(), capres.hi(), false);
        itRes = res.insert(itRes, aux);
      }  
    }

    // "After" intersection
    if(i1.hi() > capres.hi()){
      IntervalImp1 aux = IntervalImp1(capres.hi() + 1, 1, i1.hi(), false);
      IntervalImp1 right = i1.cap(aux);
      itRes = res.insert(itRes, right);
    }

    return res; 
  }

  /// @brief Return a new interval with only the min element
  IntervalImp1 min(){
    return IntervalImp1(lo(), 1, lo(), false);
  }
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

  bool empty(){
    return inters.empty();
  }

  CT<IntervalImp> ints(){
    return inters;
  }

  MultiInterImp1 cap(MultiInterImp1 &mi2){
    CT<IntervalImp> res;
    intImpIt itres = res.begin();

    intImpIt it1 = inters.begin();
    intImpIt it2 = mi2.inters().begin();
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

      for(int j = 0; j < inters().size(); j++){
        int i1 = std::distance(inters().begin(), it1), i11 = std::distance(inters().begin(), it11);

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
  MultiInterAbs(){};
  MultiInterAbs(MultiInterImp d){
    multiInterImp = d;
  };

  bool empty(){
    return multiInterImp.empty();
  };

  CT<IntervalImp> ints(){
    return multiInterImp.inters();
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

  bool empty(){
    return aset.empty();  
  }

  MultiInterImp as(){
    return aset;
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

  bool empty(){
    return aset.empty();  
  }

  MultiInterImp as(){
    return aset.as();
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
  setType sets;
 
  SetImp1(){};
  SetImp1(setType ss){
    sets = ss;
  };

  bool empty(){
    typename setType::iterator it = sets.begin();

    if(sets.empty())
      return true;

    while(it != sets.end())
      if(!(*it).empty())
        return false;
  
    return true;
  }

  setType asets(){
    return sets;
  }

  SetImp1 addAtomSet(const ASetImp &aset2){
    typename setType::iterator itsets = sets.begin();

    sets.insert(itsets, aset2);
    return SetImp1(sets);
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
    typename setType::iterator it1 = sets.begin();

    while(it1 != sets.end()){
      typename setType::iterator it2 = set2.sets.begin();

      while(it2 != set2.sets.end()){
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
      typename setType::iterator it1 = sets.begin();

      while(it1 != sets.end()){
        setType asets;
        it1 = insert(it1, *it1);
        typename setType::iterator itasets = asets.begin();
 
        typename setType::iterator it2 = setCap.begin();
        while(it2 != setCap.end()){
          setType newSet;
          typename setType::iterator itnew = newSet.begin();

          while(itasets != asets.end()){
            itnew = newSet.insert(itnew, *itasets.diff(*it2));

            ++itasets;
          }

          asets = newSet;

          ++it2;
        }

        res.addAtomSets(asets);

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

  CT<ASetImp> asets(){
    return set.asets();
  }

  SetImp addAtomSet(const ASetImp &aset2){
    return SetAbs(set.addAtomSet(aset2)); 
  }

  SetImp addAtomSets(const CT<ASetImp> &sets2){
    return SetAbs(set.addAtomSets(sets2));
  }

  SetImp cap(const SetImp &set2){
    return SetAbs(set.cap(set2));
  }

  SetImp diff(const SetImp &set2){
    return SetAbs(set.diff(set2)); 
  }

  SetImp cup(const SetImp &set2){
    return SetAbs(set.cup(set2));
  }

  private:
  SetImp set;
};

template <template<typename T, typename = std::allocator<T>> class CT,
          typename>
struct LExprImp1{
  float m;
  float h;

  LExprImp1(){};
  LExprImp1(float mm, float hh){
    m = mm;
    h = hh;
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

    cerr << "The expression is no injective";
    return LExprImp1();
  }
};

// TODO: abstract class!

template <template<typename T, typename = std::allocator<T>> class CT,
          typename LExprImp, typename IntervalImp>
struct LMIntImp1{
  IntervalImp dom;
  LExprImp expr;

  LMIntImp1(){};
  LMIntImp1(IntervalImp d, LExprImp e){
    dom = d;
    expr = e;
  }

  IntervalImp image(){
    int newLo = expr.m() * dom.lo() + expr.h();
    int newStep = expr.m() * dom.step();
    int newHi = expr.m() * dom.hi() + expr.h();

    IntervalImp res(newLo, newStep, newHi, false);
    return res;    
  }

  LMIntImp1 compose(LMIntImp1 &lm2){
    IntervalImp newDom;
    LExprImp newExpr = expr.compose(lm2.expr);

    if(lm2.m() == 0){
      if(dom.isIn(lm2.h()))
        newDom = lm2.dom;
    }

    else{
      LMIntImp1 invLMInt(lm2.image(), lm2.invExpr()); 
      newDom = invLMInt.dom();
    }

    return LMIntImp1(newDom, newExpr);
  }

  LMIntImp1 miniInv(){
    IntervalImp newDom;
    LExprImp newExpr;

    if(expr.m == 0){
      newDom = image();
      newExpr = LExprImp(0, dom.lo()); 
    }

    else{
      newDom = image();
      newExpr = expr.invExpr();
    }

    return LMIntImp1(newDom, newExpr);
  }
};

// TODO: abstract class!

/*
  MultiInterImp applyExprMI(){
    int mss = ms.size();
    int hss = hs.size();

    if(mss > 0 && hss > 0 && mss == hss){
      CT<IntervalImp> dints = dom.ints();
      if(mss > dints.size()){
        int newLo, newStep, newHi;

        typename CT<IntervalImp>::iterator it = dints.begin();
        typename CTFloat::iterator itm = ms.begin();
        typename CTFloat::iterator ith = hs.begin();

        CT<IntervalImp> res;
        typename CT<IntervalImp>::iterator itres = res.begin(); 
     
        if(*itm == (int) ((*itm) * (*it.step))){ 
          while(it != dints.end()){ 
            newLo = (*it).lo * (*itm) + (*ith);
            newHi = (*it).hi * (*itm) + (*ith); 

            if(*itm > 0)
              newStep = (*it).step * (*itm); 

            else
              newStep = 1;

            IntervalImp intRes(newLo, newStep, newHi, false);
            itres = res.insert(itRes, intRes);
            ++itres;
 
            ++it;
          }

          return MultiInterImp(res);
        }

        else
          cerr << "Linear map slope is not compatible";
      }
    }

    cerr << "Expression dimension should be larger than domain dimension";
    return MultiInterImp();
  }

  IntervalImp domCompInt(IntervalImp &int2){
  }*/

  /*
  MultiInterImp domComp(LinearMapImp1 &lm2){
    CT<IntervalImp> mi1 = dom.ints();
    typename CT<IntervalImp>::iterator itdom1 = mi1.begin();
    CT<IntervalImp> mi2 = lm2.dom.ints();
    typename CT<IntervalImp>::iterator itdom2 = mi2.begin();

    typename CTFloat::iterator itms1 = ms.begin();
    typename CTFloat::iterator iths1 = hs.begin();
    typename CTFloat::iterator itms2 = lm2.ms.begin();
    typename CTFloat::iterator iths2 = lm2.hs.begin();
    
    MultiInterImp im = applyExprMI();
    typename CT<IntervalImp>::iterator itim = im.ints().begin();

    CT<IntervalImp> res;
    typename CT<IntervalImp>::iterator itres = res.begin();

    while(itdom2 != mi2.end()){
      if(*itms2 == 0){
        // Check if it is in lm1 domain
        if(*itdom1.isIn(*iths2)){
          itres = res.insert(itres, *itdom2);
          ++itres;
        }
      }
  
      else{
        IntervalImp aux = (*itim).cap(*itdom1);
        
      }

      ++itdom1;
      ++itdom2;
    }
  }

  LinearMapImp1 compose(const LinearMapImp1 &lm2){
    typename CTFloat::iterator it1ms = ms.begin(); 
    typename CTFloat::iterator it1hs = hs.begin(); 
    typename CTFloat::iterator it2ms = lm2.ms.begin(); 
    typename CTFloat::iterator it2hs = lm2.hs.begin(); 
    CTFloat msres;
    CTFloat hsres;
    typename CTFloat::iterator itResms = msres.begin(); 
    typename CTFloat::iterator itReshs = hsres.begin(); 

    int len = min(ms.size(), le2.ms.size());
    for(int i = 0; i < len; i++){
      itResms = msres.insert(itResms, (*it1ms) * (*it2ms));
      itReshs = hsres.insert(itReshs, (*it1ms) * (*it2hs) + (*it1hs)); 

      ++it1ms;
      ++it1hs;
      ++it2ms;
      ++it2hs;
    }

    return LinearMapImp1(msres, hsres);
  }
};

template <template<typename T, typename = std::allocator<T>> class CT,
          typename LinearMapImp, typename SetImp, typename ASetImp, typename MultiInterImp,
          typename NumImp>
struct LinearMapAbs{
  LinearMapAbs(){};
  LinearMapAbs(LinearMapImp lm){
    lmap = lm;
  }

  NumImp slopes(){
    return lmap.slopes;
  }

  NumImp intercepts(){
    return lmap.intercepts();
  }

  MultiInterImp applyExprMI(MultiInterImp x){
    return lmap.applyExprMI(x);
  }

  ASetImp applyExprAS(ASetImp x){
    return lmap.applyExprAS(x);
  }

  SetImp applyExpr(SetImp x){
    return lmap.applyExpr(x);
  }

  LinearMapAbs compose(const LinearMapAbs &le2){
    return LinearMapAbs(lmap.compose(le2));
  }

  LinearMapAbs invLMap(){
    return LinearMapAbs(lmap.invLMap());
  }

  private:
  LinearMapImp lmap;
};
*/
/*
  ASetImp applyExprAS(ASetImp x){
    MultiInterImp mi = x.as();

    return ASetImp(applyExpr(mi));
  }

  SetImp applyExpr(SetImp x){
    CT<ASetImp> ass = x.asets();       
    typename CT<ASetImp>::iterator it = ass.begin();

    SetImp res;

    while(it != ass.end()){
     SetImp auxs(applyExprAS(*it));
     res = res.cup(auxs);

     ++it;
    }

    return res;
  } 
  */
/*
  ASetImp applyExprAS(ASetImp x){
    MultiInterImp mi = x.as();

    return ASetImp(applyExpr(mi));
  }

  SetImp applyExpr(SetImp x){
    CT<ASetImp> ass = x.asets();       
    typename CT<ASetImp>::iterator it = ass.begin();

    SetImp res;

    while(it != ass.end()){
     SetImp auxs(applyExprAS(*it));
     res = res.cup(auxs);

     ++it;
    }

    return res;
  } 
  */

/// @brief Careful! This class exploit linear maps's properties. 
 /*
template <template<typename T, typename = allocator<T>>, 
          typename LMapImp, SetImp>
struct PWLMapImp1{
  CT<SetImp> dom;
  CT<LMapImp> expr;

  PWLMapImp1(){}
  PWLMapImp1(CT<SetImp> d, CT<LMapImp> e){
    if(d.size() == e.size()){
      dom = d;
      expr = e; 
    }

    else
      cerr << "Domain and maps should have the same size";
  }

  SetImp PWLImage(SetImp x){
    SetImp res;

    typename CT<SetImp>::iterator itdom = dom.begin();
    typename CT<LMapImp>::iterator itexpr = expr.begin();

    while(itdom != dom.end()){
      SetImp auxs = x.cap(*itdom);
      res = res.cup((*itexpr).applyExpr(auxs));
      
      ++itdom;
      ++itexpr;
    }

    return res;
  }

  PWLMapImp1 compose(const PWLMapImp1 &m2){
    CT<SetImp> newDom;
    typename CT<SetImp>::iterator itnewdom = newDom.begin();   
    CT<LMapImp> newExpr;
    typename CT<LMapImp>::iterator itnewexpr = newExpr.begin();

    typename CT<SetImp>::iterator itdom1 = dom.begin();
    typename CT<SetImp>::iterator itdom2 = m2.dom.begin();
    typename CT<LMapImp>::iterator itexpr1 = expr.begin();
    typename CT<LMapImp>::iterator itexpr2 = m2.expr.begin();

    while(itexpr1 != expr.end()){
      while(itexpr2 != m2.expr.end()){
        LMapImp invExpr = m2.expr.invLMap();

        NumImp ms = invExpr.slopes;
        NumImp hs = invExpr.intercepts();
        typename NumImp::iterator itms = ms.begin()
        typename NumImp::iterator iths = hs.begin()

        NumImp msRes;
        NumImp hsRes;
        typename NumImp::iterator itmsres = msRes.begin();
        typename NumImp::iterator ithsres = hsRes.begin();

        SetImp domRes;

        // Domain calculus
        while(itms != ms.end()){
          if(*itms == 0){
            // Image of m2 is in dom of m1?
            bool auxCond = ((*itexpr2).applyExpr(*itdom2)).cap(*itdom1)  
            if(auxCond){
            }
          }

          else{
          }

          ++itms;
          ++iths;
        }

        // Expression calculus

        ++itdom;
        ++itexpr;
      }
    }

    return MapImp1(newDom, newExp);
  }
 
  setData preImage(setData set1){
    Option<T2> preSet = expr.inv(set1);
    T2 newSet1;

    if(preSet)
      newSet1 = dom.cap(preSet.get);
    else
      newSet1 = dom.min();

    return newSet1;
  }

  MapAbs mininv(){
    Option<T1> opExpr = expr.inv();
    T1 eres;
    T3 dres = image();

    if(opExpr)
      eres = opExpr.get();
    else
      eres = eres.constant(dom.min()); 
  
    return MapAbs(eres, dres);
  setData pieceWiseImage(const setData &set){
    setData res;

    if(!set.empty()){
      typename setType::iterator itExpr = expr.begin();
      typename setType::iterator itDom = dom.begin();

      while(itDom != dom.end()){
        setData capRes = (*itDom).cap(set);
        setData im = (*itExpr).image(capRes);
        res = res.cup(res, im);    
  
        itExpr++;
        itDom++;
      }

    }

    return res;
  }
  }

}; 
*/

/*-----------------------------------------------------------------------------------------------*/
/*-----------------------------------------------------------------------------------------------*/
// Concrete classes
/*-----------------------------------------------------------------------------------------------*/
/*-----------------------------------------------------------------------------------------------*/

typedef IntervalImp1<list> Interval;
typedef MultiInterAbs<list, MultiInterImp1<list, Interval>, Interval> MultiInterval;
typedef AtomSetAbs<list, AtomSetImp1<list, MultiInterval>, MultiInterval> AtomSet;
typedef SetAbs<list, SetImp1<list, AtomSet>, AtomSet> Set;

#endif
