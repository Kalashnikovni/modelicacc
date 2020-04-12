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
*   classes can be defined, but the client modules
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
template <template<typename T1, typename = allocator<T1>> class CT,
          typename Data>
struct MultiInterImp1{
  CT<Data> inters;
  typedef typename CT<Data>::iterator intImpIt;

  MultiInterImp1(){
    CT<Data> emptyRes;
    inters = emptyRes;  
  }
  MultiInterImp1(CT<Data> is){
    inters = is;
  }

  bool empty(){
    return inters.empty();
  }

  CT<Data> ints(){
    return inters;
  };

  MultiInterImp1 cap(MultiInterImp1 &mi2){
    CT<Data> res;
    intImpIt itres = res.begin();

    intImpIt it1 = ints().begin();
    intImpIt it2 = mi2.ints().begin();
    int minLength = std::min(inters.size(), mi2.inters.size());
    for(int i = 0; i < minLength; i++){
      Data capres = (*it1).cap(*it2);

      if(capres.empty()){
        CT<Data> emptyRes;
        return MultiInterImp1<CT, Data>(emptyRes);
      }
      
      else
        itres = res.insert(itres, capres);

      ++it1;
      ++it2;    
    }

    return MultiInterImp1<CT, Data>(res);
  }

  CT<MultiInterImp1<CT, Data>> diff(MultiInterImp1 &mi2){
    MultiInterImp1<CT, Data> capres = cap(mi2);

    if(capres.empty()){
      CT<MultiInterImp1<CT, Data>> emptyRes;
      return emptyRes;
    }

    CT<MultiInterImp1<CT, Data>> res1;
    typename CT<MultiInterImp1<CT, Data>>::iterator itres1 = res1.begin();

    intImpIt it1 = ints().begin();
    intImpIt itcap = capres.ints().begin();

    for(int i = 0; i < ints().size(); i++){
      CT<Data> res2;
      intImpIt itres2 = res2.begin();

      intImpIt it11 = ints().begin();

      for(int j = 0; j < ints().size(); j++){
        int i1 = std::distance(ints().begin(), it1), i11 = std::distance(ints().begin(), it11);

        if(i1 < i11)
          it11 = res2.insert(itres2, *itcap);

        else if(i1 == i11){
          CT<Data> diffres = diff(*it1, *itcap);
          typename CT<Data>::iterator taux = diffres.begin();
          while(taux != diffres.end())
            it11 = res2.insert(itres2, *taux);
        }

        else if(i1 > i11)
          it11 = res2.insert(itres2, *it11);
      }

      itres1 = res1.insert(itres1, res2);

      ++itcap;
    }

    return res1;
  }
};

template <template <typename T1, typename = allocator<T1>> class CT,
          typename DataImp1,
          typename DataImp2>
struct MultiInterAbs{
  typedef CT<DataImp1> CTData;
  CTData inters;
  typedef typename CTData::iterator intImpIt;

  MultiInterAbs(){
    CTData emptyRes;
    inters = emptyRes;
  };
  MultiInterAbs(CTData is){
    inters = is;
  };

  bool empty(){
    return multiInterImp.empty();
  };

  CT<DataImp1> ints(){
    return multiInterImp.ints();
  };

  MultiInterAbs cap(MultiInterAbs &mi2){
    return multiInterImp.cap(mi2.multiInterImp);
  };

  CT<MultiInterAbs<CT, DataImp1, DataImp2>> diff(MultiInterAbs &mi2){
    return multiInterImp.diff(mi2.multiInterImp);
  };

  private:
  DataImp2 multiInterImp;
};

template <template<typename T1, typename = allocator<T1>> class CT,
          typename DataImp>
struct AtomSetImp1{
  DataImp aset;

  AtomSetImp1(){
    DataImp emptyRes();
    aset = emptyRes;
  };
  AtomSetImp1(DataImp as){
    aset = as;
  }

  bool empty(){
    return aset.empty();  
  }

  AtomSetImp1 cap(const AtomSetImp1 &aset2){
    return AtomSetImp1(aset.cap(aset2.aset));
  }

  CT<AtomSetImp1<CT, DataImp>> diff(const AtomSetImp1 &aset2){
    CT<AtomSetImp1<CT, DataImp>> res;
    typename CT<AtomSetImp1<CT, DataImp>>::iterator itres = res.begin();

    CT<DataImp> atomicDiff = aset.diff(aset2);

    if(atomicDiff.empty()){
      CT<AtomSetImp1<CT, DataImp>> emptyRes;
      return emptyRes;
    } 

    else{
      typename CT<DataImp>::iterator it = atomicDiff.begin();
      for(int i = 0; i < atomicDiff.size(); i++)
        itres = res.insert(AtomSetImp1(*it)); 
    }

    return res;
  }
};

template <template<typename T1, typename = allocator<T1>> class CT,
          typename DataImp1,
          typename DataImp2>
struct AtomSetAbs{
  DataImp1 aset;

  AtomSetAbs(){
    DataImp1 emptyRes();
    aset = emptyRes;
  };
  AtomSetAbs(DataImp1 as){
    aset = as;
  }

  bool empty(){
    return aset.empty();  
  }

  DataImp2 cap(const DataImp2 &aset2){
    return aset.cap(aset2.aset);
  }

  CT<DataImp2> diff(const DataImp2 &aset2){
    return aset.diff(aset2.aset);
  }
};

/*
// T1 should be an iterative type, e.g: lists, arrays, etc
template <template<typename T1, typename = std::allocator<T1>> class CT,
          typename Data>
struct SetAbs{
  typedef CT<Data> setType;
 
  SetAbs(){};
  SetAbs(setType ss){
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

  SetAbs addAtomSet(const Data &aset2){
    typename setType::iterator itsets = sets.begin();

    sets.insert(itsets, aset2);
    return SetAbs(sets);
  }

  SetAbs addAtomSets(const setType &sets2){
    setType res;
    typename setType::iterator it2 = sets2.begin();

    while(it2 != sets2.end()){
      res = res.addAtomSet(*it2);
      ++it2;
    }

    return res;
  }

  SetAbs cap(const SetAbs &set2){
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

  SetAbs diff(const SetAbs &set2){
    CT<Data> emptyCT;
    SetAbs res(emptyCT);
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

  SetAbs cup(const SetAbs &set2){
    SetAbs res = *this;
    SetAbs resDiff = diff(set2);

    if(!resDiff.empty())
      res.addAtomSets(resDiff);

    return res;
  }

  private:
  setType sets;
};

template <template<typename T1, typename = std::allocator<T1>> class CT,
          typename Data>
struct LinearMapAbs{
  typedef CT<int> CTInt;
 
  CTInt ms;
  CTInt hs;

  LinearMapAbs(){};
  LinearMapAbs(CTInt m, CTInt h){
    ms = m;
    hs = h;
  }

  Data applyExprInt(Data x){
    if(ms.size() > 0 && hs.size() > 0){
      int m1 = *(ms.begin()); 
      int h1 = *(hs.begin());
      int newLo, newStep, newHi;

      if(m1 > 0){
        newLo = x.lo * m1 + h1;
        newStep = x.step * m1; 
        newHi = x.hi * m1 + h1; 
      }

      else{
        newLo = h1;
        newStep = 1;
        newHi = h1;
      }

      return Data(newLo, newStep, newHi, false);
    }
  }

  CTInt applyExpr(CTInt x){
    CTInt res();
    typename CTInt::iterator itres = res.begin();
    typename CTInt::iterator it = x.begin();
    typename CTInt::iterator itm = ms.begin();
    typename CTInt::iterator ith = hs.begin();

    if(x.size() <= ms.size() && x.size() <= hs.size())
      while(it != x.end()){
        int y = (*itm) * (*it) + (*ith);
        ++itm;
        ++ith;
        itres = res.insert(itres, y);
      
        ++it;
      }

    return res;
  }

  AtomSetAbs<CT, Data> atomImage(const AtomSetAbs<CT, Data> &atomSet){
    CT<Data> res;
    typename CT<Data>::iterator itres = res.begin();
    typename CT<Data>::iterator it = atomSet.aset.inters.begin();  

    while(it != atomSet.aset.inters.end()){
      Data newInt = applyExprInt(*it);
      itres = res.insert(itres, newInt);

      ++it;
    }

    MultiInterAbs<CT, Data> newRes(res); 
    return AtomSetAbs<CT, Data>(newRes);
  }

  SetAbs<CT, Data> image(const SetAbs<CT, Data> &set){
    SetAbs<CT, Data> res();

    if(!set.empty()){
      CT<AtomSetAbs<CT, Data>> doms = set.sets;
      typename CT<AtomSetAbs<CT, Data>>::iterator itDoms = doms.begin();

      while(itDoms != doms.end()){
        AtomSetAbs<CT, Data> atomIm = atomImage(*itDoms); 
        SetAbs<CT, Data> im(atomIm);
        res = res.cup(res, im);    
  
        itDoms++;
      }
    }

    return res;
  }

  LinearMapAbs compose(const LinearMapAbs &le2){
    typename CTInt::iterator it1ms = ms.begin(); 
    typename CTInt::iterator it1hs = hs.begin(); 
    typename CTInt::iterator it2ms = le2.ms.begin(); 
    typename CTInt::iterator it2hs = le2.hs.begin(); 
    CTInt msres;
    CTInt hsres;
    typename CTInt::iterator itResms = msres.begin(); 
    typename CTInt::iterator itReshs = hsres.begin(); 

    int len = min(ms.size(), le2.ms.size());
    for(int i = 0; i < len; i++){
      itResms = msres.insert(itResms, (*it1ms) * (*it2ms));
      itReshs = hsres.insert(itReshs, (*it1ms) * (*it2hs) + (*it1hs)); 

      ++it1ms;
      ++it1hs;
      ++it2ms;
      ++it2hs;
    }

    return LinearMapAbs(msres, hsres);
  }
};

template <typename setType, typename exprType>
struct MapAbs{
  setType dom;
  exprType expr;

  MapAbs(){}
  MapAbs(setType d, exprType e){
    if(d.size() == e.size()){
      dom = d;
      expr = e; 
    }

    else
      cerr << "Domain and maps should have the same size";
  }

  MapAbs compose(const MapAbs &m2){
    setType newDom =   
    exprType newExp = expr.compose(m2.expr);

    return MapAbs(newDom, newExp);
  }*/
 
/*
  setData preImage(setData set1){
    Option<T2> preSet = expr.inv(set1);
    T2 newSet1;

    if(preSet)
      newSet1 = dom.cap(preSet.get);
    else
      newSet1 = dom.min();

    return newSet1;
  }
  MapAbs compose(const MapAbs &m2){
    T1 eres = expr.compose(m2.expr());    
    T2 dres = m2.preImage(dom.cap(m2.image(m2.dom())));

    return MapAbs(eres, dres);
  }
*/

 /*
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
typedef MultiInterAbs<list, Interval, MultiInterImp1<list, Interval>> MultiInterval;
typedef AtomSetAbs<list, MultiInterval, MultiInterImp1<list, Interval>> AtomSet;
//typedef AtomSetAbs<list, MultiInterval> AtomSet;
//typedef SetAbs<list, AtomSet> Set;
// --> typedef SetAbs<list, AtomSetNew> Set;

//typedef LinearMapAbs<list, Interval> LinearMap;



#endif
