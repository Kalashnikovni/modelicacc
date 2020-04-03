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
*   Scalar definition has been removed since the set-based
*   approach is more general. First, abstract structures are
*   defined. Then, each abstraction is refined, restricting
*   the types of expressions used for vectorized variables,
*   but enhancing performance of the algorithms.
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

/*-----------------------------------------------------------------------------------------------*/
/*-----------------------------------------------------------------------------------------------*/
// Abstract classes
/*-----------------------------------------------------------------------------------------------*/
/*-----------------------------------------------------------------------------------------------*/

int lcm(int a, int b);

template <typename T>
struct IntervalAbs{
  IntervalAbs(){};
  IntervalAbs(bool isEmpty){
    empty = isEmpty;
  }
  IntervalAbs(T inter){
    interval = inter;
  }

  IntervalAbs cap(const IntervalAbs &inter2){
    return interval.cap(inter2.interval);
  }

  bool isEmpty(){
    return empty;
  }

  T getInterval(){
    return interval;
  }

  private:
  bool empty;
  T interval;
};

template <typename T, template<typename TI> class IterativeType>
struct MultiInterAbs{
  MultiInterAbs(){};
  MultiInterAbs(bool isEmpty){};
  MultiInterAbs(IterativeType<T> is){
    inters = is;
  }

  MultiInterAbs cap(MultiInterAbs &mi2){
    IterativeType<T> res();
    typename IterativeType<T>::iterator itres = res.begin();

    typename IterativeType<T>::iterator it1 = inters.begin();
    typename IterativeType<T>::iterator it2 = mi2.inters.begin();
    int minLength = std::min(inters.size(), mi2.inters.size());
    for(int i = 0; i < minLength; i++){
      T capres = *it1.cap(*it2);

      if(capres.empty())
        return IterativeType<T>();
      
      else{
        res.insert(itres, capres);
        ++itres;
      }

      ++it1;
      ++it2;    
    }

    return res;
  }

 /*
  MultiInterAbs diff(const IntervalAbs<T> &i1, const IntervalAbs<T> &i2){
    IntervalAbs<T> capi1i2 = i1.cap(i2);
    
    if(capi1i2.empty())
      return IterativeType<T>();

    int ninters = 0;

    if(capi1i2.interval().start() > i1.interval().start()){
      IntervalAbs<T>
    }
  }*/

  private:
  IterativeType<T> inters;
};

template <typename T>
struct AtomSetAbs{
  AtomSetAbs(){};
  AtomSetAbs(T s){
    set = s;
  }

  AtomSetAbs cap(const AtomSetAbs &aset2){
    return set.cap(aset2.set);
  }

  private:
  T set;
};

// T1 should be an iterative type, e.g: lists, arrays, etc
template <typename T, template<typename TI> class IterativeType>
struct SetAbs{
  SetAbs(){};
  SetAbs(IterativeType<T> ss){
    sets = ss;
  };

  /*
  SetAbs cap(const SetAbs &set2){
  }*/

  private:
  IterativeType<T> sets;
};

// First type is the type of expression, second the type of the domain,
// and finally the type of the image. T2 and T3 should be sets, but 
// represented intensively.
template <typename T1, typename T2, typename T3>
struct MapAbs{
  MapAbs(){};
  MapAbs(T1 expr, T2 dom); 

  T3 image(T2 set1){
    T2 newSet1 = dom.cap(set1);
   
    if(newSet1.empty()){
      T3 emptySet();
      return emptySet;
    }

    return expr.applyExpr(newSet1);
  }

  T2 preImage(T3 set1){
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

  MapAbs mininv(){
    Option<T1> opExpr = expr.inv();
    T1 eres;
    T3 dres = image();

    if(opExpr)
      eres = opExpr.get();
    else
      eres = eres.constant(dom.min()); 
  
    return MapAbs(eres, dres);
  }

  private:
  T1 expr;
  T2 dom;
};

/*-----------------------------------------------------------------------------------------------*/
/*-----------------------------------------------------------------------------------------------*/
// Concrete classes
/*-----------------------------------------------------------------------------------------------*/
/*-----------------------------------------------------------------------------------------------*/

// Intervals
struct IntervalImp{
  IntervalImp(){};
  IntervalImp(bool isEmpty);
  IntervalImp(int lo, int step, int hi, bool isEmpty);

  member_(bool, empty);
  member_(int, lo);
  member_(int, step);
  member_(int, hi);

  bool empty();
  IntervalImp cap(const IntervalImp &inter2);
  /// @brief Return a new interval with only the min element
  IntervalImp min();

  printable(IntervalImp);
};

typedef IntervalAbs<IntervalImp> Interval;

typedef MultiInterAbs<Interval, List> MultiInterval;

//typedef AtomSetAbs<AtomSetImp> AtomSet;

// Sets
//struct SetImp{};

//typedef SetAbs<SetImp> Set;

// Representation of linear functions
struct LinearFunc{
  LinearFunc(){};
  LinearFunc(int m, int h);

  member_(int, m); 
  member_(int, h); 
};

// Type of expression
/*
struct LFExpr : public LinearFunc{
  LFExpr() : LinearFunc(){};
  LFExpr(int m, int h); 

  LFRange applyExpr(LFRange set1);
  LFExpr compose(const LFExpr &e2);
  Option<LFExpr> inv();
  LFExpr constant(LFRange min);
};
*/ 

 
/*-----------------------------------------------------------------------------------------------*/
/*-----------------------------------------------------------------------------------------------*/
// Edges
/*-----------------------------------------------------------------------------------------------*/
/*-----------------------------------------------------------------------------------------------*/


/*-----------------------------------------------------------------------------------------------*/
/*-----------------------------------------------------------------------------------------------*/
// Graphs
/*-----------------------------------------------------------------------------------------------*/
/*-----------------------------------------------------------------------------------------------*/

#endif
