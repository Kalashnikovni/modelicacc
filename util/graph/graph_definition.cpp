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

/*! \file graph_definition.cpp
 *  This file implements the methods to work with set-based graphs.
 */

#include <iostream>

#include <util/debug.h>
#include <util/graph/graph_definition.h>

using namespace std;


// AtomSet
IntervalImp::IntervalImp(bool isEmpty) : empty_(isEmpty){};

IntervalImp::IntervalImp(int lo, int step, int hi, bool empty) 
  : empty_(empty), lo_(lo), step_(step){
  if(lo <= hi && !empty){
    int rem = std::fmod(lo - hi, step); 
    set_hi(hi - rem); 
  }

  else{ 
    cerr << "Wrong values for subscript" << endl;
    set_empty(true);
  }
}

member_imp(IntervalImp, bool, empty);
member_imp(IntervalImp, int, lo);
member_imp(IntervalImp, int, step);
member_imp(IntervalImp, int, hi);

bool IntervalImp::empty(){
  return empty(); 
}

IntervalImp IntervalImp::cap(const IntervalImp &inter2){
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
    return IntervalImp(true);

  if(newLo < 0)
    return IntervalImp(true);

  while(inter2.lo() > newLo)
    newLo += newStep;

  int newEnd = newLo;
  while(newEnd < minEnd)
    newEnd += newStep;

  return IntervalImp(newLo, newStep, newEnd, false);
}

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

IntervalImp IntervalImp::min(){
  return IntervalImp(lo(), 1, lo(), false);
}

std::ostream &operator<<(std::ostream &out, const IntervalImp &inter){
  out << "[" << inter.lo() << ":" << inter.step() << ":" << inter.hi() << "]";
  return out;
}

// LinearFuncs
LinearFunc::LinearFunc(int m, int h) : m_(m), h_(h){};

member_imp(LinearFunc, int, m);
member_imp(LinearFunc, int, h);

// LFExpr
/*
LFExpr::LFExpr(int m, int h){
  set_m(m);
  set_h(h);
};

AtomSet LFExpr::applyExpr(AtomSet set1){
  if(set1.empty())
    return AtomSet(true);

  return AtomSet(m() * set1.m(), m() * set1.h() + h(), 0, 1, set1.hi() + 1, false);
}

LFExpr LFExpr::compose(const LFExpr &e2){
  return LFExpr(e2.m() * m(), m() * e2.h() + h());
}

Option<LFExpr> LFExpr::inv(){
  if(m() == 0)
    return Option<LFExpr>();

  return Option<LFExpr>();
}

LFExpr constant(AtomSet min){
  return LFExpr(0, min.h());
}*/
