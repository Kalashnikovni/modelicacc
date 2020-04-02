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

// LinearFuncs
LinearFunc::LinearFunc(int m, int h) : m_(m), h_(h){};

member_imp(LinearFunc, int, m);
member_imp(LinearFunc, int, h);

// LFRange
LFRange::LFRange(bool isEmpty) : empty_(isEmpty){};

// We "normalize range" to simplify calculations
LFRange::LFRange(int m, int h, int lo, int step, int hi, bool empty) : empty_(empty){
  if(step != 0 && lo <= hi && !empty){
    int newM = m * step;
    set_m(newM);

    int newH = m * lo + h; 
    set_h(newH);

    int newHi = (hi - lo) / step;
    set_hi(newHi); 
  }

  cerr << "Wrong values for subscript" << endl;
  set_empty(true);
}

member_imp(LFRange, int, hi);
member_imp(LFRange, bool, empty);

LFRange LFRange::emptySet(){
  return LFRange();
}

bool LFRange::empty(){
  return empty(); 
}

LFRange LFRange::cap(const LFRange &set2){
  int newH = -1, newStep = lcm(m(), set2.m());
  int endHi = std::min(m() * hi() + h(), set2.m() * set2.hi() + set2.h());

  if(!empty() && !set2.empty())
    for(int i = 0; i < newStep / m(); i++){
      int res1 = m() * i + h();
      float res2 = (res1 - set2.h()) / set2.m();

      if(res2 == (int) res2){
        newH = res1;
        break;
      }
    }

  else
    return LFRange(true);

  if(newH < 0)
    return LFRange(true);

  while(set2.h() > newH)
    newH += newStep;

  int newEnd = newH;
  while(newEnd < endHi)
    newEnd += newStep;

  return LFRange(1, newH, 0, newStep, newEnd, false);
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

LFRange LFRange::min(){
  return LFRange(m(), h(), 0, 1, 1, false);
}

// LFExpr
LFExpr::LFExpr(int m, int h){
  set_m(m);
  set_h(h);
};

LFRange LFExpr::applyExpr(LFRange set1){
  if(set1.empty())
    return LFRange(true);

  return LFRange(m() * set1.m(), m() * set1.h() + h(), 0, 1, set1.hi() + 1, false);
}

LFExpr LFExpr::compose(const LFExpr &e2){
  return LFExpr(e2.m() * m(), m() * e2.h() + h());
}

Option<LFExpr> LFExpr::inv(){
  if(m() == 0)
    return Option<LFExpr>();

  return Option<LFExpr>();
}

LFExpr constant(LFRange min){
  return LFExpr(0, min.h());
}
