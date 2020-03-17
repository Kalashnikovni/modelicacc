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

// Set vertex
member_imp(SetVertex, int, card);
member_imp(SetVertex, int, index);
member_imp(SetVertex, std::string, name);

std::ostream &operator<<(std::ostream &os, const SetVertex &sv){
 os << sv.name();
 return os;
}

// LinearFunc
LinearFunc::LinearFunc(bool isEmpty) : empty_(isEmpty){};

LinearFunc::LinearFunc(int to, int slope, int intercept, bool isEmpty) :
 m_(slope), h_(intercept), empty_(isEmpty){
 if(to >= 0){
  set_hi(to);
 }

 else
  std::cerr << "The bounds of the linear func are swapped\n";
};

member_imp(LinearFunc, bool, empty);
member_imp(LinearFunc, int, hi);
member_imp(LinearFunc, int, m);
member_imp(LinearFunc, int, h);

bool LinearFunc::operator==(const LinearFunc &other) const{
 if(empty() && other.empty())
  return true;
 else if (!empty() && !other.empty())
  return (hi() == other.hi() && m() == other.m() && h() == other.h());
 
 return false;
}

LinearFunc LinearFunc::LFCap(const LinearFunc &other){
 if(!other.empty() && !empty()){
  if(other.m() == m() && other.h() == h()){
   int newEnd = std::min(hi(), other.hi());
   LinearFunc(newEnd, m(), h(), false); 
  }

  else if(other.m() == m()) return LinearFunc(true);

  else{
   int newM = m() * (other.h() - h()) / (m() - other.m()) + h(); 
   return LinearFunc(0, newM, 0, false); 
  }
 }

 return LinearFunc(true); 
}

// Set vertices implemented as linear functions
SetVertexLF::SetVertexLF(MultiDimLF vs, std::string name) : vertices_(vs){
 int n = 1;
 foreach_(LinearFunc lf, vs) n = n * (lf.hi() + 1);

 set_card(n);
 set_name(name);
};

member_imp(SetVertexLF, MultiDimLF, vertices);

bool SetVertexLF::operator==(const SetVertexLF &other) const {
 return (other.vertices() == vertices());
}
 
// Set Edge
member_imp(SetEdge, int, card);
member_imp(SetEdge, std::string, name);

std::ostream &operator<<(std::ostream &os, const SetEdge &se){
 os << se.name();
 return os;
}

// Set Edge LF
SetEdgeLF::SetEdgeLF(EdgeLF &es, std::string name, bool inv) : edges_(es), invert_(inv){
 set_name(name);

 this->card_ref() = 1;
 foreach_(MultiDimLF m, es)
  foreach_(LinearFunc lf, m)
   this->card_ref() = this->card_ref() * (lf.hi() + 1);
};

member_imp(SetEdgeLF, EdgeLF, edges);

bool SetEdgeLF::operator==(const SetEdgeLF &other) const {
 return (other.edges() == edges());
}
