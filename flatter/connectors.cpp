
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

#include <flatter/connectors.h>
#include <iostream>

using namespace std;
using namespace Modelica;
using namespace Modelica::AST;

#define NameToRef(X) Reference(Reference(), X, Option<ExpList>())
#define RefIndex(X, Y) Reference(Reference(), X, Option<ExpList>(Y))
#define PrintOpt(N) (N ? N.get() : "{}")

Connectors::Connectors(MMO_Class &c) : mmoclass_(c){}

member_imp(Connectors, SetBasedGraph, G);
member_imp(Connectors, MMO_Class, mmoclass);

/*|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||*/
/*-----------------------------------------------------------------------------------------------*/
// Debugging functions --------------------------------------------------------------------------//
/*-----------------------------------------------------------------------------------------------*/
/*|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||*/

void Connectors::debug(std::string filename){
  GraphPrinter gp(G(), 1);

  gp.printGraph(filename);
  cout << "Generated Connect Graph written to " << filename << endl;
}

/*|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||*/
/*-----------------------------------------------------------------------------------------------*/
// Create graph  --------------------------------------------------------------------------------//
/*-----------------------------------------------------------------------------------------------*/
/*|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||*/

void Connectors::createGraph(EquationList &eqs, OptExpList range){
  EquationList to_delete;

  foreach_(Equation &eq, eqs){
    if(is<Connect>(eq)){
      to_delete.push_back(eq);
      Connect co = boost::get<Connect>(eq);
      Expression eleft = co.left(), eright = co.right();
  
      Pair<Expression, ExpOptList> left = separate(eleft);
      Pair<Expression, ExpOptList> right = separate(eright);
    }
  }
}

Pair<Expression, ExpOptList> Connectors::separate(Expression e){
  Reference reference;
  UnaryOpType op = Plus;
  bool isU = false;

  if (is<UnaryOp>(e)){
    UnaryOp u = boost::get<UnaryOp>(e);
    if (is<Reference>(u.exp())) {
      isU = true;
      reference = boost::get<Reference>(u.exp());
      op = u.op();
    } else
      std::cerr << "ERROR: Deberia llegar una Reference" << std::endl;
  } 
 
  else if (is<Reference>(e))
    reference = boost::get<Reference>(e);

  Ref refs = reference.ref();
  if (refs.size() > 1) std::cerr << "ERROR: No deberia haber llamadas a miembros en connectors" << std::endl;
  RefTuple rf = refs.front();
  ExpOptList opti;
  if (get<1>(rf).size() > 0) opti = Option<ExpList>(get<1>(rf));
  Expression r = Reference(Reference(), get<0>(rf), Option<ExpList>());
  if (isU) r = UnaryOp(r, op);
  return Pair<Expression, ExpOptList>(r, opti);
}
