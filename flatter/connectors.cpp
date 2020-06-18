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

#include <iostream>
#include <string>

#include <flatter/connectors.h>

#include <util/ast_visitors/contains_expression.h>

using namespace std;
using namespace Modelica;
using namespace Modelica::AST;

#define NameToRef(X) Reference(Reference(), X, Option<ExpList>())
#define RefIndex(X, Y) Reference(Reference(), X, Option<ExpList>(Y))
#define PrintOpt(N) (N ? N.get() : "{}")

Connectors::Connectors(MMO_Class &c) : mmoclass_(c){}

member_imp(Connectors, int, vCount);
member_imp(Connectors, SBGraph, G);
member_imp(Connectors, MMO_Class, mmoclass);

/*|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||*/
/*-----------------------------------------------------------------------------------------------*/
// Debugging functions --------------------------------------------------------------------------//
/*-----------------------------------------------------------------------------------------------*/
/*|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||*/

void Connectors::debug(std::string filename){
  GraphPrinter gp(G(), -1);

  gp.printGraph(filename);
  cout << "Generated Connect Graph written to " << filename << endl;
}

/*|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||*/
/*-----------------------------------------------------------------------------------------------*/
// Create graph  --------------------------------------------------------------------------------//
/*-----------------------------------------------------------------------------------------------*/
/*|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||*/

void Connectors::createGraph(EquationList &eqs){
  EquationList to_delete;

  foreach_(Equation &eq, eqs){
    if(is<Connect>(eq)){
      to_delete.push_back(eq);
      Connect co = boost::get<Connect>(eq);
      Expression eleft = co.left(), eright = co.right();
  
      Pair<Expression, ExpOptList> left = separate(eleft);
      Pair<Expression, ExpOptList> right = separate(eright);

      cout << left << "\n";
      cout << right << "\n";

      ExpOptList range1 = get<1>(left);
      ExpOptList range2 = get<1>(right);

      if(checkRanges(range1, range2)){
        MultiInterval mi1;
        foreach_(Expression e1, range1){
          Interval l = parseExpToInter(e1);
          mi1.addInter(l);
        }

        MultiInterval mi2;
        foreach_(Expression e2, range2){
          Interval r = parseExpToInter(e2);
          mi2.addInter(r);
        }

        //updateGraph(mi1, mi2);
      }
    }
  }
}

// Get expression and range
Pair<Expression, ExpOptList> Connectors::separate(Expression e){
  Reference reference;
  UnaryOpType op = Plus;
  bool isU = false;

  if(is<UnaryOp>(e)){
    UnaryOp u = boost::get<UnaryOp>(e);
    if(is<Reference>(u.exp())){
      isU = true;
      reference = boost::get<Reference>(u.exp());
      op = u.op();
    } 
    else
      std::cerr << "ERROR: Deberia llegar una Reference" << std::endl;
  } 
 
  else if (is<Reference>(e))
    reference = boost::get<Reference>(e);

  Ref refs = reference.ref();
  if(refs.size() > 1) 
    std::cerr << "ERROR: No deberia haber llamadas a miembros en connectors" << std::endl;
  RefTuple rf = refs.front();
  ExpOptList opti;
  if(get<1>(rf).size() > 0) 
    opti = Option<ExpList>(get<1>(rf));
  Expression r = Reference(Reference(), get<0>(rf), Option<ExpList>());
  if(isU) 
    r = UnaryOp(r, op);
  return Pair<Expression, ExpOptList>(r, opti);
}

// Check if only one variable is used at each subscript
bool Connectors::checkRanges(ExpOptList range1, ExpOptList range2){
  std::vector<Name> vars = mmoclass().variables();

  if(range1 && range2){
    ExpList r1 = range1.get();
    ExpList r2 = range2.get();
 
    if(r1.size() == 0 || r2.size() == 0)
      return true;

    else if(r1.size() != r2.size()){
      cerr << "Unmatched dimensions in equation connect" << endl;
      return false;
    }

    else{
      ExpList::iterator it1 = r1.begin(), it2 = r2.begin();

      while(it1 != r1.end()){
        foreach_(Name n1, vars){
          Expression e1(n1);
          ContainsExpression co1(e1);
 
          bool cn11 = Apply(co1, *it1);
          bool cn21 = Apply(co1, *it2);

          // This loop checks that there is only one variable at each subscript
          foreach_(Name n2, vars){
            Expression e2(n2);
            ContainsExpression co2(e2);

            bool cn12 = Apply(co2, *it1);
            bool cn22 = Apply(co2, *it2);

            if(((cn11 && cn12) || (cn21 && cn22)) && (n1 != n2)){
              cerr << "Only one variable permitted at subscript";
              return false;
            }
          }
        }

        ++it1;
        ++it2;
      }
    }
  }

  return true;
}

Interval Connectors::parseExpToInter(Expression e){
  if(is<Integer>(e) || is<Boolean>(e)){
    EvalExpression evexp(mmoclass().syms());  
    int res = Apply(evexp, e);

    Interval i(res, 1, res);
    return i;
  }

  else if(is<UnaryOp>(e)){
    EvalExpression evexp(mmoclass().syms());
    UnaryOp uop = get<UnaryOp>(e);
    Expression euop = uop.exp();
    int res = Apply(evexp, euop);

    if(uop.op() == Not){
      if(res){
        Interval i(0, 1, 0);
        return i; 
      }
      else{
        Interval i(1, 1, 1);
        return i;
      }
    }
 
    else if(uop.op() == Plus){
      Interval i(res, 1, res);
      return i;
    }

    else{
      cerr << "Invalid negative subscript" << endl;
      Interval i;
      return i;
    }
  }

  else if(is<BinOp>(e)){
    BinOp bop = get<BinOp>(e);
    Expression left = bop.left(), right = bop.right();
    //EvalExpression evexp(mmoclass().syms());
    EvalExpFlatter evexpf(mmoclass().syms());
    int res = Apply(evexpf, e); 
    Interval i(res, 1, res);

    return i;
  }

  else if(is<Range>(e)){
    Range r = get<Range>(e);
    EvalExpression evexp(mmoclass().syms());
    Expression efrom = r.start();
    Expression eto = r.end();
    int from = Apply(evexp, efrom);
    int to = Apply(evexp, eto);
    int st = 1;

    if(r.step())
      st = Apply(evexp, r.step().get());

    Interval i(from, st, to);
    return i;
  }

  else if(is<SubAll>(e)){ // TODO
    Interval i;
    return i;
  }

  Interval i;
  return i;
}
