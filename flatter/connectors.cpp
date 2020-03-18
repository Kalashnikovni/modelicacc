
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

#include <boost/type_traits/remove_cv.hpp>

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

      ExpOptList range1 = get<1>(left);
      ExpOptList range2 = get<1>(right);

      if(checkRanges(range1, range2)){
        //SetVertexLF l = createVertex(get<0>(left), range1);
        //SetVertexLF r = createVertex(get<0>(right), range2);
      }
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

            if((cn11 && cn12) || (cn21 && cn22)){
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

Option<LinearFunc> Connectors::parseExpToLF(Expression e){
  if(is<Integer>(e) || is<Boolean>(e)){
    EvalExpression evexp(mmoclass().syms());  
    int res = Apply(evexp, e);
    return Option<LinearFunc>(LinearFunc(0, 0, res, false));
  }

  else if(is<UnaryOp>(e)){
    EvalExpression evexp(mmoclass().syms());
    UnaryOp uop = get<UnaryOp>(e);
    Expression euop = uop.exp();
    int res = Apply(evexp, euop);

    if(uop.op() == Not)
      if(res)
        return Option<LinearFunc>(LinearFunc(0, 0, 0, false)); 
      else
        return Option<LinearFunc>(LinearFunc(0, 0, 1, false));
 
    else if(uop.op() == Plus)
      return Option<LinearFunc>(LinearFunc(0, 0, res, false));

    else{
      cerr << "Invalid negative subscript" << endl;
      return Option<LinearFunc>();
    }
  }

  else if(is<BinOp>(e)){
    BinOp bop = get<BinOp>(e);
    Expression left = bop.left(), right = bop.right();
    EvalExpression evexp(mmoclass().syms());
    EvalExpFlatter evexpf(mmoclass().syms());
    int res = Apply(evexpf, e); 

    switch(bop.op()){
      case Or:
        return Option<LinearFunc>(LinearFunc(0, 0, res, false));
      case And:
        return Option<LinearFunc>(LinearFunc(0, 0, res, false));
      case Lower:
        return Option<LinearFunc>(LinearFunc(0, 0, res, false));
      case LowerEq:
        return Option<LinearFunc>(LinearFunc(0, 0, res, false));
      case Greater:
        return Option<LinearFunc>(LinearFunc(0, 0, res, false));
      case GreaterEq:
        return Option<LinearFunc>(LinearFunc(0, 0, res, false));
      case CompEq:
        return Option<LinearFunc>(LinearFunc(0, 0, res, false));
      case CompNotEq:
        return Option<LinearFunc>(LinearFunc(0, 0, res, false));
      case Add:
        parseAddSub(left, right, bop.op());
      case Sub:
        parseAddSub(left, right, bop.op());
      case Mult: //TODO
        return Option<LinearFunc>();
      default:
        return Option<LinearFunc>();
    }
    return Option<LinearFunc>();
  }

  else if(is<Range>(e)){
    Range r = get<Range>(e);
    EvalExpression evexp(mmoclass().syms());
    Expression efrom = r.start();
    Expression eto = r.end();
    int from = Apply(evexp, efrom);
    int to = Apply(evexp, eto);

    if(r.step()){
      int st = Apply(evexp, r.step().get());
      return Option<LinearFunc>(LinearFunc((int)((to - from) / st), st, from, false));
    }

    else
      return Option<LinearFunc>(LinearFunc(to - from, 1, 0, false));
  }

  else if(is<SubAll>(e)){ // TODO
    return Option<LinearFunc>();
  }

  else return Option<LinearFunc>();

  return Option<LinearFunc>();
}

bool Connectors::isParseNum(Expression e){
  if(is<Integer>(e))
    return true;

  else if(is<UnaryOp>(e)){
    UnaryOp uop = get<UnaryOp>(e);
    if(uop.op() == Plus || uop.op() == Minus)
      return true; 
  }

  else if(is<BinOp>(e)){
    BinOp bop = get<BinOp>(e);
    Expression l = bop.left(), r = bop.right(); 

    if(bop.op() == Add)
      return (isParseNum(l)) && (isParseNum(r));
    else if(bop.op() == Mult){ // Check that it is linear, not quadratic, etc.
      std::vector<Name> vars = mmoclass().variables(); 

      foreach_(Name n, vars){
        Expression en(n);
        ContainsExpression co(en);

        bool varl = Apply(co, l);
        bool varr = Apply(co, r);

        if(varl && varr)
          return false;
        else{
          bool res1 = isParseNum(l);
          bool res2 = isParseNum(r);
          return res1 && res2; 
        }
      }
    }
  }

  return false;
}

Option<LinearFunc> Connectors::parseAddSub(Expression l, Expression r, BinOpType bopt){ 
  if(isParseNum(l) && isParseNum(r)){
    Option<LinearFunc> lf1 = parseExpToLF(l);
    Option<LinearFunc> lf2 = parseExpToLF(r);

    if(lf1 && lf2){
      LinearFunc res1 = lf1.get(); 
      LinearFunc res2 = lf2.get();

      if(res1.hi() != res2.hi())
        cerr << "Subscripts should have the same range" << endl;
      else if(bopt == Add)
        return Option<LinearFunc>(LinearFunc(res1.hi(), res1.m() + res2.m(), 
                                             res1.h() + res2.h(), false));
      else if(bopt == Sub)
        return Option<LinearFunc>(LinearFunc(res1.hi(), res1.m() - res2.m(), 
                                             res1.h() - res2.h(), false));
    }
  }

  return Option<LinearFunc>();
}

Option<SetVertexLF> Connectors::createVertex(Expression e, ExpOptList range){
  MultiDimLF vs;

  if(range){
    ExpList rge = range.get();
    foreach_(Expression r, rge){
      Option<LinearFunc> rLF = parseExpToLF(r);
      if(rLF)
        vs.push_back(rLF.get());
      else{
        cerr << "Subscript " << r << " is not a linear expression" << endl;
        return Option<SetVertexLF>();
      }
    }
  }

  std::string s = std::to_string(vCount());
  return Option<SetVertexLF>(SetVertexLF(vs, "V" + s)); 
} 
