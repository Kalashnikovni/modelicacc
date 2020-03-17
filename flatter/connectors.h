
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
#ifndef CONNECTOR_H
#define CONNECTOR_H

#include <boost/type_traits/remove_cv.hpp>
#include <boost/variant/get.hpp>
#include <boost/variant/static_visitor.hpp>
#include <set>

#include <ast/ast_types.h>
#include <ast/element.h>
#include <ast/equation.h>
#include <ast/statement.h>
#include <flatter/class_finder.h>
#include <mmo/mmo_class.h>

#include <util/ast_visitors/eval_expression.h>
#include <util/ast_visitors/replace_expression.h>
#include <util/ast_visitors/constant_expression.h>
#include <util/type.h>
#include <util/graph/graph_definition.h>
#include <util/graph/graph_printer.h>

using namespace Modelica;
using namespace Modelica::AST;
using namespace Graph;

typedef Option<ExpList> ExpOptList; 

class Connectors{
  public:
  Connectors(MMO_Class &c);

  void debug(std::string filename);

  void createGraph(EquationList &eqs, OptExpList range);
  Pair<Expression, ExpOptList> separate(Expression e);

  private:
  member_(SetBasedGraph, G);
  member_(MMO_Class, mmoclass);
};

#endif
