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

#include <ast/ast_types.h>
#include <ast/equation.h>
#include <util/table.h>

namespace ICL = boost::icl;

/*-----------------------------------------------------------------------------------------------*/
/*-----------------------------------------------------------------------------------------------*/
// Vertices
/*-----------------------------------------------------------------------------------------------*/
/*-----------------------------------------------------------------------------------------------*/

struct SetVertex{
 SetVertex(){};

 member_(int, card); 
 /// @brief This is used for debugging purposes
 member_(int, index);
 member_(std::string, name); // Used to pretty print

 printable(SetVertex);
};

/* \struct LinearFunc
*  Array access can be thought as a function with domain in the
*  naturals from 0 to the size of the array minus 1. In this case functions
*  are restricted to only linear because intersection can be
*  done in constant time.
*
*  The domain of the function will be determined by the form of
*  the loop header, and the range by the expression that performs
*  access to the vectorized variable.
*/
struct LinearFunc{
 LinearFunc(){};
 LinearFunc(bool isEmpty);
 LinearFunc(int to, int slope, int intercept, bool isEmpty);

 member_(bool, empty);
 member_(int, hi); // Domain
 member_(int, m); // Range
 member_(int, h); // Range

 comparable(LinearFunc);
 LinearFunc LFCap(const LinearFunc &other);
};

typedef std::vector<LinearFunc> MultiDimLF;

struct SetVertexLF : public SetVertex{
 SetVertexLF(){};
 SetVertexLF(MultiDimLF vs, std::string name);

 member_(MultiDimLF, vertices);

 comparable(SetVertexLF);
};
 
/*-----------------------------------------------------------------------------------------------*/
/*-----------------------------------------------------------------------------------------------*/
// Edges
/*-----------------------------------------------------------------------------------------------*/
/*-----------------------------------------------------------------------------------------------*/

struct SetEdge{
 SetEdge(){};

 member_(int, card);
 member_(std::string, name); // Used to pretty print

 printable(SetEdge);
};

/// @brief It is a vector since a variable can occur more than once in an equation, we need order
typedef std::vector<MultiDimLF> EdgeLF;

struct SetEdgeLF : public SetEdge{
 SetEdgeLF(){};
 SetEdgeLF(EdgeLF &es, std::string name, bool inv);

 /// @brief For 1_N connections
 member_(bool, invert);
 member_(EdgeLF, edges);

 comparable(SetEdgeLF);
};

/*-----------------------------------------------------------------------------------------------*/
/*-----------------------------------------------------------------------------------------------*/
// Graphs
/*-----------------------------------------------------------------------------------------------*/
/*-----------------------------------------------------------------------------------------------*/

typedef boost::adjacency_list<boost::listS, boost::listS, boost::undirectedS, SetVertex, SetEdge>
 SetBasedGraph;
typedef SetBasedGraph::vertex_descriptor SetVertexDesc;
typedef SetBasedGraph::edge_descriptor SetEdgeDesc;

typedef boost::adjacency_list<boost::listS, boost::listS, boost::undirectedS, SetVertexLF, SetEdgeLF>
 SetBasedGraphLF;
typedef SetBasedGraphLF::vertex_descriptor SetVertexLFDesc;
typedef SetBasedGraphLF::edge_descriptor SetEdgeLFDesc;

#endif
