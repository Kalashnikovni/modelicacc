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
#include <list>
#include <iostream>
#include <fstream>
#include <sstream>

#include <boost/graph/adjacency_list.hpp>
#include <boost/icl/interval_set.hpp>
#include <boost/lexical_cast.hpp>

#include<ast/ast_types.h>
#include<util/graph/graph_definition.h>

using namespace std;
using namespace boost::icl;

#define MAKE_SPACE for(int __i=0; __i<depth; __i++) stri << " ";
#define TAB_SPACE 2
#define INSERT_TAB depth += TAB_SPACE;
#define DELETE_TAB depth -= TAB_SPACE;

namespace Graph{
 class GraphPrinter{
  public:
  GraphPrinter(const SetBasedGraph &g, const int mod) 
   : graph(g), mode(mod){};

  void printGraph(std::string name){
   stringstream stri;
   ofstream out(name.c_str());
   int depth = 0;
      
   stri << "digraph G{" << endl;
   INSERT_TAB
   MAKE_SPACE
   stri << "  ratio=\"fill\"" << endl;
   MAKE_SPACE
   stri << "  node[shape=\"ellipse\"]" << endl;
   INSERT_TAB
   MAKE_SPACE
 
   // Vertices printing
   printVertices(stri);
  
   DELETE_TAB
   DELETE_TAB
      
   INSERT_TAB
   INSERT_TAB
   stringstream colors;

   DELETE_TAB
   MAKE_SPACE
   stri << colors.str();
   DELETE_TAB
      
   INSERT_TAB
   MAKE_SPACE

   //Edge printing
   printEdges(stri);

   DELETE_TAB
   stri << "}" << endl;
   out << stri.str();
   out.close();
  };

  private:
  const SetBasedGraph &graph;
  const int mode;

  void printVertices(stringstream &stri){
   switch(mode){
    case 1: // Flatter
     break;
    case 2: // SetBased
     break;
    case 3: // Compacted
     break;
    default:
     typename boost::graph_traits<SetBasedGraph>::vertex_iterator vi, vi_end;
     for(boost::tie(vi, vi_end) = boost::vertices(graph); vi != vi_end; ++vi){
      stri << vPrinter(graph[*vi]) << " [label=\"" << vPrinter(graph[*vi]) << "\"]";
     }
   }
  };

  void printEdges(stringstream &stri){
   switch(mode){
    case 1: // Flatter
     break;
    case 2: // SetBased
     break;
    case 3: // Compacted
     break;
    default:
     typename boost::graph_traits<SetBasedGraph>::edge_iterator ei, ei_end;
     for(boost::tie(ei, ei_end) = boost::edges(graph); ei != ei_end; ++ei){
      SetVertexDesc v1 = boost::source(*ei, graph);
      SetVertexDesc v2 = boost::target(*ei, graph);
      stri << vPrinter(graph[v1]) << " -> " << vPrinter(graph[v2]); 
     }
   }
  };

  std::string vPrinter(SetVertex v){
   return v.name();
  };
 };
} // namespace Graph
