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

#ifndef VECTOR_GRAPH_DEFINITION_
#define VECTOR_GRAPH_DEFINITION_

#include <utility>
#include <set>
#include <boost/config.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/icl/discrete_interval.hpp>

#include <ast/equation.h>
#include <causalize/graph/graph_definition.h>

#include <ast/ast_types.h>
#include <boost/algorithm/string/join.hpp>
#include <boost/range/adaptor/transformed.hpp>
#include <sstream>
#include <string>

namespace ICL = boost::icl;
namespace Causalize {
  /// @brief This is the property for a vertex in the incidence graph. Nodes can be of two types: Equation or Unknow.

  struct VectorUnknown: Unknown {
    int dimension;
    std::vector<int> dimensionList;
    VectorUnknown(){};
    VectorUnknown(VarInfo varInfo, Modelica::AST::Reference var);
    void SetIndex(Modelica::AST::ExpList index);
  };

  struct VectorVertexProperty: VertexProperty {
  /// @brief The number of equations or unknowns left to causalize in this node
    int count;
    VectorUnknown unknown;
  };

  /// @brief A pair representing a usage of a variable in an equation
  typedef ICL::discrete_interval<int> Interval;
  inline Interval CreateInterval(int a, int b) {
    return ICL::discrete_interval<int>(a,b, ICL::interval_bounds::closed());
  }
  typedef std::list<Interval> IntervalList;
  typedef std::vector<Interval> IntervalVector;


  /*****************************************************************************
   ****                              Offset                                  ****
   *****************************************************************************/
  class Offset {
  public:
    inline Offset(std::vector<int> offset): offset(offset) { };
    inline Offset(): offset() { };
    inline bool operator<(const Offset& other) const { return this->offset < other.offset; };
    inline bool operator==(const Offset& other) const { return this->offset == other.offset; };
    Offset operator-() const;
  private:
      std::vector<int> offset;
  };
  /*****************************************************************************
   ****************************************************************************/



  /*****************************************************************************
   ****                               MDI                                   ****
   *****************************************************************************/
  class MDI { //Multi-Dimensional Interval

  public:
    MDI(int d, ... );
    MDI(IntervalList intervalList);
    inline MDI(IntervalVector intervals): intervals(intervals) { };
    inline int Dimension() const {return intervals.size(); }
    int Size () const;
    std::list<MDI> operator-(const MDI& other);
    std::list<MDI> Remove(const MDI& mdi, Offset offset);
    MDI ApplyUsage(std::vector<int> usage);
    bool operator<(const MDI& other) const;
    Option<MDI> operator&(const MDI& other) const;
    friend std::ostream& operator<<(std::ostream& os, const MDI mdi);

  private:
    IntervalVector intervals;
    std::list<int> usage;
    typedef IntervalVector::iterator iterator;
    typedef IntervalVector::const_iterator const_iterator;
    inline iterator begin() { return intervals.begin(); }
    inline const_iterator begin() const { return intervals.begin(); }
    inline iterator end() { return intervals.end(); }
//    MDI ApplyOffset(Offset offset);
    IntervalList Partition(Interval iA, Interval iB);
    std::list<MDI> PutHead(Interval i, std::list<MDI> mdiList);
    std::list<MDI> Filter(std::list<MDI> mdiList, MDI mdi);
    std::list<MDI> CartProd(std::list<MDI> mdiList);
    std::list<MDI> PutLists(MDI mdi, std::list<MDI> mdiList);
  };
  /*****************************************************************************
   ****************************************************************************/


  /*****************************************************************************
   ****                           INDEX PAIR                                ****
   *****************************************************************************/
  class IndexPair {
  public:
    inline IndexPair(MDI dom, MDI ran, Offset os): dom(dom), ran(ran), offset(os) { };
    inline MDI Dom() const { return dom; }
    inline MDI Ran() const { return ran; }
    inline Offset OS() const { return offset; }
    std::list<IndexPair> operator-(const IndexPair& other);
    std::list<IndexPair> RemoveUnknowns(MDI eqs);
    std::list<IndexPair> RemoveEquations(MDI eqs);
    bool operator<(const IndexPair& other) const;
    friend std::ostream& operator<<(std::ostream& os, const IndexPair& ip);
  private:
    MDI dom, ran;
    Offset offset;
    std::vector<int> usage;
  };
  /*****************************************************************************
   ****************************************************************************/

  typedef std::set<IndexPair> IndexPairSet;
  std::ostream& operator<<(std::ostream& os, const IndexPairSet& ips);

  /*****************************************************************************
   ****                              LABEL                                  ****
   *****************************************************************************/
  class Label {
  public:
    inline Label() {};
    inline Label(IndexPairSet ips): ips(ips) {};
    void RemovePairs(IndexPairSet ips);
    void RemoveUnknowns(MDI const mdi) const;
    void RemoveEquations(MDI const mdi) const;
    unsigned long int EdgeCount();
    inline bool IsEmpty() { return ips.size()==0; }
    inline const IndexPairSet & Pairs() const { return ips; }
    friend std::ostream& operator<<(std::ostream& os, const Label& label);
  private:
    IndexPairSet ips;
  };
  /*****************************************************************************
   ****************************************************************************/



  unsigned long int EdgeCount(IndexPairSet);


  /// @brief This is the definition of the Incidence graph for the vector case.
  typedef boost::adjacency_list<boost::listS, boost::listS, boost::undirectedS, VectorVertexProperty, Label> VectorCausalizationGraph;
  /// @brief This a node from the vectorized incidence graph
  typedef Causalize::VectorCausalizationGraph::vertex_descriptor VectorVertex;
  /// @brief An equation vertex is the same as a regular vertex
  typedef VectorVertex VectorEquationVertex;
  /// @brief An unknown vertex is the same as a regular vertex
  typedef VectorVertex VectorUnknownVertex;
  /// @brief This is an edge of the vectorized causalization graph
  typedef VectorCausalizationGraph::edge_descriptor VectorEdge;
  /// @brief This struct represents a set of causalized vars for the vector algorithm

  struct CausalizedVar{
    VectorUnknown unknown;
    Equation equation;
    IndexPairSet pairs;
  };
}
#endif
