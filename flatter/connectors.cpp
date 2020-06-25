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

Connectors::Connectors(MMO_Class &c) 
 : mmoclass_(c), eCount2_(0){
  SBGraph g;
  G = g;
}

member_imp(Connectors, vector<int>, vCount);
member_imp(Connectors, vector<int>, eCount1);
member_imp(Connectors, int, eCount2);
member_imp(Connectors, MMO_Class, mmoclass);

/*|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||*/
/*-----------------------------------------------------------------------------------------------*/
// Debugging functions --------------------------------------------------------------------------//
/*-----------------------------------------------------------------------------------------------*/
/*|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||*/

void Connectors::debug(std::string filename){
  GraphPrinter gp(G, -1);

  gp.printGraph(filename);
  cout << "Generated Connect Graph written to " << filename << endl;
}

/*|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||*/
/*-----------------------------------------------------------------------------------------------*/
// Create graph  --------------------------------------------------------------------------------//
/*-----------------------------------------------------------------------------------------------*/
/*|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||*/

void Connectors::solve(){
  int maxdim = 1;
  foreach_(Name n, mmoclass_.variables()){
    Option<VarInfo> ovi = mmoclass_.getVar(n);
    if(ovi){
      VarInfo vi = *ovi;
      Option<ExpList> oinds = vi.indices();
      if(oinds){
        ExpList inds = *oinds;
        NI1 aux = inds.size();
        maxdim = max(maxdim, aux);
      }
    }
  }  

  vector<NI1> aux(maxdim, 1);
  set_vCount(aux);
  set_eCount1(aux);

  createGraph(mmoclass_.equations_ref().equations_ref());

  VertexIt vi, vi_end;
  boost::tie(vi, vi_end) = boost::vertices(G);
  for(; vi != vi_end; ++vi){
    Name n = G[*vi].name;
    Set vs = G[*vi].vs_();
    cout << n << ": " << vs << "\n";
  }

  EdgeIt ei, ei_end;
  boost::tie(ei, ei_end) = boost::edges(G);
  for(; ei != ei_end; ++ei){
    Name n = G[*ei].name;
    PWLMap es1 = G[*ei].es1_();
    PWLMap es2 = G[*ei].es2_();
    cout << n << ": " << es1 << ", " << es2 << "\n";
  }

  debug("prueba.dot");

  PWLMap res = connectedComponents(G);
  cout << res << "\n";
  //generateCode();
}

void Connectors::createGraph(EquationList &eqs){
  foreach_(Equation &eq, eqs){
    if(is<Connect>(eq))
      connect(get<Connect>(eq));

    else if(is<ForEq>(eq)){
      ForEq feq = boost::get<ForEq>(eq);
      foreach_(Index ind, feq.range().indexes()){
        Name n = ind.name();
        OptExp e = ind.exp();
        if(e){
          TypePrefixes tp;
          Option<Comment> aux1;
          Expression aux2 = *e;
          ModAssign aux3(aux2);
          Modification aux4(aux3);
          Option<Modification> mod(aux4);
          ExpOptList aux5;
          VarInfo vi(tp, n, aux1, mod, aux5, false);
          mmoclass_.addVar(n, vi);
          EquationList el = feq.elements();
          createGraph(el);
          mmoclass_.rmVar(n);
        }
        else 
          cerr << "Should be defined\n";
      }
    }
  }
}

void Connectors::connect(Connect co){
  Expression eleft = co.left(), eright = co.right();
  
  Pair<Name, ExpOptList> left = separate(eleft);
  Pair<Name, ExpOptList> right = separate(eright);

  Name v1 = get<0>(left);
  Name v2 = get<0>(right);

  OrdCT<Interval> miv1 = createVertex(v1).inters_();
  OrdCT<Interval>::iterator itmiv1 = miv1.begin(); 
  OrdCT<Interval> miv2 = createVertex(v2).inters_();
  OrdCT<Interval>::iterator itmiv2 = miv2.begin();

  ExpOptList range1 = get<1>(left);
  ExpOptList range2 = get<1>(right);

  if(checkRanges(range1, range2)){
    OrdCT<Interval> mi11;
    OrdCT<Interval>::iterator itmi11 = mi11.begin();
    int dim = 0;
    const VarSymbolTable aux = mmoclass_.syms();
    EvalExpFlatter evexp(aux);
    if(range1){
      foreach_(Expression e1, range1.get()){
        if(is<SubAll>(e1)){
          Option<VarInfo> ovi = mmoclass_.getVar(v1);
          if(ovi){
            VarInfo vi = *ovi;
            ExpOptList oinds = vi.indices();
            if(oinds){
              ExpList inds = *oinds;
              ExpList::iterator itinds = inds.begin();

              for(int i = 0; i < dim; ++i)
                ++itinds;

              e1 = *itinds;
            }
          }
        }

        Interval ll = Apply(evexp, e1);
        NI1 auxlo = (*itmiv1).lo_() - 1;
        Interval l(auxlo + ll.lo_(), ll.step_(), auxlo + ll.hi_());
        ++itmiv1;
      
        if(!l.empty_()){
          itmi11 = mi11.insert(itmi11, l);
          ++dim;
        }

        else{
          OrdCT<Interval> aux;
          mi11 = aux;
          break;
        }
      }
    }

    else{
      NI1 auxlo = (*itmiv1).lo_();
      Interval l(auxlo, 1, auxlo);
      mi11.insert(mi11.begin(), l);
    }

    MultiInterval mi1(mi11);

    OrdCT<Interval> mi22;
    OrdCT<Interval>::iterator itmi22 = mi22.begin();
    dim = 0;
    if(range2){
      foreach_(Expression e2, range2.get()){
        if(is<SubAll>(e2)){
          Option<VarInfo> ovi = mmoclass_.getVar(v2);
          if(ovi){
            VarInfo vi = *ovi;
            ExpOptList oinds = vi.indices();
            if(oinds){
              ExpList inds = *oinds;
              ExpList::iterator itinds = inds.begin();

              for(int i = 0; i < dim; ++i)
                ++itinds;

              e2 = *itinds;
            }
          }
        }

        Interval rr = Apply(evexp, e2);
        NI1 auxlo = (*itmiv2).lo_() - 1;
        Interval r(auxlo + rr.lo_(), rr.step_(), auxlo + rr.hi_());
        ++itmiv2;

        if(!r.empty_()){
          itmi22 = mi22.insert(itmi22, r);
          ++dim;
        }

        else{
          OrdCT<Interval> aux;
          mi22 = aux;
          break;
        }
      }
    }

    else{
      NI1 auxlo = (*itmiv2).lo_();
      Interval r(auxlo, 1, auxlo);
      mi22.insert(mi22.begin(), r);
    }

    MultiInterval mi2(mi22);

    VertexIt vi1, vi2, vi_end1, vi_end2;
    boost::tie(vi1, vi_end1) = boost::vertices(G);
    SetVertexDesc d1 = *vi1, d2 = *vi2;

    for(; vi1 != vi_end1; ++vi1){
      boost::tie(vi2, vi_end2) = boost::vertices(G);
      Name aux1 = G[*vi1].name;

      for(; vi2 != vi_end2; ++vi2){
        Name aux2 = G[*vi2].name;
        if(aux1 == v1 && aux2 == v2){
          d1 = *vi1;
          d2 = *vi2;
          updateGraph(d1, d2, mi1, mi2);
        }
      }
    }
  }
}

// Get expression and range
Pair<Name, ExpOptList> Connectors::separate(Expression e){
  Reference reference;

  if(is<UnaryOp>(e)){
    UnaryOp u = boost::get<UnaryOp>(e);
    if(is<Reference>(u.exp()))
      reference = boost::get<Reference>(u.exp());
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
  Name r = get<0>(rf);
  return Pair<Name, ExpOptList>(r, opti);
}

MultiInterval Connectors::createVertex(Name n){
  MultiInterval mires; 

  VertexIt vi, vi_end;
  boost::tie(vi, vi_end) = boost::vertices(G);
  bool exists = false;
  for(; vi != vi_end; ++vi){
    Name aux = G[*vi].name;
    if(aux == n){
      exists = true;
      AtomSet auxas = *(G[*vi].vs_().asets_().begin()); 
      mires = auxas.aset_(); 
    }
  }
  
  if(!exists){
    Option<VarInfo> ovi = mmoclass_.getVar(n);
    if(ovi){
      VarInfo vi = *ovi;
      ExpOptList oinds = vi.indices();
      // Multi dimensional variable
      if(oinds){
        ExpList inds = *oinds;
  
        OrdCT<Interval> mi1;
        OrdCT<Interval>::iterator itmi1 = mi1.begin();
        int dim = 0;
        vector<NI1>::iterator itvc = vCount_.begin();
        vector<NI1> newvc;
        vector<NI1>::iterator itnew = newvc.begin();
   
        const VarSymbolTable auxt = mmoclass_.syms();
        EvalExpression evexp(auxt);
  
        foreach_(Expression e, inds){
          //cout << n << ": " << e << "\n";
          if(is<SubAll>(e) || is<Range>(e))
            ERROR("Ill-defined array");
   
          Real res = Apply(evexp, e);
          Interval i(*itvc, 1, *itvc + res - 1);
          //cout << i << "\n";
          if(!i.empty_()){
            itmi1 = mi1.insert(itmi1, i);
            itnew = newvc.insert(itnew, *itvc + res);
            ++itnew;

            ++dim;
            ++itvc;
          }
          else{
            OrdCT<Interval> aux;
            mi1 = aux;
            break;
          }
        }
  
        MultiInterval mi(mi1);

        if(!mi1.empty())
          set_vCount(newvc);

        AtomSet as(mi);
        Set s;
        s.addAtomSet(as);
  
        //cout << n << ": " << s << "\n";
  
        SetVertex V(n, s);
        SetVertexDesc v = boost::add_vertex(G);
  
        G[v] = V;
        mires = mi;
      }
  
      // Uni dimensional variable
      else{
        vector<NI1>::iterator itvc = vCount_.begin();
        NI1 auxvc = *itvc;
        Interval i(auxvc, 1, auxvc);
        vector<NI1> newvc;
        vector<NI1>::iterator itnew = newvc.begin();        
        itnew = newvc.insert(itnew, auxvc + 1);
        ++itnew;
        ++itvc;

        while(itvc != vCount_.end()){
          itnew = newvc.insert(itnew, *itvc);
          ++itnew;

          ++itvc;
        }
        set_vCount(newvc);

        MultiInterval mi;
        mi.addInter(i);
        AtomSet as(mi);
        Set s;
        s.addAtomSet(as);
  
        //cout << n << ": " << s << "\n";
  
        SetVertex V(n, s);
        SetVertexDesc v = boost::add_vertex(G);
 
        G[v] = V;
        mires = mi;
      }
    }
  }

  return mires;
}

// Check if only one variable is used at each subscript
bool Connectors::checkRanges(ExpOptList range1, ExpOptList range2){
  std::vector<Name> vars = mmoclass_.variables();

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

Option<SetEdgeDesc> Connectors::existsEdge(SetVertexDesc d1, SetVertexDesc d2){
  EdgeIt ei, ei_end;
  boost::tie(ei, ei_end) = boost::edges(G);

  for(; ei != ei_end; ++ei){
    SetVertexDesc v1 = boost::source(*ei, G);
    SetVertexDesc v2 = boost::target(*ei, G);
    Option<SetEdgeDesc> e(*ei);

    if(v1 == d1 && v2 == d2)
      return e;

    if(v1 == d2 && v2 == d1)
      return e;
  }

  Option<SetEdgeDesc> e;
  return e;
}

void Connectors::updateGraph(SetVertexDesc d1, SetVertexDesc d2, 
                             MultiInterval mi1, MultiInterval mi2){
  OrdCT<Interval> ints1 = mi1.inters_();
  OrdCT<Interval>::iterator itints1 = ints1.begin();
  OrdCT<Interval> ints2 = mi2.inters_();
  OrdCT<Interval>::iterator itints2 = ints2.begin();

  if(mi1.ndim_() == mi2.ndim_()){
    OrdCT<NI2> g1;
    OrdCT<NI2>::iterator itg1 = g1.begin();
    OrdCT<NI2> o1;
    OrdCT<NI2>::iterator ito1 = o1.begin();
    OrdCT<NI2> g2;
    OrdCT<NI2>::iterator itg2 = g2.begin();
    OrdCT<NI2> o2;
    OrdCT<NI2>::iterator ito2 = o2.begin();
    OrdCT<Interval> mi;
    OrdCT<Interval>::iterator itmi = mi.begin(); 

    vector<NI1>::iterator itec = eCount1_.begin(); 
    vector<NI1> newec;
    vector<NI1>::iterator itnew = newec.begin();

    while(itints1 != ints1.end()){
      NI1 sz1 = (*itints1).size();
      NI1 sz2 = (*itints2).size();

      if(sz1 == sz2 || sz1 == 1 || sz2 == 1){
        NI1 count = max(sz1, sz2);
        NI1 auxec = *itec;
        Interval i(auxec, 1, auxec + count - 1);
        itmi = mi.insert(itmi, i);
        ++itmi;

        NI2 g1i = (*itints1).step_();
        NI2 o1i = (-g1i) * auxec + (*itints1).lo_();

        NI2 g2i = (*itints2).step_();
        NI2 o2i = (-g2i) * auxec + (*itints2).lo_();

        if(sz1 == 1){
          itg1 = g1.insert(itg1, 0);
          ito1 = o1.insert(ito1, (*itints1).lo_());
        }

        else{
          itg1 = g1.insert(itg1, g1i);
          ito1 = o1.insert(ito1, o1i);
        }

        if(sz2 == 1){
          itg2 = g2.insert(itg2, 0);
          ito2 = o2.insert(ito2, (*itints2).lo_());
        }

        else{
          itg2 = g2.insert(itg2, g2i);
          ito2 = o2.insert(ito2, o2i);
        }
 
        itnew = newec.insert(itnew, auxec + count);
        ++itnew;

        ++itg1;
        ++ito1;
        ++itg2;
        ++ito2;
        ++itec;
      }
 
      else{ 
        cerr << "Incompatible connect\n"; 
        newec = eCount1_;
        break;
      }

      ++itints1;
      ++itints2;
    }

    set_eCount1(newec);

    AtomSet as(mi);
    Set s;
    s.addAtomSet(as); 

    LMap lm1(g1, o1);
    LMap lm2(g2, o2);
 
    OrdCT<Set> cts1;
    cts1.insert(cts1.end(), s);
    OrdCT<LMap> ctlm1;
    ctlm1.insert(ctlm1.end(), lm1); 
    PWLMap e1(cts1, ctlm1);

    OrdCT<Set> cts2;
    cts2.insert(cts2.end(), s);
    OrdCT<LMap> ctlm2;
    ctlm2.insert(ctlm2.end(), lm2); 
    PWLMap e2(cts2, ctlm2);

    //cout << e1 << ", " << e2 << "\n";

    Option<SetEdgeDesc> oedge = existsEdge(d1, d2);    

    if(oedge){
      SetEdgeDesc e = *oedge;
      SetEdge aux = G[e];

      PWLMap pwaux1 = aux.es1_();
      pwaux1.addLMSet(lm1, s);
      PWLMap pwaux2 = aux.es2_();
      pwaux2.addLMSet(lm2, s);

      SetEdge E(aux.name, pwaux1, pwaux2);
      G[e] = E;
    }

    else{
      string enm = "E" + to_string(eCount2_);
      SetEdge E(enm, e1, e2);
      SetEdgeDesc e;
      bool b;
      boost::tie(e, b) = boost::add_edge(d1, d2, G);
      G[e] = E;
      ++eCount2_;
    }
  }

  else
    cerr << "Incompatible connect\n";
}
