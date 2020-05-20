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

#include <boost/test/unit_test.hpp>
#include <boost/test/included/unit_test.hpp>
#include <boost/unordered_set.hpp>

#include <ast/expression.h>
#include <util/graph/graph_definition.h>

using namespace boost::unit_test;

/// @brief If sb-graphs containers implementation changes (uses vector, for example)
/// this typedef should also change.
typedef UnordCT<Interval> contInt1;
typedef OrdCT<Interval> contInt2;
typedef UnordCT<MultiInterval> contMulti;
typedef UnordCT<AtomSet> contAS;

typedef OrdCT<NI1> contNI1;
typedef OrdCT<NI2> contNI2;

//____________________________________________________________________________//

// -- Intervals --------------------------------------------------------------//

void TestIntCreation1(){
  Interval i(10, 3, 3, false);

  BOOST_CHECK(i.empty_() == true);
}

void TestIntCreation2(){
  Interval i(10, 20, 15, false);

  BOOST_CHECK(i.hi_() == 10);
}

void TestIntCreation3(){
  Interval i(10, 5, 23, false);

  BOOST_CHECK(i.hi_() == 20);
}

void TestIntCreation4(){
  Interval i(10, 1, Inf, false);

  BOOST_CHECK(i.hi_() == Inf);
}

void TestIntQuery1(){
  Interval i(10, 2, 20, false);

  BOOST_CHECK(!i.isIn(13));
}

void TestIntQuery2(){
  Interval i(10, 2, 20, false);

  BOOST_CHECK(i.isIn(18));
}

void TestIntQuery3(){
  Interval i(10, 2, 20, false);

  BOOST_CHECK(!i.isIn(100));
}

void TestIntQuery4(){
  Interval i(10, 2, 20, true);

  BOOST_CHECK(!i.isIn(16));
}

void TestIntQuery5(){
  Interval i1(10, 2, 20, false);
  Interval i2(0, 3, 25, false);

  bool b1 = i1.isIn(12);
  bool b2 = i2.isIn(12);

  BOOST_CHECK(b1 && b2);
}

void TestIntQuery6(){
  Interval i(true);

  BOOST_CHECK(!i.isIn(10));
}

// Cap should be commutative
void TestIntCap1(){
  Interval i1(10, 2, 20, false);
  Interval i2(0, 3, 25, false);

  Interval i3(i1.cap(i2));
  Interval i4(i2.cap(i1));

  BOOST_CHECK(i3 == i4); 
}

void TestIntCap2(){
  Interval i1(10, 2, 20, false);
  Interval i2(0, 3, 25, false);

  Interval i3 = i1.cap(i2);

  Interval i4(12, 6, 18, false);

  BOOST_CHECK(i3 == i4);
}

void TestIntCap3(){
  Interval i1(14, 2, 16, false);
  Interval i2(12, 3, 15, false);

  Interval i3 = i1.cap(i2);

  Interval i4(true);

  BOOST_CHECK(i3 == i4);
}

void TestIntCap4(){
  Interval i1(14, 2, 28, false);
  Interval i2(0, 1, Inf, false);

  Interval i3 = i1.cap(i2);
 
  Interval i4(14, 2, 28, false);

  IntervalImp tryDen;

  BOOST_CHECK(i3 == i4);
}

void TestIntDiff1(){
  Interval i1(0, 2, 30, false);
  Interval i2(true);
  
  contInt1 res1 = i1.diff(i2);

  contInt1 res2;
  res2.insert(i1);

  BOOST_CHECK(res1 == res2);
}

void TestIntDiff2(){
  Interval i1(0, 2, 30, false);
  Interval i2(10, 3, 40, false);

  Interval i3(0, 2, 8, false);
  Interval i4(12, 6, 24, false);
  Interval i5(14, 6, 26, false);
  Interval i6(30, 2, 30, false);

  contInt1 res2; 
  res2.insert(i3);
  res2.insert(i4);
  res2.insert(i5);
  res2.insert(i6);

  BOOST_CHECK(i1.diff(i2) == res2);
}

void TestIntDiff3(){
  Interval i1(0, 2, Inf, false);
  Interval i2(10, 3, 40, false);

  Interval i3(0, 2, 8, false);
  Interval i4(12, 6, 36, false);
  Interval i5(14, 6, 38, false);
  Interval i6(42, 2, Inf, false);

  contInt1 res2; 
  res2.insert(i3);
  res2.insert(i4);
  res2.insert(i5);
  res2.insert(i6);

  BOOST_CHECK(i1.diff(i2) == res2);
}

void TestIntDiff4(){
  Interval i1(0, 1, 10, false);
  Interval i2(true);

  contInt1 res2;
  res2.insert(i2);

  BOOST_CHECK(i1.diff(i1) == res2);
}

void TestIntMin1(){
  Interval i(10, 3, 40, false);

  NI1 res1 = i.minElem();

  BOOST_CHECK(res1 == 10);
}

// -- MultiIntervals --------------------------------------------------------------//
void TestMultiCreation1(){
  Interval i1(true);
  Interval i2(0, 2, 50, false);
  Interval i3(3, 1, 5, false);
  Interval i4(3, 8, 24, false);

  MultiInterval mi1;
  
  mi1.addInter(i1);
  mi1.addInter(i2);
  mi1.addInter(i3);
  mi1.addInter(i4);

  contInt2 res;
  contInt2::iterator it = res.begin();

  it = res.insert(it, i1);
  ++it;
  it = res.insert(it, i2);
  ++it;
  it = res.insert(it, i3);
  ++it;
  it = res.insert(it, i4);

  MultiInterval mi2(res);

  BOOST_CHECK(mi1 == mi2);
}

void TestMultiEmpty1(){
  MultiInterval mi;

  BOOST_CHECK(mi.empty());
}

void TestMultiEmpty2(){
  Interval i1(true);
  Interval i2(true);
  Interval i3(true);

  MultiInterval mi;

  mi.addInter(i1);
  mi.addInter(i2);
  mi.addInter(i3);

  BOOST_CHECK(mi.empty());
}

void TestMultiEmpty3(){
  Interval i1(true);
  Interval i2(0, 1, 10, false);
  Interval i3(true);

  MultiInterval mi;

  mi.addInter(i1);
  mi.addInter(i2);
  mi.addInter(i3);

  BOOST_CHECK(!mi.empty());
}

void TestMultiQuery1(){
  Interval i1(1, 1, 10, false);
  Interval i2(true);
  Interval i3(10, 2, 21, false);

  MultiInterval mi;
  mi.addInter(i1);
  mi.addInter(i2);
  mi.addInter(i3);

  contNI1 elem1;
  contNI1::iterator it1 = elem1.begin();
  it1 = elem1.insert(it1, 5);
  ++it1;
  it1 = elem1.insert(it1, 10);
  ++it1;
  elem1.insert(it1, 21);

  BOOST_CHECK(!mi.isIn(elem1));
}

void TestMultiQuery2(){
  Interval i1(1, 1, 10, false);
  Interval i2(10, 20, 10, false);
  Interval i3(10, 2, 21, false);

  MultiInterval mi;
  mi.addInter(i1);
  mi.addInter(i2);
  mi.addInter(i3);

  contNI1 elem1;
  contNI1::iterator it1 = elem1.begin();
  it1 = elem1.insert(it1, 5);
  ++it1;
  it1 = elem1.insert(it1, 10);
  ++it1;
  elem1.insert(it1, 21);

  BOOST_CHECK(!mi.isIn(elem1));
}

void TestMultiQuery3(){
  Interval i1(1, 1, 10, false);
  Interval i2(10, 20, 10, false);
  Interval i3(10, 2, 21, false);

  MultiInterval mi;
  mi.addInter(i1);
  mi.addInter(i2);
  mi.addInter(i3);

  contNI1 elem1;
  contNI1::iterator it1 = elem1.begin();
  it1 = elem1.insert(it1, 5);
  ++it1;
  it1 = elem1.insert(it1, 10);
  ++it1;
  elem1.insert(it1, 20);

  BOOST_CHECK(mi.isIn(elem1));
}

void TestMultiAddInter1(){
  Interval i1(0, 2, 10, false);
  MultiInterval mi1;

  mi1.addInter(i1);  

  contInt2 ints2;
  contInt2::iterator it2 = ints2.begin();
  ints2.insert(it2, i1);
  MultiInterval mi2(ints2); 

  BOOST_CHECK(mi1 == mi2);
}

void TestMultiAddInter2(){
  Interval i1(0, 2, 10, false);
  Interval i2(30, 2, 40, false);
  Interval i3(25, 1, 30, false);
  MultiInterval mi1;

  mi1.addInter(i1);  
  mi1.addInter(i2);
  mi1.addInter(i3);

  contInt2 ints2;
  contInt2::iterator it2 = ints2.end();
  it2 = ints2.insert(it2, i1);
  ++it2;
  it2 = ints2.insert(it2, i2);
  ++it2;
  ints2.insert(it2, i3);
  MultiInterval mi2(ints2); 

  BOOST_CHECK(mi1 == mi2);
}

void TestMultiCap1(){
  Interval i1(0, 2, 20, false);
  Interval i2(30, 2, 40, false);
  Interval i3(25, 1, 30, false);
  MultiInterval mi1;

  mi1.addInter(i1);
  mi1.addInter(i2);
  mi1.addInter(i3);

  MultiInterval mi2;

  MultiInterval mi3 = mi1.cap(mi2);
  MultiInterval mi4 = mi2.cap(mi1);

  BOOST_CHECK(mi2 == mi3 && mi3 == mi4);
}

void TestMultiCap2(){
  Interval i1(0, 2, 20, false);
  Interval i2(30, 2, 40, false);
  Interval i3(25, 1, 30, false);

  MultiInterval mi1;
  mi1.addInter(i1);
  mi1.addInter(i2);
  mi1.addInter(i3);

  Interval i4(5, 3, 15, false);
  Interval i5(true);
  Interval i6(27, 1, 35, false);

  MultiInterval mi2;
  mi2.addInter(i4);
  mi2.addInter(i5);
  mi2.addInter(i6);

  contInt2 aux;

  MultiInterval res1(aux);

  MultiInterval res2 = mi1.cap(mi2);
  MultiInterval res3 = mi2.cap(mi1);
  
  BOOST_CHECK(res1 == res2 && res2 == res3);
}

void TestMultiCap3(){
  Interval i1(1, 1, 10, false);

  MultiInterval mi1;
  mi1.addInter(i1);
  mi1.addInter(i1);
  mi1.addInter(i1);

  Interval i2(30, 1, 40, false);

  MultiInterval mi2;
  mi2.addInter(i1);
  mi2.addInter(i1);
  mi2.addInter(i2);

  MultiInterval mi3 = mi1.cap(mi2);

  BOOST_CHECK(mi3.empty());
}

void TestMultiDiff1(){
  Interval i1(0, 2, 20, false);
  Interval i2(30, 2, 40, false);
  Interval i3(25, 1, 30, false);

  MultiInterval mi1;
  mi1.addInter(i1);
  mi1.addInter(i2);
  mi1.addInter(i3);

  Interval i4(5, 3, 15, false);
  Interval i5(30, 2, 30, false);
  Interval i6(27, 1, 35, false);

  MultiInterval mi2;
  mi2.addInter(i4);
  mi2.addInter(i5);
  mi2.addInter(i6);

  contMulti res1 = mi1.diff(mi2);
  contMulti::iterator it1 = res1.begin();

  Interval i7(0, 2, 6, false);
 
  MultiInterval mi3;
  mi3.addInter(i7);
  mi3.addInter(i2);
  mi3.addInter(i3);

  Interval i8(10, 6, 10, false);

  MultiInterval mi4;
  mi4.addInter(i8);
  mi4.addInter(i2);
  mi4.addInter(i3);

  Interval i9(12, 6, 12, false);

  MultiInterval mi5;
  mi5.addInter(i9);
  mi5.addInter(i2);
  mi5.addInter(i3);

  Interval i10(16, 2, 20, false);

  MultiInterval mi6;
  mi6.addInter(i10);
  mi6.addInter(i2);
  mi6.addInter(i3);

  Interval i11(8, 6, 14, false);
  Interval i12(32, 2, 40, false);

  MultiInterval mi7;
  mi7.addInter(i11);
  mi7.addInter(i12);
  mi7.addInter(i3);

  Interval i13(30, 2, 30, false);
  Interval i14(25, 1, 26, false);

  MultiInterval mi8;
  mi8.addInter(i11);
  mi8.addInter(i13);
  mi8.addInter(i14);

  contMulti res2;
  res2.insert(mi3);
  res2.insert(mi4);
  res2.insert(mi5);
  res2.insert(mi6);
  res2.insert(mi7);
  res2.insert(mi8);

  BOOST_CHECK(res1 == res2);
}

void TestMultiDiff2(){
  Interval i1(0, 2, 20, false);
  Interval i2(30, 2, 40, false);
  Interval i3(25, 1, 30, false);

  MultiInterval mi1;
  mi1.addInter(i1);
  mi1.addInter(i2);
  mi1.addInter(i3);

  Interval i4(5, 3, 15, false);
  Interval i5(30, 2, 30, false);
  Interval i6(25, 1, 35, false);

  MultiInterval mi2;
  mi2.addInter(i4);
  mi2.addInter(i5);
  mi2.addInter(i6);

  contMulti res1 = mi1.diff(mi2);
  contMulti::iterator it1 = res1.begin();

  Interval i7(0, 2, 6, false);
 
  MultiInterval mi3;
  mi3.addInter(i7);
  mi3.addInter(i2);
  mi3.addInter(i3);

  Interval i8(10, 6, 10, false);

  MultiInterval mi4;
  mi4.addInter(i8);
  mi4.addInter(i2);
  mi4.addInter(i3);

  Interval i9(12, 6, 12, false);

  MultiInterval mi5;
  mi5.addInter(i9);
  mi5.addInter(i2);
  mi5.addInter(i3);

  Interval i10(16, 2, 20, false);

  MultiInterval mi6;
  mi6.addInter(i10);
  mi6.addInter(i2);
  mi6.addInter(i3);

  Interval i11(8, 6, 14, false);
  Interval i12(32, 2, 40, false);

  MultiInterval mi7;
  mi7.addInter(i11);
  mi7.addInter(i12);
  mi7.addInter(i3);

  contMulti res2;
  res2.insert(mi3);
  res2.insert(mi4);
  res2.insert(mi5);
  res2.insert(mi6);
  res2.insert(mi7);

  BOOST_CHECK(res1 == res2);
}

void TestMultiDiff3(){
  Interval i1(true);
  Interval i2(30, 2, 40, false);
  Interval i3(25, 1, 30, false);
  
  MultiInterval mi1;
  mi1.addInter(i1);
  mi1.addInter(i2);
  mi1.addInter(i3);

  Interval i4(5, 3, 14, false);
  Interval i5(true);
  Interval i6(true);

  MultiInterval mi2;
  mi2.addInter(i4);
  mi2.addInter(i5);
  mi2.addInter(i6);

  Interval i7(true);
  Interval i8(32, 6, 38, false);
  Interval i9(true);

  contMulti res1 = mi1.diff(mi2);
  contMulti::iterator it1 = res1.begin();

  BOOST_CHECK(true);
}

void TestMultiDiff4(){
  Interval i1(1, 1, 10, false);
  Interval i2(20, 3, 33, false);

  MultiInterval mi1;
  mi1.addInter(i1);
  mi1.addInter(i2);

  contMulti res1 = mi1.diff(mi1);

  BOOST_CHECK(res1.empty());
}

void TestMultiDiff5(){
  Interval i1(1, 1, 10, false);
  Interval i2(2, 2, 20, false);

  MultiInterval mi1;
  mi1.addInter(i1);
  mi1.addInter(i2);

  Interval i3(1, 1, 15, false);
  Interval i4(2, 2, 40, false);

  MultiInterval mi2;
  mi2.addInter(i3);
  mi2.addInter(i4);

  contMulti res1 = mi1.diff(mi2);

  BOOST_CHECK(res1.empty());
}

void TestMultiCrossProd1(){
  Interval i1(1, 1, 10, false);
  Interval i2(2, 2, 40, false);

  MultiInterval mi1;
  mi1.addInter(i1);
  mi1.addInter(i2);

  Interval i3(3, 3, 20, false);
  Interval i4(1, 50, Inf, false);

  MultiInterval mi2;
  mi2.addInter(i3);
  mi2.addInter(i4);

  MultiInterval res1 = mi1.crossProd(mi2);

  MultiInterval res2;
  res2.addInter(i1);
  res2.addInter(i2);
  res2.addInter(i3);
  res2.addInter(i4);

  BOOST_CHECK(res1 == res2);
}

void TestMultiMin1(){
  Interval i1(0, 1, 40, false);
  Interval i2(15, 3, 18, false);
  Interval i3(50, 2, 70, false);

  MultiInterval mi;

  mi.addInter(i1);
  mi.addInter(i2);
  mi.addInter(i3);

  contNI1 res1 = mi.minElem();

  contNI1 res2;
  contNI1::iterator it2 = res2.begin();
 
  it2 = res2.insert(it2, 0);
  ++it2;
  it2 = res2.insert(it2, 15);
  ++it2;
  res2.insert(it2, 50);

  BOOST_CHECK(res1 == res2);
}

// -- AtomicSets --------------------------------------------------------------//

void TestASetCreation1(){
  Interval i1(true);
  Interval i2(0, 2, 50, false);
  Interval i3(3, 1, 5, false);
  Interval i4(3, 8, 24, false);

  MultiInterval mi;

  mi.addInter(i1);
  mi.addInter(i2);
  mi.addInter(i3);
  mi.addInter(i4);

  AtomSet as(mi);

  BOOST_CHECK(mi == as.aset_());
}

void TestASetEmpty1(){
  MultiInterval mi;

  BOOST_CHECK(mi.empty());
}

void TestASetEmpty2(){
  Interval i1(true);
  Interval i2(0, 2, 50, false);
  Interval i3(3, 1, 5, false);
  Interval i4(3, 8, 24, false);

  MultiInterval mi;

  mi.addInter(i1);
  mi.addInter(i2);
  mi.addInter(i3);
  mi.addInter(i4);

  AtomSet as(mi);

  BOOST_CHECK(!as.empty());
}

void TestASetEmpty3(){
  Interval i1(true);
  Interval i2(true);
  Interval i3(true);

  MultiInterval mi;
  mi.addInter(i1);
  mi.addInter(i2);
  mi.addInter(i3);

  AtomSet as(mi);

  BOOST_CHECK(mi.empty());
}

void TestASetEmpty4(){
  Interval i1(true);
  Interval i2(true);
  Interval i3(1, 1, 10, false);

  MultiInterval mi;
  mi.addInter(i1);
  mi.addInter(i2);
  mi.addInter(i3);

  AtomSet as(mi);

  BOOST_CHECK(!mi.empty());
}

void TestASetCap1(){
  Interval i1(0, 2, 20, false);
  Interval i2(30, 2, 40, false);
  Interval i3(25, 1, 30, false);

  MultiInterval mi1;
  mi1.addInter(i1);
  mi1.addInter(i2);
  mi1.addInter(i3);

  AtomSet as1(mi1);

  Interval i4(5, 3, 15, false);
  Interval i5(true);
  Interval i6(27, 1, 35, false);

  MultiInterval mi2;
  mi2.addInter(i4);
  mi2.addInter(i5);
  mi2.addInter(i6);

  AtomSet as2(mi2);

  AtomSet res1 = as1.cap(as2);
  AtomSet res2 = as2.cap(as1);

  AtomSet res3;

  BOOST_CHECK(res1 == res2 && res2 == res3);
}

void TestASetDiff1(){
  Interval i1(0, 2, 20, false);
  Interval i2(30, 2, 40, false);
  Interval i3(25, 1, 30, false);

  MultiInterval mi1;
  mi1.addInter(i1);
  mi1.addInter(i2);
  mi1.addInter(i3);

  AtomSet as1(mi1);

  Interval i4(5, 3, 15, false);
  Interval i5(30, 2, 30, false);
  Interval i6(27, 1, 35, false);

  MultiInterval mi2;
  mi2.addInter(i4);
  mi2.addInter(i5);
  mi2.addInter(i6);

  AtomSet as2(mi2);

  contAS res1 = as1.diff(as2);

  Interval i7(0, 2, 6, false);

  MultiInterval mi3;
  mi3.addInter(i7);
  mi3.addInter(i2);
  mi3.addInter(i3);

  AtomSet as3(mi3);

  Interval i8(10, 6, 10, false);

  MultiInterval mi4;
  mi4.addInter(i8);
  mi4.addInter(i2);
  mi4.addInter(i3);

  AtomSet as4(mi4);

  Interval i9(12, 6, 12, false);

  MultiInterval mi5;
  mi5.addInter(i9);
  mi5.addInter(i2);
  mi5.addInter(i3);

  AtomSet as5(mi5);

  Interval i10(16, 2, 20, false);

  MultiInterval mi6;
  mi6.addInter(i10);
  mi6.addInter(i2);
  mi6.addInter(i3);

  AtomSet as6(mi6);

  Interval i11(8, 6, 14, false);
  Interval i12(32, 2, 40, false);

  MultiInterval mi7;
  mi7.addInter(i11);
  mi7.addInter(i12);
  mi7.addInter(i3);

  AtomSet as7(mi7);

  Interval i13(30, 2, 30, false);
  Interval i14(25, 1, 26, false);

  MultiInterval mi8;
  mi8.addInter(i11);
  mi8.addInter(i13);
  mi8.addInter(i14);

  AtomSet as8(mi8);

  contAS res2;
  res2.insert(as3);
  res2.insert(as4);
  res2.insert(as5);
  res2.insert(as6);
  res2.insert(as7);
  res2.insert(as8);

  BOOST_CHECK(res1 == res2);
}

void TestASetMin1(){
  Interval i1(0, 1, 40, false);
  Interval i2(15, 3, 18, false);
  Interval i3(50, 2, 70, false);

  MultiInterval mi;

  mi.addInter(i1);
  mi.addInter(i2);
  mi.addInter(i3);

  AtomSet as1(mi);

  contNI1 res1 = mi.minElem();

  contNI1 res2;
  contNI1::iterator it2 = res2.begin();
 
  it2 = res2.insert(it2, 0);
  ++it2;
  it2 = res2.insert(it2, 15);
  ++it2;
  res2.insert(it2, 50);

  BOOST_CHECK(res1 == res2);
}

// -- Sets -------------------------------------------------------------------//

void TestSetCreation1(){
  Interval i1(0, 2, 20, false);
  Interval i2(30, 2, 40, false);
  Interval i3(25, 1, 30, false);

  MultiInterval mi1;

  mi1.addInter(i1);
  mi1.addInter(i2);
  mi1.addInter(i3);

  AtomSet as1(mi1);

  Interval i4(0, 1, 10, false);

  MultiInterval mi2;

  mi2.addInter(i4);
  mi2.addInter(i4);
  mi2.addInter(i4);

  AtomSet as2(mi2);

  Set s1;
  s1.addAtomSet(as1);
  s1.addAtomSet(as2);

  contAS res2;
  contAS::iterator it2 = res2.begin();
  it2 = res2.insert(it2, as1);
  ++it2;
  res2.insert(it2, as2);

  Set s2(res2);

  BOOST_CHECK(s1 == s2);
}

void TestCompSets1(){
  Interval i1(0, 1, 10, false);
  
  MultiInterval mi1;
  mi1.addInter(i1);

  AtomSet as1(mi1);

  Set s1;
  s1.addAtomSet(as1);

  Interval i2(0, 1, 20, false);

  MultiInterval mi2;
  mi2.addInter(i2);

  AtomSet as2(mi2);

  Set s2;
  s2.addAtomSet(as2);
  
  BOOST_CHECK(!(s1 == s2));
}

void TestSetEmpty1(){
  Interval i7(0, 1, Inf, false);
  Interval i8(20, 3, 50, false);
  Interval i9(true);

  MultiInterval mi3;

  mi3.addInter(i7);
  mi3.addInter(i8);
  mi3.addInter(i9);

  AtomSet as3(mi3);

  Set s2;

  s2.addAtomSet(as3);

  BOOST_CHECK(!s2.empty());
}

void TestAddASets1(){
  Interval i1(0, 2, 20, false);
  Interval i2(30, 2, 40, false);
  Interval i3(25, 1, 30, false);
  MultiInterval mi1;

  mi1.addInter(i1);
  mi1.addInter(i2);
  mi1.addInter(i3);

  Interval i4(5, 3, 15, false);
  Interval i5(true);
  Interval i6(27, 1, 35, false);
  MultiInterval mi2;

  mi2.addInter(i4);
  mi2.addInter(i5);
  mi2.addInter(i6);

  AtomSet as1(mi1);
  AtomSet as2(mi2);

  Set s1;

  s1.addAtomSet(as1);
  s1.addAtomSet(as2);

  contAS aux;
  contAS::iterator itaux = aux.begin();
  itaux = aux.insert(itaux, as1);
  ++itaux;
  aux.insert(itaux, as2); 

  Set s2(aux);

  BOOST_CHECK(s1 == s2);
}

void TestSetCap1(){
  Set s1;
  Set s2;

  Set res1 = s1.cap(s2);
  Set res2 = s2.cap(s1);

  BOOST_CHECK(res1 == res2 && res1.empty() && res2.empty());
}

void TestSetCap2(){
  Set s1;

  Interval i1(0, 2, 20, false);
  Interval i2(30, 2, 40, false);
  Interval i3(25, 1, 30, false);

  MultiInterval mi1;
  mi1.addInter(i1);
  mi1.addInter(i2);
  mi1.addInter(i3);

  AtomSet as1(mi1);

  Set s2;

  s2.addAtomSet(as1);

  Set res1 = s1.cap(s2);
  Set res2 = s2.cap(s1);

  BOOST_CHECK(res1 == res2 && res1.empty() && res2.empty());
}

void TestSetCap3(){
  Interval i1(0, 2, 20, false);
  Interval i2(30, 2, 40, false);
  Interval i3(25, 1, 30, false);

  MultiInterval mi1;
  mi1.addInter(i1);
  mi1.addInter(i2);
  mi1.addInter(i3);

  Interval i4(5, 3, 15, false);
  Interval i5(35, 3, 40, false);
  Interval i6(27, 1, 35, false);

  MultiInterval mi2;
  mi2.addInter(i4);
  mi2.addInter(i5);
  mi2.addInter(i6);

  AtomSet as1(mi1);
  AtomSet as2(mi2);

  Set s1;
  s1.addAtomSet(as1);
  s1.addAtomSet(as2);

  Interval i7(0, 1, Inf, false);
  Interval i8(20, 3, 50, false);
  Interval i9(28, 1, 28, false);

  MultiInterval mi3;
  mi3.addInter(i7);
  mi3.addInter(i8);
  mi3.addInter(i9);

  AtomSet as3(mi3);

  Set s2;
  s2.addAtomSet(as3);

  Set res1 = s1.cap(s2);
  Set res2 = s2.cap(s1);

  Interval i10(0, 2, 20, false);
  Interval i11(32, 6, 38, false);

  MultiInterval mi4;
  mi4.addInter(i10);
  mi4.addInter(i11);
  mi4.addInter(i9);

  MultiInterval mi5;
  mi5.addInter(i4);
  mi5.addInter(i5);
  mi5.addInter(i9);

  AtomSet as4(mi4);
  AtomSet as5(mi5);

  Set res3;

  res3.addAtomSet(as4);
  res3.addAtomSet(as5);

  BOOST_CHECK(res1 == res2 && res2 == res3);
}

void TestSetDiff1(){
  Interval i1(0, 1, 10, false);
  Interval i2(0, 3, 9, false);

  MultiInterval mi1;
  mi1.addInter(i1);
  mi1.addInter(i2);

  AtomSet as1(mi1);

  Set s1;
  s1.addAtomSet(as1);

  Interval i3(0, 1, 10, false);
  Interval i4(0, 3, 9, false);

  MultiInterval mi2;
  mi2.addInter(i3);
  mi2.addInter(i4);

  AtomSet as2(mi2);

  Set s2;
  s2.addAtomSet(as2);
 
  Set res1 = s1.diff(s2); 
  Set res2;

  BOOST_CHECK(res1 == res2);
}

void TestSetMin1(){
  Interval i1(true);
  Interval i2(5, 1, 10, false);

  MultiInterval mi1;
  mi1.addInter(i1);
  mi1.addInter(i2);

  AtomSet as1(mi1);

  Interval i3(20, 20, 80, false);
  Interval i4(1, 1, 500, false);

  MultiInterval mi2;
  mi2.addInter(i3);
  mi2.addInter(i4);

  AtomSet as2(mi2);

  Interval i5(30, 5, 36, false);
  Interval i6(42, 3, 57, false);

  MultiInterval mi3;
  mi3.addInter(i5);
  mi3.addInter(i6);  

  AtomSet as3(mi3);

  Set s;
  s.addAtomSet(as1);
  s.addAtomSet(as2);
  s.addAtomSet(as3);

  contNI1 res1 = s.minElem();

  contNI1 res2;
  res2.insert(res2.end(), 20);
  res2.insert(res2.end(), 1);

  BOOST_CHECK(res1 == res2);   
}

// -- LinearMaps -------------------------------------------------------------------//

void TestLMCreation1(){
  contNI2 g1;
  contNI2::iterator itg1 = g1.begin();

  itg1 = g1.insert(itg1, 2);
  ++itg1;
  itg1 = g1.insert(itg1, 3);
  ++itg1;
  itg1 = g1.insert(itg1, 4);

  contNI2 o1;
  contNI2::iterator ito1 = o1.begin();

  ito1 = o1.insert(ito1, 7);
  ++ito1;
  o1.insert(ito1, 8);

  LMap res(g1, o1);

  BOOST_CHECK(res.empty());
}

void TestLMCompose1(){
  contNI2 g1;
  contNI2::iterator itg1 = g1.begin();
  
  itg1 = g1.insert(itg1, 5);
  ++itg1;
  itg1 = g1.insert(itg1, 10);
  ++itg1;
  itg1 = g1.insert(itg1, 3);

  contNI2 o1;
  contNI2::iterator ito1 = o1.begin();

  ito1 = o1.insert(ito1, 1);
  ++ito1;
  ito1 = o1.insert(ito1, 2);
  ++ito1;
  ito1 = o1.insert(ito1, 3);

  LMap lm1(g1, o1);

  contNI2 g2;
  contNI2::iterator itg2 = g2.begin();

  itg2 = g2.insert(itg2, 2);
  ++itg2;
  itg2 = g2.insert(itg2, 2);
  ++itg2;
  itg2 = g2.insert(itg2, 2);

  contNI2 o2;
  contNI2::iterator ito2 = o2.begin();

  ito2 = o2.insert(ito2, 3);
  ++ito2;
  ito2 = o2.insert(ito2, 3);
  ++ito2;
  ito2 = o2.insert(ito2, 3);

  LMap lm2(g2, o2);

  LMap res1 = lm1.compose(lm2);

  contNI2 g3;
  contNI2::iterator itg3 = g3.begin();

  itg3 = g3.insert(itg3, 10);
  ++itg3;
  itg3 = g3.insert(itg3, 20);
  ++itg3;
  itg3 = g3.insert(itg3, 6);

  contNI2 o3;
  contNI2::iterator ito3 = o3.begin();

  ito3 = o3.insert(ito3, 16);
  ++ito3;
  ito3 = o3.insert(ito3, 32);
  ++ito3;
  o3.insert(ito3, 12);

  LMap res2(g3, o3);

  BOOST_CHECK(res1 == res2);
}

void TestLMCompose2(){
  contNI2 g1;
  contNI2::iterator itg1 = g1.begin();
  
  itg1 = g1.insert(itg1, 5);
  ++itg1;
  itg1 = g1.insert(itg1, 10);
  ++itg1;
  itg1 = g1.insert(itg1, 3);

  contNI2 o1;
  contNI2::iterator ito1 = o1.begin();

  ito1 = o1.insert(ito1, 1);
  ++ito1;
  ito1 = o1.insert(ito1, 2);
  ++ito1;
  ito1 = o1.insert(ito1, 3);

  LMap lm1(g1, o1);

  contNI2 g2;
  contNI2::iterator itg2 = g2.begin();

  itg2 = g2.insert(itg2, 2);
  ++itg2;
  itg2 = g2.insert(itg2, 2);

  contNI2 o2;
  contNI2::iterator ito2 = o2.begin();

  ito2 = o2.insert(ito2, 3);
  ++ito2;
  ito2 = o2.insert(ito2, 3);

  LMap lm2(g2, o2);

  LMap res1 = lm1.compose(lm2);

  BOOST_CHECK(res1.empty());
}

void TestInvLMap1(){
  contNI2 g1;
  contNI2::iterator itg1 = g1.begin();
  
  itg1 = g1.insert(itg1, 5);
  ++itg1;
  itg1 = g1.insert(itg1, 10);
  ++itg1;
  itg1 = g1.insert(itg1, 3);

  contNI2 o1;
  contNI2::iterator ito1 = o1.begin();

  ito1 = o1.insert(ito1, 1);
  ++ito1;
  ito1 = o1.insert(ito1, 2);
  ++ito1;
  ito1 = o1.insert(ito1, 3);

  LMap lm1(g1, o1);
  LMap res1 = lm1.invLMap();
  
  contNI2 g2;
  contNI2::iterator itg2 = g2.begin();

  float v1 = 1.0 / 5.0;
  float v2 = 1.0 / 10.0;
  float v3 = 1.0 / 3.0; 
 
  itg2 = g2.insert(itg2, v1);
  ++itg2;
  itg2 = g2.insert(itg2, v2);
  ++itg2;
  g2.insert(itg2, v3);

  contNI2 o2;
  contNI2::iterator ito2 = o2.begin();

  ito2 = o2.insert(ito2, -v1);
  ++ito2;
  ito2 = o2.insert(ito2, -v1);
  ++ito2;
  o2.insert(ito2, -1);

  LMap res2(g2, o2);

  BOOST_CHECK(res1 == res2);
}

//____________________________________________________________________________//

test_suite *init_unit_test_suite(int, char *[]){
  framework::master_test_suite().p_name.value = "Set Based Graphs";

  framework::master_test_suite().add(BOOST_TEST_CASE(&TestIntCreation1));
  framework::master_test_suite().add(BOOST_TEST_CASE(&TestIntCreation2));
  framework::master_test_suite().add(BOOST_TEST_CASE(&TestIntCreation3));
  framework::master_test_suite().add(BOOST_TEST_CASE(&TestIntCreation4));
  framework::master_test_suite().add(BOOST_TEST_CASE(&TestIntQuery1));
  framework::master_test_suite().add(BOOST_TEST_CASE(&TestIntQuery2));
  framework::master_test_suite().add(BOOST_TEST_CASE(&TestIntQuery3));
  framework::master_test_suite().add(BOOST_TEST_CASE(&TestIntQuery4));
  framework::master_test_suite().add(BOOST_TEST_CASE(&TestIntQuery5));
  framework::master_test_suite().add(BOOST_TEST_CASE(&TestIntQuery6));
  framework::master_test_suite().add(BOOST_TEST_CASE(&TestIntCap1));
  framework::master_test_suite().add(BOOST_TEST_CASE(&TestIntCap2));
  framework::master_test_suite().add(BOOST_TEST_CASE(&TestIntCap3));
  framework::master_test_suite().add(BOOST_TEST_CASE(&TestIntCap4));
  framework::master_test_suite().add(BOOST_TEST_CASE(&TestIntDiff1));
  framework::master_test_suite().add(BOOST_TEST_CASE(&TestIntDiff2));
  framework::master_test_suite().add(BOOST_TEST_CASE(&TestIntDiff3));
  framework::master_test_suite().add(BOOST_TEST_CASE(&TestIntDiff4));
  framework::master_test_suite().add(BOOST_TEST_CASE(&TestIntMin1));

  framework::master_test_suite().add(BOOST_TEST_CASE(&TestMultiCreation1));
  framework::master_test_suite().add(BOOST_TEST_CASE(&TestMultiEmpty1));
  framework::master_test_suite().add(BOOST_TEST_CASE(&TestMultiEmpty2));
  framework::master_test_suite().add(BOOST_TEST_CASE(&TestMultiEmpty3));
  framework::master_test_suite().add(BOOST_TEST_CASE(&TestMultiQuery1));
  framework::master_test_suite().add(BOOST_TEST_CASE(&TestMultiQuery2));
  framework::master_test_suite().add(BOOST_TEST_CASE(&TestMultiQuery3));
  framework::master_test_suite().add(BOOST_TEST_CASE(&TestMultiAddInter1));
  framework::master_test_suite().add(BOOST_TEST_CASE(&TestMultiAddInter2));
  framework::master_test_suite().add(BOOST_TEST_CASE(&TestMultiCap1));
  framework::master_test_suite().add(BOOST_TEST_CASE(&TestMultiCap2));
  framework::master_test_suite().add(BOOST_TEST_CASE(&TestMultiCap3));
  framework::master_test_suite().add(BOOST_TEST_CASE(&TestMultiDiff1));
  framework::master_test_suite().add(BOOST_TEST_CASE(&TestMultiDiff2));
  framework::master_test_suite().add(BOOST_TEST_CASE(&TestMultiDiff3));
  framework::master_test_suite().add(BOOST_TEST_CASE(&TestMultiDiff4));
  framework::master_test_suite().add(BOOST_TEST_CASE(&TestMultiDiff5));
  framework::master_test_suite().add(BOOST_TEST_CASE(&TestMultiCrossProd1));
  framework::master_test_suite().add(BOOST_TEST_CASE(&TestMultiMin1));

  framework::master_test_suite().add(BOOST_TEST_CASE(&TestASetCreation1));
  framework::master_test_suite().add(BOOST_TEST_CASE(&TestASetEmpty1));
  framework::master_test_suite().add(BOOST_TEST_CASE(&TestASetEmpty2));
  framework::master_test_suite().add(BOOST_TEST_CASE(&TestASetEmpty3));
  framework::master_test_suite().add(BOOST_TEST_CASE(&TestASetCap1));
  framework::master_test_suite().add(BOOST_TEST_CASE(&TestASetDiff1));
  framework::master_test_suite().add(BOOST_TEST_CASE(&TestASetMin1));

  framework::master_test_suite().add(BOOST_TEST_CASE(&TestSetCreation1));
  framework::master_test_suite().add(BOOST_TEST_CASE(&TestCompSets1));
  framework::master_test_suite().add(BOOST_TEST_CASE(&TestSetEmpty1));
  framework::master_test_suite().add(BOOST_TEST_CASE(&TestAddASets1));
  framework::master_test_suite().add(BOOST_TEST_CASE(&TestSetCap1));
  framework::master_test_suite().add(BOOST_TEST_CASE(&TestSetCap2));
  framework::master_test_suite().add(BOOST_TEST_CASE(&TestSetCap3));
  framework::master_test_suite().add(BOOST_TEST_CASE(&TestSetDiff1));
  framework::master_test_suite().add(BOOST_TEST_CASE(&TestSetMin1));


  framework::master_test_suite().add(BOOST_TEST_CASE(&TestLMCreation1));
  framework::master_test_suite().add(BOOST_TEST_CASE(&TestLMCompose1));
  framework::master_test_suite().add(BOOST_TEST_CASE(&TestLMCompose2));
  framework::master_test_suite().add(BOOST_TEST_CASE(&TestInvLMap1));
  return 0;
}

//____________________________________________________________________________//

// EOF
