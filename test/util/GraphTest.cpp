#include <iostream>

#include <boost/test/unit_test.hpp>
#include <boost/test/included/unit_test.hpp>

#include <ast/expression.h>
#include <util/graph/graph_definition.h>

using namespace boost::unit_test;

/// @brief If sb-graphs containers implementation changes (uses vector, for example)
/// this typedef should also change.
typedef list<Interval> contInt;
typedef list<MultiInterval> contMulti;
typedef list<AtomSet> contAS;
typedef list<Set> contSet;

typedef list<NI> contNI;

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

// Cap should be commutative
void TestIntCap1(){
  Interval i1(10, 2, 20, false);
  Interval i2(0, 3, 25, false);

  Interval i3 = i1.cap(i2);
  Interval i4 = i2.cap(i1);

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

  BOOST_CHECK(i3 == i4);
}

void TestIntDiff1(){
  Interval i1(0, 2, 30, false);
  Interval i2(true);

  contInt aux;
  aux.insert(aux.begin(), i1);

  contInt res = i1.diff(i2); 

  BOOST_CHECK(res == aux);
}

void TestIntDiff2(){
  Interval i1(0, 2, 30, false);
  Interval i2(10, 3, 40, false);

  Interval i3(0, 2, 8, false);
  Interval i4(12, 6, 24, false);
  Interval i5(14, 6, 26, false);
  Interval i6(30, 2, 30, false);

  contInt res1 = i1.diff(i2);

  contInt res2; 
  contInt::iterator it = res2.begin();
  it = res2.insert(it, i3);
  ++it;
  it = res2.insert(it, i4);
  ++it;
  it = res2.insert(it, i5);
  ++it;
  it = res2.insert(it, i6);

  BOOST_CHECK(res1 == res2);
}

void TestIntDiff3(){
  Interval i1(0, 2, Inf, false);
  Interval i2(10, 3, 40, false);

  Interval i3(0, 2, 8, false);
  Interval i4(12, 6, 36, false);
  Interval i5(14, 6, 38, false);
  Interval i6(42, 2, Inf, false);

  contInt res1 = i1.diff(i2);

  contInt res2; 
  contInt::iterator it = res2.begin();
  it = res2.insert(it, i3);
  ++it;
  it = res2.insert(it, i4);
  ++it;
  it = res2.insert(it, i5);
  ++it;
  it = res2.insert(it, i6);

  //contInt::iterator

  BOOST_CHECK(res1 == res2);
}

void TestIntMin1(){
  Interval i(10, 3, 40, false);

  NI res1 = i.minElem();

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

  contInt res;
  contInt::iterator it = res.begin();

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

void TestMultiAddInter1(){
  Interval i1(0, 2, 10, false);
  MultiInterval mi1;

  mi1.addInter(i1);  

  contInt ints2;
  contInt::iterator it2 = ints2.begin();
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

  contInt ints2;
  contInt::iterator it2 = ints2.end();
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

  Interval i7(8, 6, 14, false);
  Interval i8(true);
  Interval i9(27, 1, 30, false);
  MultiInterval mi3;

  mi3.addInter(i7);
  mi3.addInter(i8);
  mi3.addInter(i9);

  MultiInterval mi4 = mi1.cap(mi2);
  MultiInterval mi5 = mi2.cap(mi1);
  
  BOOST_CHECK(mi3 == mi4 && mi4 == mi5);
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
  Interval i5(true);
  Interval i6(27, 1, 35, false);
  MultiInterval mi2;

  mi2.addInter(i4);
  mi2.addInter(i5);
  mi2.addInter(i6);

  Interval i7(8, 6, 14, false);
  Interval i8(true);
  Interval i9(27, 1, 30, false);

  MultiInterval capres;

  capres.addInter(i7);
  capres.addInter(i8);
  capres.addInter(i9);

  MultiInterval mi3;
  MultiInterval mi4;
  MultiInterval mi5;
  MultiInterval mi6;

  Interval i10(0, 2, 6, false);
  Interval i11(10, 6, 10, false);
  Interval i12(12, 6, 12, false);
  Interval i13(16, 2, 20, false);

  mi3.addInter(i10);
  mi3.addInter(i2);
  mi3.addInter(i3);
  mi4.addInter(i11);
  mi4.addInter(i2);
  mi4.addInter(i3);
  mi5.addInter(i12);
  mi5.addInter(i2);
  mi5.addInter(i3);
  mi6.addInter(i13);
  mi6.addInter(i2);
  mi6.addInter(i3);

  MultiInterval mi7;

  mi7.addInter(i7);
  mi7.addInter(i2);
  mi7.addInter(i3);

  MultiInterval mi8;

  Interval i14(25, 1, 26, false); 

  mi8.addInter(i7);
  mi8.addInter(i8);
  mi8.addInter(i14);

  contMulti res1 = mi1.diff(mi2);

  contMulti res2;
  contMulti::iterator it2 = res2.begin();
  it2 = res2.insert(it2, mi3);
  ++it2;
  it2 = res2.insert(it2, mi4);
  ++it2;
  it2 = res2.insert(it2, mi5);
  ++it2;
  it2 = res2.insert(it2, mi6);
  ++it2;
  it2 = res2.insert(it2, mi7);
  ++it2;
  it2 = res2.insert(it2, mi8);

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
  Interval i5(true);
  Interval i6(27, 1, 35, false);
  MultiInterval mi2;

  mi2.addInter(i4);
  mi2.addInter(i5);
  mi2.addInter(i6);

  Interval i7(8, 6, 14, false);
  Interval i8(true);
  Interval i9(27, 1, 30, false);

  MultiInterval mi3;
  MultiInterval mi4;

  Interval i10(5, 3, 5, false);
  Interval i11(11, 6, 11, false);

  mi3.addInter(i10);
  mi3.addInter(i5);
  mi3.addInter(i6);
  mi4.addInter(i11);
  mi4.addInter(i5);
  mi4.addInter(i6);

  MultiInterval mi5;

  mi5.addInter(i7);
  mi5.addInter(i5);
  mi5.addInter(i6);

  MultiInterval mi6;

  Interval i12(31, 1, 35, false);

  mi6.addInter(i7);
  mi6.addInter(i8);
  mi6.addInter(i12); 

  contMulti res1 = mi2.diff(mi1);

  contMulti res2;
  contMulti::iterator it2 = res2.begin();
  it2 = res2.insert(it2, mi3);
  ++it2;
  it2 = res2.insert(it2, mi4);
  ++it2;
  it2 = res2.insert(it2, mi5);
  ++it2;
  it2 = res2.insert(it2, mi6);

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

  contNI res1 = mi.minElem();

  contNI res2;
  contNI::iterator it2 = res2.begin();
 
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

void TestASetCap1(){
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

  Interval i7(8, 6, 14, false);
  Interval i8(true);
  Interval i9(27, 1, 30, false);
  MultiInterval mi3;

  mi3.addInter(i7);
  mi3.addInter(i8);
  mi3.addInter(i9);

  AtomSet as1(mi1);
  AtomSet as2(mi2);
  AtomSet as3(mi3);

  AtomSet as4 = as1.cap(as2);
  AtomSet as5 = as2.cap(as1);

  BOOST_CHECK(as3 == as4 && as4 == as5);
}

void TestASetDiff1(){
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

  Interval i7(8, 6, 14, false);
  Interval i8(true);
  Interval i9(27, 1, 30, false);

  MultiInterval capres;

  capres.addInter(i7);
  capres.addInter(i8);
  capres.addInter(i9);

  MultiInterval mi3;
  MultiInterval mi4;
  MultiInterval mi5;
  MultiInterval mi6;

  Interval i10(0, 2, 6, false);
  Interval i11(10, 6, 10, false);
  Interval i12(12, 6, 12, false);
  Interval i13(16, 2, 20, false);

  mi3.addInter(i10);
  mi3.addInter(i2);
  mi3.addInter(i3);
  mi4.addInter(i11);
  mi4.addInter(i2);
  mi4.addInter(i3);
  mi5.addInter(i12);
  mi5.addInter(i2);
  mi5.addInter(i3);
  mi6.addInter(i13);
  mi6.addInter(i2);
  mi6.addInter(i3);

  MultiInterval mi7;

  mi7.addInter(i7);
  mi7.addInter(i2);
  mi7.addInter(i3);

  MultiInterval mi8;

  Interval i14(25, 1, 26, false); 

  mi8.addInter(i7);
  mi8.addInter(i8);
  mi8.addInter(i14);

  AtomSet as1(mi1);
  AtomSet as2(mi2);
  contAS res1 = as1.diff(as2);

  contAS res2;
  contAS::iterator it2 = res2.begin();
  it2 = res2.insert(it2, AtomSet(mi3));
  ++it2;
  it2 = res2.insert(it2, AtomSet(mi4));
  ++it2;
  it2 = res2.insert(it2, AtomSet(mi5));
  ++it2;
  it2 = res2.insert(it2, AtomSet(mi6));
  ++it2;
  it2 = res2.insert(it2, AtomSet(mi7));
  ++it2;
  it2 = res2.insert(it2, AtomSet(mi8));

  BOOST_CHECK(res1 == res2);
}

void TestASetDiff2(){
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

  Interval i7(8, 6, 14, false);
  Interval i8(true);
  Interval i9(27, 1, 30, false);

  MultiInterval mi3;
  MultiInterval mi4;

  Interval i10(5, 3, 5, false);
  Interval i11(11, 6, 11, false);

  mi3.addInter(i10);
  mi3.addInter(i5);
  mi3.addInter(i6);
  mi4.addInter(i11);
  mi4.addInter(i5);
  mi4.addInter(i6);

  MultiInterval mi5;

  mi5.addInter(i7);
  mi5.addInter(i5);
  mi5.addInter(i6);

  MultiInterval mi6;

  Interval i12(31, 1, 35, false);

  mi6.addInter(i7);
  mi6.addInter(i8);
  mi6.addInter(i12); 

  AtomSet as1(mi1);
  AtomSet as2(mi2);
  contAS res1 = as2.diff(as1);

  contAS res2;
  contAS::iterator it2 = res2.begin();
  it2 = res2.insert(it2, AtomSet(mi3));
  ++it2;
  it2 = res2.insert(it2, AtomSet(mi4));
  ++it2;
  it2 = res2.insert(it2, AtomSet(mi5));
  ++it2;
  it2 = res2.insert(it2, AtomSet(mi6));

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

  contNI res1 = mi.minElem();

  contNI res2;
  contNI::iterator it2 = res2.begin();
 
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

  Interval i7(0, 1, Inf, false);
  Interval i8(20, 3, 50, false);
  Interval i9(true);

  MultiInterval mi3;

  AtomSet as3(mi3);

  Set s2;

  s2.addAtomSet(as3);

  Set res1 = s1.cap(s2);
  Set res2 = s2.cap(s1);

  Interval i10(0, 2, 20, false);
  Interval i11(32, 6, 38, false);
  Interval i12(true);

  MultiInterval mi4;

  mi4.addInter(i10);
  mi4.addInter(i11);
  mi4.addInter(i12);

  MultiInterval mi5;

  mi5.addInter(i4);
  mi5.addInter(i5);
  mi5.addInter(i5);

  AtomSet as4(mi4);
  AtomSet as5(mi5);

  Set res3;

  res3.addAtomSet(as4);
  res3.addAtomSet(as5);

  BOOST_CHECK(res1 == res2 && res2 == res3);
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
  framework::master_test_suite().add(BOOST_TEST_CASE(&TestIntCap1));
  framework::master_test_suite().add(BOOST_TEST_CASE(&TestIntCap2));
  framework::master_test_suite().add(BOOST_TEST_CASE(&TestIntCap3));
  framework::master_test_suite().add(BOOST_TEST_CASE(&TestIntCap4));
  framework::master_test_suite().add(BOOST_TEST_CASE(&TestIntDiff1));
  framework::master_test_suite().add(BOOST_TEST_CASE(&TestIntDiff2));
  framework::master_test_suite().add(BOOST_TEST_CASE(&TestIntDiff3));
  framework::master_test_suite().add(BOOST_TEST_CASE(&TestIntMin1));

  framework::master_test_suite().add(BOOST_TEST_CASE(&TestMultiCreation1));
  framework::master_test_suite().add(BOOST_TEST_CASE(&TestMultiEmpty1));
  framework::master_test_suite().add(BOOST_TEST_CASE(&TestMultiEmpty2));
  framework::master_test_suite().add(BOOST_TEST_CASE(&TestMultiEmpty3));
  framework::master_test_suite().add(BOOST_TEST_CASE(&TestMultiAddInter1));
  framework::master_test_suite().add(BOOST_TEST_CASE(&TestMultiAddInter2));
  framework::master_test_suite().add(BOOST_TEST_CASE(&TestMultiCap1));
  framework::master_test_suite().add(BOOST_TEST_CASE(&TestMultiCap2));
  framework::master_test_suite().add(BOOST_TEST_CASE(&TestMultiDiff1));
  framework::master_test_suite().add(BOOST_TEST_CASE(&TestMultiDiff2));
  framework::master_test_suite().add(BOOST_TEST_CASE(&TestMultiMin1));

  framework::master_test_suite().add(BOOST_TEST_CASE(&TestASetCreation1));
  framework::master_test_suite().add(BOOST_TEST_CASE(&TestASetEmpty1));
  framework::master_test_suite().add(BOOST_TEST_CASE(&TestASetEmpty2));
  framework::master_test_suite().add(BOOST_TEST_CASE(&TestASetEmpty3));
  framework::master_test_suite().add(BOOST_TEST_CASE(&TestASetCap1));
  framework::master_test_suite().add(BOOST_TEST_CASE(&TestASetDiff1));
  framework::master_test_suite().add(BOOST_TEST_CASE(&TestASetDiff2));
  framework::master_test_suite().add(BOOST_TEST_CASE(&TestASetMin1));

  framework::master_test_suite().add(BOOST_TEST_CASE(&TestSetCreation1));
  framework::master_test_suite().add(BOOST_TEST_CASE(&TestAddASets1));
  framework::master_test_suite().add(BOOST_TEST_CASE(&TestSetCap1));
  framework::master_test_suite().add(BOOST_TEST_CASE(&TestSetCap2));
  framework::master_test_suite().add(BOOST_TEST_CASE(&TestSetCap3));

  return 0;
}

//____________________________________________________________________________//

// EOF
