#include <iostream>

#include <boost/test/unit_test.hpp>
#include <boost/test/included/unit_test.hpp>

#include <ast/expression.h>
#include <util/graph/graph_definition.h>

using namespace boost::unit_test;

/// @brief If sb-graphs containers implementation changes (uses vector, for example)
/// this typedef should also change.
typedef list<Interval> contInt;

//____________________________________________________________________________//

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

  contInt i3 = i1.diff(i2);
  Interval res1(0, 2, 8, false);
  Interval res2(12, 6, 24, false);
  Interval res3(14, 6, 26, false);
  Interval res4(30, 2, 30, false);

  contInt res; 
  contInt::iterator it = res.begin();
  it = res.insert(it, res1);
  it = res.insert(it, res2);
  it = res.insert(it, res3);
  it = res.insert(it, res4);

  it = res.begin();
  contInt::iterator it2 = i3.begin();

  BOOST_CHECK(i3 == res);
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


/*
void TestCommInsertion(){
  Interval i1(0, 2, 20, false);
  Interval i2(5, 3, 15, false);

  MultiInterval mi1;
  mi1.addInter(i1);
  MultiInterval mi2;
  mi2.addInter(i2);

  AtomSet as1(mi1);
  AtomSet as2(mi2);

  Set s1;
  s1.addAtomSet(as1);
  std::cout << s1 << "\n";
  s1.addAtomSet(as2);
  s1.addAtomSet(as1);
  s1.addAtomSet(as2);
  s1.addAtomSet(as1);
  s1.addAtomSet(as2);
  Set s2;
  s2.addAtomSet(as2);
  s2.addAtomSet(as1);

  std::cout << s1 << "\n" << s2 << "\n";

  BOOST_CHECK(s1 == s2);
}*/

//____________________________________________________________________________//

test_suite *init_unit_test_suite(int, char *[]){
  framework::master_test_suite().p_name.value = "Set Based Graphs";

  framework::master_test_suite().add(BOOST_TEST_CASE(&TestIntCreation1));
  framework::master_test_suite().add(BOOST_TEST_CASE(&TestIntCreation2));
  framework::master_test_suite().add(BOOST_TEST_CASE(&TestIntCreation3));
  framework::master_test_suite().add(BOOST_TEST_CASE(&TestIntQuery1));
  framework::master_test_suite().add(BOOST_TEST_CASE(&TestIntQuery2));
  framework::master_test_suite().add(BOOST_TEST_CASE(&TestIntQuery3));
  framework::master_test_suite().add(BOOST_TEST_CASE(&TestIntQuery4));
  framework::master_test_suite().add(BOOST_TEST_CASE(&TestIntQuery5));
  framework::master_test_suite().add(BOOST_TEST_CASE(&TestIntCap1));
  framework::master_test_suite().add(BOOST_TEST_CASE(&TestIntCap2));
  framework::master_test_suite().add(BOOST_TEST_CASE(&TestIntCap3));
  framework::master_test_suite().add(BOOST_TEST_CASE(&TestIntDiff1));
  framework::master_test_suite().add(BOOST_TEST_CASE(&TestIntDiff2));

  framework::master_test_suite().add(BOOST_TEST_CASE(&TestMultiAddInter1));
  framework::master_test_suite().add(BOOST_TEST_CASE(&TestMultiAddInter2));
  framework::master_test_suite().add(BOOST_TEST_CASE(&TestMultiCap1));

  return 0;
}

//____________________________________________________________________________//

// EOF
