// Copyright 2021 Pierre Talbot

#include "lala/octagon.hpp"

#include <gtest/gtest.h>
#include <gtest/gtest-spi.h>
#include "abstract_testing.hpp"

#include "battery/vector.hpp"
#include "battery/shared_ptr.hpp"

#include "lala/universes/primitive_upset.hpp"
#include "lala/fixpoint.hpp"
#include "lala/interval.hpp"
#include "lala/logic/algorithm.hpp"

using namespace lala;
using namespace battery;

using F = TFormula<standard_allocator>;

using zi = local::ZInc;
using zd = local::ZDec;
using Itv = Interval<zi>;
using Oct = Octagon<Itv, standard_allocator>;

template <class L>
void test_extract(const L& oct, bool is_ua) {
  AbstractDeps<standard_allocator> deps(standard_allocator{});
  L copy1(oct, deps);
  EXPECT_EQ(oct.is_extractable(), is_ua);
  if(oct.is_extractable()) {
    oct.extract(copy1);
    EXPECT_EQ(oct.is_top(), copy1.is_top());
    EXPECT_EQ(oct.is_bot(), copy1.is_bot());
    for(int i = 0; i < oct.vars(); ++i) {
      EXPECT_EQ(oct[i], copy1[i]);
    }
  }
}

template<class L>
void refine_and_test(L& oct, int num_refine, const std::vector<Itv>& before, const std::vector<Itv>& after, bool is_ua, bool expect_changed = true) {
  EXPECT_EQ(oct.num_refinements(), num_refine);
  for(int i = 0; i < before.size(); ++i) {
    EXPECT_EQ(oct[i], before[i]) << "oct[" << i << "]";
  }
  local::BInc has_changed = GaussSeidelIteration{}.fixpoint(oct);
  printf("*********************************************\n");
  oct.print();
  EXPECT_EQ(has_changed, expect_changed);
  for(int i = 0; i < after.size(); ++i) {
    EXPECT_EQ(oct[i], after[i]) << "oct[" << i << "]";
  }
  test_extract(oct, is_ua);
}

template<class L>
void refine_and_test(L& oct, int num_refine, const std::vector<Itv>& before_after, bool is_ua = false) {
  refine_and_test(oct, num_refine, before_after, before_after, is_ua, false);
}

// x + y <= 5
TEST(OctagonTest, TemporalConstraint1) {
  Oct oct = create_and_interpret_and_tell<Oct>("var 0..10: x; var 0..10: y; constraint int_le(int_plus(x, y), 5);");
  refine_and_test(oct, 1, {Itv(0,10), Itv(0,10)}, {Itv(0,5), Itv(0,5)}, true);
}
// -x + y <= 5
TEST(OctagonTest, TemporalConstraint2) {
  Oct oct = create_and_interpret_and_tell<Oct>("var 0..10: x; var 0..10: y; constraint int_le(int_plus(int_neg(x), y), 5);");
//  refine_and_test(oct, 1, {Itv(0,10), Itv(0,10)}, {Itv(0,5), Itv(0,5)}, true);
}

// x0 in [1,3], x1 in [1,4], x0 - x1 <= 1, -x0 + x1 <= 1
TEST(OctagonTest, IctaiConstraints) {
  Oct oct = create_and_interpret_and_tell<Oct>("var 1..3: x0; var 1..4: x1; constraint int_le(int_minus(x0, x1), 1);  constraint int_le(int_plus(int_neg(x0), x1), 1);");
  oct.print();
  refine_and_test(oct, 1, {Itv(0,10), Itv(0,10)}, {Itv(0,5), Itv(0,5)}, true);
}

// x <= 5
TEST(OctagonTest, UnaryConstraint1) {
  VarEnv<standard_allocator> env;
  IDiagnostics diagnostics;
  battery::vector<battery::tuple<int, int, Itv::UB>, standard_allocator> intermediate;
  Oct oct = create_and_interpret_and_tell<Oct>("var 0..10: x;",env);
  refine_and_test(oct, 14, {Itv(0,10)}, {Itv(0,10)}, true);

  auto f = parse_flatzinc_str<standard_allocator>("constraint int_le(x, 5);");
  oct.interpret<IKind::TELL>(*f,env,intermediate,diagnostics);
  oct.tell(intermediate);
  refine_and_test(oct, 14, {Itv(0,5)}, {Itv(0,5)}, true);
}


TEST(OctagonTest, IncrementallyClosingOctagonPaper) {
  VarEnv<standard_allocator> env;
  IDiagnostics diagnostics;
  battery::vector<battery::tuple<int, int, Itv::UB>, standard_allocator> intermediate;
  Oct oct = create_and_interpret_and_tell<Oct>("var -inf..3: x0; var -inf..2: x1; constraint int_le(int_plus(x0,x1),6);constraint int_le(int_plus(x0,x1),6);constraint int_le(int_minus(int_neg(x0),x1),5); constraint int_le(int_neg(x0),3);",env);
  oct.print();
  refine_and_test(oct, 64, {Itv(0,10)}, {Itv(0,10)}, true);
  oct.print();

  // auto f = parse_flatzinc_str<standard_allocator>("constraint int_le(x, 5);");
  // oct.interpret<IKind::TELL>(*f,env,intermediate,diagnostics);
  // oct.tell(intermediate);
  // refine_and_test(oct, 14, {Itv(0,5)}, {Itv(0,5)}, true);
}





// x + y <= 5 /\ y + z <= 3
TEST(OctagonTest, TransivityTemporalConstraint) {
  Oct oct = create_and_interpret_and_tell<Oct>("var 0..10: x; var 0..10: y; var 0..10: z;\
    constraint int_le(int_plus(x, y), 5);\
    constraint int_le(int_plus(y, z), 3);");
  //refine_and_test(oct, 1, {Itv(0,10), Itv(0,10), Itv(0,10)}, {Itv(0,5), Itv (0,3), Itv(0,3)}, true);
}



TEST(OctagonTest,is_bot) {
  auto oct = Octagon<Itv, standard_allocator>::bot();
  ASSERT_TRUE(oct.is_bot());
}

TEST(OctagonTest,deinterpret) {
  VarEnv<standard_allocator> env;
  IDiagnostics diagnostics;

  Oct oct = create_and_interpret_and_tell<Oct>("var 0..10: x; var 0..10: y;\
    constraint int_le(int_plus(x, y), 5);",env);
  oct.print();
  printf("\n");
  auto f = oct.deinterpret(env);
  f.print(false);
  printf("\n");
  GaussSeidelIteration{}.fixpoint(oct);
  printf("\noct après des \n");
  oct.print();
  printf("\nf après des \n");
  auto f2 = oct.deinterpret(env);
  f2.print(false);

  VarEnv<standard_allocator> env2;
  auto oct2 = create_and_interpret_and_tell<Oct>(f,env2,diagnostics);
  EXPECT_TRUE(oct2.has_value());
  oct2.value().print();
  printf("\n");
  EXPECT_EQ(oct2->deinterpret(env2),f);
  oct2->deinterpret(env2).print(false);
  printf("\n");


  VarEnv<standard_allocator> env3;
  auto oct3 = create_and_interpret_and_tell<Oct>(f2,env3,diagnostics);
  EXPECT_EQ(oct3->deinterpret(env3),f2);
}