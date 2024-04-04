// Copyright 2021 Pierre Talbot

#include <gtest/gtest.h>
#include <gtest/gtest-spi.h>
#include "abstract_testing.hpp"

#include "battery/vector.hpp"
#include "battery/shared_ptr.hpp"

#include "lala/primitive_upset.hpp"
#include "lala/fixpoint.hpp"

using namespace lala;
using namespace battery;

using F = TFormula<standard_allocator>;

using zd = local::ZDec;
using Oct = Octagon<zd, standard_allocator>;

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
  Oct oct = create_and_interpret_and_tell<Oct>("var 0..10: x; var 0..10: y;\
    constraint int_le(int_plus(x, y), 5);");
  refine_and_test(oct, 1, {Itv(0,10), Itv(0,10)}, {Itv(0,5), Itv(0,5)}, true);
}

// x + y <= 5 /\ y + z <= 3
TEST(OctagonTest, TransivityTemporalConstraint) {
  Oct oct = create_and_interpret_and_tell<Oct>("var 0..10: x; var 0..10: y; var 0..10: z\
    constraint int_le(int_plus(x, y), 5);\
    constraint int_le(int_plus(y, z), 3);");
  refine_and_test(oct, 1, {Itv(0,10), Itv(0,10), Itv(0,10)}, {Itv(0,5), Itv(0,3), Itv(0,3)}, true);
}
