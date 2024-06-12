// Copyright 2024 Pierre Talbot

#ifndef LALA_OCTAGON_HPP
#define LALA_OCTAGON_HPP

#include "util.hpp"
#include "battery/vector.hpp"
#include "battery/unique_ptr.hpp"
#include "battery/shared_ptr.hpp"
#include "battery/allocator.hpp"

#include "lala/logic/logic.hpp"
#include "lala/universes/primitive_upset.hpp"
#include "lala/vstore.hpp"
#include "lala/abstract_deps.hpp"
#include "lala/interval.hpp"
#include "lala/logic/env.hpp"

namespace lala {
  /** Octagon is an abstract domain built on top of an abstract universe `U`. */
  template<class V, class Allocator>
  class Octagon {
  public:
    using U = typename V::UB;
    using universe_type = V;
    using local_universe = typename universe_type::local_type;
    using local_cell_type = typename U::local_type;
    using allocator_type = Allocator;
    using this_type = Octagon<V, Allocator>;
    using universe_list_type = battery::vector<U, allocator_type>;
    using dbm_type = battery::vector<universe_list_type, allocator_type>;
    using memory_type = typename universe_type::memory_type;
    template<class Alloc = allocator_type>
    using snapshot_type = battery::vector<battery::vector<U, Alloc>, Alloc>;
    template<class Alloc>
    using tell_type = battery::vector<battery::tuple<int, int, U>, Alloc>;

    template<class Alloc>
    using ask_type = tell_type<Alloc>;

    constexpr static const bool is_abstract_universe = false;
    constexpr static const bool sequential = universe_type::sequential;
    constexpr static const bool is_totally_ordered = false;
    constexpr static const bool preserve_bot = true;
    constexpr static const bool preserve_top = universe_type::preserve_top;
    constexpr static const bool preserve_join = universe_type::preserve_join;
    constexpr static const bool preserve_meet = universe_type::preserve_meet;
    constexpr static const bool injective_concretization = universe_type::injective_concretization;
    constexpr static const bool preserve_concrete_covers = universe_type::preserve_concrete_covers;
    constexpr static const char* name = "Octagon";

    template<class U2, class Alloc2>
    friend class Octagon;

  private:
    AType atype;
    dbm_type dbm;
    mutable int nbVars = 0;
    long floyd_steps = 0;
    long tight_steps = 0;
    long str_steps = 0;
    long _num_refinements = 0;
    BInc<memory_type> is_at_top;

    enum FormulaEnum {
      EXISTENTIAL, UNARY, BINARY, SET_IN, UNKNOWN
    };


    CUDA INLINE size_t dim1(size_t i) { return i / (dbm.size() * dbm.size()); }
    CUDA INLINE size_t dim2(size_t i) { return i % (dbm.size() * dbm.size()) / dbm.size(); }
    CUDA INLINE size_t dim3(size_t i) { return i % (dbm.size() * dbm.size()) % dbm.size(); }

    template<class F>
    CUDA NI FormulaEnum get_type(const F& f) const {
      if (f.is(F::E)) {
        return FormulaEnum::EXISTENTIAL;
      }
      if (f.is_binary() && is_arithmetic_comparison_operator<F>(f.sig())) {
        if (f.seq(0).is_binary()) {
          return FormulaEnum::BINARY;
        }
        return FormulaEnum::UNARY;
      }
      if (f.is_binary() && f.sig() == IN) {
        return FormulaEnum::SET_IN;
      }
      return FormulaEnum::UNKNOWN;
    }

    template<class F>
    [[nodiscard]] bool is_arithmetic_comparison_operator(const Sig& sig) const {
      return (sig == LEQ || sig == GEQ);
    }

    template<bool diagnose = false, class F>
    bool var_with_symbol(F& f, IDiagnostics diagnostics, battery::tuple<char, F>& result) const {
      if (f.is_variable()) {
        result = battery::make_tuple('+', f);
        return true;
      }
      if (f.is(F::Seq) && f.seq().size() == 1 && f.sig() == NEG) {
        auto var = f.seq(0);
        var.type_as(aty());
        result = battery::make_tuple('-', var);
        return true;
      }
      RETURN_INTERPRETATION_ERROR("There is an error when we interpret a variable formula");
    }

    bool index_for_vars(const char symbol_of_var, const AVar& avari, Sig& arith_op, const AVar& avarj,
                        battery::tuple<int, int,int,int>& result) const {
      int index1 = -1;
      int index2 = -1;
      int index3 = -1;
      int index4 = -1;
      switch (arith_op) {
        case SUB:
          if (symbol_of_var == '-') {
            index1 = avari.vid() * 2 + 1;
            index2 = avarj.vid() * 2;
            index3 = avarj.vid()*2+1;
            index4 = avari.vid()*2;
          }
          else {
            index1 = avari.vid() * 2;
            index2 = avarj.vid() * 2;
            index3 = index2+1;
            index4 = index1+1;
          }
          break;
        case ADD:
          if (symbol_of_var == '-') {
            index1 = avari.vid() * 2 + 1;
            index2 = avarj.vid() * 2 + 1;
            index3 = avarj.vid() * 2;
            index4 = avari.vid() * 2;
          }
          else {
            index1 = avari.vid() * 2;
            index2 = avarj.vid() * 2 + 1;
            index3 = index2 - 1;
            index4 = index1 + 1;
          }
          break;
        default:
          return false;
      }
      result = battery::make_tuple(index1, index2,index3,index4);
      return index1 != -1 && index2 != -1;
    }


    void index_for_var(const char symbol_of_var, const AVar& avari, battery::tuple<int, int>& result) const {
      int index1 = -1;
      int index2 = -1;
      if (symbol_of_var == '-') {
        index1 = avari.vid() * 2 + 1;
        index2 = avari.vid() * 2;
      }
      else {
        index1 = avari.vid() * 2;
        index2 = avari.vid() * 2 + 1;
      }
      assert(index1>=0 && index2>=0);
      result = battery::make_tuple(index1, index2);
    }

    void init_size() {
      floyd_steps = dbm.size() * dbm.size() * dbm.size();
      tight_steps = floyd_steps + dbm.size();
      str_steps = tight_steps + dbm.size() * dbm.size();
      _num_refinements = str_steps;
    }

    template<bool diagnose = false, class F, class Env, class Alloc2>
    CUDA NI bool interpret_tell_x_op_k(char symbol_x, const AVar& x, const Sig op, const F& constant, Env& env,
                                       tell_type<Alloc2>& tell,
                                       IDiagnostics& diagnostics) const {
      assert(!x.is_untyped());
      local_cell_type u;
      auto tmp = F::make_avar(AVar{});
      tmp.type_as(aty());
      if (constant.is(F::Z)) {
        logic_int k = constant.z() * 2;
        U::template interpret_tell<diagnose>(F::make_binary(tmp, op, F::make_z(k), aty()), env, u,
                                             diagnostics);
      }
      else {
        logic_real k = battery::make_tuple(battery::mul_down(battery::get<0>(constant.r()), 2.0),
                                           battery::mul_up(battery::get<1>(constant.r()), 2.0));
        U::template interpret_tell<diagnose>(F::make_binary(tmp, op, F::make_real(k), aty()), env, u,
                                             diagnostics);
      }

      battery::tuple<int, int> indexes;
      index_for_var(symbol_x, x, indexes);

      tell.emplace_back(battery::make_tuple(battery::get<0>(indexes), battery::get<1>(indexes), u));
      return true;
    }

    template<bool diagnose = false, class F, class Env, class Alloc2>
    CUDA NI bool interpret_tell_x_y_op_k(const F& f, const char symbol_x, const AVar& x, Sig arithm_sig, const AVar& y,
                                         const F& constant, Env& env, tell_type<Alloc2>& tell,
                                         IDiagnostics& diagnostics) const {
      battery::tuple<int, int, int, int> index_tuple;
      if (!index_for_vars(symbol_x, x, arithm_sig, y, index_tuple)) {
        RETURN_INTERPRETATION_ERROR("Octagons only handle + or - as arithmetic symbols.");
      }
      int index1 = battery::get<0>(index_tuple);
      int index2 = battery::get<1>(index_tuple);
      int index3 = battery::get<2>(index_tuple);
      int index4 = battery::get<3>(index_tuple);

      local_cell_type u;
      auto tmp = F::make_avar(AVar{});
      tmp.type_as(aty());
      local_cell_type::template interpret_tell<diagnose>(F::make_binary(tmp, LEQ, constant, aty()), env, u,
                                                         diagnostics);
      tell.emplace_back(battery::make_tuple(index1, index2, u));
      tell.emplace_back(battery::make_tuple(index3, index4, u));
      return true;
    }

  public:
    CUDA Octagon(const this_type& other)
      : atype(other.atype), dbm(other.dbm), is_at_top(other.is_at_top), nbVars(other.nbVars) {
      init_size();
    }

    /** Initialize an empty store. */
    CUDA Octagon(AType atype, const allocator_type& alloc = allocator_type())
      : atype(atype), dbm(alloc), is_at_top(false) {
      init_size();
    }

    CUDA Octagon(AType atype, size_t size, const allocator_type& alloc = allocator_type())
      : atype(atype), dbm(size, alloc), is_at_top(false) {
      init_size();
    }

    template<class R>
    CUDA Octagon(const Octagon<R, allocator_type>& other)
      : atype(other.atype), dbm(other.dbm), is_at_top(other.is_at_top), nbVars(other.nbVars) {
      init_size();
    }

    template<class R, class Alloc2>
    CUDA Octagon(const Octagon<R, Alloc2>& other, const allocator_type& alloc = allocator_type())
      : atype(other.atype), dbm(other.dbm, alloc), is_at_top(other.is_at_top), nbVars(other.nbVars) {
      init_size();
    }

    /** Copy the Octagon `other` in the current element.
     *  `deps` can be empty and is not used besides to get the allocator (since this abstract domain does not have dependencies). */
    template<class R, class Alloc2, class... Allocators>
    CUDA Octagon(const Octagon<R, Alloc2>& other, const AbstractDeps<Allocators...>& deps)
      : Octagon(other, deps.template get_allocator<allocator_type>()) {
      init_size();
    }

    CUDA Octagon(this_type&& other): atype(other.atype), dbm(std::move(other.dbm)), is_at_top(other.is_at_top),
                                     nbVars(other.nbVars) {
      init_size();
    }

    CUDA allocator_type get_allocator() const {
      return dbm.get_allocator();
    }

    CUDA AType aty() const {
      return atype;
    }

    CUDA static this_type bot(AType atype = UNTYPED,
                              const allocator_type& alloc = allocator_type()) {
      return Octagon{atype, alloc};
    }

    /** A special symbolic element representing top. */
    CUDA static this_type top(AType atype = UNTYPED,
                              const allocator_type& alloc = allocator_type()) {
      auto oct = Octagon(atype, alloc);
      oct->dbm[0][0] = -1;
      return oct;
    }

    template<class Env>
    CUDA static this_type bot(Env& env,
                              const allocator_type& alloc = allocator_type()) {
      AType atype = env.extends_abstract_dom();
      return bot(atype, alloc);
    }

    template<class Env>
    CUDA static this_type top(Env& env,
                              const allocator_type& alloc = allocator_type()) {
      AType atype = env.extends_abstract_dom();
      return top(atype, alloc);
    }


    template<bool diagnose = false, class F, class Env, class Alloc2>
    CUDA NI bool interpret_tell_existential(const F& f, Env& env, tell_type<Alloc2>& tell,
                                            IDiagnostics& diagnostics) const {
      DEBUG_PRINT("start\n");
      AVar avar;
      avar.type_as(aty());
      auto result = env.template interpret<diagnose>(f.map_atype(aty()), avar, diagnostics);
      DEBUG_PRINT("end\n");
      return result;
    }

    template<bool diagnose = false, class F, class Env, class Alloc2>
    CUDA NI bool interpret_tell_binary(const F& f, Env& env, tell_type<Alloc2>& tell,
                                       IDiagnostics& diagnostics) const {
      DEBUG_PRINT("start\n");
      auto left_operand = f.seq(0);
      auto right_operand = f.seq(1);

      battery::tuple<char, F> xi;
      auto result = var_with_symbol<diagnose, F>(left_operand.seq(0), diagnostics, xi);
      if (!result) {
        return false;
      }


      battery::tuple<char, F> xj;
      result = var_with_symbol<diagnose, F>(left_operand.seq(1), diagnostics, xj);
      if (!result) {
        return false;
      }

      auto constant = right_operand;
      auto var1 = env.variable_of(battery::get<1>(xi).lv());
      if (!var1.has_value()) {
        RETURN_INTERPRETATION_ERROR("The variable " + xi + " does not exist.");
      }
      auto var2 = env.variable_of(battery::get<1>(xj).lv());
      if (!var2.has_value()) {
        RETURN_INTERPRETATION_ERROR("The variable " + xj + " does not exist.");
      }
      if (!constant.is_constant()) {
        RETURN_INTERPRETATION_ERROR("The right operand of the formula must be a constant.");
      }
      AVar avari = var1.value().get().avar_of(atype).value();
      AVar avarj = var2.value().get().avar_of(atype).value();
      local_cell_type u;


      switch (f.sig()) {
        case LEQ:
        case LT:
          return interpret_tell_x_y_op_k(f, battery::get<0>(xi), avari, f.seq(0).sig(), avarj, constant, env, tell,
                                         diagnostics);
        case GEQ:
        case GT: {
          auto symbol_x = battery::get<0>(xi) == '-' ? '+' : '-';
          auto arithm_symb = f.seq(0).sig() == ADD ? SUB : ADD;

          if (constant.is(F::Z)) {
            logic_int k = constant.z() * -1;
            return interpret_tell_x_y_op_k<diagnose>(f, symbol_x, avari, arithm_symb, avarj, F::make_z(k), env, tell,
                                                     diagnostics);
          }
          else {
            logic_real k = battery::make_tuple(battery::mul_down(battery::get<0>(constant.r()), -1.0),
                                               battery::mul_up(battery::get<1>(constant.r()), -1.0));
            return interpret_tell_x_y_op_k<diagnose>(f, symbol_x, avari, arithm_symb, avarj, F::make_real(k), env, tell,
                                                     diagnostics);
          }
        }
        case EQ: {
          auto geq_formula = F::make_binary(f.seq(0), GEQ, constant, aty());
          return interpret_tell_x_y_op_k<diagnose>(f, battery::get<0>(xi), avari, f.seq(0).sig(), avarj, constant, env,
                                                   tell,
                                                   diagnostics) && interpret_tell_binary(
                   geq_formula, env, tell, diagnostics);
        }
        default:
          RETURN_INTERPRETATION_ERROR("Octagons only handle <= symbols at the moment.");
      }
    }

    template<bool diagnose = false, class F, class Env, class Alloc2>
    CUDA NI bool interpret_tell_in(const F& f, Env& env, tell_type<Alloc2>& tell, IDiagnostics& diagnostics) const {
      DEBUG_PRINT("interpret_tell_in\n");
      auto xi = f.seq(0);
      auto domain = f.seq(1);
      DEBUG_PRINT("%d \n", domain.is(F::S));
      auto var1 = env.variable_of(xi.lv());
      if (!var1.has_value()) {
        RETURN_INTERPRETATION_ERROR("The variable "+ xi +" does not exist.");
      }
      AVar avar1 = var1.value().get().avar_of(atype).value();
      avar1.type_as(aty());
      auto set = domain.s();
      auto lb = battery::get<0>(set[0]);
      auto ub = battery::get<1>(set[set.size() - 1]);

      interpret_tell_unary(F::make_binary(F::make_avar(avar1), GEQ, lb, aty()), env, tell, diagnostics);
      interpret_tell_unary(F::make_binary(F::make_avar(avar1), LEQ, ub, aty()), env, tell, diagnostics);

      DEBUG_PRINT("interpret_tell_in fin\n");
      return true;
    }

    template<bool diagnose = false, class F, class Env, class Alloc2>
    CUDA NI bool interpret_tell_unary(const F& f, Env& env, tell_type<Alloc2>& tell,
                                      IDiagnostics& diagnostics) const {
      DEBUG_PRINT("interpret_tell_unary\n");

      //we force the constant to the right
      auto local_f = move_constants_on_rhs(f);

      if (!is_arithmetic_comparison(local_f) || f.sig() == NEQ) {
        RETURN_INTERPRETATION_ERROR("Octagons only handle <=,>=,<,>,= symbols at the moment.");
      }
      auto left_operand = local_f.seq(0);
      left_operand.type_as(aty());
      battery::tuple<char, F> xit;
      auto result = var_with_symbol<diagnose, F>(left_operand, diagnostics, xit);
      if (!result) {
        return false;
      }
      auto xi = battery::get<1>(xit);
      auto constant = local_f.seq(1);

      AVar avar1;
      if (xi.is(F::LV)) {
        auto var1 = env.variable_of(xi.lv());
        if (!var1.has_value()) {
          RETURN_INTERPRETATION_ERROR("The variable "+ xi +" does not exist.");
        }
        avar1 = var1.value().get().avar_of(atype).value();
        avar1.type_as(aty());
      }
      else if (xi.is(F::V)) {
        avar1 = xi.v();
        avar1.type_as(aty());
      }
      else {
        assert(false);
      }

      local_cell_type u;
      switch (local_f.sig()) {
        case LEQ:
        case LT:
          return interpret_tell_x_op_k(battery::get<0>(xit), avar1, local_f.sig(), constant, env, tell, diagnostics);
        case EQ: {
          auto geq_formula = F::make_binary(f.seq(0), GEQ, constant, aty());
          return interpret_tell_x_op_k(battery::get<0>(xit), avar1, LEQ, constant, env, tell, diagnostics) &&
                 interpret_tell_unary(geq_formula, env, tell, diagnostics);
        }
        case GEQ:
        case GT: {
          if (constant.is(F::Z)) {
            logic_int k = constant.z() * -1;
            return interpret_tell_x_op_k(battery::get<0>(xit) == '-' ? '+' : '-', avar1, LEQ, F::make_z(k), env, tell,
                                         diagnostics);
          }
          else {
            logic_real k = battery::make_tuple(battery::mul_down(battery::get<0>(constant.r()), -1.0),
                                               battery::mul_up(battery::get<1>(constant.r()), -1.0));
            return interpret_tell_x_op_k(battery::get<0>(xit) == '-' ? '+' : '-', avar1, LEQ, F::make_real(k), env,
                                         tell, diagnostics);
          }
        }
        default:
          RETURN_INTERPRETATION_ERROR("Octagons only handle <=,>=,<,>,= symbols at the moment.");
      }
    }


    template<bool diagnose = false, class F, class Env, class Alloc2>
    CUDA NI bool interpret_tell(const F& f, Env& env, tell_type<Alloc2>& tell, IDiagnostics& diagnostics) const {
      DEBUG_PRINT("start interpret_tell\n");
      switch (auto ftype = get_type(f)) {
        case EXISTENTIAL: //∃x
          nbVars++;
          return interpret_tell_existential(f, env, tell, diagnostics);
        case BINARY: // ±x±y op k
          return interpret_tell_binary(f, env, tell, diagnostics);
        case UNARY: // ±x op k
          return interpret_tell_unary(f, env, tell, diagnostics);
        case SET_IN: // ±x ∈ [lb,ub]
          return interpret_tell_in(f, env, tell, diagnostics);
        default:
          RETURN_INTERPRETATION_ERROR(
            "We can't interpret the formula because it is not an existantial formula or a comparison constraint of the form (±x±y op k, ±x op k). ");
      }
      DEBUG_PRINT("end interpret_tell\n");
    }

    template<bool diagnose = false, class F, class Env, class Alloc2>
    CUDA NI bool interpret_ask(const F& f, const Env& env, ask_type<Alloc2>& ask, IDiagnostics& diagnostics) const {
      /*todo*/
      return interpret_tell(f, const_cast<Env &>(env), ask, diagnostics);
    }

    template<IKind kind, bool diagnose = false, class F, class Env, class I>
    CUDA NI bool interpret(const F& f, Env& env, I& intermediate, IDiagnostics& diagnostics) const {
      DEBUG_PRINT("start interpret octogon\n");
      if constexpr (kind == IKind::TELL) {
        return interpret_tell(f, env, intermediate, diagnostics);
      }
      else {
        return interpret_ask(f, env, intermediate, diagnostics);
      }
      DEBUG_PRINT("end interpret octogon\n");
    }

    template<class Alloc2, class Mem>
    CUDA this_type& tell(const tell_type<Alloc2>& t, BInc<Mem>& has_changed) {
      DEBUG_PRINT("start tell\n");
      DEBUG_PRINTLN("tell type size %ld\n", t.size());
      dbm.resize(nbVars * 2);
      init_size();
      for (int i = 0; i < t.size(); i++) {
        auto el = t[i];
        auto index_line = battery::get<0>(el);
        auto index_column = battery::get<1>(el);

        if (dbm[index_line].empty()) {
          dbm[index_line].resize(nbVars * 2);
        }
        dbm[index_line][index_column].tell(battery::get<2>(el), has_changed);
      }

      DEBUG_PRINT("end tell\n");
      return *this;
    }

    CUDA this_type& tell_top() {
      dbm[0][0].tell_top();
      return *this;
    }

    template<class Alloc2>
    CUDA this_type& tell(const tell_type<Alloc2>& t) {
      local::BInc has_changed;
      return tell(t, has_changed);
    }

    CUDA this_type& tell(AVar x, const universe_type& dom) {
      local::BInc has_changed;
      return tell(x, dom, has_changed);
    }

    template<class Mem>
    CUDA this_type& tell(AVar x, const universe_type& dom, BInc<Mem>& has_changed) {
      dbm[x.vid() * 2][x.vid() * 2 + 1].tell(dual<local_cell_type>(dom.lb()), has_changed);
      dbm[x.vid() * 2 + 1][x.vid() * 2].tell(dom.ub(), has_changed);
      return *this;
    }

    template<class Alloc2>
    CUDA local::BInc ask(const ask_type<Alloc2>& t) const {
      for (int i = 0; i < t.size(); i++) {
        auto el = t[i];
        auto index_line = battery::get<0>(el);
        auto index_column = battery::get<1>(el);
        auto value = battery::get<2>(el);

        if (dbm[index_line][index_column] != value) {
          return false;
        }
      }
      return true;
    }

    CUDA size_t num_refinements() const {
      return _num_refinements;
    }

    template<class Mem>
    CUDA void refine(size_t i, BInc<Mem>& has_changed) {
      //print_matrix(dbm);
      using local_flat = typename U::template flat_type<battery::local_memory>;

      // for(int k=0;k<dbm.size()/2;k++){
      //   for(int l=0;l<dbm.size();l++){
      //     for(int j=0;j<dbm.size();j++) {
      //       dbm[l][j].tell(
      //       U::template fun<ADD>(
      //         local_flat(dbm[l][k]),
      //         local_flat(dbm[k][l])),
      //       has_changed);
      //     }
      //   }
      //
      //   for(int l=0;l<dbm.size();l++) {
      //     for(int j=0;j<dbm.size();j++) {
      //
      //       size_t index_i_bar = l ^ 1;
      //       size_t index_j_bar = j ^ 1;
      //
      //       auto add = local_flat(U::template fun<ADD>(
      //         local_flat(dbm[l][index_i_bar]),
      //         local_flat(dbm[index_j_bar][j])));
      //       dbm[l][j].tell(U::template fun<FDIV>(add, local_flat(2)),has_changed);
      //     }
      //   }
      // }



      if (i < floyd_steps) {
        size_t k = dim1(i);
        size_t ii = dim2(i);
        size_t j = dim3(i);
        dbm[ii][j].tell(
          U::template fun<ADD>(
            local_flat(dbm[ii][k]),
            local_flat(dbm[k][j])),
          has_changed);
        if (ii == j) {
          is_at_top.tell(local::BInc(dbm[ii][j] < 0), has_changed);
        }
        if(i==floyd_steps-1 && !is_top()) {
          for(int l=0;l<dbm.size();l++) {
            dbm[l][l].tell(U(0),has_changed);
          }
        }
      }
      else if (i < tight_steps) {
        size_t index = i % dbm.size();
        size_t index_bar = index ^ 1;
        auto div2 = local_flat(U::template fun<FDIV>(local_flat(dbm[index][index_bar]), local_flat(2)));
        auto value = U::template fun<MUL>(div2, local_flat(2));
        dbm[index][index_bar].tell(value, has_changed);
      }
      else {
        size_t index_i = ((i - tight_steps) / dbm.size());
        size_t index_i_bar = index_i ^ 1;
        size_t index_j = ((i - tight_steps) % dbm.size());
        size_t index_j_bar = index_j ^ 1;

        auto add = local_flat(U::template fun<ADD>(
          local_flat(dbm[index_i][index_i_bar]),
          local_flat(dbm[index_j_bar][index_j])));
        dbm[index_i][index_j].tell(U::template fun<FDIV>(add, local_flat(2)),has_changed);
      }
      //print_matrix(dbm);
    }

    /** `true` if the underlying abstract element is top, `false` otherwise. */
    CUDA local::BInc is_top() const {
      return is_at_top;
    }

    /** `true` if the underlying abstract element is bot and there is no refinement function, `false` otherwise. */
    CUDA local::BDec is_bot() const {
      if (is_at_top) {
        return false;
      }
      for (int i = 0; i < dbm.size(); i++) {
        auto list = dbm[i];
        for (int j = 0; j < list.size(); j++) {
          auto el = list[j];
          if (!el.is_bot()) {
            return false;
          }
        }
      }
      return true;
    }

    CUDA universe_type operator[](int x) const {
      using local_flat = typename U::template flat_type<battery::local_memory>;
      auto result = universe_type(typename universe_type::LB(dbm[x * 2+1][x * 2]),
                                  typename universe_type::UB(dbm[x * 2][x * 2+1]));

      auto lb = result.lb().value();
      auto ub = result.ub().value();

      auto div2_lb = local_flat(U::template fun<FDIV>(local_flat(lb), local_flat(2)));
      auto div2_ub = local_flat(U::template fun<FDIV>(local_flat(ub), local_flat(2)));


      if (div2_lb < 0) {
        div2_lb = U::template fun<MUL>(local_flat(div2_lb), local_flat(-1));
      }

      if (div2_ub < 0) {
        div2_ub = U::template fun<MUL>(local_flat(div2_ub), local_flat(-1));
      }

      return universe_type(typename universe_type::LB(div2_lb),
                           typename universe_type::UB(div2_ub));
    }

    CUDA universe_type project(AVar x) const {
      return (*this)[x.vid()];
    }

    CUDA size_t vars() const {
      if (dbm.empty()) {
        return 0;
      }
      return dbm.size() / 2;
    }

    template<class Alloc2 = allocator_type>
    CUDA snapshot_type<Alloc2> snapshot(const Alloc2& alloc = Alloc2()) const {
      return snapshot_type<Alloc2>(dbm, alloc);
    }

    template<class Alloc2>
    CUDA void restore(const snapshot_type<Alloc2>& snap) {
      is_at_top.dtell_bot();
      for (int i = 0; i < snap.size(); ++i) {
        for (int j = 0; j < snap[i].size(); ++j) {
          dbm[i][j].dtell(snap[i][j]);
          if (snap[i][j].value() < 0) {
            is_at_top.tell(local::BInc(true));
          }
        }
      }
    }

    /** An abstract element is extractable when it is not equal to top, the refinement is at a fixpoint and the underlying abstract elements are extractable. */
    template<class ExtractionStrategy = NonAtomicExtraction>
    CUDA bool is_extractable(const ExtractionStrategy& strategy = ExtractionStrategy()) const {
      return true;
    }

    /** Extract the current element into `ua`.
     * \pre `is_extractable()` must be `true`. */
    template<class B>
    CUDA void extract(B& ua) const {
      if ((void *) &ua != (void *) this) {
        ua.dbm = dbm;
      }
    }


    CUDA void print() const {
      print_matrix(dbm);
    }


    template<class Env>
    CUDA NI TFormula<typename Env::allocator_type> deinterpret(const Env& env) const {
      using F = TFormula<typename Env::allocator_type>;
      typename F::Sequence seq{env.get_allocator()};
      using local_flat = typename U::template flat_type<battery::local_memory>;

      if (is_top()) {
        return TFormula<typename Env::allocator_type>::make_false();
      }

      for (int i = 0; i < vars(); i++) {
        AVar v(aty(), i);
        seq.push_back(F::make_exists(aty(), env.name_of(v), env.sort_of(v)));

        auto p = project(v);

        auto lb = p.lb().value();
        auto ub = p.ub().value();


        auto var = F::make_lvar(aty(), env.name_of(v));
        seq.push_back(F::make_binary(var, GEQ, F::make_z(lb), aty()));
        seq.push_back(F::make_binary(var, LEQ, F::make_z(ub), aty()));
      }


      for (int i = 0; i < dbm.size(); i++) {
        for (int j = 0; j < dbm[i].size(); j++) {
          if (i != j && !dbm[i][j].is_bot()) {
            AVar v1(aty(), i % 2 != 0 ? (i - 1) / 2 : i / 2);
            AVar v2(aty(), j % 2 != 0 ? (j - 1) / 2 : j / 2);
            if (v1 == v2) {
              continue;
            }

            Sig arithmSymb = SUB;
            if (i % 2 == 0) {
              arithmSymb = ADD;
            }
            auto var1 = F::make_lvar(aty(), env.name_of(v1));
            auto var2 = F::make_lvar(aty(), env.name_of(v2));

            if (j % 2 != 0) {
              seq.push_back(F::make_binary(F::make_binary(var2, arithmSymb, var1, aty()), LEQ,
                                           F::make_z(dbm[i][j]), aty()));
            }
            else {
              seq.push_back(F::make_binary(
                F::make_binary(F::make_unary(NEG, var2, aty()), arithmSymb, var1, aty()), LEQ,
                F::make_z(dbm[j][i]), aty()));
            }
          }
        }
      }

      return F::make_nary(AND, std::move(seq), aty());
    }

    template<class I, class Env>
    CUDA NI TFormula<typename Env::allocator_type> deinterpret(const I& intermediate, const Env& env) const {
      using F = TFormula<typename Env::allocator_type>;
      return F::make_false();
    }
  };
}

#endif
