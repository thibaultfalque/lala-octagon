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
#include "lala/abstract_deps.hpp"
#include "lala/vstore.hpp"
#include "lala/interval.hpp"
#include "lala/logic/env.hpp"

namespace lala {
    /** Octagon is an abstract domain built on top of an abstract universe `U`. */
    template<class U, class Allocator>
    class Octagon {
    public:
        using universe_type = Interval<U>;
        using local_universe_type = typename universe_type::local_type;
        using local_cell_type = typename U::local_type;
        using allocator_type = Allocator;
        using this_type = Octagon<U, allocator_type>;
        using universe_list_type = battery::vector<U>;
        using dbm_type = battery::vector<universe_list_type>;

        template<class Alloc>
        struct snapshot_type {
        };


        constexpr static const bool is_abstract_universe = false;
        constexpr static const bool sequential = universe_type::sequential;
        constexpr static const bool is_totally_ordered = false;
        constexpr static const bool preserve_bot = true;
        constexpr static const bool preserve_top = universe_type::preserve_top;
        constexpr static const bool preserve_join = universe_type::preserve_join;
        constexpr static const bool preserve_meet = universe_type::preserve_meet;
        constexpr static const bool injective_concretization = universe_type::injective_concretization;
        constexpr static const bool preserve_concrete_covers = universe_type::preserve_concrete_covers;
        constexpr static const char *name = "Octagon";

        template<class U2, class Alloc2>
        friend class Octagon;

    private:
        AType atype;
        dbm_type dbm;
        mutable int nbVars = 0;
        long floyd_steps;
        long tight_steps;
        long str_steps;
        long _num_refinements;

        enum FormulaEnum {
            EXISTENTIAL, UNARY, BINARY, SET_IN, UNKNOWN
        };


        CUDA INLINE size_t dim1(size_t i) { return i / (dbm.size() * dbm.size()); }
        CUDA INLINE size_t dim2(size_t i) { return i % (dbm.size() * dbm.size()) / dbm.size(); }
        CUDA INLINE size_t dim3(size_t i) { return i % (dbm.size() * dbm.size()) % dbm.size(); }

        template<class F>
        CUDA NI FormulaEnum get_type(const F &f) const {
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
        [[nodiscard]] bool is_arithmetic_comparison_operator(const Sig &sig) const {
            return (sig == LEQ || sig == GEQ);
        }

        template<bool diagnose = false, class F>
        bool var_with_symbol(F &f, IDiagnostics diagnostics, battery::tuple<char, F> &result) const {
            if (f.is_variable()) {
                result = battery::make_tuple('+', f);
                return true;
            }
            if (f.is(F::Seq) && f.seq().size() == 1 && f.sig() == NEG) {
                result = battery::make_tuple('-', f.seq(0));
                return true;
            }
            RETURN_INTERPRETATION_ERROR("There is an error when we interpret a variable formula");
        }

        void init_size() {
            floyd_steps = dbm.size() * dbm.size() * dbm.size();
            tight_steps = floyd_steps + dbm.size();
            str_steps = tight_steps+dbm.size()*dbm.size();
            _num_refinements =  str_steps;
        }

    public:
        template<class Alloc>
        using tell_type = battery::vector<battery::tuple<int, int, U>, Alloc>;

        template<class Alloc>
        using ask_type = battery::vector<battery::tuple<int, int, U>, Alloc>; /*TODO just for compile*/


        CUDA Octagon(const this_type &other)
            : atype(other.atype), dbm(other.dbm) {
            init_size();
        }

        /** Initialize an empty store. */
        CUDA Octagon(AType atype, const allocator_type &alloc = allocator_type())
            : atype(atype), dbm(alloc) {
            init_size();
        }

        CUDA Octagon(AType atype, size_t size, const allocator_type &alloc = allocator_type())
            : atype(atype), dbm(size, alloc) {
            init_size();
        }

        template<class R>
        CUDA Octagon(const Octagon<R, allocator_type> &other)
            : atype(other.atype), dbm(other.dbm) {
            init_size();
        }

        template<class R, class Alloc2>
        CUDA Octagon(const Octagon<R, Alloc2> &other, const allocator_type &alloc = allocator_type())
            : atype(other.atype), dbm(other.dbm, alloc) {
            init_size();
        }

        /** Copy the Octagon `other` in the current element.
         *  `deps` can be empty and is not used besides to get the allocator (since this abstract domain does not have dependencies). */
        template<class R, class Alloc2, class... Allocators>
        CUDA Octagon(const Octagon<R, Alloc2> &other, const AbstractDeps<Allocators...> &deps)
            : Octagon(other, deps.template get_allocator<allocator_type>()) {
            init_size();
        }

        CUDA Octagon(this_type &&other): atype(other.atype), dbm(std::move(other.dbm)) {
            init_size();
        }

        CUDA allocator_type get_allocator() const {
            return dbm.get_allocator();
        }

        CUDA AType aty() const {
            return atype;
        }

        CUDA static this_type bot(AType atype = UNTYPED,
                                  const allocator_type &alloc = allocator_type()) {
            return Octagon{atype, alloc};
        }

        /** A special symbolic element representing top. */
        CUDA static this_type top(AType atype = UNTYPED,
                                  const allocator_type &alloc = allocator_type()) {
        }

        template<class Env>
        CUDA static this_type bot(Env &env,
                                  const allocator_type &alloc = allocator_type()) {
            AType atype = env.extends_abstract_dom();
            return bot(atype, alloc);
        }

        template<class Env>
        CUDA static this_type top(Env &env,
                                  const allocator_type &alloc = allocator_type()) {
        }


        template<bool diagnose = false, class F, class Env, class Alloc2>
        CUDA NI bool interpret_tell_existential(const F &f, Env &env, tell_type<Alloc2> &tell,
                                                IDiagnostics &diagnostics) const {
            DEBUG_PRINT("start\n");
            AVar avar;
            auto result = env.template interpret<diagnose>(f.map_atype(aty()), avar, diagnostics);
            DEBUG_PRINT("end\n");
            return result;
        }

        template<bool diagnose = false, class F, class Env, class Alloc2>
        CUDA NI bool interpret_tell_binary(const F &f, Env &env, tell_type<Alloc2> &tell,
                                           IDiagnostics &diagnostics) const {
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
            thrust::optional<const Variable<Allocator> &> var1 = env.variable_of(battery::get<1>(xi).lv());
            if (!var1.has_value()) {
                RETURN_INTERPRETATION_ERROR("The variable "+xi+" not exists.");
            }
            thrust::optional<const Variable<Allocator> &> var2 = env.variable_of(battery::get<1>(xj).lv());
            if (!var2.has_value()) {
                RETURN_INTERPRETATION_ERROR("The variable "+xj+" not exists.");
            }
            if (!constant.is_constant()) {
                RETURN_INTERPRETATION_ERROR("The right operand of the formula must be a constant.");
            }
            AVar avari = var1.value().avar_of(atype).value();
            AVar avarj = var2.value().avar_of(atype).value();
            local_cell_type u;
            int index1 = -1;
            int index2 = -1;
            auto tmp = F::make_avar(AVar{});

            local_cell_type::template interpret_tell<diagnose>(F::make_binary(tmp, f.sig(), constant), env, u,
                                                               diagnostics);

            switch (f.sig()) {
                case LEQ:
                    switch (f.seq(0).sig()) {
                        case SUB:
                            if (battery::get<0>(xi) == '-') {
                                index1 = avari.vid() * 2 + 1;
                                index2 = avarj.vid() * 2;
                            } else {
                                index1 = avari.vid() * 2;
                                index2 = avarj.vid() * 2;
                            }
                            break;
                        case ADD:
                            if (battery::get<0>(xi) == '-') {
                                index1 = avari.vid() * 2 + 1;
                                index2 = avarj.vid() * 2 + 1;
                            } else {
                                index1 = avari.vid() * 2;
                                index2 = avarj.vid() * 2 + 1;
                            }
                            break;
                        default:
                            RETURN_INTERPRETATION_ERROR("Octagons only handle + and - arithmitic symbols.");
                    }
                    break;
                default:
                    RETURN_INTERPRETATION_ERROR("Octagons only handle <= symbols at the moment.");
            }
            tell.emplace_back(battery::make_tuple(index1, index2, u));
            DEBUG_PRINT("end\n");
            return true;
        }

        template<bool diagnose = false, class F, class Env, class Alloc2>
        CUDA NI bool interpret_tell_in(const F &f, Env &env, tell_type<Alloc2> &tell, IDiagnostics &diagnostics) const {
            DEBUG_PRINT("interpret_tell_in\n");
            f.print();
            printf("\n");
            auto xi = f.seq(0);
            auto domain = f.seq(1);
            DEBUG_PRINT("%d \n", domain.is(F::S));
            thrust::optional<const Variable<Allocator> &> var1 = env.variable_of(xi.lv());
            if (!var1.has_value()) {
                RETURN_INTERPRETATION_ERROR("The variable "+xi+" not exists.");
            }
            AVar avar1 = var1.value().avar_of(atype).value();
            auto set = domain.s();
            auto lb = battery::get<0>(set[0]);
            auto ub = battery::get<1>(set[set.size() - 1]);

            interpret_tell_unary(F::make_binary(F::make_avar(avar1), GEQ, lb), env, tell, diagnostics);
            interpret_tell_unary(F::make_binary(F::make_avar(avar1), LEQ, ub), env, tell, diagnostics);

            DEBUG_PRINT("interpret_tell_in fin\n");
            return true;
        }

        template<bool diagnose = false, class F, class Env, class Alloc2>
        CUDA NI bool interpret_tell_unary(const F &f, Env &env, tell_type<Alloc2> &tell,
                                          IDiagnostics &diagnostics) const {
            DEBUG_PRINT("interpret_tell_unary\n");
            auto xi = f.seq(0);

            auto constant = f.seq(1);

            AVar avar1;
            if (xi.is(F::LV)) {
                thrust::optional<const Variable<Allocator> &> var1 = env.variable_of(xi.lv());
                if (!var1.has_value()) {
                    RETURN_INTERPRETATION_ERROR("The variable "+xi+" not exists.");
                }
                avar1 = var1.value().avar_of(atype).value();
            } else if (xi.is(F::V)) {
                avar1 = xi.v();
            }
            if (!constant.is_constant()) {
                RETURN_INTERPRETATION_ERROR("The right operand of the formula must be a constant.");
            }

            local_cell_type u;
            auto tmp = F::make_avar(AVar{});
            switch (f.sig()) {
                case LEQ:
                    if (constant.is(F::Z)) {
                        logic_int k = constant.z() * 2;
                        U::template interpret_tell<diagnose>(F::make_binary(tmp, f.sig(), F::make_z(k)), env, u,
                                                             diagnostics);
                    } else {
                        logic_real k = battery::make_tuple(battery::get<0>(constant.r()) * 2.0,
                                                           battery::get<1>(constant.r()) * 2.0);
                        U::template interpret_tell<diagnose>(F::make_binary(tmp, f.sig(), F::make_real(k)), env, u,
                                                             diagnostics);
                    }
                    tell.emplace_back(battery::make_tuple(avar1.vid() * 2, (avar1.vid() * 2) + 1, u));
                    break;
                case GEQ:
                    if (constant.is(F::Z)) {
                        logic_int k = constant.z() * -2;
                        U::template interpret_tell<diagnose>(F::make_binary(tmp, LEQ, F::make_z(k)), env, u,
                                                             diagnostics);
                    } else {
                        //logic_real k = constant.r()*-2.0;
                        logic_real k = battery::make_tuple(battery::get<0>(constant.r()) * -2.0,
                                                           battery::get<1>(constant.r()) * -2.0);
                        U::template interpret_tell<diagnose>(F::make_binary(tmp, LEQ, F::make_real(k)), env, u,
                                                             diagnostics);
                    }
                    tell.emplace_back(battery::make_tuple(avar1.vid() * 2 + 1, avar1.vid() * 2, u));
                    break;
                default:
                    RETURN_INTERPRETATION_ERROR("Octagons only handle <= and >= symbols at the moment.");
            }
            return true;
        }


        template<bool diagnose = false, class F, class Env, class Alloc2>
        CUDA NI bool interpret_tell(const F &f, Env &env, tell_type<Alloc2> &tell, IDiagnostics &diagnostics) const {
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
        CUDA NI bool interpret_ask(const F &f, const Env &env, ask_type<Alloc2> &ask, IDiagnostics &diagnostics) const {
            /*todo*/
            return true;
        }

        template<IKind kind, bool diagnose = false, class F, class Env, class I>
        CUDA NI bool interpret(const F &f, Env &env, I &intermediate, IDiagnostics &diagnostics) const {
            DEBUG_PRINT("start interpret octogon\n");
            if constexpr (kind == IKind::TELL) {
                return interpret_tell(f, env, intermediate, diagnostics);
            } else {
                return interpret_ask(f, env, intermediate, diagnostics);
            }
            DEBUG_PRINT("end interpret octogon\n");
        }

        template<class Alloc2, class Mem>
        CUDA this_type &tell(const tell_type<Alloc2> &t, BInc<Mem> &has_changed) {
            DEBUG_PRINT("start tell\n");
            DEBUG_PRINTLN("tell type size %ld\n", t.size());
            dbm.resize(nbVars * 2);
            for (int i = 0; i < t.size(); i++) {
                auto el = t[i];
                auto index_line = battery::get<1>(el);
                auto index_column = battery::get<0>(el);

                if (dbm[index_line].empty())
                    dbm[index_line].resize(nbVars * 2);
                dbm[index_line][index_column].tell(battery::get<2>(el), has_changed);

                printf("%d %d \n", index_line, index_column);
                battery::get<2>(el).print();
                printf("\n");

                print_matrix(dbm);
            }

            DEBUG_PRINT("end tell\n");
            return *this;
        }

        template<class Alloc2>
        CUDA this_type &tell(const tell_type<Alloc2> &t) {
            local::BInc has_changed;
            return tell(t, has_changed);
        }

        CUDA this_type &tell(AVar x, const universe_type &dom) {
            local::BInc has_changed;
            return tell(x, dom, has_changed);
        }

        template<class Mem>
        CUDA this_type &tell(AVar x, const universe_type &dom, BInc<Mem> &has_changed) {
            dbm[x.vid() * 2][x.vid() * 2 + 1].tell(dom.lb(), has_changed);
            dbm[x.vid() * 2 + 1][x.vid() * 2].tell(dom.ub(), has_changed);
            return *this;
        }

        template<class Alloc2>
        CUDA local::BInc ask(const ask_type<Alloc2> &t) const {
            return false;
        }

        CUDA size_t num_refinements() const {
            return _num_refinements;
        }

        template<class Mem>
        CUDA void refine(size_t i, BInc<Mem> &has_changed) {
            using local_flat = typename U::template flat_type<battery::local_memory>;
            if(i<floyd_steps) {
                dbm[dim2(i)][dim3(i)].tell(
                    U::template fun<ADD>(
                        local_flat(dbm[dim2(i)][dim1(i)]),
                        local_flat(dbm[dim1(i)][dim3(i)])),
                    has_changed);
            }else if(i<tight_steps) {
                size_t index = i % dbm.size();
                size_t index_bar = i ^ 1;
                auto div2 = local_flat(U::template fun<FDIV>(local_flat(dbm[index][index_bar]),local_flat(2)));
                auto value = U::template fun<MUL>(div2,local_flat(2));
                dbm[index][index_bar].tell(value, has_changed);
            }else{
                auto n_2 = dbm.size()*dbm.size();
                size_t index_i = ((i - tight_steps) % n_2) / dbm.size();
                size_t index_i_bar = index_i ^ 1 ;
                size_t index_j = ((i - tight_steps) / n_2) % dbm.size();
                size_t index_j_bar = index_j ^ 1;

                auto add = local_flat(U::template fun<ADD>(
                    local_flat(dbm[index_i][index_i_bar]),
                    local_flat(dbm[index_j_bar][index_j])));
                dbm[index_i][index_j].tell(U::template fun<FDIV>(add,local_flat(2)));

            }
        }

        /** `true` if the underlying abstract element is top, `false` otherwise. */
        CUDA local::BInc is_top() const {
            for (int i = 0; i < dbm.size(); i++) {
                auto list = dbm[i];
                for (int j = 0; j < list.size(); j++) {
                    auto el = list[j];
                    if (!el.is_top()) {
                        return false;
                    }
                }
            }
            return true;
        }

        /** `true` if the underlying abstract element is bot and there is no refinement function, `false` otherwise. */
        CUDA local::BDec is_bot() const {
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
            return universe_type(typename universe_type::LB(dbm[x * 2][x * 2 + 1]),
                                 typename universe_type::UB(dbm[x * 2 + 1][x * 2]));
        }

        CUDA universe_type project(AVar x) const {
            return this[x.vid()];
        }


        CUDA size_t vars() const {
            if (dbm.empty()) {
                return 0;
            }
            return dbm.size() / 2;
        }

        template<class Alloc2 = allocator_type>
        CUDA snapshot_type<Alloc2> snapshot(const Alloc2 &alloc = Alloc2()) const {
        }

        template<class Alloc2>
        CUDA void restore(const snapshot_type<Alloc2> &snap) {
        }

        /** An abstract element is extractable when it is not equal to top, the refinement is at a fixpoint and the underlying abstract elements are extractable. */
        template<class ExtractionStrategy = NonAtomicExtraction>
        CUDA bool is_extractable(const ExtractionStrategy &strategy = ExtractionStrategy()) const {
            return true;
        }

        /** Extract the current element into `ua`.
         * \pre `is_extractable()` must be `true`. */
        template<class B>
        CUDA void extract(B &ua) const {
        }

        template<class Env>
        CUDA NI TFormula<typename Env::allocator_type> deinterpret(const Env &env) const {
        }
    };
}

#endif
