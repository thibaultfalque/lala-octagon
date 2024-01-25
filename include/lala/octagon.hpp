// Copyright 2024 Pierre Talbot

#ifndef LALA_OCTAGON_HPP
#define LALA_OCTAGON_HPP

#include "battery/vector.hpp"
#include "battery/unique_ptr.hpp"
#include "battery/shared_ptr.hpp"
#include "battery/allocator.hpp"

#include "lala/logic/logic.hpp"
#include "lala/universes/primitive_upset.hpp"
#include "lala/abstract_deps.hpp"

namespace lala {

/** Octagon is an abstract domain built on top of an abstract universe `U`. */
template <class U, class Allocator>
class Octagon {
public:
  using universe_type = U;
  using local_universe_type = typename universe_type::local_type;
  using allocator_type = Allocator;
  using this_type = Octagon<universe_type, allocator_type>;

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

  template <class U2, class Alloc2>
  friend class Octagon;

private:
  AType atype;

public:
  template <class Alloc>
  using tell_type = /* todo */void;

  template <class Alloc>
  using ask_type = /* todo */void;

  /** TODO constructors */

  CUDA allocator_type get_allocator() const {}
  CUDA AType aty() const {}

  CUDA static this_type bot(AType atype = UNTYPED,
    const allocator_type& alloc = allocator_type())
  {}

  /** A special symbolic element representing top. */
  CUDA static this_type top(AType atype = UNTYPED,
    const allocator_type& alloc = allocator_type())
  {}

  template <class Env>
  CUDA static this_type bot(Env& env,
    const allocator_type& alloc = allocator_type())
  {}

  template <class Env>
  CUDA static this_type top(Env& env,
    const allocator_type& alloc = allocator_type())
  {}

public:
  template <bool diagnose = false, class F, class Env, class Alloc2>
  CUDA NI bool interpret_tell(const F& f, Env& env, tell_type<Alloc2>& tell, IDiagnostics& diagnostics) const {}

  template <bool diagnose = false, class F, class Env, class Alloc2>
  CUDA NI bool interpret_ask(const F& f, const Env& env, ask_type<Alloc2>& ask, IDiagnostics& diagnostics) const {}

  template <IKind kind, bool diagnose = false, class F, class Env, class I>
  CUDA NI bool interpret(const F& f, Env& env, I& intermediate, IDiagnostics& diagnostics) const {}

  template <class Alloc2, class Mem>
  CUDA this_type& tell(const tell_type<Alloc2>& t, BInc<Mem>& has_changed) {
    return *this;
  }

  template <class Alloc2>
  CUDA this_type& tell(const tell_type<Alloc2>& t) {
    local::BInc has_changed;
    return tell(t, has_changed);
  }

  CUDA this_type& tell(AVar x, const universe_type& dom) {
    return *this;
  }

  template <class Mem>
  CUDA this_type& tell(AVar x, const universe_type& dom, BInc<Mem>& has_changed) {
    return *this;
  }

  template <class Alloc2>
  CUDA local::BInc ask(const ask_type<Alloc2>& t) const {
    return false;
  }

  CUDA size_t num_refinements() const {
    return 0;
  }

  template <class Mem>
  CUDA void refine(size_t i, BInc<Mem>& has_changed) {
  }

  /** `true` if the underlying abstract element is top, `false` otherwise. */
  CUDA local::BInc is_top() const {
    return false;
  }

  /** `true` if the underlying abstract element is bot and there is no refinement function, `false` otherwise. */
  CUDA local::BDec is_bot() const {
    return false;
  }

  CUDA const universe_type& operator[](int x) const {
  }

  CUDA const universe_type& project(AVar x) const {
  }

  CUDA size_t vars() const {
  }

  template <class Alloc2 = allocator_type>
  CUDA snapshot_type<Alloc2> snapshot(const Alloc2& alloc = Alloc2()) const {
  }

  template <class Alloc2>
  CUDA void restore(const snapshot_type<Alloc2>& snap) {
  }

  /** An abstract element is extractable when it is not equal to top, the refinement is at a fixpoint and the underlying abstract elements are extractable. */
  template <class ExtractionStrategy = NonAtomicExtraction>
  CUDA bool is_extractable(const ExtractionStrategy& strategy = ExtractionStrategy()) const {
  }

  /** Extract the current element into `ua`.
   * \pre `is_extractable()` must be `true`. */
  template <class B>
  CUDA void extract(B& ua) const {
  }

  template<class Env>
  CUDA NI TFormula<typename Env::allocator_type> deinterpret(const Env& env) const {
  }
};

}

#endif
