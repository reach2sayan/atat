# ATAT Modernization Plan

## Current State

ATAT already configures `CMAKE_CXX_STANDARD 17` in [CMakeLists.txt](C:\atat\CMakeLists.txt), but much of the codebase still follows pre-standard-library design patterns:

- Custom containers in [include/array.h](C:\atat\include\array.h), [include/linklist.h](C:\atat\include\linklist.h), and [include/arraylist.h](C:\atat\include\arraylist.h)
- A custom string wrapper in [include/stringo.h](C:\atat\include\stringo.h)
- Macro-heavy utilities and global aliases in [include/misc.h](C:\atat\include\misc.h) and [include/machdep.h](C:\atat\include\machdep.h)
- Widespread manual ownership with raw `new` / `delete`
- Legacy stream utilities such as `strstream` still used in files like [src/strainder.c++](C:\atat\src\strainder.c++) and [src/viewgce.c++](C:\atat\src\viewgce.c++)

The biggest modernization opportunity is not the language standard itself; it is reducing custom infrastructure and manual memory management.

## Goals

1. Preserve numerical behavior and file formats.
2. Replace custom containers with STL equivalents where practical.
3. Introduce modern ownership semantics (`std::unique_ptr`, values, references, spans/views where appropriate).
4. Use Boost only where it clearly improves maintainability over bespoke code.
5. Keep the migration incremental so existing executables can keep building during the transition.

## Priority Refactor Targets

### 1. Replace custom string handling

Current target:

- [include/stringo.h](C:\atat\include\stringo.h)
- [src/stringo.c++](C:\atat\src\stringo.c++)

Desired replacements:

- `AutoString` -> `std::string`
- C string manipulation -> `std::string`, `std::string_view`, `std::ostringstream`

Why first:

- Low algorithmic risk
- Removes manual buffer allocation
- Unlocks simpler parsing and logging code across the project

### 2. Replace custom linked list ownership patterns

Current target:

- [include/linklist.h](C:\atat\include\linklist.h)
- Call sites in [include/refine.h](C:\atat\include\refine.h), [include/mrefine.h](C:\atat\include\mrefine.h), [src/refine.c++](C:\atat\src\refine.c++), [src/tlambda.c++](C:\atat\src\tlambda.c++)

Desired replacements:

- `LinkedList<T>` storing owning raw pointers -> `std::list<T>`, `std::vector<T>`, or `std::deque<T>` depending on access pattern
- If polymorphism is required: `std::list<std::unique_ptr<T>>`

Why second:

- Highest source of ownership ambiguity
- Converts many `new` / `delete detach(...)` patterns into RAII

### 3. Replace `Array<T>` / `Array2d<T>` with STL-backed storage

Current target:

- [include/array.h](C:\atat\include\array.h)
- Consumers throughout `src/`

Desired replacements:

- One-dimensional dynamic buffers -> `std::vector<T>`
- Fixed-size math tuples -> `std::array<T, N>` where dimensions are static
- Multi-dimensional storage -> `std::vector<T>` plus indexing helpers, or Eigen where matrix semantics already exist

Migration advice:

- Do this in layers. First make `Array` a thin wrapper over `std::vector<T>` internally if needed.
- Avoid a flag day rewrite across the whole project.

### 4. Remove unsafe macros and aliases

Current target:

- [include/misc.h](C:\atat\include\misc.h)

Desired replacements:

- `MAX`, `MIN`, `ABS` -> `std::max`, `std::min`, `std::abs`
- `swap(T*, T*)` -> `std::swap`
- `typedef` aliases -> `using`
- memory macros -> STL algorithms or standard library facilities

Why now:

- Macro collisions and hidden side effects make later refactors harder

### 5. Replace deprecated/legacy stream headers

Current target:

- [src/strainder.c++](C:\atat\src\strainder.c++)
- [src/viewgce.c++](C:\atat\src\viewgce.c++)

Desired replacements:

- `strstream`, `ostrstream`, `istrstream` -> `std::stringstream`, `std::ostringstream`, `std::istringstream`

## Boost Guidance

Boost should be introduced selectively, not as a blanket replacement for STL.

Recommended Boost candidates:

- `Boost.Program_options` for complex command-line parsing if current parsing utilities become a maintenance bottleneck
- `Boost.Filesystem` only if `std::filesystem` is not sufficient for cross-platform behavior already needed by ATAT
- `Boost.MultiArray` or `Boost.uBLAS` only if a specific multidimensional container problem is not already better served by Eigen or `std::vector`

Avoid introducing Boost for:

- Smart pointers
- Strings
- Containers already covered well by the STL
- Features already handled by Eigen in numerical code

## Suggested Migration Phases

### Phase 1: Safety and tooling

- Turn on modernization diagnostics
- Remove deprecated headers like `strstream`
- Replace `NULL` with `nullptr`
- Convert `typedef` to `using` in touched files

### Phase 2: Strings and utilities

- Replace `AutoString`
- Replace common macros in [include/misc.h](C:\atat\include\misc.h)
- Introduce helper functions in namespaces instead of global macros

### Phase 3: Ownership cleanup

- Migrate owning raw pointers to `std::unique_ptr`
- Replace linked-list-based ownership patterns with standard containers

### Phase 4: Container migration

- Move high-value modules from `Array<T>` to `std::vector<T>`
- Keep compatibility shims where needed to limit churn

### Phase 5: API cleanup

- Use `const std::string&`, `std::span`, iterators, and standard algorithms
- Replace hand-written loops with STL algorithms when they improve clarity

## Recommended First Slice

The safest first implementation slice is:

1. Replace `strstream` usage with stringstream equivalents.
2. Modernize [include/stringo.h](C:\atat\include\stringo.h) away from manual `char*`.
3. Add RAII to a single ownership-heavy subsystem, likely `refine` or `tlambda`.

This gives a real improvement without forcing an immediate rewrite of every `Array<T>` call site.

## Success Criteria

- No new owning raw pointers introduced in touched code
- New code prefers STL containers and algorithms
- Boost is added only with a concrete use case
- Each phase leaves the project buildable
