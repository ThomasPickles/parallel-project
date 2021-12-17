# Parallel and Distributed Algorithms and Programs (APPD)
Programming Project
LAB6 (code generation for functions), APPD 2021-22

# Author

Thomas Pickles

# Contents

./mst-solution.c - source code file
./tests/data-*.txt -
./tests/out-*.txt -


# Howto

`make` to build executable
`make run` to run code on a basic example
`make tests` to run each algorithm against expected output for a range of system sizes and processor numbers.
`make performance` to evaluate algorithm performance against a range of criteria.

# Dependencies

mpi
smpi

# Test design

All provided tests pass.

# Design choices

new_offset() function generates a negative offset relative to some base.
The course notes suggest using a positive offset relative to stack pointer
to save registers.  To avoid coding another function to generate a positive
offset I have set the memory locations
on the stack relative to frame pointer rather than stack pointer.

# Known bugs and limitations

gcc works with 32bits integer, when riscv uses 64bits integer
Differences in output when result > 2^31.  Accordingly, when calculating 15! using
a recursive function, MiniC produces the correct result where gcc creates the
wrong result

Compatibility with gcc is not thoroughly tested.  I have provided compatibility
tests with gcc in both directions, but these are for very simple functions, so
this test coverage is currently very basic.  With more time, I would focus on
extending this test suite.

# Checklists

A check ([X]) means that the feature is implemented
and *tested* with appropriate test cases.

## Parser

- [x] Function definition
- [x] Function declaration
- [x] Function call

## Typer

- [x] Function declaration
- [x] Function definition
- [x] Function call
- [x] Function return

## Code generation

- [x] Function return
- [x] Callee-saved registers
- [x] Function call
- [x] Getting the result of a function call
- [x] Caller-saved registers
- [x] Increase the size of the stack for callee/caller-saved registers
- [x] Temporaries for giving arguments to a function call
- [x] Temporaries for retriving arguments at the beginning of a function


