# Connected Information 
This project provides a way to compute the approximate connected information and entropy of any number of random distributions.

## Installation
To install a version from git simply run the commands below in Julia REPL.
> **Note**
> Only symbols after `>` should be pasted to command line.

```shell
julia> ]
(@v1.8) pkg> dev https://github.com/ann-ib/ConnectedInformation.jl.git
```

Start using the package in Julia code:
```julia
using ConnectedInformation
```

## How to run tests
After installing the package it is possible to run tests with the `test` command:
```shell
(@v1.8) pkg> test ConnectedInformation
```

## How to update the package
It is possible to update the project with recent changes with the command:
```shell
(@v1.8) pkg> up --preserve=direct ConnectedInformation
```

## License
Connected Information project is released under the [MIT License](LICENSE.txt).
