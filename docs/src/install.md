# OpenSOAP Installation Guide
## New to Julia?
If you are new to the Julia programming language, follow these steps to install the language:
1. Install Julia for your OS [from this page](https://julialang.org/downloads/).
   a. There's a good list of [language tutorials here](https://julialang.org/learning/tutorials/).
2. Try running the Julia Read/Evaluate/Print Loop (REPL) from your terminal like this:
```bash
julia
```

## Setting up this project
Find a folder on your computer you want to store this work in. Navigate there in the terminal. When you've arrived, do:
```bash
git clone --recursive https://github.com/thanasipantazides/OpenSOAP.git
cd OpenSOAP
```

From here, you should start Julia using the `OpenSOAP` project environment:
```bash
julia --project=.
```
Then, at the Julia REPL, press `]` to enter the package manager, and activate the local environment:
```julia
julia> ]
pkg> activate .
```
Now install all required packages:
```julia
pkg> instantiate
```
This may take a few minutes. When complete, press backspace/delete to return to the main Julia REPL.

Once complete, try running the test simulation and displaying plots with:
```julia
julia> include("test/plot.jl")
```

## Building documentation
Navigate to the docs folder, and start the local Julia environment:
```bash
cd OpenSOAP/docs
julia --project=.
```

Then, at the Julia REPL, press `]` to enter the package manager, and activate the local environment:
```julia
julia> ]
pkg> activate .
```

Include the `OpenSOAP` package as a dependency:
```julia
pkg> dev ..
```
Press backspace/delete to exit the package manager and return to the main Julia REPL.
Finally, build the documentation like this:
```julia
julia> include("make.jl")
```

The PDF document will appear in `OpenSOAP/docs/build/`.
