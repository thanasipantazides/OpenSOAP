# OpenSOAP
The Spacecraft Operations Analysis Package provides some tools for performing system-level analyses of spacecraft mission operations. It could be useful for investigating power and data budgets, and playing with different ConOps designs. 

## Installation
See the installation guide included in the [docs package](docs/src/install.md).

## Documentation
Documentation is included in the [docs folder](docs/). 

## Usage
After [installing](docs/src/install.md), navigate to the OpenSOAP directory, and start Julia:
```bash
cd OpenSOAP/
julia --project=.
```

Setup for development (with the Revise.jl package), then compile the [`test/plot.jl`](test/plot.jl) script. After it finishes compiling (may take a while the first time, will be much faster after that), run `plot_main()` to see the GUI:
```julia
julia> using Revise
julia> include("test/plot.jl")
...
julia> plot_main()
```

When you run this, the dynamics will take 10s of seconds to integrate before plots are displayed.

The display is driven by the [Makie plotting library](https://docs.makie.org/stable/).

After running, you can navigate the 3D plot:
- left-click + drag to orbit camera
- right-click to pan
- scroll wheel to zoom

And detail 2D plots:
- scroll to zoom
- scroll + `X` key to zoom x-axis only
- scroll + `Y` key to zoom y-axis only
- ctrl + left-click on plot to reset view

There is a scrollbar along the bottom to advance time. The neighboring arrow buttons increment/decrement timestep.

### An example display:
![mission simulation](assets/example_screenshot.png?raw=true)