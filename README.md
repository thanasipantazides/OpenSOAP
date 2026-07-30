# OpenSOAP

Spacecraft Operations Analysis Package

## What it is

OpenSOAP is a simulator for spacecraft in LEO, with a bent for early mission design, CONOPs exploration, and technical budget analysis. When you run it, it looks like this:

![Screenshot of the simulator, showing a satellite going around Earth with some diagnostic plots off to the side.](assets/soap_screenshot.png)

The following dynamical properties are simulated:

- orbit
- attitude
- onboard power
- onboard data
- onboard operating mode (e.g. charging, science, etc.)

None of the modeled dynamics are state-of-the-art fidelity, but they are good for first-order designs.

## To run it

You will need an installation of [Julia](https://julialang.org/).

Clone the repository into a nice place on your computer, go into the directory, and launch the Julia REPL:

```console
$ git clone https://github.com/thanasipantazides/OpenSOAP.git
$ cd OpenSOAP
$ julia --project=.
```

From the REPL, hit the `]` key to go into the package manager. Your prompt should change, and you should run `instantiate` to install required packages.

```julia-repl
julia> ]
(OpenSOAP) pkg> instantiate
```

You can hit backspace to get back to the `julia>` prompt from the package prompt.

To launch the interactive orbit viewer, use the example called `both.jl`:

```julia-repl
julia> include("examples/both.jl")
```

This will probably precompile for a while. This happens the first time you run something from a new REPL session, but subsequent runs are much faster. Once that's done, run the main function:

```julia-repl
julia> main()
```

You should see a display like the one further up this page. Press the `h` key for a help menu.

To exit, press `q` in the viewer.
