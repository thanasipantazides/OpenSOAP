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

You should see a display like the one further up this page. Press the `h` key for a help menu:

### Actions in the monitor

| Key          | Action                           |
| ------------ | -------------------------------- |
| `h`          | show help                        |
| `v`          | verbose                          |
| `q`          | quit                             |
| `e`          | change earth texture             |
| `f`          | show/hide Body frame             |
| `l`          | show/hide labels                 |
| `p`          | change camera projection         |
| `k`          | perturb attitude (check control) |
| `spacebar`   | play/pause                       |
| `left/right` | slower/faster                    |

To exit, press `q` in the viewer. Maybe also ctrl-C in the REPL.

## Details

A `mode` is an operational state of the spacecraft. The power consumption is determined by `mode`. Some modes have higher priority than others. A mode may have `constraint`s which must be satisfied to run it. A mode may have one or more `target`s to look at. The rules are:

- A mode can run if _any_ of its `target`s are visible.
- A mode can run if _all_ of its `constraint`s are satisfied.

A `target` is a position or reference direction for the spacecraft to look at. Maybe it is a groundstation. Maybe it is the sun. Maybe it is something else. Maybe it comes with some metadata.

A `constraint` is a restriction of `mode` based on spacecraft position. This is used to enforce things like `only

As it stands now, the waters are muddy: I use constrain the usage of the B-field pointing mode, but there are no explicit constraints on groundstation communication, that info is rolled into the groundstation visibility cone. So the interface is not very pure. But I hope to change that.

## Roadmap

I would like to add the following features, in the following priority order:

- [x] Finish REPL interface (propagate all REPL changes to vars into monitor).
- [ ] Documentation.
  - [ ] Usage guide.
  - [ ] Theory: state dynamics, control, onboard linear models. Visibility and eclipse.
  - [ ] Custom extension of JSON + providers interface.
- [ ] Schema validation of config JSON file.
- [ ] On [Makie 0.25 release](https://github.com/MakieOrg/Makie.jl/pull/5484), a GUI with tabs for editing config.
  - Ideally, this is just a skin on the existing REPL interface.
  - I could build something like this in Qt or Tachikoma as is, but Makie would be a nice solution.
- [ ] TLE + SGP4 propagation (with `SatelliteToolbox`).
- [ ] Fancier position dynamics (SRP, drag).
- [ ] Fancier attitude dynamics (SRP, drag, gravity gradient). Hand-in-hand with some fancier (linprog) control allocation.
- [ ] Some sort of import for [ADBSat](https://github.com/nhcrisp/ADBSat) coefficient tables would be very cool.
