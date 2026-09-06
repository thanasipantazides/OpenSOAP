# Overview

A complex system is an amorphous blob in a higher dimension. All we can do is project it into 3D or 2D space from a few different angles and hope to gain some understanding of its true nature.

## Architecture

There are three parts of OpenSOAP.

1. A dynamics simulation. This takes an initial condition and pushes it forward in time.
1. A visualizer. This takes results from the simulation and displays them.
1. You, the user, who interacts with (1) and (2). <3.

These parts talk to each other using [sockets](https://docs.julialang.org/en/v1/manual/networking-and-streams/#A-simple-TCP-example). The simulation produces small `struct`s that pass over the network for the visualizer to display. 

All the important `struct`s in OpenSOAP are subtypes of `NetworkMessage`:

```
NetworkMessage
 ├─ AbstractConfig
 │   ├─ AbstractTargetConfig
 │   │   ├─ EarthConfig
 │   │   ├─ GroundConfig
 │   │   ├─ MagneticFieldConfig
 │   │   └─ SunConfig
 │   ├─ ModeConfig
 │   ├─ SatelliteConfig
 │   └─ SimConfig
 ├─ AbstractConstraint
 │   └─ OpenSOAP.LLAConstraint
 ├─ AbstractState
 │   ├─ AbstractTarget
 │   │   ├─ EarthState
 │   │   ├─ GroundState
 │   │   ├─ OpenSOAP.MagneticFieldState
 │   │   └─ SunState
 │   ├─ AttitudeState
 │   ├─ PositionState
 │   └─ SatelliteState
 └─ ControlMessage
     ├─ AskMessage
     ├─ PerturbationMessage
     ├─ PlayMessage
     ├─ QuitMessage
     └─ RateMessage
```

!!! tip
    I generate this type tree using the handy [TypeTree.jl](https://github.com/cnaak/TypeTree.jl) package:
    ```julia-repl
    julia> import TypeTree
    julia> using OpenSOAP
    julia> print(TypeTree.tt(NetworkMessage)...)
    ```

There are approximately four types of `NetworkMessage`: configs, states, constraints, and control messages. The control messages are for the REPL or visualizer to control how the simulation runs, and we'll ignore those for now. Of the remainder:

- configs (`<:AbstractConfig`) define time-fixed properties of the simulation—things like the display name of a groundstation, or how efficient solar panels are. These are things that do not change when the simulation moves forward in time. Of course, as an interactive tool, the user can change these dynamically while the simulation runs.
- states (`<:AbstractState`) define properties of the simulation that _do_ change in time—things like the position and attitude of the spacecraft, or the transform from the Earth-fixed frame to the inertial frame.
- constraints (`<:AbstractConstraint`) are used to restrict when certain operating modes can be active. 

!!! note "On terminology"
    _States_ and _modes_ as distinct things. _State_ is used in the [dynamical systems sense](https://en.wikipedia.org/wiki/State-space_representation#State_variables)—the state of a system fully describes the system at a given moment in time, and informs the forward evolution of the system. A _mode_ captures how the spacecraft is operating at a moment in time with respect to onboard systems. Usually there are dedicated modes for e.g. charging, downlinking data, detumbling, doing some mission-specific operations, etc.