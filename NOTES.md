# OpenSOAP development notes

## To do

- [ ] ~~swap underlying socket to UDP so there is no locking of sim based on TCP~~. Note: because of the current Julia implementation of sockets, there is no Unix domain datagram socket available. Shucks.
  - Another note: visually, you actually do want the TCP ACK to hold up transfer of more data. This forces the simulation to run in lockstep with the display. There are ways around this, but they are more work.
- [x] add PD control on $\mathrm{SO}(3)$ directly.
- [ ] ingest complex config. Or maybe a callback-ability for state transition events. or something.
- [ ] Make it easier to launch and exit.
- [ ] later: hot reload of config. I.e. a watcher for config file that reloads the spacecraft and target config (but not initial conditions).
- [x] Add labeling control to monitor.
- [ ] Add a `initialize_targets()` function or similar, which just calls `step!()` on each target at the initial time and spacecraft position.
- [ ] In `simulate.jl`, maybe put all configs in a big dict (by ID), ditto for all states. Maybe two dicts, one by ID and another by type. And just pass the dict around instead a bunch of parameters.
  - [ ] Or put _everything_ into one really big dict? Since the IDs are unique? Less args to throw around everywhere. To do this:
    - [ ] Everywhere `modes::Vector{ModeConfig}` appears, replace with `modes::IDDict{ModeConfig}` and update local logic accordingly.
    - [x] Everywhere `(target_states, target_configs, constraints, modes)` appear together, merge into one `IDDict`. Update local logic to also find based on type, if needed.
      - [x] `io.jl`
    - [ ] Add a cache object that gets passed around internally in the simulation steps. `cache.sat == sat.id`, `cache.groundstations == filter(p->p.id && p isa GroundState, sim_env)`. And so on.
- [x] Make all structs mutable and default-constructible via `@kwdef`. This enables the deserializer to allocate, then memcpy into them. ~~Use internal `const` fields for IDs.~~

### Complex config

The new config looks like this at the top level:

```jsonc
{
  "simulation": {},
  "spacecraft": {
    "modes": {},
    "targets": [],
    "constraints": [],
  },
}
```

A `mode` defines an operating mode and references (by name) `targets` and `constraints`. A `target` provides a reference direction or position vector for the spacecraft to look at, and maybe some additional metadata. A `constraint` limits when a `mode` can be used.

Items needed to implement this:

- [x] A type `MagneticState` or similar. Holds almost no info aside from B-field dir.
- [x] Maybe: A new abstract `ModeConstraint` or similar, which is concretized as `LLAConstraint` for this application.
  - [x] It gets an ID, and a `ModeConfig` can have a field `constraint_ids::Vector{IDType}` (similar to how it now has `target_ids::Vector{IDType}` which refer to the targets for a mode.)
  - How is this different from just being another `TargetState`? Ans: it doesn't have a `direction`. Just some criteria.
- [x] A type `LLAMaskConfig` or similar. Include a flag (in JSON too) for whether to point nadir always along field, zenith along field, along field S, or along field N.
- [x] For now, only `MagneticState` needs to be serializable.
- [x] Add a JSON loader for `MagneticState`. Should be easy. If the target has a `"mask"` key, make a `LLAMaskConfig` that is linked to `MagneticState` by ID.
- [x] Add a JSON loader for `LLAMaskConfig`.
- [x] Add a `step!` function for a `MagneticState`. Just query the IGRF and flip the vector nadir if necessary.
- [x] Add a `set_visibility!` function for a `MagneticState`.

Could imagine other kinds of `*MaskConfig`s, maybe for ECI position, or dynamic state properties.

### [`simulate.jl`](src/simulate.jl) refactor

The simulation has a _lot_ of functions that all take some permutations of the same few arguments. I need to boil this down into some nicer interfaces.

Thoughts, list read top-to-bottom goes from the bottom of the code up:

- `set_visibility!()` functions should dispatch on some `<:AbstractTarget`. They (union) need to know spacecraft position, their own `<:AbstractConfig`, and time. They don't care about `ModeConfig` it is.
- `check_constraint()` functions should dispatch on `<:AbstractConstraint`. They (union) need to know spacecraft state, and time.
- Have a high-level wrapper to look at `SatelliteState`, check `mode`, and check all constraints for that mode (union of constraints).
- Have a high-level wrapper to look at `SatelliteState`, check `mode`, and check visibility of any targets for that mode (intersection of visibilities).
- `step!()` functions should dispatch on `<:AbstractTarget`.
  - Could move the `set_visibility!()` _outside_ of `step!()` to avoid binding their signatures together, but there is something uncomfortable about having a bunch of necessary sequential mutating calls necessary to push state forward. Order matters.
- What I see emerging: there are some high-level functions on `SatelliteState` that do **specific** things. There are some low-level functions on `<:AbstractTarget` and `<:AbstractConstraint` that do **generic** things.
  - The specific functions have specialized signatures, and the generic ones should dispatch with generic signatures.
- Consider: an `<:AbstractState` is oblivious to its configuration data.
- How about this for the signatures:

```julia
const IDDict{T} = Dict{IDType, T}

function high_level!(
    sat::SatelliteState,
    sat_config::SatelliteConfig,
    target_states::IDDict{AbstractTarget},
    target_configs::IDDict{AbstractConfig},
    constraints::IDDict{AbstractConstraint},
    modes::IDDict{<:ModeConfig},
    dt::Float64,
    time::Dates.DateTime;
    kwargs...)
end

# rather than target_state, target_config be abstract here,
# they should specialize for dispatch.
function low_level!(
    target_state::AbstractTarget,
    target_config::AbstractConfig,
    sat::SatelliteState,
    sat_config::SatelliteConfig,
    constraints::IDDict{<:AbstractConstraint},
    modes::IDDict{<:ModeConfig},
    dt::Float64,
    time::Dates.DateTime;
    kwargs...)
end

# convention for low-level functions: first args should be the dispatching args. the rest can be all the context/config data.
```

(I don't want to just wrap everything simulation-related in a giant struct and pass that around. That loses the multiple dispatch-ness and creates a sad interface. But some blobbing is necessary.)

Items to do:

- [x] Should pass around mode tables (`Dict{IDType, ModeConfig}`), target state tables (`Dict{IDType, <:AbstractState}`), target config tables (`Dict{IDType, <:AbstractConfig}`), and constraint tables (`Dict{IDType, <:AbstractConstraint}`).
  - Instead of `Vector`s of the same.
  - If a vector is passed, I just need to filter it to find IDType. If a table is passed, I can query by ID directly.
  - There is a small amount of filtering/finding/mapping by type (e.g. to find the `SunState`). That can still be done by iterating on the `Dict` values.
- [x] Rewrite all the signatures in `simulate.jl` based on the `high_level()` and `low_level()` prototypes defined above.
- [x] Rewrite (and reconsider) `mode_table`. It currently is a matrix of targets (rows) by modes (cols) indicating mode dependency on targets. This is complicated by the introduction of constraints. It's not part of a clean function signature. If we kept it, we would instead multiple by `[visibilities; constraints_satisfied]` boolean vector.
- [ ] Remove `params` from signature. There's a `Dict` in `SimConfig` that is meant to hold stuff like this.
- [ ] Revise `step_satellite` to `step!()`

### REPL interactivity

Plan is to use existing network ser/des stuff to pass `<:AbstractConfig` and `<:AbstractState`s to `core` from REPL context. Something like this:

![REPL interactivity plan](assets/repl-plan.pdf)

Things to do:

- [x] Make all `<:AbstractConfig` serializable. Maybe switch over to built in serialization if needed.
  - [x] In, e.g., `ModeConfig`, `Vector{IDType}` is **not** `isbits()`. Need to use a StaticArray instead, likely. I think this is fine—size of the `ModeConfig.target_ids` list is fixed after initialization.
  - This turned out to be very involved. Basically need to rewrite the ser/des functions, but they end up being more generic which is nice. But some side effects:
    - [x] Stop sending a `(count, <:NetworkMessage)` tuple as the default State message from core to monitor.
    - [ ] If a dedicated `count` message is needed, create a new dedicated message type for that. (Currently, `count` is unused in monitor).
    - [x] Add a default constructor for all `<:NetworkMessage` mutable types. This is required for serialization of mutable structs.
- [x] Add a REPL to CORE IP/port and sockets. REPL process is server, CORE is client.
- [x] Set up the socket in CORE and add an async listening process on it. Use `Threads.@spawn`.
  - [x] Deserialize what comes off the socket and conditionally update some process-wide `IDDict`s if the contents look good.
- [x] REPL interface needs:
  - [x] Some way to identify the shared objects in a bag. Maybe a way to mark certain ones as stale.
  - [x] An async socket interface to CORE, with an `update()` function that writes (either all, or just modified) objects to the interface.

### Notes from serializability via reinterpret:

- `Vector{IDType}` is not serializable (it is not `isbits`).
- `StaticArrays.SVector{N, UInt16}` is `isbits`. But `N` shouldn't be much bigger than 100.
- `StaticArrays.SizedVector{N}` not `isbits`.
