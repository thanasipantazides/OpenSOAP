# OpenSOAP

Spacecraft Operations Analysis Package

## To do:
- [ ] swap underlying socket to UDP so there is no locking of sim based on TCP. Note: because of the current Julia implementation of sockets, there is no  
- [ ] add PD control on $\mathrm{SO}(3)$ directly.
- [ ] ingest complex config. Or maybe a callback-ability for state transition events. or something.
- [ ] later: hot reload of config. I.e. a watcher for config file that reloads the spacecraft and target config (but not initial conditions). 
- [ ] 