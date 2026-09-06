````@raw html
---
# https://vitepress.dev/reference/default-theme-home-page
layout: home

hero:
  name: "OpenSOAP"
  text: ""
  tagline: CONOPs simulation
  image:
    src: logo.svg
    alt: OpenSOAP
  actions:
    - theme: brand
      text: Get Started
      link: https://google.com
    - theme: alt
      text: View on Github
      link: https://github.com/thanasipantazides/OpenSOAP

  features:
      
---
````

This is a package for simulating satellites in LEO doing the following:

- going around in an orbit
- spinning around/controlling their attitude
- pointing at different targets of interest (groundstations, the Sun, etc.)
- consuming power and producing data, at rates driven by their current activites

This sort of simulation is useful for exploring CONOPs, closing technical budgets during early design, sizing subsystems, and performing trade studies.

## Alternatives and comparison

| Tool          | Comparison                                                                                                                                                        |
| ------------- |:---------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| GMAT          | Better for orbit propagation and orbit + transfer planning                                                                                                        |
| STK           | Probably does all this better than OpenSOAP. Probably costs a lot of money. Probably wants to talk to MATLAB for you to plug in, which also costs a lot of money. |
| FreeFlyer     | Probably does all this better than OpenSOAP. Also probably expensive, but less so than STK                                                                        |
| 42            | Specializes in orbit modeling, multibody and formation flying mechanics                                                                                           |
| GENEOS        | Closed outside of EU users                                                                                                                                        |
| Orekit        | Primarily targets orbit propagation                                                                                                                               |
| AMMOS PlanDev | Operational-focused, **very** scalable                                                                                                                            |
| Basilisk      | Primarily for high fidelity orbit and ADCS modeling                                                                                                               |
