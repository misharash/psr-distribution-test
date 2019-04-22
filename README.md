# Pulsar distribution test
Checks how the new correction to BGI model affects the results of Arzamasskiy, Beskin & Pirov (2017) paper

* mc &mdash; Monte-Carlo simulation: creates initial pulsar randomly according to predefined distribution functions, evolves them in time and adds new ones according to birth functions.
* pde &mdash; numerical solution of kinetic equation as PDE using staggered leapfrog method starting from P=0 where N dP/dt -> 0.