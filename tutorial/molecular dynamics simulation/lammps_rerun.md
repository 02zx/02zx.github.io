#The input for a 2 component LJ nemd solid-liquid system
units lj

dimension    3
boundary     p p p

lattice fcc 0.8
region simcell block 0 5 0 5 0 5
create_box 1 simcell  

mass         1 1.0
pair_style   lj/cut 3.2
pair_coeff   1 1 1.0 1.0 3.2

neighbor      0.3 bin

thermo 100


compute myrdf all rdf 300

fix frdf all ave/time 200 100 20000  c_myrdf file lj.rdf mode vector

rerun dump.atom first 0 last 20000 every 200 dump x y z purge yes add yes replace no format native
