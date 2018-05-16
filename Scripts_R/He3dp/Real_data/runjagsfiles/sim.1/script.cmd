model in "model.txt"
data in "data.txt"
compile, nchains(3)
parameters in "inits1.txt", chain(1)
parameters in "inits2.txt", chain(2)
parameters in "inits3.txt", chain(3)
initialize
adapt 100
update 100
monitor e1, thin(1)
monitor gin, thin(1)
monitor gout, thin(1)
monitor ue, thin(1)
monitor tau, thin(1)
monitor ri, thin(1)
monitor rf, thin(1)
monitor RSS, thin(1)
monitor mux0, thin(1)
monitor mux1, thin(1)
monitor mux2, thin(1)
monitor scale, thin(1)
monitor DeltaM, thin(1)
monitor e1_2, thin(1)
monitor gin_2, thin(1)
monitor gout_2, thin(1)
monitor tau_2, thin(1)
monitor ri_2, thin(1)
monitor rf_2, thin(1)
update 100
parameters to "out1.Rdump", chain(1)
parameters to "out2.Rdump", chain(2)
parameters to "out3.Rdump", chain(3)
coda *, stem(sim.1/CODA)
samplers to sim.1/samplers.csv
update 0
model clear
exit
