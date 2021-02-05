# tunedSNN
> A two-layer spiking neural network consisting of discretely tuned, disjoint inputs and a layer of recurrently connected neurons containing assemblies that reinforce the input tuning. 

Layer 1 network activity has two-dimensional shared variance by design, where the two dimensions correspond to 2 different tuned inputs. Recurrent layer network activity has high (>2) dimensional shared variance when recurrent clustering is strong.

All codes are for the manuscript: Rager, D. M., Khanna, S., Smith, M., & Doiron, B. "Assembly structure expands the dimension of shared variability
in cortical networks" (2020, pending) More information about this research can be found at [cortical.network](https://cortical.network/proj_spikingNets.html) 

Julia (v1.5) codes simulate the spiking neural network. Matlab (R2019b) codes analyze the network activity simulated by the Julia codes.

## Useage

`runsim.jl`:
1. initalizes multiple Julia processes (one per physical computer core) 
2. initializes network parameters `simParams`(simulation length in ms),`sysSize` (number of input, recurrent layer E, and recurrent layer I neurons),`connProbs` (network connection probabilities for each celltype pair),`taus`(membrane and synaptic time constants in ms),`v4OU` (strength and timescale of OU process that correlates input layer spiking activity) 
3. calls all network simulation functions.
<br><br><br>

Call
```julia
times,ns,times0,ns0,weights,bias,connStrength = simTwoPopHemiInputUnpack_NoCoupleInit(simParams,sysSize,connProbs,taus,v4OU)
```
to simulate 1 trial (default 12 seconds) of activity from a network with tuned, disjoint inputs and no recurrent connections.
<br><br><br>

```julia
times,ns,times0,ns0,weights,voltageOverTime,bias,connStrength = simTwoPopHemiInputUnpack_WeakCoupleInit(simParams,sysSize,connProbs,taus,v4OU)
```
simulates 1 trial of activity from a network with tuned, disjoint inputs and weak recurrent connections.
<br><br><br>


```julia
times,ns,times0,ns0,weights,synInputPerNeuronOverTime,bias,connStrength = simTwoPopHemiInputUnpack_StrongRecSymmClusters(simParams,sysSize,connProbs,taus,v4OU)
```
simulates 1 trial of activity from a network with tuned, disjoint inputs and strong (uniform or clustered) recurrent connections. Within the `simTwoPopHemiInputUnpack_StrongRecSymmClusters` method, you can change R, the clustering strength of the recurrent architecture. R = 1.25 corresponds to a network in which recurrent layer neurons receiving the same tuned input connect with 1.25x the strength as recurrent layer neurons receiving opposite tuned inputs. R is the 2nd argument in the `ConnStrength` constructor. Upcoming versions of the code base will make R accessible as an argument to `simTwoPopHemiInputUnpack_StrongRecSymmClusters`.
<br><br><br>

```julia 
runSimsFrozenGraph(numSims,simParams,sysSize,connProbs,connStrength,taus,v4OU,bias,weights)
```
simulates `numSims` trials of activity for a frozen network architecture that was initialized with one of the above functions.
<br><br><br>

```julia
runSimsFrozenGraphLinResp(numSims,simParams,sysSize,connProbs,connStrength,taus,v4OU,bias,weights)
```
mirrors the behavior of `runSimsFrozenGraph(numSims,simParams,sysSize,connProbs,connStrength,taus,v4OU,bias,weights)` plus runs an additional 20 trials of activity for a network with the same, frozen input structure and no recurrent connectivity. This is useful for the paper's linear response calculations.   














## Dependencies

- Julia (v 1.2+)
- [Hwloc](https://github.com/JuliaParallel/Hwloc.jl)
- [PyPlot](https://github.com/JuliaPy/PyPlot.jl)
- [JLD2](https://github.com/JuliaIO/JLD2.jl)
- [HDF5](https://github.com/JuliaIO/HDF5.jl)
- [CSV](https://juliadata.github.io/CSV.jl/stable/index.html)
- [Parameters](https://github.com/mauro3/Parameters.jl)
- [BenchmarkTools](https://github.com/JuliaCI/BenchmarkTools.jl)

##

- MATLAB (R2019b)
- [Statistics and Machine Learning Toolbox](https://www.mathworks.com/help/stats/index.html?s_tid=CRUX_lftnav)
- [Shape Language Modeling](https://www.mathworks.com/matlabcentral/fileexchange/24443-slm-shape-language-modeling)



## Contact
For bugs, comments, concerns: use the Github issue tracker.

Author: [Danielle Rager](https://cortical.network), danielle [at] cortical.network
