# tunedSNN
> A two-layer spiking neural network consisting of disjoint, tuned inputs and a layer of recurrently connected neurons containing assemblies that reinforce the input tuning. 

Layer 1 network activity has two-dimensional shared variance by design, where the two dimensions correspond to 2 different tuned inputs. Recurrent layer network activity has high (>2) dimensional shared variance when recurrent clustering is strong.

All codes are for the manuscript: Rager, D. M., Khanna, S., Smith, M., & Doiron, B. "Assembly structure expands the dimension of shared variability
in cortical networks" (2020, pending) More information about this research can be found at [cortical.network](https://cortical.network/proj_spikingNets.html) 

Julia (v1.5) codes simulate the spiking neural network. Matlab (R2019b) codes analyze the network activity simulated by the Julia codes.

## Useage



## Dependencies

- Julia (v 1.2+)
- [Hwloc](https://github.com/JuliaParallel/Hwloc.jl)
- [PyPlot](https://github.com/JuliaPy/PyPlot.jl)
- [JLD2](https://github.com/JuliaIO/JLD2.jl)
- [HDF5](https://github.com/JuliaIO/HDF5.jl)
- [CSV](https://juliadata.github.io/CSV.jl/stable/index.html)
- [Parameters](https://github.com/mauro3/Parameters.jl)
- [BenchmarkTools](https://github.com/JuliaCI/BenchmarkTools.jl)
<pre>


</pre>

- MATLAB (R2019b)
- [Statistics and Machine Learning Toolbox](https://www.mathworks.com/help/stats/index.html?s_tid=CRUX_lftnav)
- [Shape Language Modeling] (https://www.mathworks.com/matlabcentral/fileexchange/24443-slm-shape-language-modeling)



## Contact
For bugs, comments, concerns: use the Github issue tracker.

Author: [Danielle Rager](https://cortical.network), danielle [at] cortical.netwoek
