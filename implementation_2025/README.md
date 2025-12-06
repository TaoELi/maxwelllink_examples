MaxwellLink Examples
----------------------

This project stores the input, data, and post-processing files of the following publication:

- X. Ji, A. F. Bocanegra Vargas, G. Meng, and T. E. Li. **MaxwellLink: A Unified Framework for Self-Consistent Light-Matter Simulations.** *Submitted* 2025.

### 1. Reproducing the figures

Please go to [plotting/](./plotting/) and run the corresponding Python files named by *plot\*.py*.

The Python files read data locally and generate the figures exactly the same as those in the paper.

### 2. Rerunning the simulations

Please go to other folders here, such as [meep_energy_transfer_optimized/](./meep_energy_transfer_optimized/), and then check and run the bash scripts named as *run\*.sh*

**Caution #1:** [MaxwellLink v0.2.0](https://github.com/TaoELi/MaxwellLink/releases/tag/v0.2.0) was used to perform this research. There is no guarantee that a newer version of MaxwellLink may exactly generate the same results, as the input variables may change in later updates. Please always use [MaxwellLink v0.2.0](https://github.com/TaoELi/MaxwellLink/releases/tag/v0.2.0) to reproduce this research.

**Caution #2:** The [Anvil HPC system in Purdue University](https://www.rcac.purdue.edu/knowledge/anvil) was used to submit the jobs. U.S. researchers may get access to this HPC system through the [NSF ACCESS](https://access-ci.org/) program. The input SLURM parameters should be adjusted accordingly if you use other HPC systems to rerun the simulation.






