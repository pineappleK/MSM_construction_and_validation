# MSMs construction and validation

## Data

The files within the `data` directory consist of toy simulation topologies and trajectories, in addition to the corresponding toy tIC data, which are essential for constructing MSMs.

## Setup Environment

We set up the environment using [Anaconda](https://docs.anaconda.com/anaconda/install/index.html).

**Requirements**

```
python==3.7.3
pyemma==2.5.7
```

## Running MSMs construction and validation on your own MD simulations

To utilize custom simulation topologies and trajectories, replace the existing files in the `data` directory with your own files. Subsequently, modify the relevant parameters within the `MSM_construction_and_validation.py` script to align with your specific data inputs.

And you are ready to run

```
python MSM_construction_and_validation.py
```

The anticipated outputs are the following:

- Implied timescale test (50_microstate_timescale.png)
  
- Chapman–Kolmogorov test (ck-test.png)
  
- Macrostate position (5-multiple.png)
  
- Transition time between states (5-time.csv)
  
- Representative structures (5-macro-1-0.003.pdb and so on)