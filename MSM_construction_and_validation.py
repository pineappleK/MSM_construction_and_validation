import os
import mdtraj as md
import pyemma
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from itertools import product
from tqdm import tqdm, trange


# Make result directory
os.makedirs('./results', exist_ok=True)

# 1. Load data
data = np.loadtxt('./data/tIC_input_toy.txt').astype(float)

# 2. Implied timescale test for data
k = 100 # In the toy model, k is low; it is recommended to use 300 for your simulations
fig, ax = plt.subplots(1, 1, figsize=(12, 4))
cluster = pyemma.coordinates.cluster_kmeans(data, k=k, max_iter=100, stride=1, fixed_seed=True)
dtrajs_concatenated = np.concatenate(cluster.dtrajs)
its = pyemma.msm.its(cluster.dtrajs, lags=[i for i in range(1, 40)], nits=4, errors='bayes')
pyemma.plots.plot_implied_timescales(its)
plt.tight_layout()
plt.savefig(f'./results/{k}_microstate_timescale.png', dpi=300)
plt.close()

# 3. Build the MSM
msm = pyemma.msm.estimate_markov_model(cluster.dtrajs, lag=30, dt_traj='200 ps')

# 4. Divide the landscape to 5 PCCA macrostates
for macro_number in [5]:
    prop_list = []
    nstates = macro_number
    msm.pcca(nstates)
    
    # To obtain the proportion of each macrostate
    print('The proportion for macrostate {} is'.format(macro_number))

    for i, s in enumerate(msm.metastable_sets):
        print('x_{} = {:f}'.format(i + 1, msm.pi[s].sum()))

    fig, axes = plt.subplots(1, nstates, figsize=(25, 4), sharex=True, sharey=True)
    try:
        for i, ax in enumerate(axes.flat):
            pyemma.plots.plot_contour(
                *data.T,
                msm.metastable_distributions[i][dtrajs_concatenated],
                ax=ax,
                cmap='afmhot_r',
                mask=True,
                cbar_label='metastable distribution {}'.format(i + 1))
            ax.set_xlabel('tIC 1')
        axes[0].set_ylabel('tIC 2')
        fig.tight_layout()
        plt.savefig(f'./results/{nstates}-multiple.png', dpi=300)
        plt.close()
    except:
        print('ax error happens!')
        
    # Calculate the transition between macrostates by MFPT
    mfpt = np.zeros((nstates, nstates))
    for i, j in product(range(nstates), repeat=2):
        mfpt[i, j] = msm.mfpt(
            msm.metastable_sets[i],
            msm.metastable_sets[j])
    # In unit of us
    print('MFPT / us')
    pd.DataFrame(np.round(mfpt, decimals=2) * 0.2 / 1000, index=range(1, nstates + 1), columns=range(1, nstates + 1))
    data1 = pd.DataFrame(np.round(mfpt, decimals=2) * 0.2 / 1000, index=range(1, nstates + 1),
                        columns=range(1, nstates + 1))
    data1.to_csv(f'./results/{nstates}-time.csv')


# 5. Extract representative structures (with mdtraj)
d_thr = 0.003 # Consider a proper point number in your simulations
pdb = './data/r5_toy.pdb'
parm = './data/r5_toy.parm7'
traj = './data/r5_toy.nc'

for number, i in enumerate(msm.metastable_sets):
    index_list = []
    for cents in [i]:
        b = cluster.clustercenters[cents]
        for bs in b:
            dist_arr = np.sqrt(np.sum(np.square(data - bs), axis=1))
            for index, nums in enumerate(dist_arr):
                if nums < d_thr:
                    index_list.append(index)
    
    # load trajectory via mdtraj
    t0 = md.load(pdb)
    t_cluster = md.load_frame(traj, index_list[0], top=parm)
    for ix in index_list[1:]:
        t_join = md.load_frame(traj, ix, top=parm)
        t_cluster = md.join([t_cluster, t_join])
    atom_indices = [a.index for a in t_cluster.topology.atoms]
    t_cluster = t_cluster.superpose(t0, atom_indices=atom_indices)

    distances = np.empty((t_cluster.n_frames, t_cluster.n_frames))
    for num, i in enumerate(range(t_cluster.n_frames)):
        distances[i] = md.rmsd(t_cluster, t_cluster, i, atom_indices=atom_indices)

    for beta in [1]:
        index = np.exp(-beta * distances / distances.std()).sum(axis=1).argmax()
        centroid = t_cluster[index]
        print(f'Macro {number + 1}: frame {index_list[index]}')
        centroid.save(f'./results/{nstates}-macro-{number + 1}-{d_thr}.pdb')

# 6. Chapman–Kolmogorov test
bayesian_msm = pyemma.msm.bayesian_markov_model(cluster.dtrajs, lag=30, dt_traj='200 ps', conf=0.95)
for i in [5]:
    pyemma.plots.plot_cktest(bayesian_msm.cktest(i), units='ps')
plt.savefig('./results/ck-test.png', dpi=300)
plt.close()