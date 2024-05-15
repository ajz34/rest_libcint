# +
from pyscf import gto
import numpy as np
import itertools

np.set_printoptions(16)
# -

mol = gto.Mole(atom="""
    O 0 0 0
    O 0 1 1
    H 0 2 0
    H 2 0 2
""", basis="6-31G").build()

nshl = mol.nbas
nao = mol.nao
nao, nshl

# ## int3c2e

shls_slice = [[7, 14], [0, 10], [3, 12]]

ao_loc = mol.ao_loc[list(itertools.chain(*shls_slice))].reshape(-1, 2)
[list(l) for l in ao_loc]

out = np.zeros((nao, nao, nao))
out_view = out[ao_loc[0][0]:ao_loc[0][1], ao_loc[1][0]:ao_loc[1][1], ao_loc[2][0]:ao_loc[2][1]]
out_view[:] = mol.intor("int3c2e", shls_slice=list(itertools.chain(*shls_slice)))

mol.intor("int3c2e", shls_slice=list(itertools.chain(*shls_slice)))[0, 0]

mol.intor("int3c2e").sum()







print(out.sum())

print((out * np.linspace(0, 10, out.size).reshape(out.shape)).sum())

mol.ao_loc[shls_slice[0][0]:shls_slice[0][1]]

arr_i = []
for i in range(*shls_slice[0]):
    arr_j = []
    for j in range(*shls_slice[1]):
        arr_k = []
        for k in range(*shls_slice[2]):
            arr_k.append(
                out[
                    mol.ao_loc[i]:mol.ao_loc[i+1],
                    mol.ao_loc[j]:mol.ao_loc[j+1],
                    mol.ao_loc[k]:mol.ao_loc[k+1],
                ].sum()
            )
        arr_j.append(arr_k)
    arr_i.append(arr_j)

# +
# np.array(arr_i).flatten() - block_val
# -

np.array(arr_i).flatten()[:10]

(block_val == 0).sum()

(np.abs(np.array(arr_i)) < 1e-14).sum()

(np.abs(np.array(block_val)) < 1e-14).sum()













# ## Molecule information

mol._atm

mol._bas

mol._env






