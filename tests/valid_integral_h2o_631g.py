# +
from pyscf import gto
import numpy as np
import itertools

np.set_printoptions(16)
# -

mol = gto.Mole(atom="O; H 1 0.94; H 1 0.94 2 104.5", basis="6-31G").build()

nshl = mol.nbas
nao = mol.nao

# ## int3c2e

shls_slice = [[1, 4], [0, 5], [3, 9]]

ao_loc = mol.ao_loc[list(itertools.chain(*shls_slice))].reshape(-1, 2)

out = np.zeros((nao, nao, nao))
out_view = out[ao_loc[0][0]:ao_loc[0][1], ao_loc[1][0]:ao_loc[1][1], ao_loc[2][0]:ao_loc[2][1]]
out_view[:] = mol.intor("int3c2e", shls_slice=list(itertools.chain(*shls_slice)))

print(out.sum())

print((out * np.linspace(0, 10, out.size).reshape(out.shape)).sum())





# ## Molecule information

mol._atm

mol._bas

mol._env
