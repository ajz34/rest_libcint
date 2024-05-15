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
slcs = [slice(loc[0], loc[1]) for loc in ao_loc]
out_view = out[slcs[0], slcs[1], slcs[2]]
out_view[:] = mol.intor("int3c2e", shls_slice=list(itertools.chain(*shls_slice)))

print(out.sum())
print((out * np.linspace(0, 10, out.size).reshape(out.shape)).sum())

# ## int3c2e_ip1 (inplace)

shls_slice = [[2, 6], [4, 14], [0, 8]]

ao_loc = mol.ao_loc[list(itertools.chain(*shls_slice))].reshape(-1, 2)
[list(l) for l in ao_loc]

out = np.zeros((nao, nao, nao, 3))
slcs = [slice(loc[0], loc[1]) for loc in ao_loc]
out_view = out[slcs[0], slcs[1], slcs[2]]
out_view[:] = mol.intor("int3c2e_ip1", shls_slice=list(itertools.chain(*shls_slice))).transpose(1, 2, 3, 0)

print(out.sum())
print((out * np.linspace(0, 10, out.size).reshape(out.shape)).sum())

# ## int3c2e_ip1 (outplace)

shls_slice = [[2, 6], [4, 14], [0, 8]]
out = mol.intor("int3c2e_ip1", shls_slice=list(itertools.chain(*shls_slice))).transpose(1, 2, 3, 0)

print(out.sum())
print((out * np.linspace(0, 10, out.size).reshape(out.shape)).sum())

out = mol.intor("int3c2e_ip1").transpose(1, 2, 3, 0)

print(out.sum())
print((out * np.linspace(0, 10, out.size).reshape(out.shape)).sum())

# ## Molecule information

mol._atm

mol._bas

mol._env
