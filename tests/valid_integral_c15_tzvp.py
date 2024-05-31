# +
from pyscf import gto, lib
import numpy as np
import itertools

np.set_printoptions(16)
# -

lib.num_threads(32)

mol = gto.Mole(atom="""
C          0.89735       -0.00113        0.08382
C          2.41558       -0.00973        0.08398
C          2.95403       -1.39223        0.44060
C          4.47840       -1.40490        0.43927
C          5.00472       -2.79081        0.79800
C          6.52868       -2.80677        0.79924
C          7.04947       -4.19283        1.15712
C          8.57098       -4.20695        1.16283
C          9.09164       -5.59214        1.51883
C         10.61579       -5.60738        1.52100
C         11.13683       -6.99323        1.87886
C         12.66149       -7.01081        1.88281
C         13.18196       -8.39866        2.24144
C         14.70681       -8.42700        2.24831
C         15.22338       -9.80961        2.60521
H          0.52316        0.99486       -0.17283
H          0.50536       -0.27018        1.06991
H          0.50315       -0.71346       -0.64803
H          2.78364        0.73060        0.80387
H          2.78126        0.28908       -0.90547
H          2.57994       -2.13111       -0.27881
H          2.58434       -1.68881        1.42971
H          4.85726       -0.66903        1.15930
H          4.85232       -1.11324       -0.54954
H          4.62816       -3.52793        0.07814
H          4.63163       -3.08460        1.78691
H          6.90615       -2.07036        1.51945
H          6.90302       -2.51453       -0.18957
H          6.67226       -4.92766        0.43715
H          6.67377       -4.48532        2.14643
H          8.94855       -3.47382        1.88410
H          8.94922       -3.91225        0.17366
H          8.71508       -6.32967        0.79715
H          8.71691       -5.88744        2.50652
H         10.99100       -4.87238        2.24175
H         10.98995       -5.31412        0.53343
H         10.75962       -7.72809        1.15673
H         10.76022       -7.28664        2.86631
H         13.03725       -6.27521        2.60442
H         13.03678       -6.71689        0.89495
H         12.80104       -9.13172        1.51944
H         12.80163       -8.68997        3.22825
H         15.08995       -7.69771        2.97148
H         15.08954       -8.13915        1.26219
H         14.87740      -10.11040        3.59946
H         16.31717       -9.81702        2.60682
H         14.87685      -10.55420        1.88142
""", basis="def2-TZVP").build()

nshl = mol.nbas
nao = mol.nao
nao, nshl

# ## test_int3c2e_s1_full

# %%time
out = mol.intor("int3c2e")

out_c = out.transpose(2, 1, 0)
assert out_c.flags.c_contiguous
scale = np.linspace(-1, 1, out.size)
print(out_c.sum())
print((out_c.flatten() * scale).sum())

# ## int3c2e_s1_slice

# %%time
out = mol.intor("int3c2e", shls_slice=(7, 275, 12, 129, 5, 264))

out_c = out.transpose(2, 1, 0)
assert out_c.flags.c_contiguous
scale = np.linspace(-1, 1, out.size)
print(out_c.sum())
print((out_c.flatten() * scale).sum())

# ## int3c2e_ip1_s1_slice

# %%time
out = mol.intor("int3c2e_ip1", shls_slice=(7, 275, 12, 129, 5, 264))

out_c = out.transpose(0, 3, 2, 1)
assert out_c.flags.c_contiguous
scale = np.linspace(-1, 1, out.size)
print(out_c.sum())
print((out_c.flatten() * scale).sum())

# ## int2c2e_ip1_s1_slice

# %%time
out = mol.intor("int2c2e_ip1", shls_slice=(7, 275, 5, 264))

out_c = out.transpose(0, 2, 1)
assert out_c.flags.c_contiguous
scale = np.linspace(-1, 1, out.size)
print(out_c.sum())
print((out_c.flatten() * scale).sum())

# ## int2c2e_ip1ip2_s1_slice

# %%time
out = mol.intor("int2e_ip1ip2", shls_slice=(7, 30, 5, 52, 58, 89, 125, 156))

out_c = out.transpose(0, 4, 3, 2, 1).copy()
assert out_c.flags.c_contiguous
scale = np.linspace(-1, 1, out.size)
print(out_c.sum())
print((out_c.flatten() * scale).sum())

# ## int3c2e_s2ij_full

# %%time
out = mol.intor("int3c2e", aosym="s2ij")

out_c = out.transpose(1, 0)
assert out_c.flags.c_contiguous
scale = np.linspace(-1, 1, out.size)
print(out_c.sum())
print((out_c.flatten() * scale).sum())

# ## int3c2e_ip2_s2ij_slice

# %%time
out = mol.intor("int3c2e_ip2", shls_slice=(0, 275, 0, 275, 7, 264), aosym="s2ij")

out_c = out.transpose(0, 2, 1)
assert out_c.flags.c_contiguous
scale = np.linspace(-1, 1, out.size)
print(out_c.sum())
print((out_c.flatten() * scale).sum())

# **Note** In pyscf, it seems that `s2ij` in `getints3c`, starting shell should always zero, otherwise it may not seems to be correct.

# ## int2e_ip2_s2ij_slice

# %%time
out = mol.intor("int2e_ip2", shls_slice=(10, 50, 10, 50, 127, 168, 215, 272), aosym="s2ij")

out_c = out.transpose(0, 3, 2, 1).copy()
assert out_c.flags.c_contiguous
scale = np.linspace(-1, 1, out.size)
print(out_c.sum())
print((out_c.flatten() * scale).sum())

# ## int2c2e_s2ij_slice

# %%time
out = mol.intor("int2c2e", shls_slice=(10, 283, 10, 283), aosym="s2ij")

out_c = out[np.tril_indices(out.shape[0])]
assert out_c.flags.c_contiguous
scale = np.linspace(-1, 1, out_c.size)
print(out_c.sum())
print((out_c.flatten() * scale).sum())

# ## Molecule information

np.set_printoptions(threshold=1e15)

mol._atm

mol._bas

mol._env
