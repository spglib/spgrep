from spgrep import get_spacegroup_irreps
from spgrep.representation import get_character

# Rutile structure (https://materialsproject.org/materials/mp-2657/)
# P4_2/mnm (No. 136)
a = 4.603
c = 2.969
x_4f = 0.3046
lattice = [
    [a, 0, 0],
    [0, a, 0],
    [0, 0, c],
]
positions = [
    [0, 0, 0],  # Ti(2a)
    [0.5, 0.5, 0.5],  # Ti(2a)
    [x_4f, x_4f, 0],  # O(4f)
    [-x_4f, -x_4f, 0],  # O(4f)
    [-x_4f + 0.5, x_4f + 0.5, 0.5],  # O(4f)
    [x_4f + 0.5, -x_4f + 0.5, 0.5],  # O(4f)
]
numbers = [0, 0, 1, 1, 1, 1]

kpoint = [0.5, 0, 0]  # X point
irreps, rotations, translations, mapping_little_group = get_spacegroup_irreps(
    lattice, positions, numbers, kpoint
)

# Symmetry operations by spglib
assert len(rotations) == 16
assert len(translations) == 16

# At X point, the little co-group is isomorphic to mmm (order=8)
assert len(mapping_little_group) == 8
print(mapping_little_group)  # [ 0,  1,  4,  5,  8,  9, 12, 13]

# Two two-dimensional irreps
for irrep in irreps:
    print(get_character(irrep))
# [2.+0.j 0.+0.j 0.+0.j 2.+0.j 0.+0.j 0.+0.j 0.+0.j 0.+0.j]
# [2.+0.j 0.+0.j 0.+0.j -2.+0.j 0.+0.j 0.+0.j 0.+0.j 0.+0.j]
