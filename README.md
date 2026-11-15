# Common ancestors at a prescribed distance in a random recursive tree

Reproducibility materials for **The last common ancestor of two vertices at a prescribed distance in a random recursive tree**, by Marek Gałązka and Hanna Wdowicka.

A new vertex chooses a parent uniformly among earlier vertices. If two uniformly chosen distinct vertices are exactly a specified number of edges apart, when was their last common ancestor born? The scripts check the paper's finite joint probabilities and reproduce its two figures.

## Verify the finite identities

```sh
python3 checks/verify_exact.py --max-n 9
```

Only the Python standard library is required. The script enumerates every parent sequence, counts both orientations of every pair, and compares the counts with exact rational probabilities for the distance, ancestor label, first arm length and ancestor depth. A mismatch raises an assertion. At size nine it checks 40,320 trees and 2,903,040 ordered tree–pair combinations. The [recorded output](results/enumeration.txt) includes the seven-vertex example in the introduction. Among 10,224 distance-two pairs across all 720 trees, 6,496 have the root as their last common ancestor, giving 406/639.

## Reproduce the figures

The figures were generated with Python 3.12.14, NumPy 2.3.5 and Matplotlib 3.10.8. From the repository root, run the following commands.

```sh
python3 -m pip install -r requirements.txt
python3 checks/draw_tree_example.py
python3 checks/make_figure.py
```

- [Figure 1](figures/two_equal_paths.pdf) shows two paths of length two in one labelled tree. Each path and its shared ancestor are derived from the same parent list.
- [Figure 2](figures/conditional_laws.pdf) compares finite-tree probabilities and limiting curves. Panel (b) plots the inclusive upper tail at each birth label. The value at zero is one, and the drop immediately to its right is the conditional root probability. The finite values are obtained by convolving independent Bernoulli record indicators.

The scripts write the vector PDFs into `figures/`. Neither check requires external data.
