"""
Authors:
    Kee-Myoung Nam

Last updated:
    6/19/2024
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns

rng = np.random.default_rng(1234567890)

#######################################################################
def cut_from_end(lengths):
    """
    Cut a polymer of the given length from the end, into one monomer and 
    the remaining polymer.

    A polymer must have length >= 2 to be recognized.
    """
    # If there are no polymers of length >= 2, return as is
    if all(x < 2 for x in lengths):
        return lengths

    # Choose a polymer of length >= 2
    i = None
    length = 0
    while length < 2:
        i = rng.choice(len(lengths))
        length = lengths[i]

    # Cut the polymer and update lengths
    cut_lengths = 1, length - 1
    lengths = lengths[:i] + lengths[i+1:]
    lengths.append(cut_lengths[0])
    lengths.append(cut_lengths[1])
    return lengths

#######################################################################
def recognize_monomer(lengths):
    """
    Cut a polymer to the right of a randomly recognized monomer within
    the polymer.

    A polymer must have length >= 1 to be recognized. (If a polymer has 
    length == 1, then no cut is made.)
    """
    # If there are no polymers of length >= 1, return as is
    if all(x < 1 for x in lengths):
        return lengths

    # Choose a polymer of length >= 1
    i = None
    length = 0
    while length < 1:
        i = rng.choice(len(lengths))
        length = lengths[i]

    # Is the chosen polymer a monomer?
    if length == 1:
        return lengths

    # Otherwise, cut the polymer at a randomly chosen site and update lengths
    choice = rng.choice(length)
    if choice == length - 1:    # If the last monomer was chosen, no cut
        return lengths
    else:                       # Otherwise, cut
        cut_lengths = choice + 1, length - choice - 1
        lengths = lengths[:i] + lengths[i+1:]
        lengths.append(cut_lengths[0])
        lengths.append(cut_lengths[1])
        return lengths

#######################################################################
def recognize_dimer_cut_middle(lengths):
    """
    Cut a polymer in the middle of a randomly recognized dimer within 
    the polymer.

    A polymer must have length >= 2 to be recognized. 
    """
    # If there are no polymers of length >= 2, return as is
    if all(x < 2 for x in lengths):
        return lengths

    # Choose a polymer of length >= 2
    i = None
    length = 0
    while length < 2:
        i = rng.choice(len(lengths))
        length = lengths[i]

    # Cut the polymer at a randomly chosen site and update lengths
    choice = rng.choice(length - 1)
    cut_lengths = choice + 1, length - choice - 1
    lengths = lengths[:i] + lengths[i+1:]
    lengths.append(cut_lengths[0])
    lengths.append(cut_lengths[1])
    return lengths

#######################################################################
def recognize_dimer_cut_right(lengths):
    """
    Cut a polymer to the right of a randomly recognized dimer within the 
    polymer.

    A polymer must have length >= 2 to be recognized. (If a polymer has 
    length == 2, then no cut is made.)
    """
    # If there are no polymers of length >= 2, return as is
    if all(x < 2 for x in lengths):
        return lengths

    # Choose a polymer of length >= 2
    i = None
    length = 0
    while length < 2:
        i = rng.choice(len(lengths))
        length = lengths[i]

    # Cut the polymer at a randomly chosen site and update lengths
    choice = rng.choice(length - 1)
    if choice == length - 2:    # If the last dimer was chosen, no cut
        return lengths
    else:                       # Otherwise, cut
        cut_lengths = choice + 2, length - choice - 2
        lengths = lengths[:i] + lengths[i+1:]
        lengths.append(cut_lengths[0])
        lengths.append(cut_lengths[1])
        return lengths

#######################################################################
def recognize_trimer_cut_middle(lengths):
    """
    Cut a polymer in the middle of a randomly recognized trimer within 
    the polymer (between the second and third monomers of the trimer).

    A polymer must have length >= 3 to be recognized.
    """
    # If there are no polymers of length >= 3, return as is
    if all(x < 3 for x in lengths):
        return lengths

    # Choose a polymer of length >= 3
    i = None
    length = 0
    while length < 3:
        i = rng.choice(len(lengths))
        length = lengths[i]

    # Cut the polymer at a randomly chosen site and update lengths
    choice = rng.choice(length - 2)
    cut_lengths = choice + 2, length - choice - 2
    lengths = lengths[:i] + lengths[i+1:]
    lengths.append(cut_lengths[0])
    lengths.append(cut_lengths[1])
    return lengths

#######################################################################
total = 20

# Cut from end
lengths = [total]
print(lengths)
while any(x >= 2 for x in lengths):
    lengths = cut_from_end(lengths)
    print(lengths)

# Recognize dimer and cut in the middle
lengths = [total]
print(lengths)
while any(x >= 2 for x in lengths):
    lengths = recognize_dimer_cut_middle(lengths)
    print(lengths)

# Recognize dimer and cut to the right
lengths = [total]
print(lengths)
while any(x >= 3 for x in lengths):
    lengths = recognize_dimer_cut_right(lengths)
    print(lengths)

# Recognize trimer and cut in the middle
lengths = [total]
print(lengths)
while any(x >= 3 for x in lengths):
    lengths = recognize_trimer_cut_middle(lengths)
    print(lengths)

#######################################################################
fig, axes = plt.subplots(nrows=2, ncols=2, figsize=(9, 6))
total = 1000
colors = sns.color_palette()

# ------------------------------------------------------------------- #
# Cut from end
lengths = [total]
n_monomers = []
n_dimers = []
n_molecules = []
while any(x >= 2 for x in lengths):
    lengths = cut_from_end(lengths)
    n_monomers.append(sum(1 for x in lengths if x == 1))
    n_dimers.append(sum(1 for x in lengths if x == 2))
    n_molecules.append(len(lengths))
# Plot simulation 
axes[0, 0].plot(list(range(len(n_monomers))), n_monomers, c=colors[0])
axes[0, 0].plot(list(range(len(n_dimers))), n_dimers, c=colors[1])
axes[0, 0].plot(list(range(len(n_molecules))), n_molecules, c=colors[2])
axes[0, 0].set_title('Cut from end')
# Write simulation to file 
with open('cut_end.csv', 'w') as f:
    f.write('t,n_monomers,n_dimers,n_molecules\n')
    for i, (x, y, z) in enumerate(zip(n_monomers, n_dimers, n_molecules)):
        f.write('{},{},{},{}\n'.format(i + 1, x, y, z))

# ------------------------------------------------------------------- #
# Recognize dimer and cut in the middle
nsim = 10
for n in range(nsim):
    lengths = [total]
    n_monomers = []
    n_dimers = []
    n_molecules = []
    while any(x >= 2 for x in lengths):
        lengths = recognize_dimer_cut_middle(lengths)
        n_monomers.append(sum(1 for x in lengths if x == 1))
        n_dimers.append(sum(1 for x in lengths if x == 2))
        n_molecules.append(len(lengths))
    # Plot simulation
    axes[0, 1].plot(list(range(len(n_monomers))), n_monomers, c=colors[0])
    axes[0, 1].plot(list(range(len(n_dimers))), n_dimers, c=colors[1])
    axes[0, 1].plot(list(range(len(n_molecules))), n_molecules, c=colors[2])
    # Write simulation to file 
    with open('cut_dimer_mid_{}.csv'.format(n), 'w') as f:
        f.write('t,n_monomers,n_dimers,n_molecules\n')
        for i, (x, y, z) in enumerate(zip(n_monomers, n_dimers, n_molecules)):
            f.write('{},{},{},{}\n'.format(i + 1, x, y, z))
axes[0, 1].set_title('Bind dimers, cut in middle')

# ------------------------------------------------------------------- #
# Recognize dimer and cut to the right
for n in range(nsim):
    lengths = [total]
    n_monomers = []
    n_dimers = []
    n_molecules = []
    while any(x >= 3 for x in lengths):
        lengths = recognize_dimer_cut_right(lengths)
        n_monomers.append(sum(1 for x in lengths if x == 1))
        n_dimers.append(sum(1 for x in lengths if x == 2))
        n_molecules.append(len(lengths))
    # Plot simulation
    axes[1, 0].plot(list(range(len(n_monomers))), n_monomers, c=colors[0])
    axes[1, 0].plot(list(range(len(n_dimers))), n_dimers, c=colors[1])
    axes[1, 0].plot(list(range(len(n_molecules))), n_molecules, c=colors[2])
    # Write simulation to file
    with open('cut_dimer_right_{}.csv'.format(n), 'w') as f:
        f.write('t,n_monomers,n_dimers,n_molecules\n')
        for i, (x, y, z) in enumerate(zip(n_monomers, n_dimers, n_molecules)):
            f.write('{},{},{},{}\n'.format(i + 1, x, y, z))
axes[1, 0].set_title('Bind dimers, cut to the right')

# ------------------------------------------------------------------- #
# Recognize trimer and cut in the middle
for n in range(nsim):
    lengths = [total]
    n_monomers = []
    n_dimers = []
    n_molecules = []
    while any(x >= 3 for x in lengths):
        lengths = recognize_trimer_cut_middle(lengths)
        n_monomers.append(sum(1 for x in lengths if x == 1))
        n_dimers.append(sum(1 for x in lengths if x == 2))
        n_molecules.append(len(lengths))
    # Plot simulation
    axes[1, 1].plot(list(range(len(n_monomers))), n_monomers, c=colors[0])
    axes[1, 1].plot(list(range(len(n_dimers))), n_dimers, c=colors[1])
    axes[1, 1].plot(list(range(len(n_molecules))), n_molecules, c=colors[2])
    # Write simulation to file
    with open('cut_trimer_mid_{}.csv'.format(n), 'w') as f:
        f.write('t,n_monomers,n_dimers,n_molecules\n')
        for i, (x, y, z) in enumerate(zip(n_monomers, n_dimers, n_molecules)):
            f.write('{},{},{},{}\n'.format(i + 1, x, y, z))
axes[1, 1].set_title('Bind trimers, cut after 2nd monomer')

for i in range(2):
    for j in range(2):
        axes[i, j].set_xlabel('Number of binding events')
for i in range(2):
    for j in range(2):
        axes[i, j].set_ylabel(r'Number of $k$-mers')
axes[0, 0].legend(['$k = 1$', '$k = 2$', 'All $k$'])

plt.tight_layout()
plt.savefig('bind_vs_time.pdf', transparent=True)

#######################################################################
rng = np.random.default_rng(1234567890)
fig, axes = plt.subplots(nrows=2, ncols=1, figsize=(5, 8))

def parse_distribution(distribution):
    """
    Return a 7-bin "histogram" of the distribution, indicating how many 
    entries are (1) equal to 1, (2) equal to 2, (3) equal to 3, (4) equal
    to 4, (5) between 5 and 10, (6) between 11 and 100, and (7) greater 
    than 100.
    """
    return [
        sum(x == 1 for x in distribution) / len(distribution),
        sum(x == 2 for x in distribution) / len(distribution),
        sum(x == 3 for x in distribution) / len(distribution),
        sum(x == 4 for x in distribution) / len(distribution),
        sum(x > 4 and x <= 10 for x in distribution) / len(distribution),
        sum(x > 10 and x <= 100 for x in distribution) / len(distribution),
        sum(x > 100 for x in distribution) / len(distribution)
    ]

# ------------------------------------------------------------------- #
# Recognize dimer and cut to the right
lengths = [total]
i = 0
distributions = []
multiple = 1000    # Timescale for identifying snapshots of polymer length distribution
distributions.append(lengths)
while any(x >= 3 for x in lengths):
    lengths = recognize_dimer_cut_right(lengths)
    i += 1
    if i % multiple == 0:
        distributions.append(lengths)
distributions.append(lengths)

with open('cut_dimer_right_dist.csv', 'w') as f:
    f.write('t,1,2,3,4,5_to_10,11_to_100,101_to_1000\n')
    for i in range(len(distributions)):
        bars = parse_distribution(distributions[i])
        time = multiple * i
        f.write('{},'.format(time) + ','.join(['{:.10f}'.format(x) for x in bars]) + '\n')

width = 0.11
for i in range(len(distributions)):
    offset = width * i
    axes[0].bar(
        [j + offset for j in range(7)],
        parse_distribution(distributions[i]),
        width, color=colors[i]
    )
axes[0].set_xticks([i - width / 2 + width * len(distributions) / 2 for i in range(7)])
axes[0].set_xticklabels([
    r'$1$', r'$2$', r'$3$', r'$4$',
    r'$> 4,$' + '\n' + r'$\leq 10$',
    r'$> 10,$' + '\n' + r'$\leq 100$',
    r'$> 100,$' + '\n' + r'$\leq 1000$'
])
patches = []
patches.append(mpatches.Patch(color=colors[0], label='Initial'))
for i in range(len(distributions) - 2):
    patches.append(
        mpatches.Patch(color=colors[i+1], label='$t = {}$'.format(multiple * (i + 1)))
    )
patches.append(mpatches.Patch(color=colors[len(distributions)-1], label='Final'))
axes[0].legend(handles=patches, fontsize=9)
axes[0].set_title('Bind dimers, cut to the right')

# Recognize trimer and cut in the middle
lengths = [total]
i = 0
distributions = []
multiple = 100    # Timescale for identifying snapshots of polymer length distribution
distributions.append(lengths)
while any(x >= 3 for x in lengths):
    lengths = recognize_trimer_cut_middle(lengths)
    i += 1
    if i % multiple == 0:
        distributions.append(lengths)
distributions.append(lengths)

with open('cut_trimer_mid_dist.csv', 'w') as f:
    f.write('t,1,2,3,4,5_to_10,11_to_100,101_to_1000\n')
    for i in range(len(distributions)):
        bars = parse_distribution(distributions[i])
        time = multiple * i
        f.write('{},'.format(time) + ','.join(['{:.10f}'.format(x) for x in bars]) + '\n')

width = 0.11
for i in range(len(distributions)):
    offset = width * i
    axes[1].bar(
        [j + offset for j in range(7)],
        parse_distribution(distributions[i]),
        width, color=colors[i]
    )
axes[1].set_xticks([i - width / 2 + width * len(distributions) / 2 for i in range(7)])
axes[1].set_xticklabels([
    r'$1$', r'$2$', r'$3$', r'$4$',
    r'$> 4,$' + '\n' + r'$\leq 10$',
    r'$> 10,$' + '\n' + r'$\leq 100$',
    r'$> 100,$' + '\n' + r'$\leq 1000$'
])
patches = []
patches.append(mpatches.Patch(color=colors[0], label='Initial'))
for i in range(len(distributions) - 2):
    patches.append(
        mpatches.Patch(color=colors[i+1], label='$t = {}$'.format(multiple * (i + 1)))
    )
patches.append(mpatches.Patch(color=colors[len(distributions)-1], label='Final'))
axes[1].legend(handles=patches, fontsize=9)
axes[1].set_title('Bind trimers, cut after 2nd monomer')

axes[0].set_xlabel('Polymer length')
axes[1].set_xlabel('Polymer length')
axes[0].set_ylabel('Proportion')
axes[1].set_ylabel('Proportion')

plt.tight_layout()
plt.savefig('distributions.pdf', transparent=True)

