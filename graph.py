import pandas as pd
import matplotlib.pyplot as plt

# Read the CSV files
areas = pd.read_csv("areas.csv", header=None).iloc[:, 0].tolist()
maxEdges = pd.read_csv("maxEdges.csv", header=None).iloc[:, 0].tolist()
minEdges = pd.read_csv("minEdges.csv", header=None).iloc[:, 0].tolist()

# Plot histograms
plt.figure(figsize=(12, 4))

plt.subplot(1, 3, 1)
plt.title('Areas')
plt.hist(areas, bins=10000, edgecolor='black')
plt.xlabel('Area')
plt.ylabel('Frequency')

plt.subplot(1, 3, 2)
plt.title('Max Edges')
plt.hist(maxEdges, bins=10000, edgecolor='black')
plt.xlabel('Max Edge')
plt.ylabel('Frequency')

plt.subplot(1, 3, 3)
plt.title('Min Edges')
plt.hist(minEdges, bins=10000, edgecolor='black')
plt.xlabel('Min Edge')
plt.ylabel('Frequency')

plt.tight_layout()
plt.show()
