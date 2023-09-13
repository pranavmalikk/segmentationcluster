import matplotlib.pyplot as plt

def read_values_from_file(filename):
    with open(filename, 'r') as f:
        return [float(line.strip()) for line in f]

if __name__ == "__main__":
    all_stitj = read_values_from_file('all_stitj_non_normalize.txt')

    plt.hist(all_stitj, bins=100, edgecolor='black')
    plt.legend(loc='upper right')
    plt.title('Distribution of Stitj values')
    plt.xlabel('Stitj')
    plt.ylabel('Frequency')
    plt.show()