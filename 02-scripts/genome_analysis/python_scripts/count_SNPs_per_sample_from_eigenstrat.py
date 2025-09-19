# Set your file path and number of samples
geno_file = "pileupcaller.single.geno"
num_samples = 16  # Change this to your actual number of samples

# Initialize counters
diff_counts = [0] * num_samples
called_counts = [0] * num_samples
total_positions = 0

# Read geno file and count per sample
with open(geno_file, "r") as f:
    for line in f:
        line = line.strip()
        total_positions += 1
        for i in range(num_samples):
            genotype = line[i]
            if genotype != "9":
                called_counts[i] += 1
                if genotype in ("1", "2"):
                    diff_counts[i] += 1

# Output results
print("Sample\tDifferent_from_REF\tCalled_sites\tTotal_positions\tFraction_different")
for i in range(num_samples):
    diff = diff_counts[i]
    called = called_counts[i]
    frac = diff / called if called > 0 else 0
    print(f"Sample_{i+1}\t{diff}\t{called}\t{total_positions}\t{frac:.4f}")

