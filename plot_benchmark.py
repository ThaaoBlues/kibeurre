import matplotlib.pyplot as plt
from pandas import read_csv

# Load your data
encryption_data = read_csv("encryption_benchmark.csv", sep=",")
decryption_data = read_csv("decryption_benchmark.csv", sep=",")

plt.figure(figsize=(10, 6))

# Plot original lines
plt.plot(
    encryption_data["input_length"],
    encryption_data["time_taken"],
    label="Encryption",
)
plt.plot(
    decryption_data["input_length"],
    decryption_data["time_taken"],
    label="Decryption",
)


# --- Step Detection Logic ---
def mark_steps(data, label_prefix, color):
    # 1. Calculate the difference between consecutive time values
    time_diff = data["time_taken"].diff()

    # 2. Define a threshold for what constitutes a "step"
    threshold = 1e-4
    max_jump = 1e-3

    # add maximum constraint to avoid taking cpu spikes into account

    # 3. Find rows where the jump is greater than the threshold
    step_indices = data[(time_diff > threshold) & (time_diff < max_jump)]  # Adjust the max time constraint as needed

    # 4. Plot vertical lines at these indices
    for idx, row in step_indices.iterrows():
        plt.axvline(
            x=row["input_length"],
            color=color,
            linestyle="--",
            alpha=0.6,
        )
        # Optional: Label the exact X-value on the plot
        plt.text(
            row["input_length"],
            row["time_taken"],
            f"{int(row['input_length'])}",
            rotation=90,
            verticalalignment="bottom",
            fontsize=9,
            color=color,
        )

    # Print to console so you have the exact numbers
    print(
        f"{label_prefix} steps occur at input lengths: {step_indices['input_length'].tolist()}"
    )


# Detect and mark steps for both datasets
mark_steps(encryption_data, "Encryption", "blue")
mark_steps(decryption_data, "Decryption", "orange")
# ----------------------------

plt.xlabel("Input Length")
plt.ylabel("Time (seconds)")
plt.title("Benchmark Results with Step Detection")
plt.legend()
plt.tight_layout()
plt.savefig("benchmark_results.png")
plt.show()