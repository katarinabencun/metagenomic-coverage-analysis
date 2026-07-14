import os
import matplotlib.pyplot as plt


def build_full_coverage_array(buckets, genome_length, bucket_size=1000):
    num_buckets = (genome_length + bucket_size - 1) // bucket_size
    xs = list(range(num_buckets))
    ys = [buckets.get(x, 0) for x in xs]
    return xs, ys

def plot_coverage_profile_kbp(
    coverage,
    genome_lengths,
    bucket_size=1000,
    output_prefix="sim",
    output_dir="results/profile_kbp"
):
    os.makedirs(output_dir, exist_ok=True)

    for genome_id, buckets in coverage.items():
        if genome_id not in genome_lengths:
            print(f"Upozorenje: genome_length nije pronađen za {genome_id}")
            continue

        xs, ys = build_full_coverage_array(
            buckets,
            genome_lengths[genome_id],
            bucket_size=bucket_size
        )

        # pretvaranje bucket indexa u poziciju na genomu u kbp
        xs_kbp = [x * bucket_size / 1000 for x in xs]

        plt.figure(figsize=(15, 5))
        plt.fill_between(xs_kbp, ys, alpha=0.35)
        plt.plot(xs_kbp, ys, linewidth=1)

        plt.title(f"Profil pokrivenosti: {genome_id} ({output_prefix})")
        plt.xlabel("Pozicija u genomu [kbp]")
        plt.ylabel("Prosječna pokrivenost po bazi")
        plt.ylim(bottom=0)
        plt.grid(axis="y", linestyle="--", alpha=0.6)
        plt.tight_layout()

        plt.savefig(f"{output_dir}/{genome_id}_{output_prefix}_profile_kbp.png", dpi=150)
        plt.close()


def plot_coverage_comparison_stacked_kbp(
    coverage_sim,
    coverage_initial,
    coverage_final,
    genome_lengths,
    bucket_size=1000,
    output_dir="results/comparison"
):
    os.makedirs(output_dir, exist_ok=True)

    all_genomes = sorted(
        set(genome_lengths.keys())
        | set(coverage_sim.keys())
        | set(coverage_initial.keys())
        | set(coverage_final.keys())
    )

    for genome_id in all_genomes:
        if genome_id not in genome_lengths:
            print(f"Upozorenje: genome_length nije pronađen za {genome_id}")
            continue

        genome_length = genome_lengths[genome_id]
        num_buckets = (genome_length + bucket_size - 1) // bucket_size

        xs = list(range(num_buckets))
        xs_kbp = [x * bucket_size / 1000 for x in xs]

        sim_buckets = coverage_sim.get(genome_id, {})
        initial_buckets = coverage_initial.get(genome_id, {})
        final_buckets = coverage_final.get(genome_id, {})

        ys_sim = [sim_buckets.get(x, 0) for x in xs]
        ys_initial = [initial_buckets.get(x, 0) for x in xs]
        ys_final = [final_buckets.get(x, 0) for x in xs]

        ymax = max(
            max(ys_sim) if ys_sim else 0,
            max(ys_initial) if ys_initial else 0,
            max(ys_final) if ys_final else 0
        )

        if ymax == 0:
            ymax = 1
        else:
            ymax *= 1.05

        fig, axes = plt.subplots(3, 1, figsize=(15, 10), sharex=True)

        # 1. Simulator
        axes[0].fill_between(xs_kbp, ys_sim, alpha=0.35)
        axes[0].plot(xs_kbp, ys_sim, linewidth=1)
        axes[0].set_title("Simulator")
        axes[0].set_ylabel("Pokrivenost")
        axes[0].set_ylim(0, ymax)
        axes[0].grid(axis="y", linestyle="--", alpha=0.6)

        # 2. Inicijalna dodjela
        axes[1].fill_between(xs_kbp, ys_initial, alpha=0.35)
        axes[1].plot(xs_kbp, ys_initial, linewidth=1)
        axes[1].set_title("Inicijalna dodjela")
        axes[1].set_ylabel("Pokrivenost")
        axes[1].set_ylim(0, ymax)
        axes[1].grid(axis="y", linestyle="--", alpha=0.6)

        # 3. Finalna preraspodjela
        axes[2].fill_between(xs_kbp, ys_final, alpha=0.35)
        axes[2].plot(xs_kbp, ys_final, linewidth=1)
        axes[2].set_title("Finalna preraspodjela")
        axes[2].set_xlabel("Pozicija u genomu [kbp]")
        axes[2].set_ylabel("Pokrivenost")
        axes[2].set_ylim(0, ymax)
        axes[2].grid(axis="y", linestyle="--", alpha=0.6)

        fig.suptitle(f"Usporedba pokrivenosti: {genome_id}", fontsize=14)
        plt.tight_layout()
        plt.subplots_adjust(top=0.92)

        output_path = os.path.join(output_dir, f"{genome_id}_usporedba.png")
        plt.savefig(output_path, dpi=150)
        plt.close()