import os


def write_statistics(
    sim_counts,
    paf_counts,
    sim_stats,
    paf_stats,
    bucket_size=1000,
    output_dir=None,
    filename="statistics.txt"
):
    """
    Sprema osnovnu statistiku:
    - simulator / FASTQ
    - mapiranje / inicijalna dodjela

    Ovo je osnovni izvještaj koji ide u:
        statistika/osnovno/statistics.txt
    """
    if output_dir is None:
        output_dir = f"results/bucket{bucket_size}/statistika/osnovno"

    os.makedirs(output_dir, exist_ok=True)

    output_path = os.path.join(output_dir, filename)

    all_genomes = sorted(
        set(sim_counts)
        | set(paf_counts)
        | set(sim_stats)
        | set(paf_stats)
    )

    with open(output_path, "w") as f:
        f.write("STATISTIKA POKRIVENOSTI I RASPODJELE OCITANJA\n")
        f.write("=" * 70 + "\n\n")
        f.write(f"Bucket size: {bucket_size} bp\n\n")

        for genome_id in all_genomes:
            f.write(f"Genom: {genome_id}\n")
            f.write("-" * 70 + "\n")

            sim_count = sim_counts.get(genome_id, 0)
            paf_count = paf_counts.get(genome_id, 0)

            f.write("Broj ocitanja:\n")
            f.write(f"        Simulator / FASTQ            : {sim_count:.2f}\n")
            f.write(f"  Mapiranje / inicijalna dodjela     : {paf_count:.2f}\n")
            f.write(f"         Razlika PAF-SIM             : {paf_count - sim_count:.2f}\n\n")

            if genome_id in sim_stats:
                s = sim_stats[genome_id]
                f.write("Coverage statistika - simulator:\n")
                f.write(f"  broj bucket-a          : {s['num_buckets']}\n")
                f.write(f"  mean                   : {s['mean']:.4f}\n")
                f.write(f"  std                    : {s['std']:.4f}\n")
                f.write(f"  cv                     : {s['cv']:.4f}\n")
                f.write(f"  min                    : {s['min']:.4f}\n")
                f.write(f"  max                    : {s['max']:.4f}\n")
                f.write(f"  range max-min          : {s['range']:.4f}\n")
                f.write(f"  zero buckets           : {s['zero_buckets']}\n")
                f.write(f"  zero fraction          : {s['zero_fraction']:.4f}\n")
                f.write(f"  high outliers          : {s['high_outliers']}\n")
                f.write(f"  low outliers           : {s['low_outliers']}\n")
                f.write(f"  high outlier fraction  : {s['high_outlier_fraction']:.4f}\n")
                f.write(f"  low outlier fraction   : {s['low_outlier_fraction']:.4f}\n")
                f.write(f"  score cv+zero_fraction : {s['score_cv_zero']:.4f}\n\n")

            if genome_id in paf_stats:
                s = paf_stats[genome_id]
                f.write("Coverage statistika - mapiranje / inicijalna dodjela:\n")
                f.write(f"  broj bucket-a          : {s['num_buckets']}\n")
                f.write(f"  mean                   : {s['mean']:.4f}\n")
                f.write(f"  std                    : {s['std']:.4f}\n")
                f.write(f"  cv                     : {s['cv']:.4f}\n")
                f.write(f"  min                    : {s['min']:.4f}\n")
                f.write(f"  max                    : {s['max']:.4f}\n")
                f.write(f"  range max-min          : {s['range']:.4f}\n")
                f.write(f"  zero buckets           : {s['zero_buckets']}\n")
                f.write(f"  zero fraction          : {s['zero_fraction']:.4f}\n")
                f.write(f"  high outliers          : {s['high_outliers']}\n")
                f.write(f"  low outliers           : {s['low_outliers']}\n")
                f.write(f"  high outlier fraction  : {s['high_outlier_fraction']:.4f}\n")
                f.write(f"  low outlier fraction   : {s['low_outlier_fraction']:.4f}\n")
                f.write(f"  score cv+zero_fraction : {s['score_cv_zero']:.4f}\n\n")

            if genome_id in sim_stats and genome_id in paf_stats:
                sim = sim_stats[genome_id]
                paf = paf_stats[genome_id]

                f.write("Usporedba mapiranje / inicijalna dodjela - simulator:\n")
                f.write(f"  mean razlika                   : {paf['mean'] - sim['mean']:.4f}\n")
                f.write(f"  std razlika                    : {paf['std'] - sim['std']:.4f}\n")
                f.write(f"  cv razlika                     : {paf['cv'] - sim['cv']:.4f}\n")
                f.write(f"  range razlika                  : {paf['range'] - sim['range']:.4f}\n")
                f.write(f"  zero buckets razlika           : {paf['zero_buckets'] - sim['zero_buckets']}\n")
                f.write(f"  zero fraction razlika          : {paf['zero_fraction'] - sim['zero_fraction']:.4f}\n")
                f.write(f"  high outliers razlika          : {paf['high_outliers'] - sim['high_outliers']}\n")
                f.write(f"  low outliers razlika           : {paf['low_outliers'] - sim['low_outliers']}\n")
                f.write(
                    "  score cv+zero_fraction razlika : "
                    f"{paf['score_cv_zero'] - sim['score_cv_zero']:.4f}\n\n"
                )

            f.write("\n")

    return output_path