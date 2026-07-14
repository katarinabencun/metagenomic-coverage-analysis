import numpy as np

def coverage_stats(coverage, genome_lengths, bucket_size=1000):
    """
    Računa osnovne statistike pokrivenosti po genomu.

    U coverageu se mogu nalaziti samo bucket-i koji imaju neku pokrivenost,
    ali ovdje se eksplicitno uzimaju i nulti bucket-i do pune duljine genoma.
    """
    stats = {}

    for genome_id, buckets in coverage.items():
        if genome_id not in genome_lengths:
            continue

        num_buckets = (genome_lengths[genome_id] + bucket_size - 1) // bucket_size

        values = np.array(
            [buckets.get(i, 0) for i in range(num_buckets)],
            dtype=float
        )

        mean = float(np.mean(values))
        std = float(np.std(values))
        cv = float(std / mean) if mean > 0 else 0.0

        min_value = float(np.min(values))
        max_value = float(np.max(values))
        coverage_range = max_value - min_value

        zero_buckets = int(np.sum(values == 0))
        zero_fraction = float(zero_buckets / num_buckets) if num_buckets > 0 else 0.0

        high_threshold = mean + 2 * std
        low_threshold = max(0.0, mean - 2 * std)

        high_outliers = int(np.sum(values > high_threshold))
        low_outliers = int(np.sum(values < low_threshold))

        high_outlier_fraction = (
            float(high_outliers / num_buckets)
            if num_buckets > 0
            else 0.0
        )

        low_outlier_fraction = (
            float(low_outliers / num_buckets)
            if num_buckets > 0
            else 0.0
        )

        # Jednostavna kompozitna mjera:
        # manji score znači ravnomjerniji i potpuniji coverage.
        score_cv_zero = cv + zero_fraction

        stats[genome_id] = {
            "num_buckets": num_buckets,
            "mean": mean,
            "std": std,
            "cv": cv,
            "min": min_value,
            "max": max_value,
            "range": coverage_range,
            "zero_buckets": zero_buckets,
            "zero_fraction": zero_fraction,
            "high_outliers": high_outliers,
            "low_outliers": low_outliers,
            "high_outlier_fraction": high_outlier_fraction,
            "low_outlier_fraction": low_outlier_fraction,
            "score_cv_zero": score_cv_zero,
        }

    return stats


def full_coverage_values_for_genome(
    coverage,
    genome_id,
    genome_lengths,
    bucket_size=1000
):
    """
    Vraća puni niz coverage vrijednosti za jedan genom,
    uključujući bucket-e s nulom.
    """
    num_buckets = (genome_lengths[genome_id] + bucket_size - 1) // bucket_size

    return np.array(
        [
            coverage.get(genome_id, {}).get(i, 0.0)
            for i in range(num_buckets)
        ],
        dtype=float
    )


def max_zero_run_fraction_from_values(values):
    """
    Računa udio najduljeg kontinuiranog niza nultih bucket-a.

    Primjer:
        0.25 znači da najdulja rupa pokriva 25% bucket-a genoma.
    """
    if len(values) == 0:
        return 0.0

    max_run = 0
    current_run = 0

    for value in values:
        if value <= 1e-12:
            current_run += 1
            max_run = max(max_run, current_run)
        else:
            current_run = 0

    return max_run / len(values)