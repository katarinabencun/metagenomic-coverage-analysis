import os
import numpy as np

from statistika.coverage_metrics import (
    full_coverage_values_for_genome,
    max_zero_run_fraction_from_values,
)


def false_genome_coverage_stats(
    coverage_sim,
    coverage_initial,
    coverage_final,
    genome_lengths,
    bucket_size=1000
):
    """
    Računa statistiku za genome koji nisu prisutni u simulatoru.

    Za njih je idealni simulator coverage jednak 0 po cijelom genomu.

    Za lažne genome gledamo:
    - mean coverage
    - RMSE od nule
    - max coverage
    - zero fraction
    - max zero run fraction
    """
    false_stats = {}

    simulated_genomes = set(coverage_sim.keys())

    false_genomes = [
        genome_id
        for genome_id in genome_lengths
        if genome_id not in simulated_genomes
    ]

    for genome_id in sorted(false_genomes):
        initial_values = full_coverage_values_for_genome(
            coverage_initial,
            genome_id,
            genome_lengths,
            bucket_size
        )

        final_values = full_coverage_values_for_genome(
            coverage_final,
            genome_id,
            genome_lengths,
            bucket_size
        )

        if len(initial_values) == 0:
            continue

        initial_mean = float(np.mean(initial_values))
        final_mean = float(np.mean(final_values))

        initial_rmse_from_zero = float(np.sqrt(np.mean(initial_values ** 2)))
        final_rmse_from_zero = float(np.sqrt(np.mean(final_values ** 2)))

        initial_max = float(np.max(initial_values))
        final_max = float(np.max(final_values))

        initial_zero_fraction = float(np.mean(initial_values <= 1e-12))
        final_zero_fraction = float(np.mean(final_values <= 1e-12))

        initial_max_zero_run_fraction = max_zero_run_fraction_from_values(initial_values)
        final_max_zero_run_fraction = max_zero_run_fraction_from_values(final_values)

        false_stats[genome_id] = {
            "num_buckets": len(initial_values),

            "initial_mean": initial_mean,
            "final_mean": final_mean,
            "mean_change": final_mean - initial_mean,

            "initial_rmse_from_zero": initial_rmse_from_zero,
            "final_rmse_from_zero": final_rmse_from_zero,
            "rmse_from_zero_change": final_rmse_from_zero - initial_rmse_from_zero,

            "initial_max": initial_max,
            "final_max": final_max,
            "max_change": final_max - initial_max,

            "initial_zero_fraction": initial_zero_fraction,
            "final_zero_fraction": final_zero_fraction,
            "zero_fraction_change": final_zero_fraction - initial_zero_fraction,

            "initial_max_zero_run_fraction": initial_max_zero_run_fraction,
            "final_max_zero_run_fraction": final_max_zero_run_fraction,
            "max_zero_run_fraction_change": (
                final_max_zero_run_fraction - initial_max_zero_run_fraction
            ),
        }

    return false_stats


def write_false_genome_coverage_stats(
    coverage_sim,
    coverage_initial,
    coverage_final,
    genome_lengths,
    bucket_size=1000,
    output_dir=None,
    filename="false_genome_coverage_stats.txt"
):
    """
    Sprema statistiku coveragea za lažne genome.

    Lažni genomi su oni koji postoje u referentnoj bazi,
    ali nemaju očitanja u simulatoru.

    Za njih je idealno:
    - mean coverage što bliže 0
    - RMSE od nule što bliže 0
    - max coverage što manji
    - zero fraction što veći
    """
    if output_dir is None:
        output_dir = f"results/bucket{bucket_size}/statistika/dodatno"

    os.makedirs(output_dir, exist_ok=True)

    output_path = os.path.join(output_dir, filename)

    stats = false_genome_coverage_stats(
        coverage_sim=coverage_sim,
        coverage_initial=coverage_initial,
        coverage_final=coverage_final,
        genome_lengths=genome_lengths,
        bucket_size=bucket_size
    )

    total_initial_mean = 0.0
    total_final_mean = 0.0
    total_initial_rmse = 0.0
    total_final_rmse = 0.0
    total_initial_max = 0.0
    total_final_max = 0.0

    genome_count = 0

    with open(output_path, "w") as f:
        f.write("STATISTIKA COVERAGEA ZA LAZNE GENOME\n")
        f.write("=" * 80 + "\n\n")
        f.write(f"Bucket size: {bucket_size} bp\n\n")

        if not stats:
            f.write("Nema laznih genoma za analizu.\n")
            return output_path

        for genome_id, s in sorted(stats.items()):
            genome_count += 1

            total_initial_mean += s["initial_mean"]
            total_final_mean += s["final_mean"]

            total_initial_rmse += s["initial_rmse_from_zero"]
            total_final_rmse += s["final_rmse_from_zero"]

            total_initial_max += s["initial_max"]
            total_final_max += s["final_max"]

            f.write(f"Genom: {genome_id}\n")
            f.write("-" * 80 + "\n")
            f.write(f"Broj bucket-a: {s['num_buckets']}\n\n")

            f.write("MEAN COVERAGE - idealno sto blize 0\n")
            f.write(f"  initial mean : {s['initial_mean']:.4f}\n")
            f.write(f"  final mean   : {s['final_mean']:.4f}\n")
            f.write(f"  promjena     : {s['mean_change']:.4f}\n")

            if s["mean_change"] < 0:
                f.write("  zakljucak    : final ima manje prosjecnog lazno dodijeljenog coveragea\n\n")
            elif s["mean_change"] > 0:
                f.write("  zakljucak    : final ima vise prosjecnog lazno dodijeljenog coveragea\n\n")
            else:
                f.write("  zakljucak    : nema promjene mean coveragea\n\n")

            f.write("RMSE OD NULE - idealno sto blize 0\n")
            f.write(f"  initial RMSE : {s['initial_rmse_from_zero']:.4f}\n")
            f.write(f"  final RMSE   : {s['final_rmse_from_zero']:.4f}\n")
            f.write(f"  promjena     : {s['rmse_from_zero_change']:.4f}\n")

            if s["rmse_from_zero_change"] < 0:
                f.write("  zakljucak    : final je blizi idealnoj nuli\n\n")
            elif s["rmse_from_zero_change"] > 0:
                f.write("  zakljucak    : final je dalji od idealne nule\n\n")
            else:
                f.write("  zakljucak    : nema promjene RMSE od nule\n\n")

            f.write("MAX COVERAGE - idealno sto manje\n")
            f.write(f"  initial max  : {s['initial_max']:.4f}\n")
            f.write(f"  final max    : {s['final_max']:.4f}\n")
            f.write(f"  promjena     : {s['max_change']:.4f}\n\n")

            f.write("ZERO FRACTION - za lazne genome idealno sto vece\n")
            f.write(f"  initial zero fraction : {s['initial_zero_fraction']:.4f}\n")
            f.write(f"  final zero fraction   : {s['final_zero_fraction']:.4f}\n")
            f.write(f"  promjena              : {s['zero_fraction_change']:.4f}\n\n")

            f.write("MAX ZERO RUN FRACTION - najduza kontinuirana rupa\n")
            f.write(
                "  initial max zero run fraction : "
                f"{s['initial_max_zero_run_fraction']:.4f}\n"
            )
            f.write(
                "  final max zero run fraction   : "
                f"{s['final_max_zero_run_fraction']:.4f}\n"
            )
            f.write(
                "  promjena                      : "
                f"{s['max_zero_run_fraction_change']:.4f}\n\n"
            )

            f.write("\n")

        f.write("=" * 80 + "\n")
        f.write("UKUPNI SAZETAK ZA LAZNE GENOME\n")
        f.write("=" * 80 + "\n\n")

        avg_initial_mean = total_initial_mean / genome_count
        avg_final_mean = total_final_mean / genome_count

        avg_initial_rmse = total_initial_rmse / genome_count
        avg_final_rmse = total_final_rmse / genome_count

        avg_initial_max = total_initial_max / genome_count
        avg_final_max = total_final_max / genome_count

        f.write("Prosjek po laznim genomima:\n")
        f.write(f"  initial mean coverage prosjek : {avg_initial_mean:.4f}\n")
        f.write(f"  final mean coverage prosjek   : {avg_final_mean:.4f}\n")
        f.write(f"  mean promjena                 : {avg_final_mean - avg_initial_mean:.4f}\n\n")

        f.write(f"  initial RMSE od nule prosjek  : {avg_initial_rmse:.4f}\n")
        f.write(f"  final RMSE od nule prosjek    : {avg_final_rmse:.4f}\n")
        f.write(f"  RMSE promjena                 : {avg_final_rmse - avg_initial_rmse:.4f}\n\n")

        f.write(f"  initial max coverage prosjek  : {avg_initial_max:.4f}\n")
        f.write(f"  final max coverage prosjek    : {avg_final_max:.4f}\n")
        f.write(f"  max promjena                  : {avg_final_max - avg_initial_max:.4f}\n\n")

        if avg_final_mean < avg_initial_mean:
            f.write("Zakljucak po mean coverageu: final smanjuje prosjecni coverage na laznim genomima.\n")
        elif avg_final_mean > avg_initial_mean:
            f.write("Zakljucak po mean coverageu: final povecava prosjecni coverage na laznim genomima.\n")
        else:
            f.write("Zakljucak po mean coverageu: nema promjene.\n")

        if avg_final_rmse < avg_initial_rmse:
            f.write("Zakljucak po RMSE od nule: final je blizi idealnoj nuli na laznim genomima.\n")
        elif avg_final_rmse > avg_initial_rmse:
            f.write("Zakljucak po RMSE od nule: final je dalji od idealne nule na laznim genomima.\n")
        else:
            f.write("Zakljucak po RMSE od nule: nema promjene.\n")

    return output_path