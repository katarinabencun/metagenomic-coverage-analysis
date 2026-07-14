import os
import numpy as np


def coverage_distance_to_simulator(
    coverage_sim,
    coverage_initial,
    coverage_final,
    genome_lengths,
    bucket_size=1000
):
    """
    Računa udaljenost initial i final coverage grafa od simulator coverage grafa.

    Mjere:
    - MAE = mean absolute error po bucketima
    - RMSE = root mean squared error po bucketima
    - max_abs_error = najveća apsolutna razlika po bucketu

    Računa se samo za genome koji postoje u simulator coverageu,
    jer samo za njih imamo stvarnu referencu.
    """
    distances = {}

    for genome_id, sim_buckets in coverage_sim.items():
        if genome_id not in genome_lengths:
            continue

        num_buckets = (genome_lengths[genome_id] + bucket_size - 1) // bucket_size

        sim_values = np.array(
            [
                coverage_sim.get(genome_id, {}).get(i, 0.0)
                for i in range(num_buckets)
            ],
            dtype=float
        )

        initial_values = np.array(
            [
                coverage_initial.get(genome_id, {}).get(i, 0.0)
                for i in range(num_buckets)
            ],
            dtype=float
        )

        final_values = np.array(
            [
                coverage_final.get(genome_id, {}).get(i, 0.0)
                for i in range(num_buckets)
            ],
            dtype=float
        )

        initial_diff = initial_values - sim_values
        final_diff = final_values - sim_values

        initial_mae = float(np.mean(np.abs(initial_diff)))
        final_mae = float(np.mean(np.abs(final_diff)))

        initial_rmse = float(np.sqrt(np.mean(initial_diff ** 2)))
        final_rmse = float(np.sqrt(np.mean(final_diff ** 2)))

        initial_max_abs_error = float(np.max(np.abs(initial_diff)))
        final_max_abs_error = float(np.max(np.abs(final_diff)))

        distances[genome_id] = {
            "num_buckets": num_buckets,

            "initial_mae": initial_mae,
            "final_mae": final_mae,
            "mae_change": final_mae - initial_mae,

            "initial_rmse": initial_rmse,
            "final_rmse": final_rmse,
            "rmse_change": final_rmse - initial_rmse,

            "initial_max_abs_error": initial_max_abs_error,
            "final_max_abs_error": final_max_abs_error,
            "max_abs_error_change": final_max_abs_error - initial_max_abs_error,
        }

    return distances


def write_coverage_distance_to_simulator(
    coverage_sim,
    coverage_initial,
    coverage_final,
    genome_lengths,
    bucket_size=1000,
    output_dir=None,
    filename="coverage_distance_to_simulator.txt"
):
    """
    Sprema usporedbu oblika coverage grafa:
    simulator vs initial
    simulator vs final

    Manji MAE/RMSE znači da je graf bliži simulatoru.
    """
    if output_dir is None:
        output_dir = f"results/bucket{bucket_size}/statistika/dodatno"

    os.makedirs(output_dir, exist_ok=True)

    output_path = os.path.join(output_dir, filename)

    distances = coverage_distance_to_simulator(
        coverage_sim=coverage_sim,
        coverage_initial=coverage_initial,
        coverage_final=coverage_final,
        genome_lengths=genome_lengths,
        bucket_size=bucket_size
    )

    total_initial_mae = 0.0
    total_final_mae = 0.0
    total_initial_rmse = 0.0
    total_final_rmse = 0.0
    genome_count = 0

    with open(output_path, "w") as f:
        f.write("UDALJENOST COVERAGE GRAFA OD SIMULATORA PO BUCKETIMA\n")
        f.write("=" * 80 + "\n\n")
        f.write(f"Bucket size: {bucket_size} bp\n\n")

        for genome_id, d in sorted(distances.items()):
            genome_count += 1

            total_initial_mae += d["initial_mae"]
            total_final_mae += d["final_mae"]
            total_initial_rmse += d["initial_rmse"]
            total_final_rmse += d["final_rmse"]

            f.write(f"Genom: {genome_id}\n")
            f.write("-" * 80 + "\n")
            f.write(f"Broj bucket-a: {d['num_buckets']}\n\n")

            f.write("MAE - mean absolute error\n")
            f.write(f"  initial vs simulator : {d['initial_mae']:.4f}\n")
            f.write(f"  final vs simulator   : {d['final_mae']:.4f}\n")
            f.write(f"  promjena             : {d['mae_change']:.4f}\n")

            if d["mae_change"] < 0:
                f.write("  zakljucak            : final graf je blizi simulatoru po MAE\n\n")
            elif d["mae_change"] > 0:
                f.write("  zakljucak            : final graf je dalji od simulatora po MAE\n\n")
            else:
                f.write("  zakljucak            : nema promjene po MAE\n\n")

            f.write("RMSE - root mean squared error\n")
            f.write(f"  initial vs simulator : {d['initial_rmse']:.4f}\n")
            f.write(f"  final vs simulator   : {d['final_rmse']:.4f}\n")
            f.write(f"  promjena             : {d['rmse_change']:.4f}\n")

            if d["rmse_change"] < 0:
                f.write("  zakljucak            : final graf je blizi simulatoru po RMSE\n\n")
            elif d["rmse_change"] > 0:
                f.write("  zakljucak            : final graf je dalji od simulatora po RMSE\n\n")
            else:
                f.write("  zakljucak            : nema promjene po RMSE\n\n")

            f.write("Najveca apsolutna razlika po bucketu\n")
            f.write(f"  initial vs simulator : {d['initial_max_abs_error']:.4f}\n")
            f.write(f"  final vs simulator   : {d['final_max_abs_error']:.4f}\n")
            f.write(f"  promjena             : {d['max_abs_error_change']:.4f}\n\n")

            f.write("\n")

        f.write("=" * 80 + "\n")
        f.write("UKUPNI SAZETAK ZA GENOME KOJI POSTOJE U SIMULATORU\n")
        f.write("=" * 80 + "\n\n")

        if genome_count > 0:
            avg_initial_mae = total_initial_mae / genome_count
            avg_final_mae = total_final_mae / genome_count

            avg_initial_rmse = total_initial_rmse / genome_count
            avg_final_rmse = total_final_rmse / genome_count

            f.write("Prosjek po simuliranim genomima:\n")
            f.write(f"  initial MAE prosjek : {avg_initial_mae:.4f}\n")
            f.write(f"  final MAE prosjek   : {avg_final_mae:.4f}\n")
            f.write(f"  MAE promjena        : {avg_final_mae - avg_initial_mae:.4f}\n\n")

            f.write(f"  initial RMSE prosjek: {avg_initial_rmse:.4f}\n")
            f.write(f"  final RMSE prosjek  : {avg_final_rmse:.4f}\n")
            f.write(f"  RMSE promjena       : {avg_final_rmse - avg_initial_rmse:.4f}\n\n")

            if avg_final_mae < avg_initial_mae:
                f.write("Zakljucak po MAE: finalni grafovi su u prosjeku blizi simulatoru.\n")
            elif avg_final_mae > avg_initial_mae:
                f.write("Zakljucak po MAE: finalni grafovi su u prosjeku dalji od simulatora.\n")
            else:
                f.write("Zakljucak po MAE: nema promjene.\n")

            if avg_final_rmse < avg_initial_rmse:
                f.write("Zakljucak po RMSE: finalni grafovi su u prosjeku blizi simulatoru.\n")
            elif avg_final_rmse > avg_initial_rmse:
                f.write("Zakljucak po RMSE: finalni grafovi su u prosjeku dalji od simulatora.\n")
            else:
                f.write("Zakljucak po RMSE: nema promjene.\n")
        else:
            f.write("Nema simuliranih genoma za usporedbu.\n")

    return output_path