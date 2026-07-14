import os


def write_redistribution_comparison(
    sim_counts,
    initial_counts,
    final_counts,
    sim_stats,
    initial_stats,
    final_stats,
    bucket_size=1000,
    output_dir=None,
    filename="redistribution_comparison.txt"
):
    """
    Sprema usporedbu:
    - simulator / stvarno stanje
    - inicijalna dodjela
    - finalna preraspodjela nakon algoritma

    Cilj:
        vidjeti je li algoritam poboljšao raspodjelu u odnosu
        na inicijalno stanje.
    """
    if output_dir is None:
        output_dir = f"results/bucket{bucket_size}/statistika/dodatno"

    os.makedirs(output_dir, exist_ok=True)

    output_path = os.path.join(output_dir, filename)

    all_genomes = sorted(
        set(sim_counts)
        | set(initial_counts)
        | set(final_counts)
        | set(sim_stats)
        | set(initial_stats)
        | set(final_stats)
    )

    total_abs_read_error_initial = 0.0
    total_abs_read_error_final = 0.0

    total_abs_mean_error_initial = 0.0
    total_abs_mean_error_final = 0.0

    total_abs_score_error_initial = 0.0
    total_abs_score_error_final = 0.0

    genomes_with_sim_stats = 0

    with open(output_path, "w") as f:
        f.write("USPOREDBA INICIJALNE DODJELE I FINALNE PRERASPODJELE\n")
        f.write("=" * 80 + "\n\n")
        f.write(f"Bucket size: {bucket_size} bp\n\n")

        for genome_id in all_genomes:
            sim_count = sim_counts.get(genome_id, 0)
            initial_count = initial_counts.get(genome_id, 0)
            final_count = final_counts.get(genome_id, 0)

            initial_read_error = initial_count - sim_count
            final_read_error = final_count - sim_count

            total_abs_read_error_initial += abs(initial_read_error)
            total_abs_read_error_final += abs(final_read_error)

            f.write(f"Genom: {genome_id}\n")
            f.write("-" * 80 + "\n")

            f.write("BROJ OCITANJA\n")
            f.write(f"  simulator / FASTQ        : {sim_count:.2f}\n")
            f.write(f"  inicijalna dodjela       : {initial_count:.2f}\n")
            f.write(f"  finalna preraspodjela    : {final_count:.2f}\n")
            f.write(f"  initial - simulator      : {initial_read_error:.2f}\n")
            f.write(f"  final - simulator        : {final_read_error:.2f}\n")
            f.write(f"  final - initial          : {final_count - initial_count:.2f}\n\n")

            sim = sim_stats.get(genome_id)
            initial = initial_stats.get(genome_id)
            final = final_stats.get(genome_id)

            f.write("COVERAGE STATISTIKA\n")

            if sim is not None:
                f.write("  Simulator:\n")
                f.write(f"    mean                   : {sim['mean']:.4f}\n")
                f.write(f"    std                    : {sim['std']:.4f}\n")
                f.write(f"    cv                     : {sim['cv']:.4f}\n")
                f.write(f"    zero fraction          : {sim['zero_fraction']:.4f}\n")
                f.write(f"    score cv+zero_fraction : {sim['score_cv_zero']:.4f}\n")
            else:
                f.write("  Simulator: nema ocitanja / nema statistike\n")

            if initial is not None:
                f.write("  Inicijalna dodjela:\n")
                f.write(f"    mean                   : {initial['mean']:.4f}\n")
                f.write(f"    std                    : {initial['std']:.4f}\n")
                f.write(f"    cv                     : {initial['cv']:.4f}\n")
                f.write(f"    zero fraction          : {initial['zero_fraction']:.4f}\n")
                f.write(f"    score cv+zero_fraction : {initial['score_cv_zero']:.4f}\n")
            else:
                f.write("  Inicijalna dodjela: nema statistike\n")

            if final is not None:
                f.write("  Finalna preraspodjela:\n")
                f.write(f"    mean                   : {final['mean']:.4f}\n")
                f.write(f"    std                    : {final['std']:.4f}\n")
                f.write(f"    cv                     : {final['cv']:.4f}\n")
                f.write(f"    zero fraction          : {final['zero_fraction']:.4f}\n")
                f.write(f"    score cv+zero_fraction : {final['score_cv_zero']:.4f}\n")
            else:
                f.write("  Finalna preraspodjela: nema statistike\n")

            f.write("\n")

            if sim is not None and initial is not None and final is not None:
                genomes_with_sim_stats += 1

                mean_error_initial = initial["mean"] - sim["mean"]
                mean_error_final = final["mean"] - sim["mean"]

                score_error_initial = initial["score_cv_zero"] - sim["score_cv_zero"]
                score_error_final = final["score_cv_zero"] - sim["score_cv_zero"]

                total_abs_mean_error_initial += abs(mean_error_initial)
                total_abs_mean_error_final += abs(mean_error_final)

                total_abs_score_error_initial += abs(score_error_initial)
                total_abs_score_error_final += abs(score_error_final)

                f.write("USPOREDBA PREMA SIMULATORU\n")
                f.write(f"  mean greska initial      : {mean_error_initial:.4f}\n")
                f.write(f"  mean greska final        : {mean_error_final:.4f}\n")
                f.write(
                    "  mean promjena greske     : "
                    f"{abs(mean_error_final) - abs(mean_error_initial):.4f}\n"
                )

                f.write(f"  score greska initial     : {score_error_initial:.4f}\n")
                f.write(f"  score greska final       : {score_error_final:.4f}\n")
                f.write(
                    "  score promjena greske    : "
                    f"{abs(score_error_final) - abs(score_error_initial):.4f}\n\n"
                )

            if initial is not None and final is not None:
                f.write("PROMJENA FINAL - INITIAL\n")
                f.write(f"  mean promjena            : {final['mean'] - initial['mean']:.4f}\n")
                f.write(f"  std promjena             : {final['std'] - initial['std']:.4f}\n")
                f.write(f"  cv promjena              : {final['cv'] - initial['cv']:.4f}\n")
                f.write(
                    "  zero fraction promjena   : "
                    f"{final['zero_fraction'] - initial['zero_fraction']:.4f}\n"
                )
                f.write(
                    "  score promjena           : "
                    f"{final['score_cv_zero'] - initial['score_cv_zero']:.4f}\n\n"
                )

            f.write("\n")

        f.write("\n")
        f.write("=" * 80 + "\n")
        f.write("UKUPNI SAZETAK\n")
        f.write("=" * 80 + "\n\n")

        f.write("GRESKA BROJA OCITANJA PREMA SIMULATORU\n")
        f.write(f"  suma |initial - simulator| : {total_abs_read_error_initial:.2f}\n")
        f.write(f"  suma |final - simulator|   : {total_abs_read_error_final:.2f}\n")
        f.write(
            "  promjena                   : "
            f"{total_abs_read_error_final - total_abs_read_error_initial:.2f}\n\n"
        )

        if total_abs_read_error_final < total_abs_read_error_initial:
            f.write("  Zakljucak: finalna preraspodjela je poboljsala brojcanu raspodjelu ocitanja.\n\n")
        elif total_abs_read_error_final > total_abs_read_error_initial:
            f.write("  Zakljucak: finalna preraspodjela je pogorsala brojcanu raspodjelu ocitanja.\n\n")
        else:
            f.write("  Zakljucak: finalna preraspodjela nije promijenila ukupnu brojcanu gresku.\n\n")

        if genomes_with_sim_stats > 0:
            f.write("GRESKA MEAN COVERAGEA PREMA SIMULATORU\n")
            f.write(
                "  suma |initial mean - simulator mean| : "
                f"{total_abs_mean_error_initial:.4f}\n"
            )
            f.write(
                "  suma |final mean - simulator mean|   : "
                f"{total_abs_mean_error_final:.4f}\n"
            )
            f.write(
                "  promjena                            : "
                f"{total_abs_mean_error_final - total_abs_mean_error_initial:.4f}\n\n"
            )

            f.write("GRESKA SCOREA PREMA SIMULATORU\n")
            f.write(
                "  suma |initial score - simulator score| : "
                f"{total_abs_score_error_initial:.4f}\n"
            )
            f.write(
                "  suma |final score - simulator score|   : "
                f"{total_abs_score_error_final:.4f}\n"
            )
            f.write(
                "  promjena                              : "
                f"{total_abs_score_error_final - total_abs_score_error_initial:.4f}\n\n"
            )

    return output_path