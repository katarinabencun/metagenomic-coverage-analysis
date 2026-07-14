import os
import argparse

from statistika import (
    count_reads_simulator,
    coverage_stats,
    write_statistics,
    write_redistribution_comparison,
    write_coverage_distance_to_simulator,
    write_false_genome_coverage_stats,
    write_genome_truth_summary,
    write_full_assignment_evaluation,
    write_cleanup_truth_evaluation,
    organize_cleanup_comparison_images,
)

from parse_i_coverage import (
    parse_fastq,
    parse_paf_grouped,
    parse_paf_grouped_with_cigar,
    get_genome_lengths,
    compute_coverage_from_intervals,
    compute_coverage_from_assignments,
    compute_coverage_from_assignments_cigar,
)

from vizualizacija import (
    plot_coverage_profile_kbp,
    plot_coverage_comparison_stacked_kbp,
)

from inicijalna_preraspodjela import (
    initial_balanced_assignment,
    assignment_counts,
    write_initial_assignment_summary,
)

from algoritam_preraspodjele import (
    redistribute_algorithm,
    write_redistribution_summary,
    cleanup_simple,
    write_cleanup_simple_summary,
)

def compute_assignment_coverage(assignments, use_cigar, bucket_size):
    """
    Računa coverage iz dodjela, ovisno o tome koristi li se CIGAR ili ne.
    Kasnija logika pipelinea ostaje ista.
    """
    if use_cigar:
        return compute_coverage_from_assignments_cigar(
            assignments,
            bucket_size=bucket_size
        )

    return compute_coverage_from_assignments(
        assignments,
        bucket_size=bucket_size
    )

def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Analiza pokrivenosti i preraspodjela metagenomskih očitanja."
    )

    parser.add_argument(
        "--fastq",
        required=True,
        help="Putanja do FASTQ datoteke simuliranih očitanja."
    )

    parser.add_argument(
        "--fasta",
        required=True,
        help="Putanja do FASTA datoteke referentnih genoma."
    )

    parser.add_argument(
        "--paf",
        required=True,
        help="Putanja do PAF datoteke s mapiranjima."
    )

    parser.add_argument(
        "--bucket-size",
        type=int,
        default=5000,
        help="Veličina pretinca u baznim parovima. Zadano: 5000."
    )

    parser.add_argument(
        "--experiment-name",
        default="experiment",
        help="Naziv eksperimenta i izlazne mape."
    )

    parser.add_argument(
        "--without-cigar",
        action="store_true",
        help="PAF datoteka ne sadrži CIGAR zapis."
    )

    parser.add_argument(
        "--skip-redistribution",
        action="store_true",
        help="Preskače glavnu preraspodjelu."
    )

    parser.add_argument(
        "--skip-cleanup",
        action="store_true",
        help="Preskače cleanup fazu."
    )

    parser.add_argument(
        "--without-drainage",
        action="store_true",
        help="Isključuje drainage signal u cleanup fazi."
    )

    return parser.parse_args()

def main():
    args = parse_arguments()

    fastq_path = args.fastq
    fasta_path = args.fasta
    paf_path = args.paf
    bucket_size = args.bucket_size

    # Provjera postoje li zadane ulazne datoteke
    for input_path in (fastq_path, fasta_path, paf_path):
        if not os.path.isfile(input_path):
            raise FileNotFoundError(
                f"Ulazna datoteka ne postoji: {input_path}"
            )

    use_cigar = not args.without_cigar
    run_main_redistribution = not args.skip_redistribution
    run_cleanup = not args.skip_cleanup
    use_cleanup_drainage = not args.without_drainage

    results_dir = f"results/bucket{bucket_size}_{args.experiment_name}"

    statistika_dir = f"{results_dir}/statistika"
    statistika_osnovno_dir = f"{statistika_dir}/osnovno"
    statistika_dodatno_dir = f"{statistika_dir}/dodatno"
    statistika_summaries_dir = f"{statistika_dir}/summaries"

    os.makedirs(results_dir, exist_ok=True)
    os.makedirs(statistika_dir, exist_ok=True)
    os.makedirs(statistika_osnovno_dir, exist_ok=True)
    os.makedirs(statistika_dodatno_dir, exist_ok=True)
    os.makedirs(statistika_summaries_dir, exist_ok=True)

    genome_lengths = get_genome_lengths(fasta_path)

    # 1. Simulator / FASTQ
    reads_sim = list(parse_fastq(fastq_path))

    coverage_sim = compute_coverage_from_intervals(
        reads_sim,
        bucket_size=bucket_size
    )

    sim_counts = count_reads_simulator(reads_sim)

    sim_stats = coverage_stats(
        coverage_sim,
        genome_lengths,
        bucket_size=bucket_size
    )

    plot_coverage_profile_kbp(
        coverage_sim,
        genome_lengths,
        bucket_size=bucket_size,
        output_prefix=f"sim",
        output_dir=f"{results_dir}/profile_kbp"
    )

    # 2. Minimap2 / PAF best-MAPQ
    if use_cigar:
        grouped_reads = parse_paf_grouped_with_cigar(paf_path)
    else:
        grouped_reads = parse_paf_grouped(paf_path)

    seed = 42

    initial_assignments = initial_balanced_assignment(
        grouped_reads,
        seed=seed
    )

    initial_counts = assignment_counts(initial_assignments)

    initial_assignment_coverage = compute_assignment_coverage(
        initial_assignments,
        use_cigar=use_cigar,
        bucket_size=bucket_size
    )

    initial_stats = coverage_stats(
        initial_assignment_coverage,
        genome_lengths,
        bucket_size=bucket_size
    )

    plot_coverage_profile_kbp(
        initial_assignment_coverage,
        genome_lengths,
        bucket_size=bucket_size,
        output_prefix="initial_assignment",
        output_dir=f"{results_dir}/profile_kbp"
    )

    if run_main_redistribution:
        final_assignments, redistribution_summary = redistribute_algorithm(
            initial_assignments=initial_assignments,
            grouped_reads=grouped_reads,
            genome_lengths=genome_lengths,
            bucket_size=bucket_size,

            max_iterations=10,
            support_weight=1.0,
            suspicion_weight=8.0,

            unique_low_threshold=0.40,
            zero_start=0.04,
            zero_full=0.08,
            long_zero_start=0.015,
            long_zero_full=0.4,
        )
    else:
        final_assignments = dict(initial_assignments)
        redistribution_summary = {
            "algorithm": "redistribution_disabled",
            "iterations": [],
            "total_moves": 0,
        }

    if run_cleanup:
        final_assignments, cleanup_summary = cleanup_simple(
            assignments=final_assignments,
            grouped_reads=grouped_reads,
            genome_lengths=genome_lengths,
            bucket_size=bucket_size,

            suspicion_threshold=0.05,

            # osnovni signal
            zero_start=0.40,
            zero_full=0.70,
            long_zero_start=0.10,
            long_zero_full=0.25,

            # otočni signal
            top_fraction=0.05,
            top_coverage_start=0.20,
            top_coverage_full=0.50,

            spike_ratio_start=6.0,
            spike_ratio_full=12.0,

            # signal pražnjenja
            baseline_assignments=initial_assignments,
            use_drainage_signal=use_cleanup_drainage,
            drainage_rel_rmse_drop_start=0.10,
            drainage_rel_rmse_drop_full=0.20, #bilo je na 0.25 sasvim ok,

            # premještanje
            cleanup_mapq_delta=60,
            mapq_penalty_weight=0.05,
            non_primary_weight=0.05,
            alignment_quality_weight=0.5,
        )
    else:
        cleanup_summary = {
            "algorithm": "cleanup_disabled",
            "suspicious_genomes": [],
            "suspicion_details": {},
            "total_cleanup_moves": 0,
            "skipped_no_alternative": 0,
            "skipped_only_suspicious_alternatives": 0,
            "moves_from": {},
            "moves_to": {},
        }

    final_counts = assignment_counts(final_assignments)

    final_coverage = compute_assignment_coverage(
        final_assignments,
        use_cigar=use_cigar,
        bucket_size=bucket_size
    )

    final_stats = coverage_stats(
        final_coverage,
        genome_lengths,
        bucket_size=bucket_size
    )

    plot_coverage_profile_kbp(
        final_coverage,
        genome_lengths,
        bucket_size=bucket_size,
        output_prefix="final_redistribution",
        output_dir=f"{results_dir}/profile_kbp"
    )

    plot_coverage_comparison_stacked_kbp(
        coverage_sim=coverage_sim,
        coverage_initial=initial_assignment_coverage,
        coverage_final=final_coverage,
        genome_lengths=genome_lengths,
        bucket_size=bucket_size,
        output_dir=f"{results_dir}/usporedba"
    )

    if run_cleanup:
        organize_cleanup_comparison_images(
            sim_counts=sim_counts,
            genome_lengths=genome_lengths,
            cleanup_summary=cleanup_summary,
            comparison_dir=f"{results_dir}/usporedba"
        )

    # 3. Statistika
    write_genome_truth_summary(
        sim_counts=sim_counts,
        genome_lengths=genome_lengths,
        output_dir=statistika_summaries_dir,
        filename="genome_truth_summary.txt"
    )

    write_initial_assignment_summary(
        grouped_reads,
        initial_assignments,
        output_dir=statistika_summaries_dir,
        filename="summary_init_assign.txt"
    )

    write_redistribution_summary(
        redistribution_summary,
        output_dir=statistika_summaries_dir,
        filename="redistribution_summary.txt"
    )

    write_cleanup_simple_summary(
        cleanup_summary,
        output_dir=statistika_summaries_dir,
        filename="cleanup_simple_summary.txt"
    )

    write_statistics(
        sim_counts=sim_counts,
        paf_counts=initial_counts,
        sim_stats=sim_stats,
        paf_stats=initial_stats,
        bucket_size=bucket_size,
        output_dir=statistika_osnovno_dir,
        filename="statistics.txt"
    )

    write_redistribution_comparison(
        sim_counts=sim_counts,
        initial_counts=initial_counts,
        final_counts=final_counts,
        sim_stats=sim_stats,
        initial_stats=initial_stats,
        final_stats=final_stats,
        bucket_size=bucket_size,
        output_dir=statistika_dodatno_dir,
        filename="redistribution_comparison.txt"
    )

    write_coverage_distance_to_simulator(
        coverage_sim=coverage_sim,
        coverage_initial=initial_assignment_coverage,
        coverage_final=final_coverage,
        genome_lengths=genome_lengths,
        bucket_size=bucket_size,
        output_dir=statistika_dodatno_dir,
        filename="coverage_distance_to_simulator.txt"
    )

    write_false_genome_coverage_stats(
        coverage_sim=coverage_sim,
        coverage_initial=initial_assignment_coverage,
        coverage_final=final_coverage,
        genome_lengths=genome_lengths,
        bucket_size=bucket_size,
        output_dir=statistika_dodatno_dir,
        filename="false_genome_coverage_stats.txt"
    )

    if run_cleanup:
        write_cleanup_truth_evaluation(
            sim_counts=sim_counts,
            genome_lengths=genome_lengths,
            cleanup_summary=cleanup_summary,
            output_dir=statistika_dodatno_dir,
            filename="cleanup_truth_evaluation.txt"
        )

    write_full_assignment_evaluation(
        fastq_path=fastq_path,
        initial_assignments=initial_assignments,
        final_assignments=final_assignments,
        genome_lengths=genome_lengths,
        output_dir=f"{statistika_dodatno_dir}/assignment_evaluation"
    )

    print(f"Rezultati spremljeni u: {results_dir}")

if __name__ == "__main__":
    main()