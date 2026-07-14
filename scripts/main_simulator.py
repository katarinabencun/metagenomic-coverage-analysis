import os

from statistika import (
    count_reads_simulator,
    coverage_stats,
    write_statistics,
    write_redistribution_comparison,
    write_coverage_distance_to_simulator,
    write_false_genome_coverage_stats,
    write_genome_truth_summary,
    write_full_assignment_evaluation,
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


BUCKET_SIZE = 5000 # radim s 1000 i 5000
DATASET = "2b4s"   # radim s 2b4s, 5b15s i 5b15s_add
OZNAKA = "c1"      # radim s c1 i c4 za 2b4s, a za 5b15s i 5b15s_add samo s c4
USE_CIGAR = True   # True = mapping_c1_cigar.paf, False = mapping_c1.paf
EXPERIMENT_NAME = "finalno"  # "samo_preraspodjela", "preraspodjela_cleanup", "sve" = preraspodjela_cleanup_drainage

RUN_MAIN_REDISTRIBUTION = True
RUN_CLEANUP = True
USE_CLEANUP_DRAINAGE = True

RUN_NAME = f"bucket{BUCKET_SIZE}_{DATASET}_{OZNAKA}_{EXPERIMENT_NAME}"
RESULTS_DIR = f"results/{RUN_NAME}"

DATASET_DIR = f"data/{DATASET}/{DATASET}"
SAMPLE_DIR = f"{DATASET_DIR}/{OZNAKA}"

if DATASET == "2b4s":
    FASTQ_PATH = f"{SAMPLE_DIR}/2bacteria4strains_{OZNAKA}.fastq"

    if OZNAKA == "c1":
        FASTA_PATH = f"{SAMPLE_DIR}/2bacteria4strains_c1.reduced_db.fa"
    elif OZNAKA == "c4":
        # nema posebne c4 reducirane baze, koristi istu reduciranu bazu iz c1.
        FASTA_PATH = f"{DATASET_DIR}/c1/2bacteria4strains_c1.reduced_db.fa"
    else:
        raise ValueError(f"Nepoznata OZNAKA za DATASET={DATASET}: {OZNAKA}")

elif DATASET == "5b15s":
    FASTQ_PATH = f"{SAMPLE_DIR}/5bacteria15strains_{OZNAKA}.fastq.gz"
    FASTA_PATH = f"{SAMPLE_DIR}/5bacteria15strains_{OZNAKA}.reduced_db.fasta"

elif DATASET == "5b15s_add":
    FASTQ_PATH = f"{SAMPLE_DIR}/5bacteria15strains_additional_{OZNAKA}.fastq.gz"
    FASTA_PATH = f"{SAMPLE_DIR}/5bacteria15strains_additional_{OZNAKA}.reduced_db.fa"

else:
    raise ValueError(f"Nepoznat DATASET: {DATASET}")

PAF_PATH = f"minimap_output/{DATASET}/mapping_{OZNAKA}.paf"
PAF_CIGAR_PATH = f"minimap_output/{DATASET}/mapping_{OZNAKA}_cigar.paf"

ACTIVE_PAF_PATH = PAF_CIGAR_PATH if USE_CIGAR else PAF_PATH

STATISTIKA_DIR = f"{RESULTS_DIR}/statistika"
STATISTIKA_OSNOVNO_DIR = f"{STATISTIKA_DIR}/osnovno"
STATISTIKA_DODATNO_DIR = f"{STATISTIKA_DIR}/dodatno"
STATISTIKA_SUMMARYS_DIR = f"{STATISTIKA_DIR}/summaries"


def compute_assignment_coverage(assignments):
    """
    Računa coverage iz dodjela, ovisno o tome koristi li se CIGAR ili ne.
    Kasnija logika pipelinea ostaje ista.
    """
    if USE_CIGAR:
        return compute_coverage_from_assignments_cigar(
            assignments,
            bucket_size=BUCKET_SIZE
        )

    return compute_coverage_from_assignments(
        assignments,
        bucket_size=BUCKET_SIZE
    )


def main():
    os.makedirs(RESULTS_DIR, exist_ok=True)

    genome_lengths = get_genome_lengths(FASTA_PATH)

    # 1. Simulator / FASTQ
    reads_sim = list(parse_fastq(FASTQ_PATH))

    coverage_sim = compute_coverage_from_intervals(
        reads_sim,
        bucket_size=BUCKET_SIZE
    )

    sim_counts = count_reads_simulator(reads_sim)

    sim_stats = coverage_stats(
        coverage_sim,
        genome_lengths,
        bucket_size=BUCKET_SIZE
    )

    plot_coverage_profile_kbp(
        coverage_sim,
        genome_lengths,
        bucket_size=BUCKET_SIZE,
        output_prefix=f"sim",
        output_dir=f"{RESULTS_DIR}/profile_kbp"
    )

    # 2. Minimap2 / PAF best-MAPQ
    if USE_CIGAR:
        grouped_reads = parse_paf_grouped_with_cigar(ACTIVE_PAF_PATH)
    else:
        grouped_reads = parse_paf_grouped(ACTIVE_PAF_PATH)

    seed = 42

    initial_assignments = initial_balanced_assignment(
        grouped_reads,
        seed=seed
    )

    initial_counts = assignment_counts(initial_assignments)

    initial_assignment_coverage = compute_assignment_coverage(
        initial_assignments
    )

    initial_stats = coverage_stats(
        initial_assignment_coverage,
        genome_lengths,
        bucket_size=BUCKET_SIZE
    )

    plot_coverage_profile_kbp(
        initial_assignment_coverage,
        genome_lengths,
        bucket_size=BUCKET_SIZE,
        output_prefix="initial_assignment",
        output_dir=f"{RESULTS_DIR}/profile_kbp"
    )

    if RUN_MAIN_REDISTRIBUTION:
        final_assignments, redistribution_summary = redistribute_algorithm(
            initial_assignments=initial_assignments,
            grouped_reads=grouped_reads,
            genome_lengths=genome_lengths,
            bucket_size=BUCKET_SIZE,

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

    if RUN_CLEANUP:
        final_assignments, cleanup_summary = cleanup_simple(
            assignments=final_assignments,
            grouped_reads=grouped_reads,
            genome_lengths=genome_lengths,
            bucket_size=BUCKET_SIZE,

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
            use_drainage_signal=USE_CLEANUP_DRAINAGE,
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
        final_assignments
    )

    final_stats = coverage_stats(
        final_coverage,
        genome_lengths,
        bucket_size=BUCKET_SIZE
    )

    plot_coverage_profile_kbp(
        final_coverage,
        genome_lengths,
        bucket_size=BUCKET_SIZE,
        output_prefix="final_redistribution",
        output_dir=f"{RESULTS_DIR}/profile_kbp"
    )

    plot_coverage_comparison_stacked_kbp(
        coverage_sim=coverage_sim,
        coverage_initial=initial_assignment_coverage,
        coverage_final=final_coverage,
        genome_lengths=genome_lengths,
        bucket_size=BUCKET_SIZE,
        output_dir=f"{RESULTS_DIR}/usporedba"
    )

    # 3. Statistika
    write_genome_truth_summary(
        sim_counts=sim_counts,
        genome_lengths=genome_lengths,
        output_dir=f"{STATISTIKA_SUMMARYS_DIR}",
        filename="genome_truth_summary.txt"
    )

    write_initial_assignment_summary(
        grouped_reads,
        initial_assignments,
        output_dir=f"{STATISTIKA_SUMMARYS_DIR}",
        filename="summary_init_assign.txt"
    )

    write_redistribution_summary(
        redistribution_summary,
        output_dir=f"{STATISTIKA_SUMMARYS_DIR}",
        filename="redistribution_summary.txt"
    )

    write_cleanup_simple_summary(
        cleanup_summary,
        output_dir=f"{STATISTIKA_SUMMARYS_DIR}",
        filename="cleanup_simple_summary.txt"
    )

    write_statistics(
        sim_counts=sim_counts,
        paf_counts=initial_counts,
        sim_stats=sim_stats,
        paf_stats=initial_stats,
        bucket_size=BUCKET_SIZE,
        output_dir=f"{STATISTIKA_OSNOVNO_DIR}",
        filename="statistics.txt"
    )

    write_redistribution_comparison(
        sim_counts=sim_counts,
        initial_counts=initial_counts,
        final_counts=final_counts,
        sim_stats=sim_stats,
        initial_stats=initial_stats,
        final_stats=final_stats,
        bucket_size=BUCKET_SIZE,
        output_dir=f"{STATISTIKA_DODATNO_DIR}",
        filename="redistribution_comparison.txt"
    )

    write_coverage_distance_to_simulator(
        coverage_sim=coverage_sim,
        coverage_initial=initial_assignment_coverage,
        coverage_final=final_coverage,
        genome_lengths=genome_lengths,
        bucket_size=BUCKET_SIZE,
        output_dir=f"{STATISTIKA_DODATNO_DIR}",
        filename="coverage_distance_to_simulator.txt"
    )

    write_false_genome_coverage_stats(
        coverage_sim=coverage_sim,
        coverage_initial=initial_assignment_coverage,
        coverage_final=final_coverage,
        genome_lengths=genome_lengths,
        bucket_size=BUCKET_SIZE,
        output_dir=f"{STATISTIKA_DODATNO_DIR}",
        filename="false_genome_coverage_stats.txt"
    )

    write_full_assignment_evaluation(
        fastq_path=FASTQ_PATH,
        initial_assignments=initial_assignments,
        final_assignments=final_assignments,
        genome_lengths=genome_lengths,
        output_dir=f"{STATISTIKA_DODATNO_DIR}/assignment_evaluation"
    )

    print(f"Rezultati spremljeni u: {RESULTS_DIR}")

if __name__ == "__main__":
    main()