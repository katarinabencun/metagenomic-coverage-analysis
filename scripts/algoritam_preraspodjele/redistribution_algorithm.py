from copy import deepcopy
from collections import defaultdict
import math
import os

from parse_i_coverage.parse_paf import parse_cigar
from inicijalna_preraspodjela.initial_assignment import (
    get_best_alignments,
    get_candidate_genomes
)


# =============================================================================
# COVERAGE POMOCNE FUNKCIJE
# =============================================================================

def alignment_bucket_contributions(aln, bucket_size=1000):
    """
    Vraca doprinos jednog poravnanja po bucketima.

    Povratna vrijednost:
    {
        bucket_index: coverage_doprinos
    }

    Doprinos je vec normaliziran s bucket_size.
    Koristi se CIGAR, pa se u coverage broje samo M, = i X operacije.
    """
    contributions = defaultdict(float)

    ref_pos = aln["start"]
    cigar = aln["cigar"]

    for length, op in parse_cigar(cigar):
        if op in ("M", "=", "X"):
            segment_start = ref_pos
            segment_end = ref_pos + length

            start_bucket = segment_start // bucket_size
            end_bucket = (segment_end - 1) // bucket_size

            for bucket_index in range(start_bucket, end_bucket + 1):
                bucket_start = bucket_index * bucket_size
                bucket_end = bucket_start + bucket_size

                overlap_start = max(segment_start, bucket_start)
                overlap_end = min(segment_end, bucket_end)

                overlap = overlap_end - overlap_start

                if overlap > 0:
                    contributions[bucket_index] += overlap / bucket_size

            ref_pos += length

        elif op in ("D", "N"):
            # Pomice referentnu poziciju, ali ne dodaje coverage ocitanjem.
            ref_pos += length

        elif op in ("I", "S", "H", "P"):
            # Ne doprinosi coverageu referentnog genoma.
            continue

    return dict(contributions)


def add_alignment_to_coverage(coverage, aln, bucket_size=1000, weight=1.0):
    """
    Dodaje ili uklanja jedno poravnanje iz trenutnog coveragea.

    weight =  1.0 -> dodaj ocitanje
    weight = -1.0 -> makni ocitanje
    """
    genome_id = aln["genome_id"]

    if genome_id not in coverage:
        coverage[genome_id] = {}

    contributions = alignment_bucket_contributions(
        aln,
        bucket_size=bucket_size
    )

    for bucket_index, value in contributions.items():
        new_value = coverage[genome_id].get(bucket_index, 0.0) + weight * value

        if abs(new_value) < 1e-12:
            new_value = 0.0

        coverage[genome_id][bucket_index] = new_value


def build_coverage_from_assignments(assignments, bucket_size=1000):
    """
    Gradi coverage iz trenutne diskretne dodjele ocitanja.
    """
    coverage = {}

    for read_id, aln in assignments.items():
        add_alignment_to_coverage(
            coverage,
            aln,
            bucket_size=bucket_size,
            weight=1.0
        )

    return coverage


def genome_num_buckets(genome_id, genome_lengths, bucket_size=1000):
    """
    Vraca broj bucket-a za zadani genom.
    """
    return (genome_lengths[genome_id] + bucket_size - 1) // bucket_size


def genome_mean_coverage(coverage, genome_id, genome_lengths, bucket_size=1000):
    """
    Racuna srednji coverage genoma, ukljucujuci bucket-e s nulom.
    """
    if genome_id not in genome_lengths:
        return 0.0

    num_buckets = genome_num_buckets(
        genome_id,
        genome_lengths,
        bucket_size
    )

    if num_buckets == 0:
        return 0.0

    total_coverage = sum(coverage.get(genome_id, {}).values())

    return total_coverage / num_buckets


def genome_coverage_values(coverage, genome_id, genome_lengths, bucket_size=1000):
    """
    Vraca puni niz coverage vrijednosti za genom,
    ukljucujuci i bucket-e koji nemaju coverage.
    """
    if genome_id not in genome_lengths:
        return []

    num_buckets = genome_num_buckets(
        genome_id,
        genome_lengths,
        bucket_size
    )

    genome_cov = coverage.get(genome_id, {})

    return [
        genome_cov.get(bucket_index, 0.0)
        for bucket_index in range(num_buckets)
    ]


def genome_zero_fraction_from_coverage(
    coverage,
    genome_id,
    genome_lengths,
    bucket_size=1000
):
    """
    Racuna udio bucket-a s coverageom 0.
    Veci zero_fraction znaci fragmentiraniju / nepotpuniju pokrivenost.
    """
    values = genome_coverage_values(
        coverage,
        genome_id,
        genome_lengths,
        bucket_size
    )

    if not values:
        return 1.0

    zero_count = sum(1 for value in values if value <= 1e-12)

    return zero_count / len(values)


def genome_max_zero_run_fraction(
    coverage,
    genome_id,
    genome_lengths,
    bucket_size=1000,
    zero_epsilon=1e-12
):
    """
    Racuna najduzi uzastopni niz bucket-a s nulom,
    normaliziran ukupnim brojem bucket-a genoma.

    Primjer:
        0.05 znaci da je najduza rupa 5% genoma.
        0.40 znaci da postoji velika kontinuirana rupa.
    """
    values = genome_coverage_values(
        coverage,
        genome_id,
        genome_lengths,
        bucket_size
    )

    if not values:
        return 1.0

    max_run = 0
    current_run = 0

    for value in values:
        if value <= zero_epsilon:
            current_run += 1
            max_run = max(max_run, current_run)
        else:
            current_run = 0

    return max_run / len(values)

# =============================================================================
# LOKALNI SSE KRITERIJ
# =============================================================================

def local_sse_before_after_remove(
    coverage,
    genome_id,
    contributions,
    mean
):
    """
    Racuna lokalnu promjenu sume kvadratnih odstupanja ako se ocitanje makne.

    Negativan rezultat znaci da se lokalni coverage poboljsava.
    Pozitivan rezultat znaci da se lokalni coverage pogorsava.
    """
    before = 0.0
    after = 0.0

    genome_cov = coverage.get(genome_id, {})

    for bucket_index, value in contributions.items():
        old_cov = genome_cov.get(bucket_index, 0.0)
        new_cov = old_cov - value

        before += (old_cov - mean) ** 2
        after += (new_cov - mean) ** 2

    return after - before


def local_sse_before_after_add(
    coverage,
    genome_id,
    contributions,
    mean
):
    """
    Racuna lokalnu promjenu sume kvadratnih odstupanja ako se ocitanje doda.

    Negativan rezultat znaci da se lokalni coverage poboljsava.
    Pozitivan rezultat znaci da se lokalni coverage pogorsava.
    """
    before = 0.0
    after = 0.0

    genome_cov = coverage.get(genome_id, {})

    for bucket_index, value in contributions.items():
        old_cov = genome_cov.get(bucket_index, 0.0)
        new_cov = old_cov + value

        before += (old_cov - mean) ** 2
        after += (new_cov - mean) ** 2

    return after - before


# =============================================================================
# UNIQUE SUPPORT I SUMNJIVOST GENOMA
# =============================================================================

def compute_unique_support(grouped_reads):
    """
    Broji koliko svaki genom ima jednoznacno najboljih ocitanja.

    Jednoznacno najbolje ocitanje je ono za koje, nakon izdvajanja najboljeg MAPQ-a,
    postoji samo jedan kandidatski genom.

    Taj broj koristimo kao signal da je genom vjerojatno stvarno prisutan.
    """
    support = defaultdict(int)

    for read_id, alignments in grouped_reads.items():
        best_alignments = get_best_alignments(alignments)

        if not best_alignments:
            continue

        candidate_genomes = get_candidate_genomes(best_alignments)

        if len(candidate_genomes) == 1:
            genome_id = candidate_genomes[0]
            support[genome_id] += 1

    return dict(support)


def normalize_support(unique_support):
    """
    Normalizira unique support u interval [0, 1].

    Koristi se log skala da genomi s jako puno jednoznacnih ocitanja
    ne dominiraju previse nad ostalima.
    """
    if not unique_support:
        return {}

    max_support = max(unique_support.values())

    if max_support == 0:
        return {
            genome_id: 0.0
            for genome_id in unique_support
        }

    normalized = {}

    for genome_id, count in unique_support.items():
        normalized[genome_id] = math.log1p(count) / math.log1p(max_support)

    return normalized


def clamp01(value):
    return max(0.0, min(1.0, value))


def compute_suspicion_score(
    coverage,
    genome_id,
    genome_lengths,
    support,
    bucket_size=1000,
    unique_low_threshold=0.40,
    zero_start=0.10,
    zero_full=1.0,
    long_zero_start=0.01,
    long_zero_full=1.0
):
    """
    Racuna meku sumnjivost genoma.

    Sumnjivost se ne temelji samo na malom unique supportu.
    Genom postaje sumnjiv tek kada istovremeno ima:
      1. nizak normalizirani unique support
      2. los oblik coveragea, tj. puno nultih bucket-a,
         dugu kontinuiranu rupu

    Formula:
        suspicion_score = unique_part * bad_shape_part

    Vrijednost je u intervalu [0, 1].
    Veca vrijednost znaci da je genom losiji kandidat za primanje ocitanja.
    """
    normalized_support = support.get(genome_id, 0.0)

    unique_part = clamp01(
        (unique_low_threshold - normalized_support) / unique_low_threshold
    )

    zero_fraction = genome_zero_fraction_from_coverage(
        coverage,
        genome_id,
        genome_lengths,
        bucket_size
    )

    zero_part = clamp01(
        (zero_fraction - zero_start) / (zero_full - zero_start)
    )

    long_zero_fraction = genome_max_zero_run_fraction(
        coverage,
        genome_id,
        genome_lengths,
        bucket_size
    )

    long_zero_part = clamp01(
        (long_zero_fraction - long_zero_start)
        / (long_zero_full - long_zero_start)
    )
    
    bad_shape_part = max(zero_part, long_zero_part)

    return unique_part * bad_shape_part


def find_alignment_for_genome(best_alignments, genome_id):
    """
    Iz liste najboljih poravnanja bira poravnanje za zadani genom.

    Ako ih ima vise za isti genom, bira stabilno prema koordinatama.
    """
    candidates = [
        aln for aln in best_alignments
        if aln["genome_id"] == genome_id
    ]

    if not candidates:
        return None

    candidates = sorted(
        candidates,
        key=lambda aln: (
            0 if aln.get("is_primary") else 1,
            aln["start"],
            aln["end"],
            aln["mapq"],
            aln["genome_id"]
        )
    )

    return candidates[0]


# =============================================================================
# GLAVNI ALGORITAM PRERASPODJELE
# =============================================================================

def redistribute_algorithm(
    initial_assignments,
    grouped_reads,
    genome_lengths,
    bucket_size=1000,
    max_iterations=10,
    min_improvement=1e-9,
    support_weight=1.0,
    suspicion_weight=0.25,
    unique_low_threshold=0.40,
    zero_start=0.20,
    zero_full=0.60,
    long_zero_start=0.05,
    long_zero_full=0.25
):
    """
    Glavni algoritam preraspodjele ocitanja.

    Algoritam radi nad pocetnom diskretnom dodjelom.
    U svakoj iteraciji prolazi kroz ocitanja koja imaju vise jednako najboljih
    kandidatnih genoma i provjerava bi li premjestanje na drugi kandidat
    poboljsalo kriterij.

    Kriterij:

        total_delta =
            local_delta
            - support_weight * support_delta
            + suspicion_penalty

    gdje je:

        local_delta:
            promjena lokalnog SSE coveragea nakon micanja ocitanja sa starog
            genoma i dodavanja na novi genom

        support_delta:
            normalized_support(novi_genom) - normalized_support(stari_genom)

        suspicion_penalty:
            kazna za dodavanje ocitanja u genom koji ima nizak unique support
            i los oblik coveragea

    Premjestanje se prihvaca ako je total_delta negativan.
    """
    final_assignments = deepcopy(initial_assignments)

    coverage = build_coverage_from_assignments(
        final_assignments,
        bucket_size=bucket_size
    )

    unique_support = compute_unique_support(grouped_reads)
    support = normalize_support(unique_support)

    summary = {
        "algorithm": "local_sse_with_support_and_suspicion",
        "iterations": [],
        "total_moves": 0,

        "support_weight": support_weight,
        "suspicion_weight": suspicion_weight,

        "unique_low_threshold": unique_low_threshold,
        "zero_start": zero_start,
        "zero_full": zero_full,
        "long_zero_start": long_zero_start,
        "long_zero_full": long_zero_full,

        "unique_support": unique_support,
        "suspicion_scores_initial": {},
    }

    for genome_id in genome_lengths:
        summary["suspicion_scores_initial"][genome_id] = compute_suspicion_score(
            coverage=coverage,
            genome_id=genome_id,
            genome_lengths=genome_lengths,
            support=support,
            bucket_size=bucket_size,
            unique_low_threshold=unique_low_threshold,
            zero_start=zero_start,
            zero_full=zero_full,
            long_zero_start=long_zero_start,
            long_zero_full=long_zero_full,
        )

    for iteration in range(max_iterations):
        moves_in_iteration = 0
        total_delta_in_iteration = 0.0
        total_suspicion_penalty_in_iteration = 0.0

        moves_from = defaultdict(int)
        moves_to = defaultdict(int)
        moves_to_suspicious = defaultdict(int)

        for read_id in sorted(grouped_reads.keys()):
            if read_id not in final_assignments:
                continue

            alignments = grouped_reads[read_id]
            best_alignments = get_best_alignments(alignments)

            if not best_alignments:
                continue

            candidate_genomes = get_candidate_genomes(best_alignments)

            if len(candidate_genomes) <= 1:
                continue

            current_alignment = final_assignments[read_id]
            current_genome = current_alignment["genome_id"]

            if current_genome not in genome_lengths:
                continue

            current_contrib = alignment_bucket_contributions(
                current_alignment,
                bucket_size=bucket_size
            )

            current_mean = genome_mean_coverage(
                coverage,
                current_genome,
                genome_lengths,
                bucket_size=bucket_size
            )

            remove_delta = local_sse_before_after_remove(
                coverage,
                current_genome,
                current_contrib,
                current_mean
            )

            current_support = support.get(current_genome, 0.0)

            best_alignment = current_alignment
            best_total_delta = 0.0
            best_suspicion_penalty = 0.0
            best_suspicion_score = 0.0

            for candidate_genome in candidate_genomes:
                if candidate_genome == current_genome:
                    continue

                if candidate_genome not in genome_lengths:
                    continue

                candidate_alignment = find_alignment_for_genome(
                    best_alignments,
                    candidate_genome
                )

                if candidate_alignment is None:
                    continue

                candidate_contrib = alignment_bucket_contributions(
                    candidate_alignment,
                    bucket_size=bucket_size
                )

                candidate_mean = genome_mean_coverage(
                    coverage,
                    candidate_genome,
                    genome_lengths,
                    bucket_size=bucket_size
                )

                add_delta = local_sse_before_after_add(
                    coverage,
                    candidate_genome,
                    candidate_contrib,
                    candidate_mean
                )

                local_delta = remove_delta + add_delta

                candidate_support = support.get(candidate_genome, 0.0)
                support_delta = candidate_support - current_support

                suspicion_score = compute_suspicion_score(
                    coverage=coverage,
                    genome_id=candidate_genome,
                    genome_lengths=genome_lengths,
                    support=support,
                    bucket_size=bucket_size,
                    unique_low_threshold=unique_low_threshold,
                    zero_start=zero_start,
                    zero_full=zero_full,
                    long_zero_start=long_zero_start,
                    long_zero_full=long_zero_full,
                )

                suspicion_penalty = suspicion_weight * suspicion_score

                total_delta = (
                    local_delta
                    - support_weight * support_delta
                    + suspicion_penalty
                )

                if total_delta < best_total_delta - min_improvement:
                    best_total_delta = total_delta
                    best_alignment = candidate_alignment
                    best_suspicion_penalty = suspicion_penalty
                    best_suspicion_score = suspicion_score

            if best_alignment is current_alignment:
                continue

            new_genome = best_alignment["genome_id"]

            add_alignment_to_coverage(
                coverage,
                current_alignment,
                bucket_size=bucket_size,
                weight=-1.0
            )

            add_alignment_to_coverage(
                coverage,
                best_alignment,
                bucket_size=bucket_size,
                weight=1.0
            )

            final_assignments[read_id] = best_alignment

            moves_in_iteration += 1
            total_delta_in_iteration += best_total_delta
            total_suspicion_penalty_in_iteration += best_suspicion_penalty
            summary["total_moves"] += 1

            moves_from[current_genome] += 1
            moves_to[new_genome] += 1

            if best_suspicion_score > 0:
                moves_to_suspicious[new_genome] += 1

        summary["iterations"].append({
            "iteration": iteration + 1,
            "moves": moves_in_iteration,
            "total_delta": total_delta_in_iteration,
            "total_suspicion_penalty": total_suspicion_penalty_in_iteration,
            "moves_from": dict(moves_from),
            "moves_to": dict(moves_to),
            "moves_to_suspicious": dict(moves_to_suspicious),
        })

        if moves_in_iteration == 0:
            break

    return final_assignments, summary


# =============================================================================
# ZAPIS SAZETKA
# =============================================================================

def write_redistribution_summary(
    summary,
    output_dir,
    filename="redistribution_summary.txt"
):
    """
    Sprema sazetak algoritma preraspodjele u tekstualnu datoteku.
    """
    os.makedirs(output_dir, exist_ok=True)

    output_path = os.path.join(output_dir, filename)

    with open(output_path, "w") as f:
        f.write("SAZETAK ALGORITMA PRERASPODJELE\n")
        f.write("=" * 70 + "\n\n")

        if "algorithm" in summary:
            f.write(f"Algoritam: {summary['algorithm']}\n")

        f.write(f"Ukupan broj premjestanja: {summary.get('total_moves', 0)}\n")

        if "support_weight" in summary:
            f.write(f"Support weight: {summary['support_weight']}\n")

        if "suspicion_weight" in summary:
            f.write(f"Suspicion weight: {summary['suspicion_weight']}\n")

        f.write("\n")

        if "unique_support" in summary:
            f.write("UNIQUE SUPPORT PO GENOMU\n")
            f.write("-" * 70 + "\n")

            for genome_id, count in sorted(summary["unique_support"].items()):
                f.write(f"{genome_id}: {count}\n")

            f.write("\n")

        if "suspicion_scores_initial" in summary:
            f.write("POCETNI SUSPICION SCORE PO GENOMU\n")
            f.write("-" * 70 + "\n")

            for genome_id, score in sorted(summary["suspicion_scores_initial"].items()):
                f.write(f"{genome_id}: {score:.6f}\n")

            f.write("\n")

        f.write("PARAMETRI SUMNJIVOSTI\n")
        f.write("-" * 70 + "\n")
        f.write(f"unique_low_threshold : {summary.get('unique_low_threshold')}\n")
        f.write(f"zero_start           : {summary.get('zero_start')}\n")
        f.write(f"zero_full            : {summary.get('zero_full')}\n")
        f.write(f"long_zero_start      : {summary.get('long_zero_start')}\n")
        f.write(f"long_zero_full       : {summary.get('long_zero_full')}\n")

        f.write("\n")

        f.write("ITERACIJE\n")
        f.write("-" * 70 + "\n")

        for iteration_info in summary.get("iterations", []):
            f.write(f"Iteracija {iteration_info.get('iteration')}:\n")
            f.write(f"  broj premjestanja: {iteration_info.get('moves', 0)}\n")

            if "total_delta" in iteration_info:
                f.write(f"  ukupni delta: {iteration_info['total_delta']:.6f}\n")

            if "total_suspicion_penalty" in iteration_info:
                f.write(
                    "  ukupna suspicion kazna: "
                    f"{iteration_info['total_suspicion_penalty']:.6f}\n"
                )

            moves_to_suspicious = iteration_info.get("moves_to_suspicious", {})

            if moves_to_suspicious:
                f.write("  Premjestanja U sumnjive genome:\n")

                for genome_id, count in sorted(moves_to_suspicious.items()):
                    f.write(f"    {genome_id}: {count}\n")

            moves_from = iteration_info.get("moves_from", {})

            if moves_from:
                f.write("  Premjestanja IZ genoma:\n")

                for genome_id, count in sorted(moves_from.items()):
                    f.write(f"    {genome_id}: {count}\n")

            moves_to = iteration_info.get("moves_to", {})

            if moves_to:
                f.write("  Premjestanja U genome:\n")

                for genome_id, count in sorted(moves_to.items()):
                    f.write(f"    {genome_id}: {count}\n")

            f.write("\n")

    return output_path