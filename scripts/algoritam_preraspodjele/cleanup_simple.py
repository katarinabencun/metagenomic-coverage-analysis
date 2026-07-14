from copy import deepcopy
from collections import defaultdict
import math
import os

from algoritam_preraspodjele.redistribution_algorithm import (
    build_coverage_from_assignments,
    add_alignment_to_coverage,
    alignment_bucket_contributions,
    local_sse_before_after_remove,
    local_sse_before_after_add,
    genome_mean_coverage,
    genome_coverage_values,
    genome_zero_fraction_from_coverage,
    genome_max_zero_run_fraction,
)

# =============================================================================
# POMOĆNE FUNKCIJE
# =============================================================================

def clamp01(x):
    return max(0.0, min(1.0, x))


def safe_part(value, start, full):
    """
    Linearno pretvara vrijednost u interval [0, 1].

    Ako je value manji od start, rezultat je 0.
    Ako je value veći od full, rezultat je 1.
    Između start i full raste linearno.
    """
    if abs(full - start) <= 1e-12:
        return 1.0 if value >= full else 0.0

    return clamp01((value - start) / (full - start))


def genome_top_coverage_fraction(
    coverage,
    genome_id,
    genome_lengths,
    bucket_size=1000,
    top_fraction=0.05
):
    """
    Računa koliki dio ukupne pokrivenosti nose najpokriveniji bucketi.

    Primjer:
    top_fraction = 0.05 znači da se uzima najviših 5% bucket-a
    i računa se koliki dio ukupnog coveragea dolazi iz njih.

    Velika vrijednost znači da je coverage koncentriran u malom dijelu genoma,
    što upućuje na otok pokrivenosti.
    """
    values = genome_coverage_values(
        coverage,
        genome_id,
        genome_lengths,
        bucket_size
    )

    if not values:
        return 0.0

    total_coverage = sum(values)

    if total_coverage <= 1e-12:
        return 0.0

    number_of_top_buckets = max(
        1,
        math.ceil(len(values) * top_fraction)
    )

    sorted_values = sorted(values, reverse=True)
    top_sum = sum(sorted_values[:number_of_top_buckets])

    return top_sum / total_coverage


def genome_spike_ratio_from_coverage(
    coverage,
    genome_id,
    genome_lengths,
    bucket_size=1000
):
    """
    Računa omjer najveće pokrivenosti jednog bucketa i prosječne pokrivenosti genoma.

    Velika vrijednost znači da genom ima lokalni šiljak koji je višestruko
    veći od prosjeka coverage profila.
    """
    values = genome_coverage_values(
        coverage,
        genome_id,
        genome_lengths,
        bucket_size
    )

    if not values:
        return 0.0

    mean = sum(values) / len(values)

    if mean <= 1e-12:
        return 0.0

    return max(values) / mean


def genome_rmse_from_zero(
    coverage,
    genome_id,
    genome_lengths,
    bucket_size=1000
):
    """
    RMSE coveragea od idealne nule.

    U ovoj verziji current RMSE ne ulazi u cleanup score.
    Ostaje samo informativna dijagnostika u izvještaju.
    """
    values = genome_coverage_values(
        coverage,
        genome_id,
        genome_lengths,
        bucket_size
    )

    if not values:
        return 0.0

    return math.sqrt(
        sum(v ** 2 for v in values) / len(values)
    )


def alignment_quality_penalty(candidate_alignment, best_alignment, weight=0.001):
    """
    Mala kazna za kandidata koji je lošiji od najboljeg dostupnog kandidata
    prema PAF signalima.

    Koristi se samo pri odabiru alternative za očitanje, ne pri određivanju
    sumnjivosti genoma.

    Signali:
    - AS: veće je bolje
    - s1: veće je bolje
    - cm: veće je bolje
    - NM: manje je bolje
    """
    penalty = 0.0

    candidate_as = candidate_alignment.get("AS")
    best_as = best_alignment.get("AS")

    if candidate_as is not None and best_as is not None:
        penalty += max(0, best_as - candidate_as) * weight

    candidate_s1 = candidate_alignment.get("s1")
    best_s1 = best_alignment.get("s1")

    if candidate_s1 is not None and best_s1 is not None:
        penalty += max(0, best_s1 - candidate_s1) * weight

    candidate_cm = candidate_alignment.get("cm")
    best_cm = best_alignment.get("cm")

    if candidate_cm is not None and best_cm is not None:
        penalty += max(0, best_cm - candidate_cm) * weight

    candidate_nm = candidate_alignment.get("NM")
    best_nm = best_alignment.get("NM")

    if candidate_nm is not None and best_nm is not None:
        penalty += max(0, candidate_nm - best_nm) * weight

    return penalty


# =============================================================================
# DETEKCIJA SUMNJIVIH GENOMA
# =============================================================================

def compute_cleanup_simple_details(
    coverage,
    grouped_reads,
    genome_lengths,
    bucket_size=1000,

    suspicion_threshold=0.05,

    # osnovni signal
    zero_start=0.40,
    zero_full=0.70,
    long_zero_start=0.10,
    long_zero_full=0.25,

    # island signal
    top_fraction=0.05,
    top_coverage_start=0.20,
    top_coverage_full=0.50,

    # spike signal
    spike_ratio_start=6.0,
    spike_ratio_full=10.0,

    # drainage signal
    baseline_coverage=None,
    use_drainage_signal=False,
    drainage_rel_rmse_drop_start=0.10,
    drainage_rel_rmse_drop_full=0.25,
):
    """
    Jednostavna cleanup detekcija.

    Score ima tri grane:

    1. Base score:
        loš oblik coveragea
       prema zero / long zero signalima

    2. Island score:
        coverage je koncentriran
       u malom broju najpokrivenijih bucket-a.

    3. Spike score:
        jedan ili nekoliko lokalnih
       šiljaka koji su višestruko veći od prosječnog coveragea.

    4. Drainage score:
       koristi se samo ako base/island/spike već nisu označili genom.
       Hvata genome koji su se relativno ispraznili tijekom glavne
       preraspodjele.

    U cleanup_score ne ulaze count, mean, MAPQ, AS, NM ni current RMSE.
    MAPQ/AS/NM koriste se tek pri odabiru alternative za konkretno očitanje.
    """


    details = {}

    for genome_id in genome_lengths:

        # ------------------------------------------------------------
        # 1. Rupe u coverage profilu
        # ------------------------------------------------------------
        zero_fraction = genome_zero_fraction_from_coverage(
            coverage,
            genome_id,
            genome_lengths,
            bucket_size
        )

        zero_part = safe_part(
            zero_fraction,
            zero_start,
            zero_full
        )

        long_zero_fraction = genome_max_zero_run_fraction(
            coverage,
            genome_id,
            genome_lengths,
            bucket_size
        )

        long_zero_part = safe_part(
            long_zero_fraction,
            long_zero_start,
            long_zero_full
        )

        bad_shape_part = max(
            zero_part,
            long_zero_part
        )

        base_score = (
            bad_shape_part
        )

        # ------------------------------------------------------------
        # 2. Island / koncentracija coveragea
        # ------------------------------------------------------------
        top_coverage_fraction = genome_top_coverage_fraction(
            coverage,
            genome_id,
            genome_lengths,
            bucket_size,
            top_fraction=top_fraction
        )

        top_coverage_part = safe_part(
            top_coverage_fraction,
            top_coverage_start,
            top_coverage_full
        )

        island_score = (
             top_coverage_part
        )

        # ------------------------------------------------------------
        # 3. Spike / lokalni šiljci
        # ------------------------------------------------------------
        spike_ratio = genome_spike_ratio_from_coverage(
            coverage,
            genome_id,
            genome_lengths,
            bucket_size
        )

        spike_ratio_part = safe_part(
            spike_ratio,
            spike_ratio_start,
            spike_ratio_full
        )

        spike_score = (
             spike_ratio_part
        )

        primary_score = max(
            base_score,
            island_score,
            spike_score
        )

        # ------------------------------------------------------------
        # 4. Drainage / pražnjenje
        # ------------------------------------------------------------
        drainage_score = 0.0

        baseline_rmse_from_zero = 0.0
        current_rmse_from_zero = 0.0
        rmse_drop_from_baseline = 0.0
        relative_rmse_drop_from_baseline = 0.0
        drainage_rel_drop_part = 0.0

        if use_drainage_signal and baseline_coverage is not None:
            baseline_rmse_from_zero = genome_rmse_from_zero(
                baseline_coverage,
                genome_id,
                genome_lengths,
                bucket_size
            )

            current_rmse_from_zero = genome_rmse_from_zero(
                coverage,
                genome_id,
                genome_lengths,
                bucket_size
            )

            rmse_drop_from_baseline = (
                baseline_rmse_from_zero - current_rmse_from_zero
            )

            if baseline_rmse_from_zero > 1e-12:
                relative_rmse_drop_from_baseline = (
                    rmse_drop_from_baseline / baseline_rmse_from_zero
                )
            else:
                relative_rmse_drop_from_baseline = 0.0

            # Negativan pad znači da se genom nije ispraznio.
            # To ne smije povećavati cleanup score.
            relative_rmse_drop_from_baseline = max(
                0.0,
                relative_rmse_drop_from_baseline
            )

            drainage_rel_drop_part = safe_part(
                relative_rmse_drop_from_baseline,
                drainage_rel_rmse_drop_start,
                drainage_rel_rmse_drop_full
            )

            # Drainage je sekundarni signal:
            # koristi se samo za genome koje base/island/spike nisu već označili.
            if primary_score < suspicion_threshold:
                drainage_score = drainage_rel_drop_part

        cleanup_score = max(
            primary_score,
            drainage_score
        )

        details[genome_id] = {
            "cleanup_score": cleanup_score,
            "primary_score": primary_score,

            "zero_fraction": zero_fraction,
            "zero_part": zero_part,

            "long_zero_fraction": long_zero_fraction,
            "long_zero_part": long_zero_part,

            "bad_shape_part": bad_shape_part,
            "base_score": base_score,

            "top_fraction": top_fraction,
            "top_coverage_fraction": top_coverage_fraction,
            "top_coverage_part": top_coverage_part,
            "island_score": island_score,

            "spike_ratio": spike_ratio,
            "spike_ratio_part": spike_ratio_part,
            "spike_score": spike_score,

            "baseline_rmse_from_zero": baseline_rmse_from_zero,
            "current_rmse_from_zero": current_rmse_from_zero,
            "rmse_drop_from_baseline": rmse_drop_from_baseline,
            "relative_rmse_drop_from_baseline": relative_rmse_drop_from_baseline,

            "drainage_rel_drop_part": drainage_rel_drop_part,
            "drainage_score": drainage_score,
        }

    return details


# =============================================================================
# CLEANUP PREMJEŠTANJE
# =============================================================================

def cleanup_simple(
    assignments,
    grouped_reads,
    genome_lengths,
    bucket_size=1000,

    suspicion_threshold=0.05,

    zero_start=0.40,
    zero_full=0.70,
    long_zero_start=0.10,
    long_zero_full=0.25,

    top_fraction=0.05,
    top_coverage_start=0.20,
    top_coverage_full=0.50,

    spike_ratio_start=6.0,
    spike_ratio_full=10.0,

    baseline_assignments=None,
    use_drainage_signal=True,
    drainage_rel_rmse_drop_start=0.10,
    drainage_rel_rmse_drop_full=0.20,

    cleanup_mapq_delta=60,
    mapq_penalty_weight=0.05,
    non_primary_weight=0.05,
    alignment_quality_weight=0.5,
):
    """
    Jednostavna cleanup faza.

    1. Izračuna sumnjivost genoma iz coverage profila.
    2. Označi genome čiji je cleanup_score >= suspicion_threshold.
    3. Za očitanja trenutno dodijeljena sumnjivim genomima traži alternativu.
    4. Alternativa mora imati MAPQ unutar cleanup_mapq_delta od najboljeg MAPQ-a.
    5. Očitanje se ne premješta u drugi sumnjivi genom.
    6. Ako postoji barem jedna dopuštena nesumnjiva alternativa, očitanje se premješta
       na najbolju alternativu prema lokalnom SSE učinku i kvaliteti poravnanja.

    Cleanup uvijek premješta dostupna očitanja sa sumnjivih genoma.
    Nema dodatnog praga koji bi nakon odabira kandidata blokirao premještanje.
    """

    cleaned_assignments = deepcopy(assignments)

    coverage = build_coverage_from_assignments(
        cleaned_assignments,
        bucket_size=bucket_size
    )

    if baseline_assignments is not None:
        baseline_coverage = build_coverage_from_assignments(
            baseline_assignments,
            bucket_size=bucket_size
        )
    else:
        baseline_coverage = None

    details = compute_cleanup_simple_details(
        coverage=coverage,
        grouped_reads=grouped_reads,
        genome_lengths=genome_lengths,
        bucket_size=bucket_size,

        suspicion_threshold=suspicion_threshold,

        zero_start=zero_start,
        zero_full=zero_full,
        long_zero_start=long_zero_start,
        long_zero_full=long_zero_full,

        top_fraction=top_fraction,
        top_coverage_start=top_coverage_start,
        top_coverage_full=top_coverage_full,

        spike_ratio_start=spike_ratio_start,
        spike_ratio_full=spike_ratio_full,

        baseline_coverage=baseline_coverage,
        use_drainage_signal=use_drainage_signal,
        drainage_rel_rmse_drop_start=drainage_rel_rmse_drop_start,
        drainage_rel_rmse_drop_full=drainage_rel_rmse_drop_full,
    )

    suspicious_genomes = {
        genome_id
        for genome_id, d in details.items()
        if d["cleanup_score"] >= suspicion_threshold
    }

    summary = {
        "algorithm": "cleanup_simple",
        "suspicion_threshold": suspicion_threshold,

        "zero_start": zero_start,
        "zero_full": zero_full,
        "long_zero_start": long_zero_start,
        "long_zero_full": long_zero_full,

        "top_fraction": top_fraction,
        "top_coverage_start": top_coverage_start,
        "top_coverage_full": top_coverage_full,

        "spike_ratio_start": spike_ratio_start,
        "spike_ratio_full": spike_ratio_full,

        "use_drainage_signal": use_drainage_signal,
        "drainage_rel_rmse_drop_start": drainage_rel_rmse_drop_start,
        "drainage_rel_rmse_drop_full": drainage_rel_rmse_drop_full,

        "cleanup_mapq_delta": cleanup_mapq_delta,
        "mapq_penalty_weight": mapq_penalty_weight,
        "non_primary_weight": non_primary_weight,
        "alignment_quality_weight": alignment_quality_weight,

        "suspicious_genomes": sorted(suspicious_genomes),
        "suspicion_details": details,

        "total_cleanup_moves": 0,
        "skipped_no_alternative": 0,
        "skipped_only_suspicious_alternatives": 0,

        "moves_from": defaultdict(int),
        "moves_to": defaultdict(int),
    }

    for read_id in sorted(grouped_reads.keys()):
        if read_id not in cleaned_assignments:
            continue

        current_alignment = cleaned_assignments[read_id]
        current_genome = current_alignment["genome_id"]

        if current_genome not in suspicious_genomes:
            continue

        alignments = grouped_reads[read_id]

        if not alignments:
            continue

        best_mapq = max(aln["mapq"] for aln in alignments)

        candidate_alignments = [
            aln for aln in alignments
            if aln["mapq"] >= best_mapq - cleanup_mapq_delta
        ]

        candidate_alignments = [
            aln for aln in candidate_alignments
            if aln["genome_id"] != current_genome
        ]

        if not candidate_alignments:
            summary["skipped_no_alternative"] += 1
            continue

        non_suspicious_candidates = [
            aln for aln in candidate_alignments
            if aln["genome_id"] not in suspicious_genomes
            and aln["genome_id"] in genome_lengths
        ]

        if not non_suspicious_candidates:
            summary["skipped_only_suspicious_alternatives"] += 1
            continue

        best_quality_alignment = max(
            non_suspicious_candidates,
            key=lambda aln: (
                aln.get("AS") if aln.get("AS") is not None else -1,
                aln.get("s1") if aln.get("s1") is not None else -1,
                aln.get("cm") if aln.get("cm") is not None else -1,
                -(aln.get("NM") if aln.get("NM") is not None else 10**9),
            )
        )

        current_contrib = alignment_bucket_contributions(
            current_alignment,
            bucket_size=bucket_size
        )

        current_mean = genome_mean_coverage(
            coverage,
            current_genome,
            genome_lengths,
            bucket_size
        )

        remove_delta = local_sse_before_after_remove(
            coverage,
            current_genome,
            current_contrib,
            current_mean
        )

        best_candidate = None
        best_score = None

        for candidate_alignment in non_suspicious_candidates:
            candidate_genome = candidate_alignment["genome_id"]

            candidate_contrib = alignment_bucket_contributions(
                candidate_alignment,
                bucket_size=bucket_size
            )

            candidate_mean = genome_mean_coverage(
                coverage,
                candidate_genome,
                genome_lengths,
                bucket_size
            )

            add_delta = local_sse_before_after_add(
                coverage,
                candidate_genome,
                candidate_contrib,
                candidate_mean
            )

            mapq_penalty = mapq_penalty_weight * (
                best_mapq - candidate_alignment["mapq"]
            )

            non_primary_penalty = (
                0.0 if candidate_alignment.get("is_primary") else non_primary_weight
            )

            quality_penalty = alignment_quality_penalty(
                candidate_alignment=candidate_alignment,
                best_alignment=best_quality_alignment,
                weight=alignment_quality_weight
            )

            score = (
                remove_delta
                + add_delta
                + mapq_penalty
                + non_primary_penalty
                + quality_penalty
            )

            if best_score is None or score < best_score:
                best_score = score
                best_candidate = candidate_alignment

        if best_candidate is None:
            summary["skipped_no_alternative"] += 1
            continue

        new_genome = best_candidate["genome_id"]

        add_alignment_to_coverage(
            coverage,
            current_alignment,
            bucket_size=bucket_size,
            weight=-1.0
        )

        add_alignment_to_coverage(
            coverage,
            best_candidate,
            bucket_size=bucket_size,
            weight=1.0
        )

        cleaned_assignments[read_id] = best_candidate

        summary["total_cleanup_moves"] += 1
        summary["moves_from"][current_genome] += 1
        summary["moves_to"][new_genome] += 1

    summary["moves_from"] = dict(summary["moves_from"])
    summary["moves_to"] = dict(summary["moves_to"])

    return cleaned_assignments, summary

# =============================================================================
# IZVJEŠTAJ
# =============================================================================

def write_cleanup_simple_summary(
    summary,
    output_dir,
    filename="cleanup_simple_summary.txt"
):
    os.makedirs(output_dir, exist_ok=True)

    output_path = os.path.join(output_dir, filename)

    with open(output_path, "w", encoding="utf-8") as f:
        f.write("SAZETAK CLEANUP SIMPLE FAZE\n")
        f.write("=" * 80 + "\n\n")

        f.write(f"Algoritam: {summary.get('algorithm')}\n")
        f.write(f"Suspicion threshold: {summary.get('suspicion_threshold')}\n")
        f.write(f"Cleanup MAPQ delta: {summary.get('cleanup_mapq_delta')}\n")
        f.write(f"MAPQ penalty weight: {summary.get('mapq_penalty_weight')}\n")
        f.write(f"Non-primary weight: {summary.get('non_primary_weight')}\n")
        f.write(f"Alignment quality weight: {summary.get('alignment_quality_weight')}\n")
        f.write(f"Ukupan broj cleanup premjestanja: {summary.get('total_cleanup_moves', 0)}\n")
        f.write(f"Preskoceno - nema alternative: {summary.get('skipped_no_alternative', 0)}\n")
        f.write(
            "Preskoceno - postoje samo sumnjive alternative: "
            f"{summary.get('skipped_only_suspicious_alternatives', 0)}\n\n"
        )

        f.write("PARAMETRI DETEKCIJE\n")
        f.write("-" * 80 + "\n")
        f.write(f"zero start                        : {summary.get('zero_start')}\n")
        f.write(f"zero full                         : {summary.get('zero_full')}\n")
        f.write(f"long zero start                   : {summary.get('long_zero_start')}\n")
        f.write(f"long zero full                    : {summary.get('long_zero_full')}\n")
        f.write(f"top fraction                      : {summary.get('top_fraction')}\n")
        f.write(f"top coverage start                : {summary.get('top_coverage_start')}\n")
        f.write(f"top coverage full                 : {summary.get('top_coverage_full')}\n")
        f.write(f"spike ratio start                 : {summary.get('spike_ratio_start')}\n")
        f.write(f"spike ratio full                  : {summary.get('spike_ratio_full')}\n")
        f.write(f"use drainage signal               : {summary.get('use_drainage_signal')}\n")
        f.write(f"drainage rel rmse drop start      : {summary.get('drainage_rel_rmse_drop_start')}\n")
        f.write(f"drainage rel rmse drop full       : {summary.get('drainage_rel_rmse_drop_full')}\n\n")

        f.write("SUMNJIVI GENOMI ZA CLEANUP\n")
        f.write("-" * 80 + "\n")

        suspicious = summary.get("suspicious_genomes", [])

        if not suspicious:
            f.write("Nema genoma oznacenih za cleanup.\n\n")
        else:
            for genome_id in suspicious:
                f.write(f"{genome_id}\n")
            f.write("\n")

        f.write("DETALJI SUMNJIVOSTI PO GENOMU\n")
        f.write("-" * 80 + "\n")

        details = summary.get("suspicion_details", {})

        for genome_id, d in sorted(details.items()):
            f.write(f"Genom: {genome_id}\n")
            f.write(f"  cleanup score             : {d['cleanup_score']:.6f}\n")
            f.write(f"  primary score             : {d['primary_score']:.6f}\n\n")

            f.write("  Oblik coveragea:\n")
            f.write(f"    zero fraction            : {d['zero_fraction']:.6f}\n")
            f.write(f"    zero part                : {d['zero_part']:.6f}\n")
            f.write(f"    long zero fraction       : {d['long_zero_fraction']:.6f}\n")
            f.write(f"    long zero part           : {d['long_zero_part']:.6f}\n")
            f.write(f"    bad shape part           : {d['bad_shape_part']:.6f}\n")
            f.write(f"    base score               : {d['base_score']:.6f}\n\n")

            f.write("  Island / koncentracija coveragea:\n")
            f.write(f"    top fraction             : {d['top_fraction']:.6f}\n")
            f.write(f"    top coverage fraction    : {d['top_coverage_fraction']:.6f}\n")
            f.write(f"    top coverage part        : {d['top_coverage_part']:.6f}\n")
            f.write(f"    island score             : {d['island_score']:.6f}\n\n")

            f.write("  Spike / lokalni siljci:\n")
            f.write(f"    spike ratio              : {d['spike_ratio']:.6f}\n")
            f.write(f"    spike ratio part         : {d['spike_ratio_part']:.6f}\n")
            f.write(f"    spike score              : {d['spike_score']:.6f}\n\n")

            f.write("  Drainage / praznjenje:\n")
            f.write(f"    baseline RMSE od nule    : {d['baseline_rmse_from_zero']:.6f}\n")
            f.write(f"    current RMSE od nule     : {d['current_rmse_from_zero']:.6f}\n")
            f.write(f"    RMSE drop                : {d['rmse_drop_from_baseline']:.6f}\n")
            f.write(f"    relative RMSE drop       : {d['relative_rmse_drop_from_baseline']:.6f}\n")
            f.write(f"    rel drop part            : {d['drainage_rel_drop_part']:.6f}\n")
            f.write(f"    drainage score           : {d['drainage_score']:.6f}\n\n")

        f.write("PREMJESTANJA IZ GENOMA\n")
        f.write("-" * 80 + "\n")
        for genome_id, move_count in sorted(summary.get("moves_from", {}).items()):
            f.write(f"{genome_id}: {move_count}\n")
        f.write("\n")

        f.write("PREMJESTANJA U GENOME\n")
        f.write("-" * 80 + "\n")
        for genome_id, move_count in sorted(summary.get("moves_to", {}).items()):
            f.write(f"{genome_id}: {move_count}\n")
        f.write("\n")

    return output_path
