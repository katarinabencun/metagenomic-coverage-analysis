import os
import random
from collections import defaultdict

def get_best_alignments(alignments):
    if not alignments:
        return []

    best_mapq = max(aln["mapq"] for aln in alignments)

    return [
        aln for aln in alignments
        if aln["mapq"] == best_mapq
    ]


def get_candidate_genomes(best_alignments):
    return tuple(sorted(set(aln["genome_id"] for aln in best_alignments)))


def group_reads_by_candidates(grouped_reads):
    groups = defaultdict(list)
    unique_assignments = {}

    for read_id, alignments in grouped_reads.items():
        best_alignments = get_best_alignments(alignments)

        if not best_alignments:
            continue

        candidate_genomes = get_candidate_genomes(best_alignments)

        if len(candidate_genomes) == 1:
            unique_assignments[read_id] = find_alignment_for_genome(
                best_alignments,
                candidate_genomes[0]
            )
        else:
            groups[candidate_genomes].append({
                "read_id": read_id,
                "best_alignments": best_alignments
            })

    return unique_assignments, groups


def find_alignment_for_genome(best_alignments, genome_id):
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
            aln["genome_id"],
            aln["start"],
            aln["end"],
            aln["mapq"]
        )
    )

    return candidates[0]


def balanced_pseudorandom_group_assignment(groups, seed=42):
    assignments = {}

    for group_index, (candidate_genomes, read_items) in enumerate(sorted(groups.items())):
        # Poseban RNG za svaku grupu.
        # Time se dobiva ponovljivost, ali grupe ne ovise jedna o drugoj.
        group_seed = f"{seed}_{group_index}_{'_'.join(candidate_genomes)}"
        rng = random.Random(group_seed)

        # Stabilni početni poredak prije miješanja
        read_items = sorted(read_items, key=lambda item: item["read_id"])
        candidate_genomes = list(candidate_genomes)

        # Pseudo-slučajno, ali ponovljivo miješanje
        rng.shuffle(read_items)
        rng.shuffle(candidate_genomes)

        num_reads = len(read_items)
        num_candidates = len(candidate_genomes)

        base_quota = num_reads // num_candidates
        remainder = num_reads % num_candidates

        quotas = {}

        for i, genome_id in enumerate(candidate_genomes):
            quotas[genome_id] = base_quota

            if i < remainder:
                quotas[genome_id] += 1

        read_pos = 0

        for genome_id in candidate_genomes:
            quota = quotas[genome_id]

            for _ in range(quota):
                item = read_items[read_pos]
                read_id = item["read_id"]
                best_alignments = item["best_alignments"]

                selected_alignment = find_alignment_for_genome(
                    best_alignments,
                    genome_id
                )

                if selected_alignment is not None:
                    assignments[read_id] = selected_alignment

                read_pos += 1

    return assignments


def initial_balanced_assignment(grouped_reads, seed=42):
    unique_assignments, groups = group_reads_by_candidates(grouped_reads)

    group_assignments = balanced_pseudorandom_group_assignment(
        groups,
        seed=seed
    )

    assignments = {}
    assignments.update(unique_assignments)
    assignments.update(group_assignments)

    return assignments


def assignment_counts(assignments):
    counts = defaultdict(int)

    for read_id, aln in assignments.items():
        counts[aln["genome_id"]] += 1

    return dict(counts)


def group_summary(grouped_reads):
    _, groups = group_reads_by_candidates(grouped_reads)

    summary = {}

    for candidate_genomes, read_items in groups.items():
        summary[candidate_genomes] = len(read_items)

    return summary

def assignment_reduction_summary(grouped_reads, assignments):
    """
    Računa sažetak prije i poslije početne diskretne dodjele.
    """
    total_reads = len(grouped_reads)
    total_alignments_before = 0
    total_best_alignments = 0
    uniquely_mapped_reads = 0
    multi_best_reads = 0
    ignored_weaker_alignments = 0

    for read_id, alignments in grouped_reads.items():
        if not alignments:
            continue

        total_alignments_before += len(alignments)

        best_alignments = get_best_alignments(alignments)
        total_best_alignments += len(best_alignments)

        ignored_weaker_alignments += len(alignments) - len(best_alignments)

        candidate_genomes = get_candidate_genomes(best_alignments)

        if len(candidate_genomes) == 1:
            uniquely_mapped_reads += 1
        else:
            multi_best_reads += 1

    assigned_reads = len(assignments)
    removed_by_discrete_choice = total_best_alignments - assigned_reads
    total_removed_alignments = total_alignments_before - assigned_reads

    return {
        "total_reads": total_reads,
        "total_alignments_before": total_alignments_before,
        "total_best_alignments": total_best_alignments,
        "assigned_reads_after": assigned_reads,
        "uniquely_mapped_reads": uniquely_mapped_reads,
        "multi_best_reads": multi_best_reads,
        "ignored_weaker_alignments": ignored_weaker_alignments,
        "removed_by_discrete_choice": removed_by_discrete_choice,
        "total_removed_alignments": total_removed_alignments,
    }


def write_initial_assignment_summary(
    grouped_reads,
    assignments,
    output_dir,
    filename="summary_init_assign.txt"
):
    """
    Sprema sažetak početne diskretne dodjele u .txt datoteku.
    """
    os.makedirs(output_dir, exist_ok=True)

    summary = assignment_reduction_summary(grouped_reads, assignments)
    counts = assignment_counts(assignments)
    groups = group_summary(grouped_reads)

    output_path = os.path.join(output_dir, filename)

    with open(output_path, "w") as f:
        f.write("SAZETAK POCETNE GRUPNE URAVNOTEZENE DODJELE\n")
        f.write("=" * 70 + "\n\n")

        f.write("OPCI SAZETAK\n")
        f.write("-" * 70 + "\n")
        f.write(f"Ukupan broj ocitanja u PAF-u: {summary['total_reads']}\n")
        f.write(f"Ukupan broj mapiranja prije dodjele: {summary['total_alignments_before']}\n")
        f.write(
            "Broj najboljih MAPQ mapiranja nakon odbacivanja slabijih: "
            f"{summary['total_best_alignments']}\n"
        )
        f.write(f"Broj konacno dodijeljenih ocitanja: {summary['assigned_reads_after']}\n\n")

        f.write("STRUKTURA OCITANJA\n")
        f.write("-" * 70 + "\n")
        f.write(f"Jednoznacna ocitanja: {summary['uniquely_mapped_reads']}\n")
        f.write(f"Visestruko jednako dobra ocitanja: {summary['multi_best_reads']}\n\n")

        f.write("MAKNUTO / IGNORIRANO\n")
        f.write("-" * 70 + "\n")
        f.write(
            "Slabija mapiranja ignorirana zbog losijeg MAPQ-a: "
            f"{summary['ignored_weaker_alignments']}\n"
        )
        f.write(
            "Visak najboljih mapiranja maknut diskretnim izborom jednog kandidata: "
            f"{summary['removed_by_discrete_choice']}\n"
        )
        f.write(f"Ukupno maknutih mapiranja: {summary['total_removed_alignments']}\n\n")

        if summary["total_alignments_before"] > 0:
            kept_percent = 100 * summary["assigned_reads_after"] / summary["total_alignments_before"]
            removed_percent = 100 * summary["total_removed_alignments"] / summary["total_alignments_before"]

            f.write("OMJERI\n")
            f.write("-" * 70 + "\n")
            f.write(f"Zadrzano mapiranja nakon dodjele: {kept_percent:.2f}%\n")
            f.write(f"Maknuto mapiranja nakon dodjele: {removed_percent:.2f}%\n\n")

        f.write("BROJ OCITANJA PO GENOMU NAKON POCETNE DODJELE\n")
        f.write("-" * 70 + "\n")
        for genome_id, count in sorted(counts.items()):
            f.write(f"{genome_id}: {count}\n")

        f.write("\nGRUPE VISESTRUKO JEDNAKO DOBRO MAPIRANIH OCITANJA\n")
        f.write("-" * 70 + "\n")
        for candidate_genomes, count in sorted(groups.items(), key=lambda x: (-x[1], x[0])):
            f.write(f"{candidate_genomes}: {count}\n")

    return output_path
