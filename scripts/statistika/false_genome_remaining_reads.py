import os
from collections import Counter, defaultdict
from statistika.assignment_evaluation import parse_fastq_origin_sets
from statistika.chimera_types import parse_fastq_chimera_types


def write_false_genome_remaining_reads(fastq_path, final_assignments, sim_counts,
                                       genome_lengths, cleanup_summary, output_dir):
    """Preostala ocitanja na laznim genomima, podrijetlo i vrsta kimere."""
    os.makedirs(output_dir, exist_ok=True)
    false_genomes = {g for g in genome_lengths if sim_counts.get(g, 0) == 0}
    suspicious = set(cleanup_summary.get("suspicious_genomes", []))
    origins = parse_fastq_origin_sets(fastq_path)
    chimera_types = parse_fastq_chimera_types(fastq_path)
    grouped = defaultdict(Counter)

    for read_id, aln in final_assignments.items():
        genome_id = aln["genome_id"]
        if genome_id not in false_genomes:
            continue
        if read_id not in origins:
            raise ValueError(f"Ocitanje {read_id} nije pronadjeno u FASTQ-u")
        genome_origins = origins[read_id]
        kind = chimera_types[read_id]
        label = " + ".join(genome_origins) if genome_origins else "NO_GENOMIC_ORIGIN"
        if kind == "chimera_w/junk":
            label += " + junk_seq"
        elif kind == "chimera_w/random":
            label += " + random_seq"
        # Grupiranje po (podrijetlo, vrsta) cuva oznake i ako imaju jedan izvor.
        grouped[genome_id][(label, kind)] += 1

    for filename, genomes in (
        ("false_genomes_undetected.txt", sorted(false_genomes - suspicious)),
        ("false_genomes_detected_remaining.txt", sorted(false_genomes & suspicious)),
    ):
        with open(os.path.join(output_dir, filename), "w", encoding="utf-8") as out:
            if not genomes:
                out.write("Nema genoma u ovoj skupini.\n")
            for g in genomes:
                counts = grouped[g]
                total = sum(counts.values())
                out.write(f"FALSE GENOME: {g}\n")
                out.write(f"Assigned reads after cleanup: {total}\n\n")
                out.write("True origin according to simulator:\n")
                for (label, kind), n in sorted(counts.items(), key=lambda x: (-x[1], x[0][0], x[0][1] or "")):
                    suffix = f" ({kind})" if kind else ""
                    out.write(f"    {label}: {n}{suffix}\n")
                if not counts:
                    out.write("    No remaining reads.\n")
                out.write(f"\nTOTAL: {total}\n" + "-" * 60 + "\n\n")
