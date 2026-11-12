import os
from collections import Counter, defaultdict
from statistika.assignment_evaluation import parse_fastq_origin_sets
from statistika.chimera_types import parse_fastq_chimera_types
from statistika.ani_lookup import format_ani_comparison


def write_false_genome_remaining_reads(fastq_path, final_assignments, sim_counts,
                                       genome_lengths, cleanup_summary, output_dir,
                                       ani_pairs=None):
    """Preostala ocitanja na laznim genomima, podrijetlo i vrsta kimere."""
    os.makedirs(output_dir, exist_ok=True)
    ani_pairs = ani_pairs if ani_pairs is not None else {}
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
        # Grupiranje po (skup genomskih izvora, vrsta) zadrzava izvorne ID-jeve
        # za ANI/AF usporedbe bez parsiranja tekstualne oznake.
        grouped[genome_id][(genome_origins, kind)] += 1

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
                for (genome_origins, kind), n in sorted(
                    counts.items(), key=lambda x: (-x[1], x[0][0], x[0][1] or "")
                ):
                    label = " + ".join(genome_origins) if genome_origins else "NO_GENOMIC_ORIGIN"
                    if kind == "chimera_w/junk":
                        label += " + junk_seq"
                    elif kind == "chimera_w/random":
                        label += " + random_seq"
                    suffix = f" ({kind})" if kind else ""
                    out.write(f"    {label}: {n}{suffix}\n")

                    for origin in genome_origins:
                        comparison = format_ani_comparison(
                            ani_pairs, origin, g, "AF_origin", "AF_false"
                        )
                        if len(genome_origins) == 1:
                            out.write(f"        {comparison}\n")
                        else:
                            out.write(f"        {origin}: {comparison}\n")
                    if not genome_origins:
                        out.write("        ANI: NA | AF_origin: NA | AF_false: NA\n")
                if not counts:
                    out.write("    No remaining reads.\n")
                out.write(f"\nTOTAL: {total}\n" + "-" * 60 + "\n\n")
