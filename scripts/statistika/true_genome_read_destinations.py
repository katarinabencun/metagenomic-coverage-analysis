import os
from collections import Counter, defaultdict
from statistika.assignment_evaluation import parse_fastq_origin_sets
from statistika.chimera_types import parse_fastq_chimera_types, format_chimera_counts


def write_true_genome_read_destinations(fastq_path, final_assignments, output_dir):
    """Odredista ocitanja po stvarnom genomu, uz oznake vrsta kimera.

    Ocitanje s dva izvorna genoma pojavljuje se u odjeljku oba genoma.
    """
    os.makedirs(output_dir, exist_ok=True)
    origins = parse_fastq_origin_sets(fastq_path)
    chimera_types = parse_fastq_chimera_types(fastq_path)
    destinations = defaultdict(Counter)
    chimera_counts = defaultdict(lambda: defaultdict(Counter))
    no_origin = Counter()
    no_origin_chimeras = defaultdict(Counter)

    unknown = set(final_assignments) - set(origins)
    if unknown:
        raise ValueError(f"Dodijeljeni read_id nije u FASTQ-u: {sorted(unknown)[0]}")

    for read_id, genomes in origins.items():
        aln = final_assignments.get(read_id)
        destination = (aln.get("genome_id") if aln else None) or "UNASSIGNED"
        kind = chimera_types[read_id]
        if not genomes:
            no_origin[destination] += 1
            if kind:
                no_origin_chimeras[destination][kind] += 1
        for true_genome in genomes:
            destinations[true_genome][destination] += 1
            if kind:
                chimera_counts[true_genome][destination][kind] += 1

    output_path = os.path.join(output_dir, "true_genomes_read_destinations.txt")
    with open(output_path, "w", encoding="utf-8") as out:
        out.write("NOTE: Chimeric reads are counted under EACH origin genome.\n")
        out.write("An assignment to ANY origin genome is considered correct.\n")
        out.write("Totals across genome sections are not unique read counts.\n\n")
        for true_genome in sorted(destinations):
            counts = destinations[true_genome]
            total = sum(counts.values())
            out.write(f"TRUE GENOME: {true_genome}\n\nSimulator reads: {total}\n\n")
            out.write("Final assignments:\n")
            n = counts.get(true_genome, 0)
            out.write(f"    {true_genome}: {n}{format_chimera_counts(chimera_counts[true_genome][true_genome])}\n")
            for genome, n in sorted(counts.items(), key=lambda x: (-x[1], x[0])):
                if genome not in (true_genome, "UNASSIGNED"):
                    out.write(f"    {genome}: {n}{format_chimera_counts(chimera_counts[true_genome][genome])}\n")
            n = counts.get("UNASSIGNED", 0)
            out.write(f"    UNASSIGNED: {n}{format_chimera_counts(chimera_counts[true_genome]['UNASSIGNED'])}\n")
            out.write("\n" + "-" * 60 + "\n\n")
        out.write("READS WITHOUT GENOMIC ORIGIN (random_seq/junk_seq):\n")
        out.write(f"Total: {sum(no_origin.values())}\n")
        for genome, n in sorted(no_origin.items(), key=lambda x: (-x[1], x[0])):
            out.write(f"    {genome}: {n}{format_chimera_counts(no_origin_chimeras[genome])}\n")
    return output_path
