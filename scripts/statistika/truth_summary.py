import os


def write_genome_truth_summary(
    sim_counts,
    genome_lengths,
    output_dir=None,
    filename="genome_truth_summary.txt"
):
    """
    Sprema pregled stvarno prisutnih i lažnih genoma.

    Pravi genomi:
        postoje u simulatoru, tj. imaju barem jedno očitanje u FASTQ-u.

    Lažni genomi:
        postoje u referentnoj bazi,
        ali nemaju očitanja u simulatoru.
    """
    if output_dir is None:
        output_dir = "results/statistika/summarys"

    os.makedirs(output_dir, exist_ok=True)

    output_path = os.path.join(output_dir, filename)

    true_genomes = sorted(
        genome_id
        for genome_id, count in sim_counts.items()
        if count > 0
    )

    false_genomes = sorted(
        genome_id
        for genome_id in genome_lengths
        if genome_id not in true_genomes
    )

    with open(output_path, "w") as f:
        f.write("PREGLED PRAVIH I LAZNIH GENOMA\n")
        f.write("=" * 80 + "\n\n")

        f.write("Definicija:\n")
        f.write(
            "  Pravi genomi su oni iz kojih su ocitanja stvarno generirana "
            "u simulatoru.\n"
        )
        f.write(
            "  Lazni genomi su oni koji postoje u referentnoj bazi, "
            "ali nemaju ocitanja u simulatoru.\n\n"
        )

        f.write(f"Ukupan broj referentnih genoma: {len(genome_lengths)}\n")
        f.write(f"Broj pravih genoma: {len(true_genomes)}\n")
        f.write(f"Broj laznih genoma: {len(false_genomes)}\n\n")

        f.write("PRAVI GENOMI / IMAJU OCITANJA U SIMULATORU\n")
        f.write("-" * 80 + "\n")

        for genome_id in true_genomes:
            f.write(f"{genome_id}: {sim_counts.get(genome_id, 0)} ocitanja\n")

        f.write("\n")

        f.write("LAZNI GENOMI / NEMAJU OCITANJA U SIMULATORU\n")
        f.write("-" * 80 + "\n")

        for genome_id in false_genomes:
            f.write(f"{genome_id}: 0 ocitanja\n")

    return output_path