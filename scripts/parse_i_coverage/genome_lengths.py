def get_genome_lengths(fasta_path):
    genome_lengths = {}
    current_id = None
    seq_parts = []

    with open(fasta_path) as f:
        for line in f:
            line = line.strip()

            if not line:
                continue

            if line.startswith(">"):
                if current_id is not None:
                    genome_lengths[current_id] = len("".join(seq_parts))

                header = line[1:]

                # ako header izgleda kao kraken:taxid|123|NZ_CP085939.1
                if "|" in header:
                    current_id = header.split("|")[-1].split()[0]
                else:
                    current_id = header.split()[0]

                seq_parts = []
            else:
                seq_parts.append(line)

    if current_id is not None:
        genome_lengths[current_id] = len("".join(seq_parts))

    return genome_lengths