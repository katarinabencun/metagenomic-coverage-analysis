from collections import defaultdict


def count_reads_simulator(reads):
    """
    Broji koliko je očitanja simulator stvarno generirao iz svakog genoma.

    Ulaz:
        reads: iterable zapisa oblika:
            (genome_id, start, end)

    Izlaz:
        {
            genome_id: broj_ocitanja
        }
    """
    counts = defaultdict(int)

    for genome_id, start, end in reads:
        counts[genome_id] += 1

    return dict(counts)