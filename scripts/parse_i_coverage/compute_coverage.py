from parse_i_coverage.parse_paf import parse_cigar

def add_interval_fractional(coverage, genome_id, start, end, bucket_size=1000, weight=1.0):
    start_bucket = start // bucket_size
    end_bucket = (end - 1) // bucket_size

    for b in range(start_bucket, end_bucket + 1):
        bucket_start = b * bucket_size
        bucket_end = bucket_start + bucket_size

        overlap_start = max(start, bucket_start)
        overlap_end = min(end, bucket_end)

        overlap = overlap_end - overlap_start

        if overlap > 0:
            coverage[genome_id][b] = coverage[genome_id].get(b, 0) + overlap * weight 

def compute_coverage_from_assignments(assignments, bucket_size=1000):
    """
    Računa coverage iz diskretne dodjele očitanja.

    assignments očekuje oblik:
    {
        read_id: {
            "genome_id": ...,
            "start": ...,
            "end": ...,
            "mapq": ...
        }
    }

    Svako očitanje doprinosi samo jednom genomu.
    """
    coverage = {}

    for read_id, aln in assignments.items():
        genome_id = aln["genome_id"]
        start = aln["start"]
        end = aln["end"]

        if genome_id not in coverage:
            coverage[genome_id] = {}

        start_bucket = start // bucket_size
        end_bucket = (end - 1) // bucket_size

        for b in range(start_bucket, end_bucket + 1):
            bucket_start = b * bucket_size
            bucket_end = bucket_start + bucket_size

            overlap_start = max(start, bucket_start)
            overlap_end = min(end, bucket_end)

            overlap = overlap_end - overlap_start

            if overlap > 0:
                coverage[genome_id][b] = coverage[genome_id].get(b, 0) + overlap

    for genome_id in coverage:
        for b in coverage[genome_id]:
            coverage[genome_id][b] /= bucket_size

    return coverage

def compute_coverage_from_assignments_cigar(assignments, bucket_size=1000):
    coverage = {}

    for read_id, aln in assignments.items():
        genome_id = aln["genome_id"]
        ref_pos = aln["start"]
        cigar = aln["cigar"]

        if genome_id not in coverage:
            coverage[genome_id] = {}

        for length, op in parse_cigar(cigar):
            if op in ("M", "=", "X"):
                segment_start = ref_pos
                segment_end = ref_pos + length

                add_interval_fractional(
                    coverage,
                    genome_id,
                    segment_start,
                    segment_end,
                    bucket_size=bucket_size,
                    weight=1.0
                )

                ref_pos += length

            elif op in ("D", "N"):
                ref_pos += length

            elif op in ("I", "S", "H", "P"):
                continue

    for genome_id in coverage:
        for b in coverage[genome_id]:
            coverage[genome_id][b] /= bucket_size

    return coverage

#originalna za stvarnu (simulator) raspodejlu
def compute_coverage_from_intervals(reads, bucket_size=1000):
    coverage = {}

    for genome_id, start, end in reads:

        if genome_id not in coverage:
            coverage[genome_id] = {}

        start_bucket = start // bucket_size
        end_bucket = (end - 1) // bucket_size

        for b in range(start_bucket, end_bucket + 1):

            bucket_start = b * bucket_size
            bucket_end = bucket_start + bucket_size

            # overlap
            overlap_start = max(start, bucket_start)
            overlap_end = min(end, bucket_end)

            overlap = overlap_end - overlap_start

            if overlap > 0:
                coverage[genome_id][b] = coverage[genome_id].get(b, 0) + overlap

    # normalizacija - dijeljenje s bucket_size da dobijem decimalni broj
    for genome_id in coverage:
        for b in coverage[genome_id]:
            coverage[genome_id][b] /= bucket_size

    return coverage