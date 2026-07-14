import re

def normalize_genome_id(header):
    if "|" in header:
        return header.split("|")[-1].split()[0]

    return header.split()[0]


def parse_cigar(cigar):
    return [(int(length), op) for length, op in re.findall(r"(\d+)([MIDNSHP=X])", cigar)]


def parse_paf_grouped_with_cigar(file_path):
    reads = {}

    with open(file_path) as f:
        for line in f:
            fields = line.strip().split()

            if len(fields) < 12:
                continue

            cigar = None
            tp = None
            as_score = None
            s1_score = None
            s2_score = None
            cm_score = None
            nm = None
            de = None

            for field in fields[12:]:
                if field.startswith("cg:Z:"):
                    cigar = field[5:]
                elif field.startswith("tp:A:"):
                    tp = field[5:]
                elif field.startswith("AS:i:"):
                    as_score = int(field[5:])
                elif field.startswith("s1:i:"):
                    s1_score = int(field[5:])
                elif field.startswith("s2:i:"):
                    s2_score = int(field[5:])
                elif field.startswith("cm:i:"):
                    cm_score = int(field[5:])
                elif field.startswith("NM:i:"):
                    nm = int(field[5:])
                elif field.startswith("de:f:"):
                    de = float(field[5:])

            if cigar is None:
                continue

            try:
                read_id = fields[0]
                query_length = int(fields[1])
                query_start = int(fields[2])
                query_end = int(fields[3])
                genome_id = normalize_genome_id(fields[5])
                start = int(fields[7])
                end = int(fields[8])
                mapq = int(fields[11])
            except (IndexError, ValueError):
                continue

            if read_id not in reads:
                reads[read_id] = []

            reads[read_id].append({
                "genome_id": genome_id,
                "start": start,
                "end": end,
                "mapq": mapq,
                "cigar": cigar,
                "tp": tp,
                "is_primary": tp == "P",

                "query_length": query_length,
                "query_start": query_start,
                "query_end": query_end,

                "AS": as_score,
                "s1": s1_score,
                "s2": s2_score,
                "cm": cm_score,
                "NM": nm,
                "de": de,
            })

    return reads


def parse_paf_grouped(file_path):
    reads = {}

    with open(file_path) as f:
        for line in f:
            fields = line.strip().split()

            if len(fields) < 12:
                continue

            try:
                read_id = fields[0]
                query_length = int(fields[1])
                query_start = int(fields[2])
                query_end = int(fields[3])
                genome_id = normalize_genome_id(fields[5])
                start_pos = int(fields[7])
                end_pos = int(fields[8])
                mapq = int(fields[11])
            except (IndexError, ValueError):
                continue

            tp = None
            for field in fields[12:]:
                if field.startswith("tp:A:"):
                    tp = field[5:]
                    break

            if read_id not in reads:
                reads[read_id] = []

            reads[read_id].append({
                "genome_id": genome_id,
                "start": start_pos,
                "end": end_pos,
                "mapq": mapq,
                "tp": tp,
                "is_primary": tp == "P",
                "query_length": query_length,
                "query_start": query_start,
                "query_end": query_end,
                "raw_line": line.strip()
            })

    return reads
