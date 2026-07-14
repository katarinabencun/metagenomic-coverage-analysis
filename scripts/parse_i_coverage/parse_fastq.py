import gzip


def open_text_maybe_gzip(file_path):
    if file_path.endswith(".gz"):
        return gzip.open(file_path, "rt")
    return open(file_path, "r")

def parse_fastq(file_path):
    with open_text_maybe_gzip(file_path) as f:
        while True:
            header = f.readline().strip()
            seq = f.readline()
            plus = f.readline()
            qual = f.readline()

            if not header:
                break

            parts = header.split()
            if len(parts) < 2:
                continue

            try:
                genome_part = parts[1].split("|")[-1]
                genome_id = genome_part.split(",")[0]
                position_range = genome_part.split(",")[2]
                start_pos = int(position_range.split("-")[0])
                end_pos = int(position_range.split("-")[1])
            except:
                continue

            yield genome_id, start_pos, end_pos