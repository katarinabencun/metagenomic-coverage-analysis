"""Vrsta kimere iz FASTQ zaglavlja; ne mijenja evaluacijsku istinu."""

from statistika.assignment_evaluation import open_text_maybe_gzip


CHIMERA_TYPES = ("chimera_w/genom", "chimera_w/junk", "chimera_w/random")


def parse_fastq_chimera_types(fastq_path):
    """read_id -> oznaka vrste kimere ili None za obicno ocitanje.

    Za ovaj format simulatora: `chimera` odvaja dva segmenta.
    Nepoznate kombinacije prijavljujemo, ne oznacavamo pogresno.
    """
    result = {}
    with open_text_maybe_gzip(fastq_path) as fastq:
        while True:
            header = fastq.readline()
            if not header:
                break
            if not (fastq.readline() and fastq.readline() and fastq.readline()):
                raise ValueError(f"Nepotpun FASTQ zapis: {header.strip()}")
            if not header.startswith("@"):
                raise ValueError(f"Neispravno FASTQ zaglavlje: {header.strip()}")
            parts = header.split()
            read_id = parts[0][1:]
            if read_id in result:
                raise ValueError(f"Duplicirani read_id: {read_id}")
            if "chimera" not in parts[1:]:
                result[read_id] = None
                continue
            segments = [t for t in parts[1:] if t.startswith("kraken:taxid|")]
            junk = "junk_seq" in parts[1:]
            random = "random_seq" in parts[1:]
            if len(segments) >= 2 and not junk and not random:
                kind = "chimera_w/genom"
            elif len(segments) == 1 and junk and not random:
                kind = "chimera_w/junk"
            elif len(segments) == 1 and random and not junk:
                kind = "chimera_w/random"
            else:
                raise ValueError(f"Nepoznata kombinacija dijelova kimere: {header.strip()}")
            result[read_id] = kind
    return result


def format_chimera_counts(counts):
    """Npr. ' (3 chimera_w/genom, 1 chimera_w/junk)' ili ''."""
    fields = [f"{counts[k]} {k}" for k in CHIMERA_TYPES if counts.get(k, 0)]
    return " (" + ", ".join(fields) + ")" if fields else ""
