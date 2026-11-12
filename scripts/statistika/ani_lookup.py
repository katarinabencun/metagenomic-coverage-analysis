"""Ucitavanje ANI/AF tablice i ispis usporedbi u dodatnim izvjestajima.

ANI, AF_1 i AF_2 iz data/ani/<dataset>/ani_pairs.tsv su postoci.
AF_1 pripada genome_1, a AF_2 pripada genome_2, neovisno o
redoslijedu u kojem genomima pristupamo u izvjestajima.
"""

import csv
import math
from pathlib import Path


DATASETS = ("5b15s_add", "5b15s", "2b4s")
PROJECT_ROOT = Path(__file__).resolve().parents[2]
REQUIRED_COLUMNS = ("genome_1", "genome_2", "ANI", "AF_1", "AF_2")


def resolve_ani_table_path(fasta_path, explicit_path=None):
    """Automatski bira tablicu za referentni skup ili koristi --ani-table.

    Prepoznaje se oznaka skupa u putanji FASTA datoteke, npr.
    data/5b15s/5b15s/c4/....fasta -> data/ani/5b15s/ani_pairs.tsv.
    Ako skup nije prepoznat, vraca None (izvjestaji koriste NA).
    """
    if explicit_path is not None:
        return Path(explicit_path).expanduser().resolve()

    parts = Path(fasta_path).expanduser().resolve().parts
    dataset = next((name for name in DATASETS if name in parts), None)
    if dataset is None:
        return None
    return PROJECT_ROOT / "data" / "ani" / dataset / "ani_pairs.tsv"


def _percentage(value, name, line_no):
    try:
        number = float(value)
    except (TypeError, ValueError):
        raise ValueError(f"ANI tablica, redak {line_no}: neispravan {name}: {value!r}") from None
    if not math.isfinite(number) or not 0 <= number <= 100:
        raise ValueError(f"ANI tablica, redak {line_no}: {name} izvan raspona 0-100: {value!r}")
    return number


def load_ani_pairs(path, *, required=False):
    """(sortirani genom_1, genom_2) -> (ANI, AF_sortirani_1, AF_sortirani_2).

    Nedostajuca tablica je dopustena samo kod automatskog pronalazenja,
    jer se tada dogovoreno ispisuje NA umjesto lazne vrijednosti 0.
    """
    if path is None:
        print("[ANI] Skup nije prepoznat iz FASTA putanje: ANI/AF u izvjestajima = NA."
              " Po potrebi zadaj --ani-table.")
        return {}

    path = Path(path)
    if not path.is_file():
        if required:
            raise FileNotFoundError(f"Zadana ANI tablica ne postoji: {path}")
        print(f"[ANI] Tablica nije pronadena: {path}. ANI/AF u izvjestajima = NA.")
        return {}

    pairs = {}
    with path.open("r", encoding="utf-8-sig", newline="") as file:
        reader = csv.DictReader(file, delimiter="\t")
        if not reader.fieldnames or not set(REQUIRED_COLUMNS).issubset(reader.fieldnames):
            raise ValueError(
                f"ANI tablica nema stupce {', '.join(REQUIRED_COLUMNS)}: {path}"
            )
        for line_no, row in enumerate(reader, start=2):
            if None in row:
                raise ValueError(f"ANI tablica, redak {line_no}: visak stupaca")
            first = (row["genome_1"] or "").strip()
            second = (row["genome_2"] or "").strip()
            if not first or not second or first == second:
                raise ValueError(f"ANI tablica, redak {line_no}: neispravan par genoma")
            ani = _percentage(row["ANI"], "ANI", line_no)
            af1 = _percentage(row["AF_1"], "AF_1", line_no)
            af2 = _percentage(row["AF_2"], "AF_2", line_no)
            if first > second:
                first, second = second, first
                af1, af2 = af2, af1
            key = (first, second)
            if key in pairs:
                raise ValueError(f"ANI tablica: duplikat para {first}, {second}")
            pairs[key] = (ani, af1, af2)

    print(f"[ANI] Ucitano {len(pairs)} usporedbi iz {path}.")
    return pairs


def format_ani_comparison(pairs, first_genome, second_genome, first_af_label, second_af_label):
    """Formatira ANI i AF svake strane u smjeru izvjestaja; NA ako par nedostaje."""
    key = tuple(sorted((first_genome, second_genome)))
    record = pairs.get(key) if first_genome != second_genome else None
    if record is None:
        return f"ANI: NA | {first_af_label}: NA | {second_af_label}: NA"
    ani, af1, af2 = record
    af_first, af_second = (af1, af2) if first_genome == key[0] else (af2, af1)
    return (
        f"ANI: {ani:.2f}% | {first_af_label}: {af_first:.2f}%"
        f" | {second_af_label}: {af_second:.2f}%"
    )
