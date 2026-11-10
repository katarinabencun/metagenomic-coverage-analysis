#!/usr/bin/env python3
"""Izracun ANI/AF vrijednosti svih parova referentnih genoma pomocu skanija.

Fiksna naredba (ostali parametri ostaju na zadanim vrijednostima skanija):
    skani triangle -i -E REFERENCE.fa -o PRIVREMENI_IZLAZ.tsv

Jedini trajni izlaz: data/ani/<dataset>/ani_pairs.tsv
Stupci: genome_1, genome_2, ANI, AF_1, AF_2 (ANI i AF u postocima).

Za ponovljivost zabiljeziti verziju instaliranog skanija u projektnoj
README dokumentaciji; skripta pri pokretanju ispisuje verziju ako je dostupna.
Nema izmjena ostatka projekta niti automatskog pokretanja u eksperimentima.
"""

import argparse
import csv
import math
import subprocess
import sys
import tempfile
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[2]
DATASETS = ("2b4s", "5b15s", "5b15s_add")
OUT_COLUMNS = ("genome_1", "genome_2", "ANI", "AF_1", "AF_2")
REQUIRED_COLUMNS = (
    "ANI", "Align_fraction_ref", "Align_fraction_query", "Ref_name", "Query_name"
)

def normalize_genome_id(header: str) -> str:
    header = header.strip()
    if "|" in header:
        return header.split("|")[-1].split()[0]
    return header.split()[0]

def read_fasta_ids(fasta: Path) -> set[str]:
    """Vraca identifikatore iz FASTA zaglavlja (prvi token nakon '>')."""
    ids = set()
    count = 0
    with fasta.open("r", encoding="utf-8-sig") as stream:
        for line in stream:
            if not line.startswith(">"):
                continue
            if not line[1:].strip():
                raise ValueError("Prazno FASTA zaglavlje.")
            accession = normalize_genome_id(line[1:])
            if accession in ids:
                raise ValueError(f"Duplicirani FASTA identifikator: {accession}")
            ids.add(accession)
            count += 1
    if count < 2:
        raise ValueError("FASTA mora sadrzavati barem dva razlicita genoma.")
    return ids


def percentage(raw: str, name: str, line_no: int) -> str:
    """Provjerava broj u postocima, ali ne mijenja njegovu preciznost."""
    try:
        number = float(raw)
    except (TypeError, ValueError):
        raise ValueError(f"Redak {line_no}: neispravan {name}: {raw!r}") from None
    if not math.isfinite(number) or not 0 <= number <= 100:
        raise ValueError(f"Redak {line_no}: {name} mora biti u rasponu 0–100: {raw!r}")
    return raw


def normalize_skani_output(raw_path: Path, output_path: Path, fasta_ids: set[str]) -> int:
    """Pretvara skani edge-list u dogovorenih pet stupaca.

    Ako se alfabetni poredak parova mijenja, AF vrijednosti mijenjaju mjesta.
    Neuspjesne usporedbe se ne dodaju, ne upisuju se nule.
    """
    pairs: dict[tuple[str, str], tuple[str, str, str]] = {}
    with raw_path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames:
            raise ValueError("Skani nije vratio TSV zaglavlje.")
        missing = set(REQUIRED_COLUMNS) - set(reader.fieldnames)
        if missing:
            raise ValueError("Neocekivani skani izlaz: nedostaju stupci " + ", ".join(sorted(missing)))
        for line_no, row in enumerate(reader, start=2):
            if None in row:
                raise ValueError(f"Redak {line_no}: visak TSV polja.")
            if not (row["Ref_name"] or "").strip() or not (row["Query_name"] or "").strip():
                raise ValueError(f"Redak {line_no}: nedostaje identifikator genoma.")

            first = normalize_genome_id(row["Ref_name"])
            second = normalize_genome_id(row["Query_name"])
            if first not in fasta_ids or second not in fasta_ids:
                raise ValueError(
                    f"Redak {line_no}: skani identifikatori nisu u FASTA datoteci: {first}, {second}"
                )
            if first == second:
                continue  # Ne pohranjujemo samousporedbe.
            ani = percentage(row["ANI"], "ANI", line_no)
            af_ref = percentage(row["Align_fraction_ref"], "AF_ref", line_no)
            af_query = percentage(row["Align_fraction_query"], "AF_query", line_no)
            if first > second:
                first, second = second, first
                af_ref, af_query = af_query, af_ref
            pair = (first, second)
            if pair in pairs:
                raise ValueError(f"Duplicirana usporedba genoma: {first}, {second}")
            pairs[pair] = (ani, af_ref, af_query)

    # Pisanje u privremenu datoteku sprecava ostecivanje postojece tablice.
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(
        mode="w", encoding="utf-8", newline="", prefix=".ani_", suffix=".tsv",
        dir=output_path.parent, delete=False,
    ) as stream:
        temp_output = Path(stream.name)
        try:
            writer = csv.writer(stream, delimiter="\t", lineterminator="\n")
            writer.writerow(OUT_COLUMNS)
            for (first, second), (ani, af1, af2) in sorted(pairs.items()):
                writer.writerow((first, second, ani, af1, af2))
        except Exception:
            temp_output.unlink(missing_ok=True)
            raise
    temp_output.replace(output_path)
    return len(pairs)


def main() -> int:
    parser = argparse.ArgumentParser(description="Izracun svih dostupnih ANI/AF parova za jedan referentni skup.")
    parser.add_argument("--fasta", required=True, type=Path, help="Multi-FASTA s referentnim genomima.")
    parser.add_argument("--dataset", required=True, choices=DATASETS, help="Referentni skup, ne oznaka uzorka c1/c4.")
    parser.add_argument("--force", action="store_true", help="Dopusti zamjenu postojece ANI tablice.")
    args = parser.parse_args()

    fasta = args.fasta.resolve()
    if not fasta.is_file():
        parser.error(f"FASTA ne postoji: {fasta}")
    output = PROJECT_ROOT / "data" / "ani" / args.dataset / "ani_pairs.tsv"
    if output.exists() and not args.force:
        parser.error(f"Tablica vec postoji: {output}. Za zamjenu koristi --force.")

    try:
        fasta_ids = read_fasta_ids(fasta)
        try:
            version = subprocess.run(["skani", "--version"], capture_output=True, text=True, check=False)
        except FileNotFoundError:
            parser.error("skani nije pronaden u PATH-u. Instaliraj skani prije pokretanja.")
        print("skani verzija:", (version.stdout or version.stderr).strip() or "nedostupna")
        print(f"FASTA: {fasta} ({len(fasta_ids)} genoma)")
        with tempfile.TemporaryDirectory(prefix="skani_ani_") as temporary:
            raw_path = Path(temporary) / "skani_raw.tsv"
            command = ["skani", "triangle", "-i", "-E", str(fasta), "-o", str(raw_path)]
            print("Naredba:", " ".join(command))
            subprocess.run(command, check=True)
            if not raw_path.is_file():
                raise RuntimeError("Skani nije stvorio ocekivanu izlaznu datoteku.")
            count = normalize_skani_output(raw_path, output, fasta_ids)
        print(f"Spremljeno: {output}")
        print(f"Uspjesni parovi: {count} / {len(fasta_ids) * (len(fasta_ids) - 1) // 2}")
    except (ValueError, OSError, RuntimeError, subprocess.CalledProcessError) as error:
        print(f"Greska: {error}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
