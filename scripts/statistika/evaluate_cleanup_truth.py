import os
import re

# nije pripremljeno za koristenje u mainu nad bilo kojim setovima

#postoji na dnu zakomentiran main ako treba pokrenuti kao zasebnu datoteku

# ------------------------------------------------------------
# POSTAVKE
# ------------------------------------------------------------

# BASE_DIR = "results/mentor_sastanak/bucket5000_5b15s_add_c4_sve/statistika"

# GENOME_TRUTH_PATH = f"{BASE_DIR}/summaries/genome_truth_summary.txt"
# CLEANUP_SUMMARY_PATH = f"{BASE_DIR}/summaries/cleanup_simple_summary.txt"
# OUTPUT_PATH = f"{BASE_DIR}/dodatno/cleanup_truth_evaluation.txt"

# ------------------------------------------------------------
# PARSIRANJE GENOME TRUTH DATOTEKE
# ------------------------------------------------------------

def parse_genome_truth(path):
    """
    Čita genome_truth_summary.txt i vraća:
    {
        genome_id: "true" ili "false"
    }

    true  = genom ima očitanja u simulatoru
    false = genom nema očitanja u simulatoru
    """

    truth = {}

    current_section = None

    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()

            if line.startswith("PRAVI GENOMI"):
                current_section = "true"
                continue

            if line.startswith("LAZNI GENOMI"):
                current_section = "false"
                continue

            if current_section not in ("true", "false"):
                continue

            # očekivani format:
            # NC_002662.1: 1596 ocitanja
            match = re.match(r"^([A-Za-z0-9_.]+):\s+(\d+)\s+ocitanja", line)

            if match:
                genome_id = match.group(1)
                truth[genome_id] = current_section

    return truth


# ------------------------------------------------------------
# PARSIRANJE CLEANUP SUMMARY DATOTEKE
# ------------------------------------------------------------

def parse_cleanup_suspicious_genomes(path):
    """
    Čita cleanup_summary.txt i vraća skup genoma koje je cleanup označio
    kao sumnjive.

    Gleda sekciju:
    SUMNJIVI GENOMI ZA CLEANUP
    """

    suspicious = set()

    in_suspicious_section = False

    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()

            if line.startswith("SUMNJIVI GENOMI ZA CLEANUP"):
                in_suspicious_section = True
                continue

            if in_suspicious_section:
                if line.startswith("DETALJI SUMNJIVOSTI PO GENOMU"):
                    break

                if not line:
                    continue

                if set(line) == {"-"}:
                    continue

                if line.startswith("Nema genoma"):
                    continue

                suspicious.add(line)

    return suspicious


# ------------------------------------------------------------
# EVALUACIJA
# ------------------------------------------------------------

def evaluate_cleanup_against_truth(truth, suspicious):
    """
    Uspoređuje cleanup oznake sa stvarnim stanjem iz simulatora.

    Interpretacija:
    - cleanup označi lažni genom  -> true positive
    - cleanup označi pravi genom  -> false positive
    - cleanup ne označi lažni     -> false negative
    - cleanup ne označi pravi     -> true negative
    """

    all_genomes = set(truth.keys())

    true_genomes = {
        genome_id
        for genome_id, label in truth.items()
        if label == "true"
    }

    false_genomes = {
        genome_id
        for genome_id, label in truth.items()
        if label == "false"
    }

    true_positives = suspicious & false_genomes
    false_positives = suspicious & true_genomes
    false_negatives = false_genomes - suspicious
    true_negatives = true_genomes - suspicious

    unknown_suspicious = suspicious - all_genomes

    precision = (
        len(true_positives) / len(suspicious)
        if len(suspicious) > 0
        else 0.0
    )

    recall = (
        len(true_positives) / len(false_genomes)
        if len(false_genomes) > 0
        else 0.0
    )

    f1 = (
        2 * precision * recall / (precision + recall)
        if precision + recall > 0
        else 0.0
    )

    return {
        "total_genomes": len(all_genomes),
        "true_genomes_count": len(true_genomes),
        "false_genomes_count": len(false_genomes),

        "suspicious_count": len(suspicious),

        "true_positives": sorted(true_positives),
        "false_positives": sorted(false_positives),
        "false_negatives": sorted(false_negatives),
        "true_negatives": sorted(true_negatives),
        "unknown_suspicious": sorted(unknown_suspicious),

        "precision": precision,
        "recall": recall,
        "f1": f1,
    }


# ------------------------------------------------------------
# ISPIS
# ------------------------------------------------------------

def write_evaluation_report(result, output_path):
    os.makedirs(os.path.dirname(output_path), exist_ok=True)

    with open(output_path, "w", encoding="utf-8") as f:
        f.write("EVALUACIJA CLEANUP FAZE PREMA GENOME TRUTH DATOTECI\n")
        f.write("=" * 80 + "\n\n")

        f.write("DEFINICIJE\n")
        f.write("-" * 80 + "\n")
        f.write("Pravi genom  = genom ima ocitanja u simulatoru.\n")
        f.write("Lazni genom  = genom nema ocitanja u simulatoru.\n")
        f.write("Cleanup pozitivan = genom je oznacen kao sumnjiv za cleanup.\n\n")

        f.write("SAZETAK\n")
        f.write("-" * 80 + "\n")
        f.write(f"Ukupan broj genoma                : {result['total_genomes']}\n")
        f.write(f"Broj pravih genoma                : {result['true_genomes_count']}\n")
        f.write(f"Broj laznih genoma                : {result['false_genomes_count']}\n")
        f.write(f"Broj genoma oznacenih za cleanup  : {result['suspicious_count']}\n\n")

        f.write("MATRICA ODLUKE\n")
        f.write("-" * 80 + "\n")
        f.write(
            "True positive  - cleanup je oznacio lazni genom       : "
            f"{len(result['true_positives'])}\n"
        )
        f.write(
            "False positive - cleanup je oznacio pravi genom       : "
            f"{len(result['false_positives'])}\n"
        )
        f.write(
            "False negative - lazni genom nije oznacen za cleanup  : "
            f"{len(result['false_negatives'])}\n"
        )
        f.write(
            "True negative  - pravi genom nije oznacen             : "
            f"{len(result['true_negatives'])}\n\n"
        )

        f.write("MJERE\n")
        f.write("-" * 80 + "\n")
        f.write(f"Precision : {result['precision']:.4f}\n")
        f.write(f"Recall    : {result['recall']:.4f}\n")
        f.write(f"F1        : {result['f1']:.4f}\n\n")

        f.write("OBJASNJENJE MJERA\n")
        f.write("-" * 80 + "\n")
        f.write("Precision govori koliki udio oznacenih genoma su stvarno lazni.\n")
        f.write("Recall govori koliki udio svih laznih genoma je cleanup uspio uhvatiti.\n")
        f.write("F1 je zajednicka mjera precisiona i recalla.\n\n")

        f.write("TOCNO PREPOZNATI LAZNI GENOMI\n")
        f.write("-" * 80 + "\n")
        if result["true_positives"]:
            for genome_id in result["true_positives"]:
                f.write(f"{genome_id}: prepoznat kao LAZAN / stvarno je LAZAN\n")
        else:
            f.write("Nema tocno prepoznatih laznih genoma.\n")
        f.write("\n")

        f.write("KRIVO PREPOZNATI GENOMI - FALSE POSITIVE\n")
        f.write("-" * 80 + "\n")
        f.write("Ovo su genomi koje je cleanup oznacio kao sumnjive, ali su u simulatoru pravi.\n\n")

        if result["false_positives"]:
            for genome_id in result["false_positives"]:
                f.write(
                    f"{genome_id}: cleanup ga je prepoznao kao LAZAN/SUMNJIV, "
                    f"ali trebao bi biti PRAVI\n"
                )
        else:
            f.write("Nema false positive gresaka. Cleanup nije oznacio nijedan pravi genom.\n")
        f.write("\n")

        f.write("PROMASENI LAZNI GENOMI - FALSE NEGATIVE\n")
        f.write("-" * 80 + "\n")
        f.write("Ovo su lažni genomi koje cleanup nije oznacio kao sumnjive.\n\n")

        if result["false_negatives"]:
            for genome_id in result["false_negatives"]:
                f.write(
                    f"{genome_id}: cleanup ga je ostavio kao NESUMNJIV, "
                    f"ali trebao bi biti LAZAN\n"
                )
        else:
            f.write("Nema false negative gresaka. Cleanup je oznacio sve lazne genome.\n")
        f.write("\n")

        if result["unknown_suspicious"]:
            f.write("GENOMI OZNAČENI U CLEANUPU, ALI IH NEMA U TRUTH DATOTECI\n")
            f.write("-" * 80 + "\n")
            for genome_id in result["unknown_suspicious"]:
                f.write(f"{genome_id}\n")
            f.write("\n")

    return output_path


# ------------------------------------------------------------
# MAIN
# ------------------------------------------------------------

# def main():
#     truth = parse_genome_truth(GENOME_TRUTH_PATH)
#     suspicious = parse_cleanup_suspicious_genomes(CLEANUP_SUMMARY_PATH)

#     result = evaluate_cleanup_against_truth(
#         truth=truth,
#         suspicious=suspicious
#     )

#     output_path = write_evaluation_report(
#         result=result,
#         output_path=OUTPUT_PATH
#     )

#     print(f"Evaluacija spremljena u: {output_path}")
#     print()
#     print("Kratki sazetak:")
#     print(f"  oznaceni za cleanup : {result['suspicious_count']}")
#     print(f"  true positives      : {len(result['true_positives'])}")
#     print(f"  false positives     : {len(result['false_positives'])}")
#     print(f"  false negatives     : {len(result['false_negatives'])}")
#     print(f"  precision           : {result['precision']:.4f}")
#     print(f"  recall              : {result['recall']:.4f}")
#     print(f"  F1                  : {result['f1']:.4f}")


# if __name__ == "__main__":
#     main()