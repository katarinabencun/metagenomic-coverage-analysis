import os
import shutil

def organize_cleanup_comparison_images(
    sim_counts,
    genome_lengths,
    cleanup_summary,
    comparison_dir
):
    """
    Kopira postojeće usporedne grafove u podmape prema rezultatu
    evaluacije cleanup faze.
    """
    truth = {
        genome_id: (
            "true"
            if sim_counts.get(genome_id, 0) > 0
            else "false"
        )
        for genome_id in genome_lengths
    }

    suspicious = set(
        cleanup_summary.get("suspicious_genomes", [])
    )

    result = evaluate_cleanup_against_truth(
        truth=truth,
        suspicious=suspicious
    )

    cleanup_dir = os.path.join(comparison_dir, "cleanup")

    categories = {
        "true_positive": result["true_positives"],
        "false_positive": result["false_positives"],
        "false_negative": result["false_negatives"],
        "true_negative": result["true_negatives"],
    }

    for category, genome_ids in categories.items():
        category_dir = os.path.join(cleanup_dir, category)
        os.makedirs(category_dir, exist_ok=True)

        # Uklanjanje starih kopija ako se eksperiment ponovno pokrene
        for filename in os.listdir(category_dir):
            if filename.endswith("_usporedba.png"):
                os.remove(os.path.join(category_dir, filename))

        for genome_id in genome_ids:
            filename = f"{genome_id}_usporedba.png"
            source_path = os.path.join(comparison_dir, filename)
            destination_path = os.path.join(category_dir, filename)

            if os.path.isfile(source_path):
                shutil.copy2(source_path, destination_path)
            else:
                print(
                    "Upozorenje: usporedni graf nije pronađen za "
                    f"{genome_id}: {source_path}"
                )

    return cleanup_dir

def evaluate_cleanup_against_truth(truth, suspicious):
    """
    Uspoređuje genome označene kao sumnjive sa stvarnim stanjem
    poznatim iz simulatora.

    Interpretacija:
    - cleanup označi lažni genom  -> true positive
    - cleanup označi pravi genom -> false positive
    - cleanup ne označi lažni    -> false negative
    - cleanup ne označi pravi    -> true negative
    """
    all_genomes = set(truth)

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
        if suspicious
        else 0.0
    )

    recall = (
        len(true_positives) / len(false_genomes)
        if false_genomes
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


def write_evaluation_report(result, output_path):
    """
    Sprema rezultate evaluacije cleanup faze u tekstualnu datoteku.
    """
    output_dir = os.path.dirname(output_path)

    if output_dir:
        os.makedirs(output_dir, exist_ok=True)

    with open(output_path, "w", encoding="utf-8") as f:
        f.write("EVALUACIJA CLEANUP FAZE PREMA SIMULATORU\n")
        f.write("=" * 80 + "\n\n")

        f.write("DEFINICIJE\n")
        f.write("-" * 80 + "\n")
        f.write("Pravi genom = genom ima ocitanja u simulatoru.\n")
        f.write("Lazni genom = genom nema ocitanja u simulatoru.\n")
        f.write(
            "Cleanup pozitivan = genom je oznacen kao sumnjiv "
            "za cleanup.\n\n"
        )

        f.write("SAZETAK\n")
        f.write("-" * 80 + "\n")
        f.write(
            f"Ukupan broj genoma                : "
            f"{result['total_genomes']}\n"
        )
        f.write(
            f"Broj pravih genoma                : "
            f"{result['true_genomes_count']}\n"
        )
        f.write(
            f"Broj laznih genoma                : "
            f"{result['false_genomes_count']}\n"
        )
        f.write(
            f"Broj genoma oznacenih za cleanup  : "
            f"{result['suspicious_count']}\n\n"
        )

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
        f.write(
            "Precision govori koliki su udio oznacenih genoma "
            "stvarno lazni genomi.\n"
        )
        f.write(
            "Recall govori koliki je udio svih laznih genoma "
            "cleanup uspio prepoznati.\n"
        )
        f.write(
            "F1 je zajednicka mjera precisiona i recalla.\n\n"
        )

        f.write("TOCNO PREPOZNATI LAZNI GENOMI\n")
        f.write("-" * 80 + "\n")

        if result["true_positives"]:
            for genome_id in result["true_positives"]:
                f.write(
                    f"{genome_id}: oznacen kao sumnjiv i stvarno je lazan\n"
                )
        else:
            f.write("Nema tocno prepoznatih laznih genoma.\n")

        f.write("\n")

        f.write("KRIVO OZNACENI PRAVI GENOMI - FALSE POSITIVE\n")
        f.write("-" * 80 + "\n")

        if result["false_positives"]:
            for genome_id in result["false_positives"]:
                f.write(
                    f"{genome_id}: oznacen kao sumnjiv, "
                    f"ali je u simulatoru pravi genom\n"
                )
        else:
            f.write(
                "Nema false positive pogresaka. "
                "Cleanup nije oznacio nijedan pravi genom.\n"
            )

        f.write("\n")

        f.write("NEPREPOZNATI LAZNI GENOMI - FALSE NEGATIVE\n")
        f.write("-" * 80 + "\n")

        if result["false_negatives"]:
            for genome_id in result["false_negatives"]:
                f.write(
                    f"{genome_id}: nije oznacen kao sumnjiv, "
                    f"ali je u simulatoru lazan genom\n"
                )
        else:
            f.write(
                "Nema false negative pogresaka. "
                "Cleanup je oznacio sve lazne genome.\n"
            )

        f.write("\n")

        f.write("TOCNO NEOZNACENI PRAVI GENOMI - TRUE NEGATIVE\n")
        f.write("-" * 80 + "\n")

        if result["true_negatives"]:
            for genome_id in result["true_negatives"]:
                f.write(
                    f"{genome_id}: nije oznacen kao sumnjiv "
                    f"i stvarno je pravi genom\n"
                )
        else:
            f.write("Nema true negative genoma.\n")

        f.write("\n")

        if result["unknown_suspicious"]:
            f.write(
                "SUMNJIVI GENOMI KOJI NISU PRONADJENI "
                "U REFERENTNOJ BAZI\n"
            )
            f.write("-" * 80 + "\n")

            for genome_id in result["unknown_suspicious"]:
                f.write(f"{genome_id}\n")

            f.write("\n")

    return output_path


def write_cleanup_truth_evaluation(
    sim_counts,
    genome_lengths,
    cleanup_summary,
    output_dir,
    filename="cleanup_truth_evaluation.txt"
):
    """
    Evaluira genome označene kao sumnjive usporedbom s poznatim
    stanjem iz simulatora i sprema izvještaj.
    """
    truth = {
        genome_id: (
            "true"
            if sim_counts.get(genome_id, 0) > 0
            else "false"
        )
        for genome_id in genome_lengths
    }

    suspicious = set(
        cleanup_summary.get("suspicious_genomes", [])
    )

    result = evaluate_cleanup_against_truth(
        truth=truth,
        suspicious=suspicious
    )

    output_path = os.path.join(output_dir, filename)

    return write_evaluation_report(
        result=result,
        output_path=output_path
    )