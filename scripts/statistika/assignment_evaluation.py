import os
import gzip
from collections import defaultdict


# =============================================================================
# FASTQ / SIMULATOR TRUTH PARSIRANJE
# =============================================================================

def open_text_maybe_gzip(file_path):
    if file_path.endswith(".gz"):
        return gzip.open(file_path, "rt")
    return open(file_path, "r")


def parse_fastq_truth_with_read_id(file_path):
    """
    Parsira FASTQ simulator datoteku i vraća stvarno podrijetlo svakog očitanja.

    Izlaz:
    {
        read_id: {
            "true_genome": genome_id,
            "true_start": start,
            "true_end": end
        }
    }

    Pretpostavka:
    - read_id je prvi token u FASTQ headeru, bez početnog znaka '@'
    - simulator informacije su u drugom tokenu, kao i u postojećem parse_fastq kodu
    """

    truth = {}

    with open_text_maybe_gzip(file_path) as f:
        while True:
            header = f.readline().strip()
            seq = f.readline()
            plus = f.readline()
            qual = f.readline()

            if not header:
                break

            if not header.startswith("@"):
                continue

            parts = header.split()

            if len(parts) < 2:
                continue

            try:
                read_id = parts[0][1:]

                genome_part = parts[1].split("|")[-1]
                true_genome = genome_part.split(",")[0]

                position_range = genome_part.split(",")[2]
                true_start = int(position_range.split("-")[0])
                true_end = int(position_range.split("-")[1])

            except (IndexError, ValueError):
                continue

            truth[read_id] = {
                "true_genome": true_genome,
                "true_start": true_start,
                "true_end": true_end,
            }

    return truth


# =============================================================================
# POMOĆNE FUNKCIJE ZA DODJELE
# =============================================================================

def get_assignment_field(assignments, read_id, field_name, default="NA"):
    aln = assignments.get(read_id)

    if aln is None:
        return default

    return aln.get(field_name, default)


def get_assigned_genome(assignments, read_id):
    aln = assignments.get(read_id)

    if aln is None:
        return None

    return aln.get("genome_id")


def is_correct_assignment(true_genome, assigned_genome):
    return assigned_genome is not None and assigned_genome == true_genome


def classify_read_status(initial_correct, final_correct, initial_genome, final_genome):
    """
    Status opisuje što se dogodilo s očitanjem između inicijalne i finalne dodjele.
    """

    if initial_genome is None and final_genome is None:
        return "unassigned_both"

    if initial_genome is None and final_genome is not None:
        if final_correct:
            return "assigned_later_correct"
        return "assigned_later_wrong"

    if initial_genome is not None and final_genome is None:
        if initial_correct:
            return "lost_correct"
        return "lost_wrong"

    if initial_correct and final_correct:
        return "kept_correct"

    if (not initial_correct) and final_correct:
        return "fixed"

    if initial_correct and (not final_correct):
        return "worsened"

    return "kept_wrong"


def percent(part, total):
    if total == 0:
        return 0.0

    return 100.0 * part / total


# =============================================================================
# 1. TABLICA PO OČITANJU
# =============================================================================

def write_read_assignment_table(
    fastq_path,
    initial_assignments,
    final_assignments,
    output_dir,
    filename="read_assignment_table.tsv"
):
    """
    Sprema tablicu u kojoj se za svako očitanje vidi:
    - iz kojeg genoma stvarno dolazi prema simulatoru
    - kojem genomu je dodijeljeno inicijalno
    - kojem genomu je dodijeljeno finalno
    - je li inicijalna/finalna dodjela točna
    - status promjene
    """

    os.makedirs(output_dir, exist_ok=True)

    truth = parse_fastq_truth_with_read_id(fastq_path)

    output_path = os.path.join(output_dir, filename)

    with open(output_path, "w", encoding="utf-8") as f:
        f.write(
            "read_id\t"
            "true_genome\ttrue_start\ttrue_end\t"
            "initial_genome\tinitial_start\tinitial_end\tinitial_mapq\tinitial_tp\t"
            "final_genome\tfinal_start\tfinal_end\tfinal_mapq\tfinal_tp\t"
            "initial_correct\tfinal_correct\tstatus\n"
        )

        for read_id in sorted(truth):
            true_info = truth[read_id]
            true_genome = true_info["true_genome"]

            initial_genome = get_assigned_genome(initial_assignments, read_id)
            final_genome = get_assigned_genome(final_assignments, read_id)

            initial_correct = is_correct_assignment(true_genome, initial_genome)
            final_correct = is_correct_assignment(true_genome, final_genome)

            status = classify_read_status(
                initial_correct=initial_correct,
                final_correct=final_correct,
                initial_genome=initial_genome,
                final_genome=final_genome
            )

            f.write(
                f"{read_id}\t"
                f"{true_genome}\t"
                f"{true_info['true_start']}\t"
                f"{true_info['true_end']}\t"

                f"{initial_genome if initial_genome is not None else 'NA'}\t"
                f"{get_assignment_field(initial_assignments, read_id, 'start')}\t"
                f"{get_assignment_field(initial_assignments, read_id, 'end')}\t"
                f"{get_assignment_field(initial_assignments, read_id, 'mapq')}\t"
                f"{get_assignment_field(initial_assignments, read_id, 'tp')}\t"

                f"{final_genome if final_genome is not None else 'NA'}\t"
                f"{get_assignment_field(final_assignments, read_id, 'start')}\t"
                f"{get_assignment_field(final_assignments, read_id, 'end')}\t"
                f"{get_assignment_field(final_assignments, read_id, 'mapq')}\t"
                f"{get_assignment_field(final_assignments, read_id, 'tp')}\t"

                f"{int(initial_correct)}\t"
                f"{int(final_correct)}\t"
                f"{status}\n"
            )

    return output_path


# =============================================================================
# 2. CONFUSION MATRIX TRUE_GENOME -> ASSIGNED_GENOME
# =============================================================================

def build_assignment_confusion_matrices(
    truth,
    initial_assignments,
    final_assignments
):
    """
    Gradi dvije matrice:
    - true_genome -> initial_assigned_genome
    - true_genome -> final_assigned_genome

    Ako očitanje nije dodijeljeno, assigned_genome je "UNASSIGNED".
    """

    initial_matrix = defaultdict(int)
    final_matrix = defaultdict(int)

    for read_id, true_info in truth.items():
        true_genome = true_info["true_genome"]

        initial_genome = get_assigned_genome(initial_assignments, read_id)
        final_genome = get_assigned_genome(final_assignments, read_id)

        if initial_genome is None:
            initial_genome = "UNASSIGNED"

        if final_genome is None:
            final_genome = "UNASSIGNED"

        initial_matrix[(true_genome, initial_genome)] += 1
        final_matrix[(true_genome, final_genome)] += 1

    return initial_matrix, final_matrix


def write_assignment_confusion_matrix(
    fastq_path,
    initial_assignments,
    final_assignments,
    output_dir,
    filename="assignment_confusion_matrix.tsv"
):
    """
    Sprema TSV tablicu:
        true_genome assigned_genome initial_count final_count

    Korisno za analizu toga koji se genomi međusobno miješaju.
    """

    os.makedirs(output_dir, exist_ok=True)

    truth = parse_fastq_truth_with_read_id(fastq_path)

    initial_matrix, final_matrix = build_assignment_confusion_matrices(
        truth=truth,
        initial_assignments=initial_assignments,
        final_assignments=final_assignments
    )

    all_pairs = sorted(set(initial_matrix) | set(final_matrix))

    output_path = os.path.join(output_dir, filename)

    with open(output_path, "w", encoding="utf-8") as f:
        f.write("true_genome\tassigned_genome\tinitial_count\tfinal_count\tchange\n")

        for true_genome, assigned_genome in all_pairs:
            initial_count = initial_matrix.get((true_genome, assigned_genome), 0)
            final_count = final_matrix.get((true_genome, assigned_genome), 0)

            f.write(
                f"{true_genome}\t"
                f"{assigned_genome}\t"
                f"{initial_count}\t"
                f"{final_count}\t"
                f"{final_count - initial_count}\n"
            )

    return output_path


# =============================================================================
# 3. STATISTIČKI SAŽETAK EVALUACIJE
# =============================================================================

def evaluate_assignments_against_truth(
    truth,
    initial_assignments,
    final_assignments,
    genome_lengths
):
    """
    Evaluira inicijalnu i finalnu dodjelu prema simulator truth podacima.

    Mjere:
    - ukupna točnost po očitanju
    - broj popravljenih i pogoršanih očitanja
    - recall po pravom genomu
    - broj lažno dodijeljenih očitanja po lažnim genomima
    """

    total_truth_reads = len(truth)

    initial_assigned_count = 0
    final_assigned_count = 0

    initial_correct_count = 0
    final_correct_count = 0

    status_counts = defaultdict(int)

    true_counts_by_genome = defaultdict(int)

    initial_assigned_by_genome = defaultdict(int)
    final_assigned_by_genome = defaultdict(int)

    initial_correct_by_true_genome = defaultdict(int)
    final_correct_by_true_genome = defaultdict(int)

    initial_wrong_to_genome = defaultdict(int)
    final_wrong_to_genome = defaultdict(int)

    for read_id, true_info in truth.items():
        true_genome = true_info["true_genome"]
        true_counts_by_genome[true_genome] += 1

        initial_genome = get_assigned_genome(initial_assignments, read_id)
        final_genome = get_assigned_genome(final_assignments, read_id)

        if initial_genome is not None:
            initial_assigned_count += 1
            initial_assigned_by_genome[initial_genome] += 1

        if final_genome is not None:
            final_assigned_count += 1
            final_assigned_by_genome[final_genome] += 1

        initial_correct = is_correct_assignment(true_genome, initial_genome)
        final_correct = is_correct_assignment(true_genome, final_genome)

        if initial_correct:
            initial_correct_count += 1
            initial_correct_by_true_genome[true_genome] += 1
        elif initial_genome is not None:
            initial_wrong_to_genome[initial_genome] += 1

        if final_correct:
            final_correct_count += 1
            final_correct_by_true_genome[true_genome] += 1
        elif final_genome is not None:
            final_wrong_to_genome[final_genome] += 1

        status = classify_read_status(
            initial_correct=initial_correct,
            final_correct=final_correct,
            initial_genome=initial_genome,
            final_genome=final_genome
        )

        status_counts[status] += 1

    true_genomes = {
        genome_id
        for genome_id, count in true_counts_by_genome.items()
        if count > 0
    }

    false_genomes = {
        genome_id
        for genome_id in genome_lengths
        if genome_id not in true_genomes
    }

    per_genome = {}

    for genome_id in sorted(genome_lengths):
        simulator_count = true_counts_by_genome.get(genome_id, 0)

        initial_assigned = initial_assigned_by_genome.get(genome_id, 0)
        final_assigned = final_assigned_by_genome.get(genome_id, 0)

        initial_correct = initial_correct_by_true_genome.get(genome_id, 0)
        final_correct = final_correct_by_true_genome.get(genome_id, 0)

        initial_wrong_received = initial_wrong_to_genome.get(genome_id, 0)
        final_wrong_received = final_wrong_to_genome.get(genome_id, 0)

        if simulator_count > 0:
            initial_recall = initial_correct / simulator_count
            final_recall = final_correct / simulator_count
        else:
            initial_recall = None
            final_recall = None

        if initial_assigned > 0:
            initial_precision = initial_correct / initial_assigned
        else:
            initial_precision = None

        if final_assigned > 0:
            final_precision = final_correct / final_assigned
        else:
            final_precision = None

        per_genome[genome_id] = {
            "simulator_count": simulator_count,

            "initial_assigned": initial_assigned,
            "final_assigned": final_assigned,
            "assigned_change": final_assigned - initial_assigned,

            "initial_correct": initial_correct,
            "final_correct": final_correct,
            "correct_change": final_correct - initial_correct,

            "initial_wrong_received": initial_wrong_received,
            "final_wrong_received": final_wrong_received,
            "wrong_received_change": final_wrong_received - initial_wrong_received,

            "initial_recall": initial_recall,
            "final_recall": final_recall,
            "recall_change": (
                final_recall - initial_recall
                if initial_recall is not None and final_recall is not None
                else None
            ),

            "initial_precision": initial_precision,
            "final_precision": final_precision,
            "precision_change": (
                final_precision - initial_precision
                if initial_precision is not None and final_precision is not None
                else None
            ),

            "is_true_genome": genome_id in true_genomes,
            "is_false_genome": genome_id in false_genomes,
        }

    result = {
        "total_truth_reads": total_truth_reads,

        "initial_assigned_count": initial_assigned_count,
        "final_assigned_count": final_assigned_count,

        "initial_unassigned_count": total_truth_reads - initial_assigned_count,
        "final_unassigned_count": total_truth_reads - final_assigned_count,

        "initial_correct_count": initial_correct_count,
        "final_correct_count": final_correct_count,

        "initial_accuracy": (
            initial_correct_count / total_truth_reads
            if total_truth_reads > 0
            else 0.0
        ),
        "final_accuracy": (
            final_correct_count / total_truth_reads
            if total_truth_reads > 0
            else 0.0
        ),

        "accuracy_change": (
            (final_correct_count / total_truth_reads)
            - (initial_correct_count / total_truth_reads)
            if total_truth_reads > 0
            else 0.0
        ),

        "status_counts": dict(status_counts),

        "true_genomes": sorted(true_genomes),
        "false_genomes": sorted(false_genomes),

        "per_genome": per_genome,
    }

    return result


def format_optional_float(value, digits=4):
    if value is None:
        return "NA"

    return f"{value:.{digits}f}"


def write_assignment_evaluation_summary(
    fastq_path,
    initial_assignments,
    final_assignments,
    genome_lengths,
    output_dir,
    filename="assignment_evaluation_summary.txt"
):
    """
    Sprema čitljiv evaluacijski izvještaj:
    - ukupna točnost
    - statusi očitanja
    - per-genome evaluacija
    - sažetak za lažne genome
    """

    os.makedirs(output_dir, exist_ok=True)

    truth = parse_fastq_truth_with_read_id(fastq_path)

    result = evaluate_assignments_against_truth(
        truth=truth,
        initial_assignments=initial_assignments,
        final_assignments=final_assignments,
        genome_lengths=genome_lengths
    )

    output_path = os.path.join(output_dir, filename)

    with open(output_path, "w", encoding="utf-8") as f:
        f.write("EVALUACIJA DODJELE OCITANJA PREMA SIMULATORU\n")
        f.write("=" * 80 + "\n\n")

        f.write("DEFINICIJA\n")
        f.write("-" * 80 + "\n")
        f.write("Simulator / FASTQ daje stvarni genom iz kojeg je ocitanje generirano.\n")
        f.write("Inicijalna dodjela i finalna preraspodjela usporeduju se s tim stvarnim genomom.\n")
        f.write("Dodjela je tocna ako je assigned_genome jednak true_genome.\n\n")

        f.write("UKUPNI SAZETAK\n")
        f.write("-" * 80 + "\n")
        f.write(f"Ukupan broj ocitanja u simulatoru       : {result['total_truth_reads']}\n")
        f.write(f"Broj ocitanja s inicijalnom dodjelom    : {result['initial_assigned_count']}\n")
        f.write(f"Broj ocitanja s finalnom dodjelom       : {result['final_assigned_count']}\n")
        f.write(f"Nedodijeljeno inicijalno                : {result['initial_unassigned_count']}\n")
        f.write(f"Nedodijeljeno finalno                   : {result['final_unassigned_count']}\n\n")

        f.write("TOCNOST PO OCITANJIMA\n")
        f.write("-" * 80 + "\n")
        f.write(f"Initial correct                         : {result['initial_correct_count']}\n")
        f.write(f"Final correct                           : {result['final_correct_count']}\n")
        f.write(f"Initial accuracy                        : {100 * result['initial_accuracy']:.2f}%\n")
        f.write(f"Final accuracy                          : {100 * result['final_accuracy']:.2f}%\n")
        f.write(f"Promjena accuracyja                     : {100 * result['accuracy_change']:+.2f} postotnih bodova\n\n")

        f.write("STATUS PROMJENE PO OCITANJIMA\n")
        f.write("-" * 80 + "\n")

        status_order = [
            "kept_correct",
            "fixed",
            "worsened",
            "kept_wrong",
            "assigned_later_correct",
            "assigned_later_wrong",
            "lost_correct",
            "lost_wrong",
            "unassigned_both",
        ]

        for status in status_order:
            count = result["status_counts"].get(status, 0)
            f.write(f"{status:28s}: {count}\n")

        f.write("\n")
        f.write("OBJASNJENJE STATUSA\n")
        f.write("-" * 80 + "\n")
        f.write("kept_correct          : initial je bio tocan i final je ostao tocan\n")
        f.write("fixed                 : initial je bio kriv, final je postao tocan\n")
        f.write("worsened              : initial je bio tocan, final je postao kriv\n")
        f.write("kept_wrong            : initial je bio kriv i final je ostao kriv\n")
        f.write("assigned_later_correct: initial nije imao dodjelu, final je tocno dodijeljen\n")
        f.write("assigned_later_wrong  : initial nije imao dodjelu, final je krivo dodijeljen\n")
        f.write("lost_correct          : initial je bio tocan, final je ostao bez dodjele\n")
        f.write("lost_wrong            : initial je bio kriv, final je ostao bez dodjele\n")
        f.write("unassigned_both       : nema dodjele ni inicijalno ni finalno\n\n")

        f.write("EVALUACIJA PO GENOMU\n")
        f.write("=" * 80 + "\n\n")

        for genome_id, g in sorted(result["per_genome"].items()):
            genome_type = "PRAVI" if g["is_true_genome"] else "LAZNI"

            f.write(f"Genom: {genome_id} ({genome_type})\n")
            f.write("-" * 80 + "\n")

            f.write(f"Simulator count              : {g['simulator_count']}\n")
            f.write(f"Initial assigned             : {g['initial_assigned']}\n")
            f.write(f"Final assigned               : {g['final_assigned']}\n")
            f.write(f"Assigned change              : {g['assigned_change']:+d}\n\n")

            f.write(f"Initial correct              : {g['initial_correct']}\n")
            f.write(f"Final correct                : {g['final_correct']}\n")
            f.write(f"Correct change               : {g['correct_change']:+d}\n\n")

            f.write(f"Initial wrong received       : {g['initial_wrong_received']}\n")
            f.write(f"Final wrong received         : {g['final_wrong_received']}\n")
            f.write(f"Wrong received change        : {g['wrong_received_change']:+d}\n\n")

            f.write(f"Initial recall               : {format_optional_float(g['initial_recall'])}\n")
            f.write(f"Final recall                 : {format_optional_float(g['final_recall'])}\n")
            f.write(f"Recall change                : {format_optional_float(g['recall_change'], digits=4)}\n\n")

            f.write(f"Initial precision            : {format_optional_float(g['initial_precision'])}\n")
            f.write(f"Final precision              : {format_optional_float(g['final_precision'])}\n")
            f.write(f"Precision change             : {format_optional_float(g['precision_change'], digits=4)}\n")

            f.write("\n\n")

        f.write("=" * 80 + "\n")
        f.write("SAZETAK ZA LAZNE GENOME\n")
        f.write("=" * 80 + "\n\n")

        false_initial_total = 0
        false_final_total = 0
        false_initial_wrong_total = 0
        false_final_wrong_total = 0

        for genome_id in result["false_genomes"]:
            g = result["per_genome"][genome_id]

            false_initial_total += g["initial_assigned"]
            false_final_total += g["final_assigned"]

            false_initial_wrong_total += g["initial_wrong_received"]
            false_final_wrong_total += g["final_wrong_received"]

        f.write(f"Broj laznih genoma                         : {len(result['false_genomes'])}\n")
        f.write(f"Initial dodjele na lazne genome            : {false_initial_total}\n")
        f.write(f"Final dodjele na lazne genome              : {false_final_total}\n")
        f.write(f"Promjena dodjela na lazne genome           : {false_final_total - false_initial_total:+d}\n\n")

        f.write(f"Initial krivo primljena ocitanja na lazne  : {false_initial_wrong_total}\n")
        f.write(f"Final krivo primljena ocitanja na lazne    : {false_final_wrong_total}\n")
        f.write(f"Promjena krivo primljenih na lazne         : {false_final_wrong_total - false_initial_wrong_total:+d}\n\n")

        f.write("Napomena: Za lazne genome idealno je da finalni broj dodijeljenih ocitanja bude sto blizi nuli.\n")

    return output_path


# =============================================================================
# 4. TSV SAŽETAK PO GENOMU
# =============================================================================

def write_assignment_per_genome_table(
    fastq_path,
    initial_assignments,
    final_assignments,
    genome_lengths,
    output_dir,
    filename="assignment_per_genome_summary.tsv"
):
    """
    Sprema strojno čitljivu TSV tablicu po genomu.
    Ovo je praktično za kasniju obradu u Excelu, Pythonu ili usporedbu eksperimenata.
    """

    os.makedirs(output_dir, exist_ok=True)

    truth = parse_fastq_truth_with_read_id(fastq_path)

    result = evaluate_assignments_against_truth(
        truth=truth,
        initial_assignments=initial_assignments,
        final_assignments=final_assignments,
        genome_lengths=genome_lengths
    )

    output_path = os.path.join(output_dir, filename)

    with open(output_path, "w", encoding="utf-8") as f:
        f.write(
            "genome_id\tgenome_type\t"
            "simulator_count\t"
            "initial_assigned\tfinal_assigned\tassigned_change\t"
            "initial_correct\tfinal_correct\tcorrect_change\t"
            "initial_wrong_received\tfinal_wrong_received\twrong_received_change\t"
            "initial_recall\tfinal_recall\trecall_change\t"
            "initial_precision\tfinal_precision\tprecision_change\n"
        )

        for genome_id, g in sorted(result["per_genome"].items()):
            genome_type = "true" if g["is_true_genome"] else "false"

            f.write(
                f"{genome_id}\t"
                f"{genome_type}\t"

                f"{g['simulator_count']}\t"

                f"{g['initial_assigned']}\t"
                f"{g['final_assigned']}\t"
                f"{g['assigned_change']}\t"

                f"{g['initial_correct']}\t"
                f"{g['final_correct']}\t"
                f"{g['correct_change']}\t"

                f"{g['initial_wrong_received']}\t"
                f"{g['final_wrong_received']}\t"
                f"{g['wrong_received_change']}\t"

                f"{format_optional_float(g['initial_recall'])}\t"
                f"{format_optional_float(g['final_recall'])}\t"
                f"{format_optional_float(g['recall_change'])}\t"

                f"{format_optional_float(g['initial_precision'])}\t"
                f"{format_optional_float(g['final_precision'])}\t"
                f"{format_optional_float(g['precision_change'])}\n"
            )

    return output_path


# =============================================================================
# 5. JEDAN POZIV KOJI GENERIRA SVE EVALUACIJSKE DATOTEKE
# =============================================================================

def write_full_assignment_evaluation(
    fastq_path,
    initial_assignments,
    final_assignments,
    genome_lengths,
    output_dir
):
    """
    Generira sve evaluacijske datoteke za dodjelu očitanja.

    Datoteke:
    - read_assignment_table.tsv
    - assignment_confusion_matrix.tsv
    - assignment_evaluation_summary.txt
    - assignment_per_genome_summary.tsv
    """

    os.makedirs(output_dir, exist_ok=True)

    read_table_path = write_read_assignment_table(
        fastq_path=fastq_path,
        initial_assignments=initial_assignments,
        final_assignments=final_assignments,
        output_dir=output_dir,
        filename="read_assignment_table.tsv"
    )

    confusion_path = write_assignment_confusion_matrix(
        fastq_path=fastq_path,
        initial_assignments=initial_assignments,
        final_assignments=final_assignments,
        output_dir=output_dir,
        filename="assignment_confusion_matrix.tsv"
    )

    summary_path = write_assignment_evaluation_summary(
        fastq_path=fastq_path,
        initial_assignments=initial_assignments,
        final_assignments=final_assignments,
        genome_lengths=genome_lengths,
        output_dir=output_dir,
        filename="assignment_evaluation_summary.txt"
    )

    per_genome_path = write_assignment_per_genome_table(
        fastq_path=fastq_path,
        initial_assignments=initial_assignments,
        final_assignments=final_assignments,
        genome_lengths=genome_lengths,
        output_dir=output_dir,
        filename="assignment_per_genome_summary.tsv"
    )

    return {
        "read_assignment_table": read_table_path,
        "assignment_confusion_matrix": confusion_path,
        "assignment_evaluation_summary": summary_path,
        "assignment_per_genome_summary": per_genome_path,
    }