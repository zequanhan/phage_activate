from Bio import SeqIO, SearchIO
from pathlib import Path
import re
import sys
import os
import subprocess
import argparse

here = Path(__file__).resolve().parent

__version__ = "1.0.0"


def parse_arguments():
    parser = argparse.ArgumentParser(
        prog="prokka.py",
        description="Predict phage draft genomes in metagenomic bins.",
        # Display default values when printing help
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    # Rearrange arguments to proper groups
    # to show up better in the help
    optionalArgs = parser._action_groups.pop()
    optionalArgs.title = "Optional arguments"

    # Required arguments
    requiredArgs = parser.add_argument_group("Required arguments")
    requiredArgs.add_argument(
        "-i",
        "--input-dir",
        required=True,
        type=lambda p: Path(p).resolve(strict=True),
        dest="input_dir",
        default="",
        help="Path to a folder containing metagenomic bins in .fa or .fasta "
             "format",
    )

    requiredArgs.add_argument(
        "-o",
        "--output-dir",
        required=True,
        type=lambda p: Path(p).resolve(),
        dest="output_dir",
        default="",
        help="Location to store results in. It is created if it doesn't exist",
    )

    # Optional args
    optionalArgs.add_argument(
        "-t",
        "--threads",
        dest="threads",
        required=False,
        default=1,
        type=int,
        help="Number of CPU threads to be used by Prokka and hmmscan",
    )

    optionalArgs.add_argument(
        "-m",
        "--models-dir",
        required=False,
        dest="models_dir",
        type=lambda p: Path(p).resolve(),
        default=here / Path("models"),
        help="Path to directory where all models are stored.",
    )

    optionalArgs.add_argument(
        "-f",
        "--files-dir",
        required=False,
        dest="files_dir",
        default=here / Path("files"),
        help="Files directory provided with vHULK",
    )

    optionalArgs.add_argument(
        "--all",
        required=False,
        dest="write_all",
        action="store_true",
        help="Write predictions for all input bins/genomes, even if they "
             "were skipped (size filtered or hmmscan failed)",
    )

    optionalArgs.add_argument(
        "-v",
        "--version",
        action="version",
        version=__version__,
    )

    parser._action_groups.append(optionalArgs)

    return parser.parse_args()


def get_bin_name(fasta_p):
    """
    Strip relevant pre and suffixes from a given fasta name

    Arguments:
        fasta_p: pathlib.Path instance: The path to the fasta file
    Return:
        bin_name: str: The name of the bin as string
    """
    bin_name = fasta_p.name
    fa_suffix = fasta_p.suffix
    prokka_prefix = "prokka_results_"
    # Check if it is prefixed with prokka_results
    if bin_name.startswith(prokka_prefix):
        bin_name = fasta_p.name.replace(prokka_prefix, "")
    # Remove the .fa[a|sta] suffix
    bin_name = bin_name.replace(fa_suffix, "")
    return bin_name


def run_prokka(fasta_in, output_dir, threads):
    """
    Run prokka for a given fasta file.

    Raises CalledProcessError if command doesn't end succesfully

    Arguments:
        fasta_in: pathlib.Path instance: Path to fasta file
        output_dir: pathlib.Path instance: Path to directory where results will
            be stored
        threads: int: Number of cpus/threads for prokka to use
    Return:
        -
    """
    out_prefix = get_bin_name(fasta_in)
    genome_dir = output_dir / Path(out_prefix)
    command_line = (
        "prokka --kingdom Viruses --centre X --compliant "
        "--gcode 11 --cpus {} --force --quiet --prefix prokka_results_{} "
        "--fast --norrna --notrna --outdir {} "
        "--cdsrnaolap --noanno {}".format(
            threads, out_prefix, genome_dir, fasta_in
        )
    )
    return_code = subprocess.run(command_line, shell=True)

    return_code.check_returncode()


def prokka_main():
    args = parse_arguments()

    input_dir = args.input_dir
    output_dir = args.output_dir
    models_dir = args.models_dir
    files_dir = args.files_dir
    threads = args.threads

    list_bins = []
    for entry in input_dir.iterdir():
        if entry.is_file() and entry.suffix.lower() in [".fa", ".fasta"]:
            list_bins.append(entry)

    count_prokka = 0
    prokka_skipped = {}
    prokka_dir = output_dir / Path("prokka")

    print("**Prokka has started, this may take a while. Be patient.")

    for bin_fasta in list_bins:
        len_bin = 0
        bin_name = get_bin_name(bin_fasta)
        for record in SeqIO.parse(bin_fasta, "fasta"):
            len_bin += len(record.seq)
        if len_bin < 5000:
            print(
                "**prokka has found a genome or bin, which is too short to "
                "code proteins (< 5000 bp). As CDSs are an import feature for "
                "vHULK, we will be skipping this: " + bin_fasta.name
            )
            prokka_skipped[bin_name] = bin_fasta
            continue

        run_prokka(bin_fasta, prokka_dir, threads)

        count_prokka += 1
        if count_prokka % 10 == 0:
            print("**Done with {} genomes...".format(count_prokka))

    print("\n**PROKKA finished with no errors")
    print(
        "**{:>{width}} : {}".format("Successfully run", count_prokka, width=20)
    )


if __name__ == "__main__":
    prokka_main()
