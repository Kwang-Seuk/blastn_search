import json
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from Bio.Blast import NCBIWWW, NCBIXML
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Blast.Record import Blast
from typing import List
from tqdm import tqdm


def blastn(sequence: str) -> NCBIXML.Blast:
    result_handle = NCBIWWW.qblast("blastn", "nt", sequence)
    blast_record = NCBIXML.read(result_handle)
    return blast_record


def parse_results(
    blast_record: NCBIXML.Blast, fasta: SeqIO.SeqRecord, hit_seqs: int
) -> List[dict]:
    result_list = []
    for alignment in blast_record.alignments[:hit_seqs]:
        for hsp in alignment.hsps:
            species = " ".join(alignment.title.split(" ")[1:3])
            result_dict = {
                "sequence_id": fasta.id,
                "sequence_description": fasta.description,
                "alignment_length": alignment.length,
                "e_value": hsp.expect,
                "species": species,
            }
            result_list.append(result_dict)
    return result_list


def blastn_and_parse(file_path, hit_seqs):
    fasta_sequences = list(SeqIO.parse(open(file_path), "fasta"))
    results = []

    with ThreadPoolExecutor() as executor:
        for i in range(0, len(fasta_sequences), 30):
            step_futures = {
                executor.submit(blastn, str(fasta.seq)): fasta
                for fasta in fasta_sequences[i : i + 10]
            }
            completed = []

            while (
                len(completed) < len(step_futures) / 2
            ):  # 50% of results received
                time.sleep(3)
                completed = [f for f in step_futures if f.done()]

            step_futures.update(
                {
                    executor.submit(blastn, str(fasta.seq)): fasta
                    for fasta in fasta_sequences[i + 10 : i + 20]
                }
            )
            completed = []
            while (
                len(completed) < len(step_futures) / 2
            ):  # 50% of results received
                time.sleep(3)
                completed = [f for f in step_futures if f.done()]

            step_futures.update(
                {
                    executor.submit(blastn, str(fasta.seq)): fasta
                    for fasta in fasta_sequences[i + 20 : i + 30]
                }
            )

            for future in as_completed(step_futures):
                fasta = step_futures[future]
                try:
                    blast_record = future.result()
                    results.extend(
                        parse_results(blast_record, fasta, hit_seqs)
                    )
                except Exception as e:
                    print(f"Failed to get results for {fasta.id}: {e}")

    with open("blast_results.json", "w") as f:
        json.dump(results, f, indent=4)
