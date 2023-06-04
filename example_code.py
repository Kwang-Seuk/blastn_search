import time
from src.blastn_search import blastn_and_parse


input_fasta = "./blastn_testdata.fas"
hit_seqs = 10

start_time = time.time()
blastn_and_parse(input_fasta, hit_seqs)
end_time = time.time()
duration = end_time - start_time

print("Sequences search and data download completed.")
print(f"Total {duration} sec elapsed.")
