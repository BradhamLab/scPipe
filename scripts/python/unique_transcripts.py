import csv

def main(sam, out):
    n_transcripts = 0 
    transcripts = set()
    out_csv = open(out, 'w')
    writer = csv.writer(out_csv, delimiter=',')
    writer.writerow(['n.reads', 'unique.transcripts'])
    with open(sam, 'r') as f: 
        sam = csv.reader(f, delimiter='\t') 
        next(sam) 
        for i, row in enumerate(sam): 
            hit = row[1].replace('SN:', '').replace('_single', '') 
            if hit not in transcripts: 
                transcripts.add(hit) 
                n_transcripts += 1 
            writer.writerow([i + 1, n_transcripts])
    out_csv.close()

if __name__ == "__main__":
    try:
        snakemake
    except NameError:
        snakemake = None
    if snakemake is not None:
        main(snakemake.input['sam'], snakemake.output['csv'])
