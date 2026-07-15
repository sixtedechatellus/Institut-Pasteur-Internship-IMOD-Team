import sys
import csv

log_file = sys.argv[1]
out_file = sys.argv[2] if len(sys.argv) > 2 else "results_dude.csv"

targets = []
lines = open(log_file).readlines()

i = 0
while i < len(lines):
    line = lines[i].strip()
    if line.isdigit() and i + 5 < len(lines):
        try:
            num       = int(line)
            n_mols    = int(lines[i+1].strip())
            shape     = lines[i+2].strip()
            name      = lines[i+3].strip()
            parts5    = lines[i+4].strip().split()
            parts6    = lines[i+5].strip().split()
            if len(parts5) == 2 and len(parts6) == 2:
                targets.append({
                    "num": num,
                    "target": name,
                    "n_molecules": n_mols,
                    "n_actives": int(parts5[0]),
                    "n_decoys": int(parts5[1]),
                    "auroc": float(parts6[0]),
                    "bedroc": float(parts6[1]),
                })
                i += 6
                continue
        except (ValueError, IndexError):
            pass
    i += 1

with open(out_file, "w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=[
        "num", "target", "n_molecules",
        "n_actives", "n_decoys", "auroc", "bedroc"
    ])
    writer.writeheader()
    writer.writerows(targets)

print(f"{len(targets)} cibles extraites dans {out_file}")
