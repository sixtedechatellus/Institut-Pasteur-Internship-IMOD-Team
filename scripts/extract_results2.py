import sys
import csv

log_file = sys.argv[1]
out_file = sys.argv[2] if len(sys.argv) > 2 else "results_dude.csv"

targets = []
lines = open(log_file).read().splitlines()

i = 0
while i < len(lines):
    line = lines[i].strip()
    # repérer la ligne du shape, souvent "(1, 128)"
    if line.startswith("(") and i + 4 < len(lines):
        try:
            shape = line
            name = lines[i+1].strip()
            parts5 = lines[i+2].strip().split()
            parts6 = lines[i+3].strip().split()
            num    = int(lines[i+4].strip())

            if len(parts5) == 2 and len(parts6) == 2:
                targets.append({
                    "num": num,
                    "target": name,
                    "shape": shape,
                    "n_actives": int(parts5[0]),
                    "n_decoys": int(parts5[1]),
                    "auroc": float(parts6[0]),
                    "bedroc": float(parts6[1]),
                })
                i += 5
                continue
        except Exception:
            pass
    i += 1

# écriture du CSV
with open(out_file, "w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=[
        "num", "target", "shape", "n_actives", "n_decoys", "auroc", "bedroc"
    ])
    writer.writeheader()
    writer.writerows(targets)

print(f"{len(targets)} cibles extraites dans {out_file}")
