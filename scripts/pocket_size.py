import lmdb, pickle, glob

base = "/pasteur/helix/projects/Imod-moulin/Sixte/Drugclip/Drug-The-Whole-Genome/data/ValidationBiolip"
sizes = []

lmdbs = glob.glob(f"{base}/*/pocket.lmdb")
for f in lmdbs:
    env = lmdb.open(f, readonly=True, lock=False, subdir=False)
    with env.begin() as txn:
        for _, v in txn.cursor():
            d = pickle.loads(v)
            sizes.append(len(d["pocket_atoms"]))
    env.close()

if sizes:
    print("---- Recap ----")
    print(sizes)
    print(f"Number of pockets : {len(sizes)}")
    print(f"Average size : {sum(sizes)/len(sizes):.1f}")
    print(f"Maximum size     : {max(sizes)}")
else:
    print("No pocket.lmdb found")
