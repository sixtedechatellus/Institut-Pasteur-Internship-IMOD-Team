#!/usr/bin/env python3
import sys, json
smis = [l.strip().split()[0] for l in open(sys.argv[1]) if l.strip()]
json.dump({"smiles": smis}, open(sys.argv[2], "w"))
