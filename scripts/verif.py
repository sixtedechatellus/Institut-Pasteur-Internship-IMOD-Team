import json

f = "data/ValidationBiolip/1A4M/molecules_1A4M.json"
data = json.load(open(f))
print("len =", len(data))
print("sample =", data[:3])
