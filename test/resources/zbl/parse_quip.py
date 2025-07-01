fn = 'quip-train.out'

with open(fn, 'r') as fin:
    lines = [x[3:] for x in fin.readlines() if 'AT ' in x]
with open(fn+'.xyz', 'w') as fout:
    for line in lines:
        fout.write(line)