"""Assert recovery of known synthetic transcript proportions."""
import csv
import json
from pathlib import Path
root=Path('results/tests/salmon')
with (root/'quant/quant.sf').open() as f:rows=list(csv.DictReader(f,delimiter='\t'))
counts={r['Name']:float(r['NumReads']) for r in rows}
total=sum(counts.values())
assert total>9900,(total,counts)
for name,expected in [('synthetic_tx1',.9),('synthetic_tx2',.08),('synthetic_tx3',.02)]:
    assert abs(counts[name]/total-expected)<.03,(name,counts)
(root/'verified.json').write_text(json.dumps({'status':'PASS','data_type':'synthetic','assigned_pairs':total,'counts':counts},indent=2)+'\n')
