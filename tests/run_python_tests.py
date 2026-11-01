"""Execute input-integrity tests inside the processing container."""
import json
from pathlib import Path
import unittest
suite=unittest.defaultTestLoader.discover('tests',pattern='test_*.py')
result=unittest.TextTestRunner(verbosity=2).run(suite)
if not result.wasSuccessful():raise SystemExit(1)
out=Path('results/tests/python/verified.json');out.parent.mkdir(parents=True,exist_ok=True)
out.write_text(json.dumps({'status':'PASS','tests':result.testsRun,'data_type':'synthetic fixtures and public metadata'},indent=2)+'\n')
