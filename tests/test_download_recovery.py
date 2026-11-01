"""Bounded transfer retries preserve partial bytes and never publish bad inputs."""
import hashlib
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'workflow/scripts'))
from download import fetch, transfer

class DownloadRecovery(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.output=Path(self.temp.name)/'read.fastq.gz'
        self.part=Path(str(self.output)+'.part')
    def test_connection_reset_resumes_current_bytes(self):
        self.part.write_bytes(b'ab')
        seen=[]
        def simulate(command,**kwargs):
            seen.append(self.part.read_bytes())
            if len(seen)==1:
                with self.part.open('ab') as f:f.write(b'c')
                return subprocess.CompletedProcess(command,56,stdout='206')
            with self.part.open('ab') as f:f.write(b'def')
            return subprocess.CompletedProcess(command,0,stdout='206')
        with patch('download.subprocess.run',side_effect=simulate) as run,patch('download.time.sleep'):
            fetch('https://example.org/read.gz',self.output,hashlib.md5(b'abcdef').hexdigest(),6)
        self.assertEqual(seen,[b'ab',b'abc'])
        self.assertEqual(self.output.read_bytes(),b'abcdef')
        self.assertEqual(run.call_count,2)
        for call in run.call_args_list:
            cmd=call.args[0];self.assertEqual(cmd[cmd.index('--continue-at')+1],'-')
    def test_complete_payload_followed_by_reset_is_verified_without_retry(self):
        def interrupted_after_last_byte(command,**kwargs):
            self.part.write_bytes(b'abcdef')
            return subprocess.CompletedProcess(command,56,stdout='206')
        with patch('download.subprocess.run',side_effect=interrupted_after_last_byte) as run,patch('download.time.sleep') as sleep:
            fetch('https://example.org/read.gz',self.output,hashlib.md5(b'abcdef').hexdigest(),6)
        self.assertEqual(run.call_count,1);sleep.assert_not_called()
        self.assertEqual(self.output.read_bytes(),b'abcdef')
    def test_range_certificate_and_write_errors_stop(self):
        for code,status in [(22,'416'),(60,'000'),(23,'200')]:
            with self.subTest(code=code):
                self.part.write_bytes(b'ab')
                with patch('download.subprocess.run',return_value=subprocess.CompletedProcess([],code,stdout=status)) as run,patch('download.time.sleep') as sleep:
                    with self.assertRaises(subprocess.CalledProcessError):transfer('https://example.org/read.gz',self.part)
                self.assertEqual(run.call_count,1);sleep.assert_not_called()
                self.assertEqual(self.part.read_bytes(),b'ab')

    def test_transient_http_retries(self):
        responses=[subprocess.CompletedProcess([],22,stdout='503'),subprocess.CompletedProcess([],0,stdout='200')]
        with patch('download.subprocess.run',side_effect=responses) as run,patch('download.time.sleep'):
            self.assertEqual(transfer('https://example.org/read.gz',self.part),2)
        self.assertEqual(run.call_count,2)
    def test_permanent_error_not_retried_or_deleted(self):
        self.part.write_bytes(b'ab')
        with patch('download.subprocess.run',return_value=subprocess.CompletedProcess([],22,stdout='404')) as run,patch('download.time.sleep') as sleep:
            with self.assertRaises(subprocess.CalledProcessError):transfer('https://example.org/read.gz',self.part)
        self.assertEqual(run.call_count,1);sleep.assert_not_called()
        self.assertEqual(self.part.read_bytes(),b'ab')
    def test_retry_budget_exhausted_preserves_partial(self):
        self.part.write_bytes(b'ab')
        with patch('download.subprocess.run',return_value=subprocess.CompletedProcess([],56,stdout='206')) as run,patch('download.time.sleep') as sleep:
            with self.assertRaises(subprocess.CalledProcessError):transfer('https://example.org/read.gz',self.part)
        self.assertEqual(run.call_count,4);self.assertEqual(sleep.call_count,3)
        self.assertEqual(self.part.read_bytes(),b'ab');self.assertFalse(self.output.exists())
    def test_complete_partial_is_verified_without_transfer(self):
        self.part.write_bytes(b'abcdef')
        with patch('download.subprocess.run') as run:
            fetch('https://example.org/read.gz',self.output,hashlib.md5(b'abcdef').hexdigest(),6)
        run.assert_not_called();self.assertEqual(self.output.read_bytes(),b'abcdef')
    def test_bad_complete_partial_not_published(self):
        self.part.write_bytes(b'abcdef')
        with patch('download.subprocess.run') as run:
            with self.assertRaises(ValueError):fetch('https://example.org/read.gz',self.output,'0'*32,6)
        run.assert_not_called();self.assertFalse(self.output.exists())
        self.assertEqual(self.part.read_bytes(),b'abcdef')
    def test_oversize_partial_not_overwritten(self):
        self.part.write_bytes(b'toolong')
        with patch('download.subprocess.run') as run:
            with self.assertRaises(ValueError):fetch('https://example.org/read.gz',self.output,'0'*32,2)
        run.assert_not_called();self.assertEqual(self.part.read_bytes(),b'toolong')

if __name__=='__main__':unittest.main()
