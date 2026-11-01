"""Lock drift must stop preparation without altering an existing environment."""
import importlib.util
import os
from pathlib import Path
import shlex
import shutil
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location('workflow_environment', ROOT / 'scripts/check_workflow_environment.py')
CHECKER = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(CHECKER)
PACKAGE_A = 'https://example.org/channel/linux-64/python-3.12.1-h123_0.conda#' + 'a' * 32
PACKAGE_B = 'https://example.org/channel/noarch/snakemake-9.27.1-py_0.tar.bz2#' + 'b' * 64
EXPLICIT = '@EXPLICIT\n' + PACKAGE_A + '\n' + PACKAGE_B + '\n'


class ExplicitEnvironment(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix='workflow-lock-')
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.lock = self.root / 'locked.txt'
        self.exported = self.root / 'exported.txt'
        self.lock.write_text(EXPLICIT)
        self.exported.write_text(EXPLICIT)

    def test_identical_packages_pass(self):
        self.assertEqual(CHECKER.compare_environment(self.lock, self.exported), 2)

    def test_order_comments_and_surrounding_whitespace_do_not_change_identity(self):
        self.exported.write_text('# different export header\n\n @EXPLICIT \n\t' + PACKAGE_B + '\t\n# comment\n ' + PACKAGE_A + ' \n')
        self.assertEqual(CHECKER.compare_environment(self.lock, self.exported), 2)

    def test_changed_version_build_url_or_hash_fails(self):
        variants = {
            'version': PACKAGE_A.replace('3.12.1', '3.12.2'),
            'build': PACKAGE_A.replace('h123_0', 'h123_1'),
            'channel': PACKAGE_A.replace('/channel/', '/other/'),
            'hash': PACKAGE_A[:-1] + 'c',
            'missing_hash': PACKAGE_A.split('#')[0],
        }
        for change, entry in variants.items():
            with self.subTest(change=change):
                self.exported.write_text(EXPLICIT.replace(PACKAGE_A, entry))
                with self.assertRaisesRegex(ValueError, 'differs from the committed explicit lock'):
                    CHECKER.compare_environment(self.lock, self.exported)

    def test_missing_package_fails(self):
        self.exported.write_text('@EXPLICIT\n' + PACKAGE_A + '\n')
        with self.assertRaisesRegex(ValueError, 'Missing locked entries'):
            CHECKER.compare_environment(self.lock, self.exported)

    def test_extra_package_fails(self):
        self.exported.write_text(EXPLICIT + 'https://example.org/noarch/extra-1.0-py_0.conda\n')
        with self.assertRaisesRegex(ValueError, 'Unexpected exported entries'):
            CHECKER.compare_environment(self.lock, self.exported)

    def test_empty_or_malformed_lists_cannot_match_each_other(self):
        cases = {
            'empty': '',
            'comments_only': '# empty export\n',
            'marker_only': '@EXPLICIT\n',
            'missing_marker': PACKAGE_A + '\n',
            'yaml': 'name: workflow\ndependencies:\n  - snakemake\n',
            'non_url': '@EXPLICIT\nnot-a-package\n',
            'bad_extension': '@EXPLICIT\nhttps://example.org/package.txt\n',
            'bad_hash': '@EXPLICIT\nhttps://example.org/package.conda#bad\n',
            'empty_hash': '@EXPLICIT\nhttps://example.org/package.conda#\n',
            'invalid_port': '@EXPLICIT\nhttps://example.org:invalid/package.conda\n',
            'embedded_space': '@EXPLICIT\nhttps://example.org/a b.conda\n',
            'credentials': '@EXPLICIT\nhttps://user:secret@example.org/package.conda\n',
            'duplicate_marker': '@EXPLICIT\n@EXPLICIT\n' + PACKAGE_A + '\n',
        }
        for case, content in cases.items():
            with self.subTest(case=case):
                self.lock.write_text(content)
                self.exported.write_text(content)
                with self.assertRaises(ValueError):
                    CHECKER.compare_environment(self.lock, self.exported)

    def test_duplicate_entries_fail(self):
        self.exported.write_text(EXPLICIT + PACKAGE_A + '\n')
        with self.assertRaisesRegex(ValueError, 'duplicate package entry'):
            CHECKER.compare_environment(self.lock, self.exported)

    def test_cli_reports_drift_and_preserves_both_inputs(self):
        self.exported.write_text(EXPLICIT.replace('3.12.1', '3.12.2'))
        before = [p.read_bytes() for p in (self.lock, self.exported)]
        result = subprocess.run([sys.executable, '-B', str(ROOT / 'scripts/check_workflow_environment.py'),
                                 str(self.lock), str(self.exported)], capture_output=True, text=True, timeout=10)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('Missing locked entries', result.stderr)
        self.assertIn('Unexpected exported entries', result.stderr)
        self.assertIn('preserved', result.stderr)
        self.assertEqual([p.read_bytes() for p in (self.lock, self.exported)], before)

    def test_missing_file_returns_failure_without_traceback(self):
        self.exported.unlink()
        result = subprocess.run([sys.executable, '-B', str(ROOT / 'scripts/check_workflow_environment.py'),
                                 str(self.lock), str(self.exported)], capture_output=True, text=True, timeout=10)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('check failed', result.stderr)
        self.assertNotIn('Traceback', result.stderr)

    def test_committed_project_lock_is_supported(self):
        path = ROOT / 'config/locks/workflow-linux-64.explicit.txt'
        self.assertGreater(CHECKER.compare_environment(path, path), 0)


class PreparationControlFlow(unittest.TestCase):
    """Run the real preparation shell with external operations replaced by stubs."""
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix='workflow-preparation-')
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        for directory in ['scripts', 'config/locks', 'containers', '.workflow-env/bin', 'stub-bin']:
            (self.root / directory).mkdir(parents=True, exist_ok=True)
        for name in ['prepare_remote.sh', 'check_workflow_environment.py']:
            shutil.copy2(ROOT / 'scripts' / name, self.root / 'scripts' / name)
        self.lock = self.root / 'config/locks/workflow-linux-64.explicit.txt'
        self.lock.write_text(EXPLICIT)
        (self.root / 'config/workflow-environment.yaml').write_text('name: synthetic-workflow\n')
        (self.root / 'containers/base-image.json').write_text('{"repository":"synthetic/base","manifest_digest":"sha256:fixture"}\n')
        self.exported = self.root / 'synthetic-export.txt'
        self.exported.write_text(EXPLICIT)
        self.events = self.root / 'events.txt'
        self.executable('stub-bin/uname', 'case "$1" in -s) echo Linux;; -m) echo x86_64;; *) exit 92;; esac\n')
        self.executable('stub-bin/docker', '[ "$1" = info ] || exit 93\necho "docker info stub" >> "$TEST_EVENTS"\n')
        # Only the disk-space preflight is bypassed. The checker and JSON read use
        # the actual test interpreter, so the guard itself is not mocked.
        self.executable('stub-bin/python3',
                        'if [ "$1" = - ]; then cat >/dev/null; echo "disk space stub" >> "$TEST_EVENTS"; exit 0; fi\n'
                        + 'exec ' + shlex.quote(sys.executable) + ' -B "$@"\n')
        self.executable('stub-bin/mamba',
                        'if [ "$1" = env ] && [ "$2" = export ]; then\n'
                        '  echo "mamba export stub" >> "$TEST_EVENTS"\n'
                        '  cat "$TEST_EXPORTED"\n'
                        'else echo "unexpected mamba command" >&2; exit 94; fi\n')
        self.executable('.workflow-env/bin/snakemake',
                        '[ "$1" = --version ] || exit 95\n'
                        'echo "snakemake version stub" >> "$TEST_EVENTS"\necho 9.27.1\n')
        self.executable('containers/build.sh', 'echo "container build stub" >> "$TEST_EVENTS"\n')
        self.executable('run.sh', '[ "$1" = container-test ] || exit 96\necho "container test stub" >> "$TEST_EVENTS"\n')
        self.environment = os.environ.copy()
        self.environment.update(PATH=str(self.root / 'stub-bin') + os.pathsep + os.defpath,
                                MAMBA_BIN=str(self.root / 'stub-bin/mamba'),
                                TEST_EVENTS=str(self.events), TEST_EXPORTED=str(self.exported),
                                OMP_NUM_THREADS='2', OPENBLAS_NUM_THREADS='2', MKL_NUM_THREADS='2',
                                PYTHONDONTWRITEBYTECODE='1')

    def executable(self, relative, body):
        path = self.root / relative
        path.write_text('#!/bin/sh\nset -eu\n' + body)
        path.chmod(0o755)

    def prepare(self):
        return subprocess.run(['bash', 'scripts/prepare_remote.sh'], cwd=self.root,
                              env=self.environment, capture_output=True, text=True, timeout=10)

    def test_matching_reused_environment_reaches_downstream_stubs(self):
        result = self.prepare()
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn('matches the committed explicit lock (2 packages)', result.stdout)
        self.assertEqual(self.events.read_text().splitlines(), [
            'docker info stub', 'disk space stub', 'snakemake version stub',
            'mamba export stub', 'container build stub', 'container test stub'])

    def test_drift_stops_before_container_build_and_preserves_environment_and_export(self):
        self.exported.write_text(EXPLICIT.replace('3.12.1', '3.12.2'))
        environment_file = self.root / '.workflow-env/bin/snakemake'
        before = environment_file.read_bytes()
        result = self.prepare()
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('differs from the committed explicit lock', result.stderr)
        self.assertNotIn('container build stub', self.events.read_text())
        self.assertNotIn('container test stub', self.events.read_text())
        self.assertEqual(environment_file.read_bytes(), before)
        self.assertEqual((self.root / 'results/environment/workflow-linux-64.explicit.txt').read_bytes(),
                         self.exported.read_bytes())
        self.assertEqual(self.lock.read_text(), EXPLICIT)

    def test_malformed_existing_lock_blocks_downstream(self):
        self.lock.write_text('@EXPLICIT\n')
        result = self.prepare()
        self.assertNotEqual(result.returncode, 0)
        self.assertNotIn('container build stub', self.events.read_text())
        self.assertNotIn('container test stub', self.events.read_text())

    def test_yaml_fallback_is_explicitly_unlocked(self):
        for content in [None, '']:
            with self.subTest(lock_content=content):
                if content is None:
                    self.lock.unlink()
                else:
                    self.lock.write_text(content)
                result = self.prepare()
                self.assertEqual(result.returncode, 0, result.stderr)
                self.assertIn('unlocked and not lock-verified', result.stderr)
                self.assertNotIn('matches the committed explicit lock', result.stdout)


if __name__ == '__main__':
    unittest.main()
