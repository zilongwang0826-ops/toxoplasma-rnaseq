"""Execute a workflow step in a verified, content-addressed Docker image."""
import hashlib
import argparse
import json
import os
from pathlib import Path
import subprocess


def run(key, command, threads=2, memory_mb=8000, expected_id=None):
    project = Path(__file__).resolve().parents[2]
    lockfile = project / 'containers/images.lock.json'
    if not lockfile.exists():
        raise RuntimeError('Container lock missing. Build images after approval.')
    entry = json.loads(lockfile.read_text())['images'][key]
    image = entry['id']
    if expected_id and image != expected_id:
        raise ValueError('Image lock changed after DAG construction; restart the workflow')
    for path,digest in entry['recipe_sha256'].items():
        if hashlib.sha256((project/path).read_bytes()).hexdigest() != digest:
            raise ValueError('Container recipe changed; rebuild and validate images first')
    if not image.startswith('sha256:') or len(image) != 71:
        raise ValueError('A content-addressed Docker image ID is required')
    actual = subprocess.check_output(['docker','image','inspect','--format','{{.Id}}',image],text=True).strip()
    if actual != image:
        raise ValueError('Docker image identity mismatch')
    args = ['docker','run','--rm','--pull=never','--user',f'{os.getuid()}:{os.getgid()}',
            '--cpus',str(threads),'--memory',f'{memory_mb}m','--memory-swap',f'{memory_mb}m',
            '--security-opt','no-new-privileges','--cap-drop','ALL',
            '--mount',f'type=bind,source={project},target=/work','--workdir','/work',
            '--env','HOME=/tmp','--env',f'OMP_NUM_THREADS={threads}',
            '--env',f'OPENBLAS_NUM_THREADS={threads}','--env',f'MKL_NUM_THREADS={threads}',
            '--env','LC_ALL=C.UTF-8']
    subprocess.run(args + [image] + command,check=True)


if __name__ == '__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--image',choices=['tools','statistics'],required=True)
    p.add_argument('--expected-id')
    p.add_argument('--threads',type=int,default=2)
    p.add_argument('--memory-mb',type=int,default=8000)
    p.add_argument('command',nargs=argparse.REMAINDER)
    a=p.parse_args()
    command=a.command[1:] if a.command and a.command[0]=='--' else a.command
    if not command or a.threads<1 or a.memory_mb<256:
        p.error('Command and positive resource limits required')
    run(a.image,command,a.threads,a.memory_mb,a.expected_id)
