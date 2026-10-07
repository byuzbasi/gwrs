"""Release integrity and local dependency handling for the paper workflows."""
from pathlib import Path
from datetime import datetime, timezone
import hashlib
import json
import os
import shutil
import subprocess
import sys
import tempfile

sys.dont_write_bytecode = True
ROOT = Path(__file__).resolve().parent
THREADS = ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS',
           'BLIS_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS',
           'RCPP_PARALLEL_NUM_THREADS')
VERSIONS = {'simulation':'0.4.0','real-data':'0.4.0.9006','data':'0.4.0.9006'}

def require(value, message):
    if not value: raise RuntimeError(message)

def sha256(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda:stream.read(1024*1024),b''): h.update(block)
    return h.hexdigest()

def read_json(path):
    return json.loads(Path(path).read_text())

def atomic_json(path, value, replace=False):
    path=Path(path); path.parent.mkdir(parents=True,exist_ok=True)
    fd,temp=tempfile.mkstemp(prefix='.write-',dir=path.parent)
    try:
        with os.fdopen(fd,'w') as f:
            json.dump(value,f,indent=2,sort_keys=True,allow_nan=False); f.write('\n')
            f.flush(); os.fsync(f.fileno())
        if replace: os.replace(temp,path)
        else: os.link(temp,path)
    finally:
        if os.path.exists(temp): os.unlink(temp)

def verify():
    manifest=read_json(ROOT/'manifest-sha256.json')
    for row in manifest['files']:
        path=ROOT/row['file']
        require(path.is_file() and not path.is_symlink(),f"Missing/linked release file: {row['file']}")
        require(path.stat().st_size==row['bytes'] and sha256(path)==row['sha256'],
                f"Changed release file: {row['file']}")
    print(f"RELEASE_INTEGRITY_PASS: {len(manifest['files'])} files",flush=True)
    return manifest

def library(study):
    return ROOT/'runtime'/('gwrs-'+VERSIONS[study])/'library'

def environment(study):
    env=os.environ.copy()
    # Inherited study controls must not silently alter the frozen design.
    for key in list(env):
        if key.startswith(('GWRS_SELECTION_SIM_','GWRS_STUDY_','GWRS_REPRO_')):
            env.pop(key)
    for key in THREADS: env[key]='1'
    env['PYTHONDONTWRITEBYTECODE']='1'
    libs=[str(library(study)),str(ROOT/'runtime/dependencies')]
    # Keep installed declared dependencies available without altering a user library.
    prior=env.get('R_LIBS','')
    if prior: libs.append(prior)
    env['R_LIBS']=os.pathsep.join(libs)
    env['GWRS_REPRO_ROOT']=str(ROOT)
    return env

def rscript(args, study, **kwargs):
    require(shutil.which('Rscript'),'Rscript is not on PATH; install R first.')
    return subprocess.run(['Rscript','--vanilla',os.path.relpath(args[0], ROOT),*map(str,args[1:])],env=environment(study),
                          cwd=ROOT,check=True,**kwargs)

def setup(study, install_dependencies=False):
    verify()
    dest=library(study); deps=ROOT/'runtime/dependencies'
    deps.mkdir(parents=True,exist_ok=True)
    if install_dependencies:
        rscript([ROOT/'setup.R','install',study],study)
    rscript([ROOT/'setup.R','dependencies',study],study)
    archive=ROOT/'packages'/f'gwrs_{VERSIONS[study]}.tar.gz'
    receipt=dest.parent/'source-sha256.txt'
    if dest.exists():
        require(receipt.is_file() and receipt.read_text().strip()==sha256(archive),
                'Project library source receipt is absent or different; use a fresh extracted copy.')
        require((dest/'gwrs').is_dir(),'Incomplete private library; retain logs and use a new copy.')
    else:
        dest.mkdir(parents=True)
        log=dest.parent/'install.log'
        print('Building the study version of gwrs in its project library.',flush=True)
        with log.open('x') as out:
            subprocess.run(['R','CMD','INSTALL','--preclean','--clean','--no-multiarch',
                            f'--library={dest}',str(archive)],env=environment(study),
                           stdout=out,stderr=subprocess.STDOUT,check=True,cwd=ROOT)
        receipt.write_text(sha256(archive)+'\n')
    rscript([ROOT/'setup.R','verify',study],study)

def utc():
    return datetime.now(timezone.utc).isoformat()
