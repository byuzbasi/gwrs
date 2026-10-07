"""Portable release guard for the unchanged county estimation engine."""
from pathlib import Path
import os
import sys

APP=Path(__file__).resolve().parents[1]
ROOT=APP.parents[1]
sys.path.insert(0,str(ROOT))
import portable
from portable import atomic_json, read_json, require, sha256

def implementation_signature():
    return (APP/'config/implementation-signature.txt').read_text().strip()

def verify(require_frozen_release=True):
    portable.verify()
    cfg=read_json(APP/'config/numerics.json')
    policy=read_json(APP/'config/policy.json')
    require(cfg['n']==3107 and cfg['p']==26,'Unexpected sample dimensions')
    require(policy['full']['total_tasks']==2272,'Unexpected task count')
    binding=read_json(APP/'config/binding.json')
    for row in binding['files']:
        require(sha256(ROOT/row['file'])==row['sha256'],'Changed bound input/source: '+row['file'])
    encoded=__import__('json').dumps(binding['files'],sort_keys=True,separators=(',',':')).encode()
    require(__import__('hashlib').sha256(encoded).hexdigest()==implementation_signature(),
            'Implementation signature mismatch')

def environment():
    env=portable.environment('real-data')
    lib=portable.library('real-data')
    require((lib/'gwrs').is_dir(),'Run setup real-data first.')
    env.update(GWRS_STUDY_V2_APP=str(APP.parent/'study-v2'),
               GWRS_STUDY_V1_APP=str(APP.parent/'study-v1'),
               GWRS_STUDY_INPUT_ROOT=str(APP.parent/'study-v1/frozen/input'),
               GWRS_STUDY_R_LIB=str(lib),GWRS_STUDY_POST_THREADS='1')
    return env
