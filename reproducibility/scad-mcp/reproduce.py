#!/usr/bin/env python3
"""Run the paper workflows on macOS/Linux from any extracted directory."""
from pathlib import Path
import argparse
import csv
from datetime import datetime, timezone, timedelta
import fcntl
import json
import os
import shutil
import subprocess
import sys
import time

sys.dont_write_bytecode=True
import portable as p
ROOT=p.ROOT

def output_path(study, run_id):
    p.require(run_id and all(c.isalnum() or c in '-_' for c in run_id),
              'Run ID must contain only letters, digits, hyphens and underscores.')
    if study=='simulation': return ROOT/'simulation/results'/run_id
    return ROOT/'real-data/study-v3/results'/run_id

def sim_env(args, action):
    env=p.environment('simulation')
    env.update(GWRS_SELECTION_SIM_OUTPUT=str(output_path('simulation',args.run_id)/'study'),
               GWRS_SELECTION_SIM_MODE='production' if args.mode=='full' else 'smoke',
               GWRS_SELECTION_SIM_THREADS=str(args.threads),
               GWRS_SELECTION_SIM_CONFIRM_LARGE='YES' if args.mode=='full' else 'NO',
               GWRS_SELECTION_SIM_SEED='20260826',GWRS_SELECTION_SIM_KERNEL='gaussian',
               GWRS_SELECTION_SIM_SOLVER_TOLERANCE='1e-7',
               GWRS_SELECTION_SIM_MAX_ITERATIONS='2000',
               GWRS_SELECTION_SIM_SELECTION_TOLERANCE='1e-8',
               GWRS_SELECTION_SIM_TRUTH_TOLERANCE='1e-12',GWRS_SELECTION_SIM_GRAIN_SIZE='16',
               GWRS_SELECTION_SIM_KEEP_TASK_DATA='NO',GWRS_REPRO_ACTION=action)
    if args.max_new_tasks: env['GWRS_SELECTION_SIM_MAX_TASKS']=str(args.max_new_tasks)
    return env

def simulation(args):
    output=output_path(args.study,args.run_id)
    cmd=['Rscript','--vanilla','simulation/code/penalized-selection-study.R']
    if args.action in ('preflight','verify'):
        subprocess.run(cmd,env=sim_env(args,args.action),cwd=ROOT,check=True)
        return
    if args.action=='status':
        info=p.read_json(output/'progress.json')
        age=(datetime.now(timezone.utc)-datetime.fromisoformat(info['heartbeat_utc'])).total_seconds()
        info['heartbeat_age_seconds']=age
        info['stale_running_heartbeat']=info['state']=='running' and age>max(120,3*args.heartbeat)
        print(json.dumps(info,indent=2)); return
    if args.action=='report':
        p.require(args.mode=='full','The paper summarizer requires the complete full experiment.')
        subprocess.run(cmd,env=sim_env(args,'verify'),cwd=ROOT,check=True)
        summary=output/'summary-v1'
        if not summary.exists():
            p.rscript([ROOT/'simulation/code/summarize-penalized-selection-study.R',
                       output/'study',summary],'simulation')
        p.rscript([ROOT/'simulation/code/build-report.R',output],'simulation')
        return
    fresh=args.action in ('run','smoke')
    p.require(not fresh or not output.exists(),'Output exists. Use resume or choose a new --run-id.')
    p.require(fresh or output.is_dir(),'Cannot resume an absent run.')
    if args.mode=='full': p.require(args.confirm_full,'Full fitting requires --confirm-full.')
    output.mkdir(parents=True,exist_ok=True)
    lock=(output/'controller.lock').open('a+')
    try: fcntl.flock(lock.fileno(),fcntl.LOCK_EX|fcntl.LOCK_NB)
    except BlockingIOError: raise RuntimeError('Another controller owns this run.')
    started=time.monotonic(); started_iso=p.utc()
    status_file=output/'study/task-status.csv'
    def rows():
        if not status_file.exists(): return []
        with status_file.open() as f: return list(csv.DictReader(f))
    baseline={r['task_id'] for r in rows() if r['status']=='completed'}
    invocation=output/'invocations'/datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%f')
    invocation.mkdir(parents=True)
    last_completed=''; proc=None
    def publish(state):
        nonlocal last_completed
        task_rows=rows(); total=len(task_rows) or (1900 if args.mode=='full' else 2)
        done=[r for r in task_rows if r['status']=='completed']
        failed=sum(r['status'] in ('failed','invalid') for r in task_rows)
        progress={}
        csv_path=output/'study/progress.csv'
        if csv_path.exists():
            with csv_path.open() as f: progress=next(csv.DictReader(f),{})
        current=progress.get('current_task_id','') if state=='running' else ''
        running=int(bool(current) and current not in {r['task_id'] for r in done})
        pending=max(0,total-len(done)-failed-running)
        recent=[r for r in done if r['task_id'] not in baseline]
        if recent: last_completed=recent[-1]['task_id']
        elapsed=time.monotonic()-started
        rate=len(recent)/elapsed if elapsed else 0
        # Cost-weighted ETA is unknown until every pending cell has an observed
        # task from this invocation. No cross-cell extrapolation is presented.
        estimates={}; plan_path=output/'study/task-plan.csv'
        if plan_path.exists():
            with plan_path.open() as f: plan={r['task_id']:r for r in csv.DictReader(f)}
            for r in recent:
                cell=plan[r['task_id']].get('cell_id','')
                estimates.setdefault(cell,[]).append(float(r['elapsed_seconds']))
            remain=[r for r in task_rows if r['status']!='completed']
            calibrated=all(plan[r['task_id']].get('cell_id','') in estimates for r in remain)
            eta=sum(sum(estimates[plan[r['task_id']].get('cell_id','')])/len(estimates[plan[r['task_id']].get('cell_id','')]) for r in remain) if calibrated and remain else None
        else: eta=None
        if state=='complete': eta=0
        now=datetime.now(timezone.utc)
        value=dict(run_id=args.run_id,job_id=os.getenv('SLURM_JOB_ID',''),state=state,
                   phase=progress.get('phase','initializing'),work_unit='one validated dataset with eight method fits',
                   total=total,completed=len(done),running=running,failed=failed,pending=pending,
                   percent_complete=100*len(done)/total,elapsed_seconds=elapsed,
                   reused_completed=len(baseline),throughput_units_per_hour=3600*rate,
                   eta_seconds=eta,estimated_completion_utc=(now+timedelta(seconds=eta)).isoformat() if eta is not None else None,
                   eta_scope='numerical task stage; aggregation, verification and optional figures excluded',
                   eta_status='measured same-cell durations' if eta is not None else 'unknown: pending cell costs not calibrated in this invocation',
                   remaining_slurm_seconds=(max(0,float(os.environ['SLURM_JOB_END_TIME'])-time.time()) if os.getenv('SLURM_JOB_END_TIME','').isdigit() else None),heartbeat_utc=now.isoformat(),current_unit=current,
                   last_completed_unit=last_completed,invocation_started_utc=started_iso)
        p.atomic_json(output/'progress.json',value,replace=True)
        hist=output/'progress.tsv'; new=not hist.exists()
        with hist.open('a',newline='') as f:
            writer=csv.DictWriter(f,fieldnames=list(value),delimiter='\t')
            if new: writer.writeheader()
            writer.writerow(value); f.flush()
        print(f"{value['heartbeat_utc']} {state} {len(done)}/{total} ({value['percent_complete']:.1f}%) phase={value['phase']} ETA={eta if eta is not None else 'unknown'}",flush=True)
    try:
        with (invocation/'run.log').open('x') as log:
            proc=subprocess.Popen(cmd,env=sim_env(args,'run'),cwd=ROOT,stdout=log,stderr=subprocess.STDOUT)
            publish('running')
            while True:
                try: code=proc.wait(timeout=args.heartbeat); break
                except subprocess.TimeoutExpired: publish('running')
            if code: raise RuntimeError('Simulation failed; inspect the preserved invocation log.')
        complete=(output/'study/COMPLETED').is_file()
        if complete: subprocess.run(cmd,env=sim_env(args,'verify'),cwd=ROOT,check=True)
        publish('complete' if complete else 'incomplete')
        p.require(complete or args.max_new_tasks is not None,'Run stopped before validated completion.')
    except BaseException:
        if proc is not None and proc.poll() is None: proc.terminate(); proc.wait()
        publish('failed'); raise
    finally: lock.close()

def real_data(args):
    app=ROOT/'real-data/study-v3'; output=output_path(args.study,args.run_id)
    action='run' if args.action in ('smoke','resume') else args.action
    fresh=args.action in ('run','smoke')
    p.require(not fresh or not output.exists(),'Output exists. Use resume or a new --run-id.')
    if args.action=='resume': p.require(output.is_dir(),'Cannot resume an absent run.')
    env=p.environment('real-data')
    if args.mode=='full' and action=='run':
        p.require(args.confirm_full,'Full fitting requires --confirm-full.')
        env['GWRS_STUDY_LOCAL_FULL']='YES'; env['GWRS_STUDY_PRODUCTION']='YES'
    command=[sys.executable,str(app/'code/run.py'),action,'--mode',args.mode,
             '--output',str(output),'--workers',str(args.workers),'--heartbeat',str(args.heartbeat)]
    if action=='run': command+=['--author-run']
    if args.max_new_tasks: command+=['--max-new-tasks',str(args.max_new_tasks),'--allow-incomplete']
    if action=='report':
        p.require(args.mode=='full','County publication figures require a full real-data run.')
        subprocess.run([sys.executable,str(app/'code/run.py'),'verify','--mode','full','--output',str(output)],env=env,cwd=ROOT,check=True)
        for script,suffix in [('render-results.R','publication'),('distance-profile.R','distance')]:
            dest=output.parent/(args.run_id+'-'+suffix)
            verb='verify' if (dest/'COMPLETED').exists() else 'create'
            subprocess.run(['Rscript','--vanilla',os.path.relpath(app/'code'/script,ROOT),verb,str(output),str(dest)],env=env,cwd=ROOT,check=True)
    else: subprocess.run(command,env=env,cwd=ROOT,check=True)

def prepare_data(args):
    p.require(args.action in ('prepare','download'),'Use data prepare or data download.')
    work=ROOT/'work'/args.run_id
    if not work.exists(): shutil.copytree(ROOT/'data-preparation',work)
    if args.action=='download':
        for name in ('00_download_metadata.py','01_download_data.py'):
            subprocess.run([sys.executable,str(work/'code'/name)],check=True,cwd=work)
        return
    for name in ('00_validate_metadata.R','02_audit_sample.R','99_validate_design.R',
                 '03_prepare_spatial_design.R','05_validate_spatial_design_v2.R'):
        p.rscript([work/'code'/name],'data')
    p.rscript([ROOT/'validate-prepared-data.R',work],'data')

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('study',choices=('simulation','real-data','data'))
    parser.add_argument('action',choices=('setup','check','preflight','smoke','run','resume','status','verify','report','prepare','download'))
    parser.add_argument('--mode',choices=('smoke','full'),default='smoke')
    parser.add_argument('--run-id',default='smoke-v1')
    parser.add_argument('--workers',type=int,default=1)
    parser.add_argument('--threads',type=int,default=1)
    parser.add_argument('--heartbeat',type=float,default=30)
    parser.add_argument('--max-new-tasks',type=int)
    parser.add_argument('--confirm-full',action='store_true')
    parser.add_argument('--install-dependencies',action='store_true')
    args=parser.parse_args()
    p.require(args.workers>=1 and args.threads>=1 and args.heartbeat>=1,'Invalid resource setting.')
    p.require(args.max_new_tasks is None or args.max_new_tasks>0,'Task budget must be positive.')
    output_path('simulation',args.run_id)  # validate the identifier for all actions
    if args.action=='setup': p.setup(args.study,args.install_dependencies); return
    p.verify()
    if args.action=='check': return
    p.rscript([ROOT/'setup.R','verify',args.study],args.study)
    if args.action=='smoke': args.mode='smoke'
    if args.study=='simulation': simulation(args)
    elif args.study=='real-data': real_data(args)
    else: prepare_data(args)

if __name__=='__main__':
    try: main()
    except Exception as error:
        print('ERROR:',error,file=sys.stderr,flush=True); raise SystemExit(1)
