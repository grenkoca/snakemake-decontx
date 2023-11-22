#!/usr/bin/env python3

# submits sge jobs

import os
import sys
import math
# from itertools import izip # for py2.7 and old method
from snakemake.utils import read_job_properties

# below is a snakemake way to do this
jobscript = sys.argv[-1]
job_properties = read_job_properties(jobscript)


cluster = job_properties.get('cluster', {})

group = cluster.get('group', 'snakemake')  # job group
queue = cluster.get('queue', None)  # job queue
gpu = cluster.get('gpu', None)  # job queue
job_name = cluster.get('name', 'snakemake')  # job name
time = cluster.get('time', '10-00:00:00')  # job name

threads = job_properties.get('threads', 1)
memory = cluster.get('memory', 10)  # default to 10Gb
if memory:
    # mem is specified as total memory in Gb, but we
    # need to give it to SGE as memory per thread in Mb
    #memory = int(math.floor(float(memory * 1e9) / (threads * 1e6)))
    resources = '--mem={}g'.format(
        memory
    )

output = cluster.get('output', None)
if output:
    output = os.path.realpath(output)
    tdir = os.path.dirname(output)
    if not os.path.exists(tdir):
        os.makedirs(tdir)
error = cluster.get('error', None)
if error:
    error = os.path.realpath(error)
    tdir = os.path.dirname(error)
    if not os.path.exists(tdir):
        os.makedirs(tdir)

if cluster.get('no_sub', 0):
    cmd = '{}'.format(jobscript)
    os.system(cmd)
else:
    # build the command
    cmd = 'sbatch'
    if queue:
        cmd = '{} -q {}'.format(cmd, queue)

    # use gpu
    if gpu:
        cmd = '{} --partition=gpu --gres=gpu:{}'.format(cmd, gpu)

    if time:
        cmd = '{} --time={}'.format(cmd, time)

    cmd = '{} "{}"'.format(
        cmd,
        resources
    )
    if threads >= 1:
        cmd = '{} --cpus-per-task={}'.format(cmd, threads)
    if output:
        cmd = '{} --output {}'.format(cmd, output)
    if error:
        cmd = '{} --error {}'.format(cmd, error)
    cmd = '{} {}'.format(cmd, jobscript)

    print(cmd)

    # run the command
    os.system(cmd)
