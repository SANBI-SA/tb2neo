from __future__ import print_function
import time
import random
import shlex
import subprocess
import os
try:
    from io import StringIO
except:
    from StringIO import StringIO


WORDS = ['again', 'ages', 'almost', 'america', 'ancient', 'another',
         'antagonisms', 'apprentices', 'armies', 'arrangement', 'away',
         'background', 'before', 'between', 'bourgeois', 'bourgeoisie',
         'burgesses', 'burghers', 'camps', 'cape', 'capital', 'carried',
         'chartered', 'chinese', 'class', 'classes', 'closed', 'colonies',
         'colonisation', 'commerce', 'commodities', 'common', 'communication',
         'complicated', 'conditions', 'constant', 'contending', 'corporate',
         'course', 'demand', 'developed', 'development', 'different',
         'directly', 'discovery', 'distinct', 'division', 'done', 'down',
         'each', 'earlier', 'earliest', 'eastindian', 'either', 'element',
         'elements', 'ended', 'epoch', 'epochs', 'established', 'even', 'ever',
         'every', 'everywhere', 'exchange', 'existing', 'extended',
         'extension', 'face', 'facing', 'feature', 'feudal', 'fight', 'find',
         'first', 'forms', 'freeman', 'fresh', 'from', 'gave', 'generally',
         'giant', 'given', 'gradation', 'gradations', 'great', 'ground',
         'growing', 'guildmaster', 'guildmasters', 'guilds', 'handed', 'have',
         'hidden', 'history', 'hitherto', 'hostile', 'however', 'immense',
         'impulse', 'increase', 'increased', 'industrial', 'industry', 'into',
         'itself', 'journeyman', 'journeymen', 'kept', 'knights', 'known',
         'labour', 'land', 'large', 'leaders', 'long', 'longer', 'lord',
         'lords', 'machinery', 'manifold', 'manufacture', 'manufacturer',
         'manufacturing', 'market', 'markets', 'means', 'meantime', 'middle',
         'millionaires', 'modern', 'modes', 'monopolised', 'more',
         'navigation', 'never', 'ones', 'open', 'opened', 'opposition',
         'oppressed', 'oppression', 'oppressor', 'orders', 'other',
         'patrician', 'patricians', 'paved', 'place', 'plebeian', 'plebeians',
         'possesses', 'product', 'production', 'proletariat', 'proportion',
         'pushed', 'railways', 'rank', 'rapid', 'reacted', 'reconstitution',
         'revolutionary', 'revolutionised', 'revolutions', 'rising', 'rome',
         'rounding', 'ruin', 'ruins', 'same', 'serf', 'serfs', 'series',
         'side', 'simplified', 'single', 'slave', 'slaves', 'social',
         'society', 'splitting', 'sprang', 'sprouted', 'steam', 'stood',
         'struggle', 'struggles', 'subordinate', 'sufficed', 'system', 'taken',
         'that', 'thereby', 'therefore', 'thereupon', 'these', 'this', 'time',
         'took', 'tottering', 'towns', 'trade', 'turn', 'uninterrupted',
         'vanished', 'various', 'vassals', 'wants', 'were', 'which', 'whole',
         'with', 'word', 'workshop', 'world']


def kill_docker(proc, container_name):
    cmd_str = "docker kill " + container_name
    cmd = shlex.split(cmd_str)
    subprocess.call(cmd)
    proc.kill()
    proc.wait()


def launch_neo4j_docker(outputdir, container_name=None, use_bolt=False,
                        image_name='quay.io/thoba/neo_ie:3.1'):
    """Launch a Neo4j server in a Docker container

    outputdir:str -
    """
    if container_name is None:
        container_name = '_'.join([random.choice(WORDS), random.choice(WORDS)])
    if use_bolt:
        bolt_string = ' -e ENABLE_BOLT=true '
    else:
        bolt_string = ''
    cmd_str = (
        ("docker run --rm -v {outputdir}:/data -p 7474 -p 7687 " +
         "-e USER_UID={uid} -e USER_GID={gid} -e MONITOR_TRAFFIC=false " +
         "{bolt} --name={name} -e NEO4J_AUTH=none {image_name}").
        format(
            outputdir=outputdir,
            uid=os.getuid(),
            gid=os.getgid(),
            name=container_name,
            bolt=bolt_string,
            image_name=image_name
        )
    )
    print("STR:", cmd_str)
    cmd = shlex.split(cmd_str)
    proc = subprocess.Popen(cmd)
    time.sleep(10)
    return (proc, container_name)


def find_docker_portmapping(container_name):
    cmd_str = "docker port {}".format(container_name)
    cmd = shlex.split(cmd_str)
    output = subprocess.check_output(cmd)
    port_mapping = dict()
    for line in StringIO(output.decode("utf-8")):
        (dest_str, src_str) = line.split(' -> ')
        # print(dest_str.split('/'))
        dest_port = int(dest_str.split('/')[0])
        src_port = int(src_str.split(':')[1])
        port_mapping[dest_port] = src_port
    return(port_mapping)
