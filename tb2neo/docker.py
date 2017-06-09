from __future__ import print_function
import time
import random
import shlex
import subprocess
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


def launch_neo4j_docker(outputdir, container_name=None):
    """Launch a Neo4j server in a Docker container

    outputdir:str -
    """
    if container_name is None:
        container_name = '_'.join([random.choice(WORDS), random.choice(WORDS)])
    cmd_str = (
        "docker run --rm -v {}:/data -p 7474 -p 7687 --name={} -e NEO4J_AUTH=none neo4j:3.2.0".
        format(
            outputdir, container_name
        )
    )
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
