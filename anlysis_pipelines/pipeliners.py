# pipelining script for creation of simple bioinformat analysis pipelines.
# set __verbose_level to control detail printed to STDOUT as required.
#
# Runs commands by creating daughter processes and checking for exit status. 
#
# Nicholas Gleadall - nick.gleadall@googlemail.com

import subprocess

__verbose_level = 0

VERBOSE_LEVELS = {'DEBUG' : 0,
                  'INFO'  : 1,
                  'WARN'  : 2,
                  'ERROR' : 3,
                  'FATAL' : 4}
REV_VERBOSE_LEVELS = { 0 : 'DEBUG',
                       1 : 'INFO',
                       2 : 'WARN',
                       3 : 'ERROR',
                       4 : 'FATAL' };



__steps = list();

def set_verbose_level( new_level ):

    global __verbose_level

    new_level = new_level.upper()

    if ( new_level in VERBOSE_LEVELS ):
        __verbose_level = VERBOSE_LEVELS[ new_level ]
        print "New verbose level: " + REV_VERBOSE_LEVELS[ __verbose_level ]
    else:
        print "Unknown verbosity level: " + new_level




def verbose_print( message, level ):

    if (__verbose_level > VERBOSE_LEVELS[ level ] ):
        return

    print REV_VERBOSE_LEVELS[ __verbose_level ] + " :: " + message


def add_step(name, cmd):
    global __steps
    __steps.append([ name, cmd ]);


def run_steps():

    for step in steps:
        (name, cmd) = step
        system_call( name, cmd)


#------------------------------------------------------------
# Make system call fucntion that checks for exit status of spawned process.

def system_call( step_name, cmd ):

    verbose_print(cmd, 'DEBUG')

    try:
        subprocess.check_call(cmd, shell=True)

    except subprocess.CalledProcessError as scall:

        verbose_print("Script failed at %s stage - exit code was %s, ouput = %s" % (step_name, scall.returncode, scall.output), 'DEBUG')
        verbose_print("Script failed at %s stage - exit code was %s" % (step_name, scall.returncode), 'INFO')
        exit()
