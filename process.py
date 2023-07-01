import os
import subprocess



def run(cmd:str, chained:bool=True, stdin:str=None, cwd:str=None, shell:bool=True):
    ''' Run commands in bash shell'''
    proc = subprocess.Popen(cmd, cwd=cwd, shell=shell, stdin=stdin, stdout=subprocess.PIPE)
    if not chained:
        pid = os.getpgid(proc.pid)
        proc.communicate()
        try:
            os.waitpid(pid, 0)
        except (ProcessLookupError, ChildProcessError):
            pass
    
    return proc 


if __name__ == "__main__":
    pass
