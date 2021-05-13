#!/bin/python2
# Script that installs the necessary submodules to build with cmake

# Imports
import os

# Colors
PROGRESS = '\033[1m\033[35m'
WARNING = '\033[33m'
PASS = '\033[32m'
DEFAULT = '\033[0m'

# Submodule Class
class Submodule:
    def _init_(self, path=None, url=None, commit=None):
        self.path = path
        self.url = url
        self.commit = commit

    def put_path(self, path):
        self.path = path

    def put_url(self, url):
        self.url = url

    def put_commit(self, commit):
        self.commit = commit

    # Used for debugging
    def __str__(self):
        return 'path = {}\nurl = {}\ncommit = {}'.format(self.path, self.url, self.commit)

def color_print(message, color=DEFAULT):
    print '{}{}{}'.format(color, message, DEFAULT)

if __name__ == '__main__':
    
    # Get the submodules
    color_print('Getting submodule info...', PROGRESS)
    submodules = []
    with open('.gitmodules') as module_file:
        line = module_file.readline().strip()
        index = -1
        while line:
            if line.startswith('[submodule'):
                submodules.append(Submodule())
                index += 1
            else:
                data = line[1 + line.index('='):].strip()
                if line.startswith('path'):
                    submodules[index].put_path(data)
                elif line.startswith('url'):
                    submodules[index].put_url(data)
                elif line.startswith('branch'):
                    submodules[index].put_commit(data)
            line = module_file.readline().strip()
    color_print('DONE', PROGRESS)

    # Prepare for install
    color_print('Preparing for install...', PROGRESS)
    for submodule in submodules:
        if os.path.isdir(submodule.path):
            color_print('{} exists. Removing...'.format(submodule.path), WARNING)
            os.system('rm -rf {}'.format(submodule.path))
        else:
            color_print('{} does not exist. Good'.format(submodule.path), PASS)
    color_print('DONE', PROGRESS)
    
    # Install submodules
    color_print('Installing submodules...', PROGRESS)
    for submodule in submodules:
        os.system('git clone {} {}'.format(submodule.url, submodule.path))
    color_print('DONE', PROGRESS)

    # Checkout correct commits
    color_print('Checking out correct commits...', PROGRESS)
    workingdir = os.getcwd()
    for submodule in submodules:
        os.chdir(submodule.path)
        try:
            os.system('git checkout {}'.format(submodule.commit))
        except:
            color_print('No commit to check out', WARNING)
        os.chdir(workingdir)
        pass
    color_print('DONE', PROGRESS)
