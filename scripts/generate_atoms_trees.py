#!/usr/bin/env python3

import os
import shutil
import sys
import subprocess

data_dir = './data'
species = 'unicorn'
histories = []

def generate_atoms_tree(history):
  results_dir = 'dupstemp_{}-{}-dna'.format(history, species)
  if os.path.exists(results_dir):
    return
  command = '../scripts/sample_trees.pl {}-{}.dna {}.atoms'.format(history, species, history)
  print(command)
  try:
    subprocess.run(command, stdout=sys.stdout, stderr=sys.stderr, shell=True, check=True)
  except:
    shutil.rmtree(results_dir)


with os.scandir(data_dir) as it:
    for entry in it:
        if entry.is_dir() and entry.name.startswith('generated-'):
            histories.append(entry.name)

histories.sort()

os.chdir(data_dir)

for history in histories:
  generate_atoms_tree(history)
