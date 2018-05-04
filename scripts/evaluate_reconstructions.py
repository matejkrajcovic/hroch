#!/usr/bin/env python3

import sys
import subprocess

def frange(x, y, jump):
  while x < y:
    yield x
    x += jump

def average(l):
  sum = 0
  for el in l:
    sum += el
  return sum / len(l)

data_dir = './data'
histories = list(map(lambda h: 'generated-{}'.format(h), [
  'F4100',
  'F4101',
  'F4102',
  'F4103',
  'F4104',
  'F4105',
  'F4106',
  'F4107',
  'F4108',
  'F4109',
  'F4110',
  'F4111',
  'F4112',
  'F4113',
  'F4114',
  'F4115',
  'F4116',
  'F4117',
  'F4118',
  'F4119',
#]))
#histories = list(map(lambda h: 'generated-{}'.format(h), [
  'F6100',
  'F6101',
  'F6102',
  'F6103',
  'F6104',
  'F6105',
  'F6106',
  'F6107',
  'F6108',
  'F6109',
  'F6110',
  'F6111',
  'F6112',
  'F6113',
  'F6114',
  'F6115',
  'F6116',
  'F6117',
  'F6118',
  'F6119',
]))

def evaluate_reconstruction(history, starting_temperature, prob_previously_used_event, baseline=False, scoring='likelihood'):
  command = './hroch.bin --evaluate-reconstructions --atoms_file data/{}.atoms --trees_dir data/dupstemp_{}-unicorn-dna --reconstructions_file {} --scoring {}'.format(history, history, reconstructions_file, scoring)
  #print(command)
  results = subprocess.run(command, stdout=subprocess.PIPE, stderr=sys.stderr, shell=True, check=True, encoding='utf-8')

  (num_events, orig_events, ji) = results.stdout.strip().split(' ')
  if baseline:
    prob_previously_used_event = "baseline"
  print('{},{},{},{},{}'.format(
    history,
    starting_temperature,
    prob_previously_used_event,
    #results.stdout.strip()))
    num_events,
    ji))

def evaluate_reconstruction_annealing(history, suffix):
  reconstructions_file = './outputs/hroch_data-{}.atoms-{}.histories'.format(history, suffix)
  command = './hroch.bin --evaluate-reconstructions --atoms_file data/{}.atoms --trees_dir data/dupstemp_{}-unicorn-dna --reconstructions_file {}'.format(history, history, reconstructions_file)
  results = subprocess.run(command, stdout=subprocess.PIPE, stderr=sys.stderr, shell=True, check=True, encoding='utf-8')

  (num_events, orig_events, ji) = results.stdout.strip().split(' ')
  print('{},{},{},{},{}'.format(
    history,
    suffix,
    num_events,
    orig_events,
    ji))

""" for history in histories:
  for starting_temperature in frange(starting_prob, ending_prop, step_prob):
    for prob_previously_used_event in frange(starting_prob, ending_prop, step_prob):
      evaluate_reconstruction(history, starting_temperature, prob_previously_used_event)
    evaluate_reconstruction(history, starting_temperature, 0, baseline=True)
    print() """

for history in histories:
  evaluate_reconstruction_annealing(history, 'annealing-likelihood-100')
  evaluate_reconstruction_annealing(history, 'annealing-likelihood-prob-0.3')
  evaluate_reconstruction_annealing(history, 'annealing-likelihood-prob-0.6')
  evaluate_reconstruction_annealing(history, 'annealing-likelihood-prob-0.9')
  evaluate_reconstruction_annealing(history, 'annealing-num_events-0.04')
  evaluate_reconstruction_annealing(history, 'annealing-num_events-prob-0.3')
  evaluate_reconstruction_annealing(history, 'annealing-num_events-prob-0.6')
  evaluate_reconstruction_annealing(history, 'annealing-num_events-prob-0.9')
  evaluate_reconstruction_annealing(history, 'original')
  print()
