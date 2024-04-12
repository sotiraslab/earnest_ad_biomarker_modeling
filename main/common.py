#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 12 15:02:45 2024

@author: tom.earnest
"""

import argparse
import os

def parse():
    parser = argparse.ArgumentParser()
    parser.add_argument('--short', action='store_true')
    return parser.parse_args()

def setup_output(call_file, short=False):
    this_dir = os.path.dirname(os.path.abspath(__file__))
    outputs_dir = os.path.join(this_dir, 'outputs')
    foldername = os.path.splitext(os.path.basename(call_file))[0]
    if short:
        foldername += '_short'
    requested_dir = os.path.join(outputs_dir, foldername)
    if not os.path.isdir(requested_dir):
        os.mkdir(requested_dir)
    return requested_dir