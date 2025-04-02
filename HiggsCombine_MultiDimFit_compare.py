#!/usr/bin/env python3

import ROOT
import os
import argparse

def scan_branch(file_name, branch_name):
    """Open a ROOT file and scan the given branch from the 'limit' tree."""
    if not os.path.exists(file_name):
        print(f"Error: File {file_name} not found!")
        return None
    
    file = ROOT.TFile.Open(file_name)
    tree = file.Get("limit")
    
    if not tree:
        print(f"Error: 'limit' tree not found in {file_name}!")
        return None
    
    branch_values = []
    for entry in tree:
        if hasattr(entry, branch_name):
            branch_values.append(getattr(entry, branch_name))
        else:
            print(f"Error: Branch '{branch_name}' not found in {file_name}!")
            return None
    
    return branch_values

# Parse command-line arguments
parser = argparse.ArgumentParser(description="Compare values from two specified branches in two ROOT files.")
parser.add_argument("file1", help="First ROOT file")
parser.add_argument("file2", help="Second ROOT file")
parser.add_argument("branch1", help="First branch name to compare")
parser.add_argument("branch2", help="Second branch name to compare")
args = parser.parse_args()

# Scan values from the specified branches in both files
values1_branch1 = scan_branch(args.file1, args.branch1)
values2_branch1 = scan_branch(args.file2, args.branch1)
values1_branch2 = scan_branch(args.file1, args.branch2)
values2_branch2 = scan_branch(args.file2, args.branch2)

# Compare the outputs
success = True
if values1_branch1 is None or values2_branch1 is None or values1_branch2 is None or values2_branch2 is None:
    success = False
elif values1_branch1 != values2_branch1:
    print(f"Test unsuccessful: Values in branch '{args.branch1}' do not match.")
    success = False
elif values1_branch2 != values2_branch2:
    print(f"Test unsuccessful: Values in branch '{args.branch2}' do not match.")
    success = False

if success:
    print("nllbackend test successful")
else:
    print("nllbackend test unsuccessful")
