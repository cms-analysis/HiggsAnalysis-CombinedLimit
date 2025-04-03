#!/usr/bin/env python3

import ROOT
import os
import argparse

def check_file_extension(file, allowed_extensions):
    file_extension = os.path.splitext(file)[1].lower()
    if file_extension not in allowed_extensions:
        raise ValueError(f"File '{file}' must have one of the following extensions: {', '.join(allowed_extensions)}")


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

def scan_plain_text_files(file1, file2, output_file):
    try:
        differences_found = False
        line_number = 1  # To track line numbers

        with open(file1, 'r') as f1, open(file2, 'r') as f2, open(output_file, 'w') as output:
            # Compare lines from both files
            for line1, line2 in zip(f1, f2):
                if line1 != line2:
                    output.write(f"Line {line_number}:\nFile 1: {line1}File 2: {line2}\n")
                    differences_found = True
                line_number += 1

            # Handle remaining lines in case the files are of unequal length
            for line in f1:
                output.write(f"Line {line_number}:\nFile 1: {line}File 2: <No line>\n")
                line_number += 1
                differences_found = True

            for line in f2:
                output.write(f"Line {line_number}:\nFile 1: <No line>\nFile 2: {line}\n")
                line_number += 1
                differences_found = True

        # Explicitly check the last lines of both files
        with open(file1, 'r') as f1_last:
            last_line_f1 = None
            for last_line_f1 in f1_last:
                pass  # Read till the last line of file1

        with open(file2, 'r') as f2_last:
            last_line_f2 = None
            for last_line_f2 in f2_last:
                pass  # Read till the last line of file2

        if last_line_f1 != last_line_f2:
            with open(output_file, 'a') as output:  # Append to avoid overwriting previous results
                output.write(f"Last Line:\nFile 1: {last_line_f1 or '<No line>'}File 2: {last_line_f2 or '<No line>'}\n")
            differences_found = True

        # Output result based on whether differences were found
        if differences_found:
            print(f"Differences written to {output_file}")
        else:
            print("Exact match between files.")

    except Exception as e:
        print(f"An error occurred: {e}")

# Parse command-line arguments
parser = argparse.ArgumentParser(description="Compare values from two specified branches in two output files.")
parser.add_argument("-m", "--method", type=str, help="Comparison method: plain text--based or ROOT-based. Arguments: txt, txt")
parser.add_argument("-f1", "--file1", help="First file")
parser.add_argument("-f2", "--file2", help="Second file")
parser.add_argument("-o", "--out", help="Output file")
parser.add_argument("-b1", "--branch1", help="First branch name to compare")
parser.add_argument("-b2", "--branch2", help="Second branch name to compare")
args = parser.parse_args()

if (args.method != "txt") and (args.method != "root"):
    raise ValueError("Incorrect method")
if (args.method == "txt" or args.method == "root") and (args.file1 is None or args.file2 is None):
    raise ValueError("Specify files")
if args.method == "txt":
    allowed_extensions = [".txt", ".output", ".out"]
    check_file_extension(args.file1, allowed_extensions)
    check_file_extension(args.file2, allowed_extensions)
    check_file_extension(args.out, allowed_extensions)
    scan_plain_text_files(args.file1, args.file2, args.out)

if (args.method == "root"):
    allowed_extensions = [".root"]
    check_file_extension(args.file1, allowed_extensions)
    check_file_extension(args.file2, allowed_extensions)

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

if args.method == "txt" and (args.branch1 is not None or args.branch2 is not None):
    raise ValueError("txt method does not accept branch option")

