import sys
import os
import ROOT
from enum import Enum, auto


class ComparisonResult(Enum):
    OK = 0
    NO_LIMIT_TREE = auto()
    DIFFERENT_ENTRIES = auto()
    DIFFERENT_POIS = auto()
    VALUE_MISMATCH = auto()


def detect_keys(startpath) -> set:
    keys = set()
    dirs = os.listdir(startpath)
    for d in dirs:
        if os.path.isdir(os.path.join(startpath, d)):
            if d.endswith("_codegen"):
                key = d[: -len("_codegen")]
                keys.add(key)
            else:
                keys.add(d)

    def get_mtime(k):
        nominal_path = os.path.join(startpath, k)
        return os.path.getmtime(nominal_path) if os.path.exists(nominal_path) else 0

    return sorted(keys, key=get_mtime)


def check_codegen_counterparts(startpath, keys) -> list:
    missing_codegen = []
    for key in keys:
        codegen_dir = os.path.join(startpath, f"{key}_codegen")
        if not os.path.exists(codegen_dir):
            missing_codegen.append(key)
    return missing_codegen


def check_codegen_counterpart_files(startpath, keys, missing_codegen) -> dict[str, set]:
    discrepancies = {}
    for key in keys:
        if key in missing_codegen:
            continue

        nominal_files = set(os.listdir(os.path.join(startpath, key)))
        codegen_files = set(os.listdir(os.path.join(startpath, f"{key}_codegen")))

        # Union should equal both sets
        missing_in_codegen = nominal_files - codegen_files
        missing_in_nominal = codegen_files - nominal_files
        discrepancies[key] = {"missing_in_codegen": missing_in_codegen, "missing_in_nominal": missing_in_nominal}
    return discrepancies


def compare_file_contents(file1, file2, tol: float = 1e-3) -> ComparisonResult:
    f1 = ROOT.TFile.Open(file1)
    f2 = ROOT.TFile.Open(file2)
    try:
        # For now, just check the 'limit' tree
        keys1 = [k.GetName() for k in f1.GetListOfKeys()]
        keys2 = [k.GetName() for k in f2.GetListOfKeys()]
        if "limit" not in keys1 or "limit" not in keys2:
            if "fitDiagnosticsTest" in file1:
                return ComparisonResult.OK  # Skip fitDiagnosticsTest for now
            return ComparisonResult.NO_LIMIT_TREE

        tree1 = f1.Get("limit")
        tree2 = f2.Get("limit")
        if tree1.GetEntries() != tree2.GetEntries():
            return ComparisonResult.DIFFERENT_ENTRIES

        # Check POIs match
        pois_1 = [b.GetName() for b in tree1.GetListOfBranches() if "r_" in b.GetName() or "r" == b.GetName()]
        pois_2 = [b.GetName() for b in tree2.GetListOfBranches() if "r_" in b.GetName() or "r" == b.GetName()]
        if set(pois_1) != set(pois_2):
            return ComparisonResult.DIFFERENT_POIS

        # Check entry by entry
        for i in range(tree1.GetEntries()):
            tree1.GetEntry(i)
            tree2.GetEntry(i)
            for poi in pois_1:
                val1 = getattr(tree1, poi)
                val2 = getattr(tree2, poi)
                deltaNLL1 = tree1.deltaNLL
                deltaNLL2 = tree2.deltaNLL
                if abs(val1 - val2) > tol or abs(deltaNLL1 - deltaNLL2) > tol:
                    return ComparisonResult.VALUE_MISMATCH

        return ComparisonResult.OK
    finally:
        f1.Close()
        f2.Close()


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python compare_codegen.py <comparison_input_directory>")
        sys.exit(1)
    comparison_input_dir = sys.argv[1]
    keys = detect_keys(comparison_input_dir)
    status = {k: "OK" for k in keys}

    # ---- Comparisons ---

    # 1. Check every directory has a codegen counterpart
    missing_codegen = check_codegen_counterparts(comparison_input_dir, keys)
    for key in missing_codegen:
        status[key] = "MISSING_CODEGEN"

    # 2. Check every codegen directory has the same files as its counterpart
    discrepancies = check_codegen_counterpart_files(comparison_input_dir, keys, missing_codegen)
    for key, diff in discrepancies.items():
        if diff["missing_in_codegen"] or diff["missing_in_nominal"]:
            status[key] = "FILE_MISMATCH"

    # 3. Check file contents
    for key in keys:
        if status[key] != "OK":
            continue

        nominal_files = os.listdir(os.path.join(comparison_input_dir, key))
        codegen_files = os.listdir(os.path.join(comparison_input_dir, f"{key}_codegen"))
        if len(nominal_files) == 0 or len(codegen_files) == 0:
            status[key] = "NO_FILES"
            continue

        for f in nominal_files:
            nominal_file = os.path.join(comparison_input_dir, key, f)
            codegen_file = os.path.join(comparison_input_dir, f"{key}_codegen", f)
            res = compare_file_contents(nominal_file, codegen_file)
            if res != ComparisonResult.OK:
                status[key] = res.name

    # ---- Report ----
    print("Comparison Summary:")
    for key, stat in status.items():
        print(f"  {key}: {stat}")
