#!/usr/bin/env python3
"""Refresh HZZ4L regression expectations after an update

With BUILD_TESTS=ON:

    python3 test/updateHZZ4LReferences.py --build-dir build-lcg
    python3 test/updateHZZ4LReferences.py --build-dir build-lcg --write

The first command rebuilds and previews a diff, --write updates the
expectations block in test/testHZZ4LReadFile.cxx, rebuilds, and runs its CTest

Captures cover BR and matrix checksums, sampled widths/interpolation, HDECAY,
and spline knots, midpoints and extrapolation
"""

import argparse
import difflib
import json
import math
from pathlib import Path
import subprocess
import sys
import tempfile

BEGIN = "  // BEGIN GENERATED HZZ4L EXPECTATIONS\n"
END = "  // END GENERATED HZZ4L EXPECTATIONS\n"
MODELS = {"ggh": 5, "bw": 5, "gghNo": 2, "vbf": 5, "vbfNo": 2}
SPLINES = {"GGH": 41, "VBF8TeV": 35, "VBF7TeV": 29}


def number(value):
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ValueError("reference values must be numbers")
    if not math.isfinite(value):
        raise ValueError("refusing to record a non-finite reference value")
    return repr(float(value))


def checksum(value):
    if type(value) is not int or not 0 <= value < 2**64:
        raise ValueError("a checksum must be an unsigned 64-bit integer")
    return f"UINT64_C(0x{value:016x})"


def shape(values, dimensions, validate):
    if not dimensions:
        validate(values)
    elif not isinstance(values, list) or len(values) != dimensions[0]:
        raise ValueError(f"unexpected reference dimensions; expected {dimensions}")
    else:
        for value in values:
            shape(value, dimensions[1:], validate)


def array_rows(values, format_value, indent="    "):
    lines = []
    for row in values:
        if any(isinstance(value, list) for value in row):
            lines.append(indent + "{{")
            lines.extend(array_rows(row, format_value, indent + "  "))
            lines.append(indent + "}},")
        else:
            tokens = [format_value(value) for value in row]
            compact = indent + "{{" + ", ".join(tokens) + "}},"
            if len(compact) <= 110:
                lines.append(compact)
            else:
                lines.append(indent + "{{")
                lines.extend(indent + "  " + token + "," for token in tokens)
                lines.append(indent + "}},")
    return lines


def array_definition(type_name, name, values, format_value=number):
    return [f"  constexpr {type_name} {name} = {{{{"] + array_rows(values, format_value) + ["  }};", ""]


def render_references(data):
    if not isinstance(data, dict) or type(data.get("schema")) is not int or data["schema"] != 1:
        raise ValueError("unsupported HZZ4L reference schema")
    if not isinstance(data["matrices"], dict) or not isinstance(data["splines"], dict):
        raise ValueError("expected matrix and spline reference objects")
    if set(data["matrices"]) != set(MODELS) or set(data["splines"]) != set(SPLINES):
        raise ValueError("the reference payload does not contain every expected model")
    lines = [
        BEGIN.rstrip("\n"),
        "  // clang-format off",
        "  // Keep captured literals and C++17 aggregate initializers unchanged.",
        "  // NOLINTBEGIN(modernize-use-designated-initializers,modernize-use-std-numbers)",
    ]
    for name in ("brSMChecksum", "brSM4Checksum", "hdecayChecksum"):
        lines.append(f"  constexpr std::uint64_t {name} = {checksum(data[name])};")
    lines.append(f"  constexpr double hdecayUpperWidth = {number(data['hdecayUpperWidth'])};")
    lines.append("")
    for name in ("smWidths", "sm4Widths"):
        rows = data[name]
        shape(rows, [4, 3], number)
        lines.append(f"  constexpr std::array<WidthReference, 4> {name} = {{{{")
        for channel, mass, width in rows:
            if type(channel) is not int or channel not in (0, 1):
                raise ValueError("unexpected SM/SM4 reference channel")
            lines.append(f"    {{{channel}, {number(mass)}, {number(width)}}},")
        lines.extend(["  }};", ""])
    for model, count in MODELS.items():
        reference = data["matrices"][model]
        hashes = reference["hashes"]
        parameters = reference["parameters"]
        shape(hashes, [2, count], checksum)
        shape(parameters, [2, 2, count], number)
        lines.extend(
            array_definition(
                f"std::array<std::array<std::uint64_t, {count}>, 2>",
                model + "ExpectedHashes",
                hashes,
                checksum,
            )
        )
        lines.extend(
            array_definition(
                f"std::array<std::array<std::array<double, {count}>, 2>, 2>",
                model + "ExpectedParameters",
                parameters,
            )
        )
    for spline, count in SPLINES.items():
        values = data["splines"][spline]
        shape(values, [count, 2], number)
        if any(left[0] >= right[0] for left, right in zip(values, values[1:])):
            raise ValueError("spline sample masses must be strictly increasing")
        lines.extend(array_definition(f"std::array<std::array<double, 2>, {count}>", "spline" + spline, values))
    lines.extend(
        [
            "  // NOLINTEND(modernize-use-designated-initializers,modernize-use-std-numbers)",
            "  // clang-format on",
            END.rstrip("\n"),
        ]
    )
    return "\n".join(lines) + "\n"


def replace_references(source, references):
    if source.count(BEGIN) != 1 or source.count(END) != 1:
        raise ValueError("expected exactly one marked reference block in the HZZ4L test file")
    start = source.index(BEGIN)
    end = source.index(END) + len(END)
    if start >= end - len(END):
        raise ValueError("invalid reference block markers")
    return source[:start] + references + source[end:]


def run(command):
    print("Running: " + " ".join(map(str, command)), file=sys.stderr)
    subprocess.run(list(map(str, command)), check=True, stdout=sys.stderr)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--build-dir",
        type=Path,
        required=True,
        help="configured CMake build with BUILD_TESTS=ON",
    )
    parser.add_argument(
        "--write",
        action="store_true",
        help="replace the marked expectations, rebuild, and run CTest",
    )
    parser.add_argument(
        "--jobs",
        type=int,
        default=4,
        help="maximum build parallelism (default: 4)",
    )
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    source_path = root / "test/testHZZ4LReadFile.cxx"
    build_dir = args.build_dir.resolve()
    try:
        cache = (build_dir / "CMakeCache.txt").read_text()
        home = next(line.split("=", 1)[1] for line in cache.splitlines() if line.startswith("CMAKE_HOME_DIRECTORY:"))
        if Path(home).resolve() != root:
            raise ValueError("the build directory belongs to a different source checkout")
        original = source_path.read_text()
        # Validate the markers before building or attempting any update.
        replace_references(original, BEGIN + END)
        build_command = [
            "cmake",
            "--build",
            build_dir,
            "--target",
            "hzz4l-readfile",
            "--parallel",
            args.jobs,
        ]
        run(build_command)
        executable = build_dir / "test/hzz4l-readfile"
        if not executable.is_file():
            raise ValueError(f"cannot find {executable}; use a configured build with BUILD_TESTS=ON")
        with tempfile.TemporaryDirectory(prefix="hzz4l-references-") as temporary:
            output = Path(temporary) / "references.json"
            run([executable, "--dump-hzz4l-references", output])
            reference_data = json.loads(output.read_text())
        updated = replace_references(original, render_references(reference_data))
        if source_path.read_text() != original:
            raise ValueError("the test file changed during capture; refusing to overwrite concurrent edits")
        if original == updated:
            print("HZZ4L expectations are up to date.")
        else:
            sys.stdout.writelines(
                difflib.unified_diff(
                    original.splitlines(keepends=True),
                    updated.splitlines(keepends=True),
                    fromfile="a/test/testHZZ4LReadFile.cxx",
                    tofile="b/test/testHZZ4LReadFile.cxx",
                )
            )
            if args.write:
                source_path.write_text(updated)
                print("Updated the marked HZZ4L expectations.", file=sys.stderr)
            else:
                print("Preview only; use --write to apply this diff.", file=sys.stderr)
        if args.write:
            run(build_command)
            run(
                [
                    "ctest",
                    "--test-dir",
                    build_dir / "test",
                    "-R",
                    "^gtest-hzz4l-readfile$",
                    "--output-on-failure",
                ]
            )
        return 0
    except (
        OSError,
        ValueError,
        KeyError,
        TypeError,
        StopIteration,
        subprocess.CalledProcessError,
    ) as error:
        print(f"HZZ4L reference update failed: {error}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
