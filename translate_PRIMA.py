#!/bin/python
import os
import re
import shutil
import subprocess
from pathlib import Path


def trim_file(file_path):
    """Removes leading and trailing empty lines and whitespaces from a file."""
    try:
        content = file_path.read_text(encoding="utf-8", errors="ignore")
        # .strip() removes all leading/trailing whitespace and newlines
        # We append a single "\n" at the end as per standard POSIX
        # file conventions
        file_path.write_text(content.strip() + "\n", encoding="utf-8")
    except Exception as e:
        print(f"Warning: Could not trim {file_path.name}. Details: {e}")


def run_preprocessor(compiler, flags, src_file, target_file, include_dir):
    """Helper function to run the compiler for preprocessing."""
    cmd = (
        [compiler]
        + flags
        + [str(src_file), f"-I{include_dir}", "-o", str(target_file)]
    )
    return subprocess.run(cmd, check=True)


def preprocess_sources(src_dir, src_files, output_dir):
    preprocessed_files = []
    for file_abs_path in src_files:
        file_rel_path = Path(os.path.relpath(file_abs_path, start=src_dir))

        full_path = src_dir / file_rel_path
        target_path = output_dir / file_rel_path
        preprocessed_files.append(target_path)

        if not full_path.exists():
            print(f"Skipping (not found): {file_rel_path}")
            continue

        # Create target subdirectories automatically
        target_path.parent.mkdir(parents=True, exist_ok=True)

        suffix = full_path.suffix

        # .F90 -> Preprocess with gfortran
        if suffix == ".F90":
            run_preprocessor(
                "gfortran",
                ["-E", "-cpp", "-P"],
                full_path,
                target_path,
                src_dir,
            )
            trim_file(target_path)
        # .f90 or other -> Copy directly
        else:
            shutil.copy2(full_path, target_path)

    return preprocessed_files


def extract_source_files(src_dir):
    """Extracts source file paths from CMake."""
    cmake_file = src_dir / "CMakeLists.txt"

    if not cmake_file.exists():
        print(f"Error: Cannot find {cmake_file}")
        return

    # Parse CMakeLists.txt for source files
    content = cmake_file.read_text(encoding="utf-8", errors="ignore")
    pattern = re.compile(r"add_library\s*\((.*?)\)", re.DOTALL | re.IGNORECASE)
    matches = pattern.findall(content)

    if not matches:
        print("Error: No 'add_library' statement found in CMakeLists.txt.")
        return

    # Extract file tokens and filter out keywords/variables
    raw_tokens = matches[0].split()
    keywords = {"SHARED", "STATIC", "INTERFACE", "MODULE", "OBJECT", "ALIAS"}
    return [
        src_dir / Path(t.strip("'\""))
        for t in raw_tokens[1:]
        if t.upper() not in keywords and not t.startswith("$")
    ]


def translate_sources(
    src_dir, preprocessed_files, output_dir, translator_exec="4ft2pm"
):
    project_dir = Path("/home/yuchen/Syncthing/graduationThesis/4ft2pm")
    cmd = [
        "cabal",
        "run",
        "4ft2pm",
        "--",
        "-r",
        str(src_dir),
        "-o",
        str(output_dir),
        "--create-setup-m",
        "--try-full-consistency",
    ]
    return subprocess.run(
        cmd,
        cwd=project_dir,
    )


if __name__ == "__main__":
    current_dir = Path(__file__).parent
    src_dir = (current_dir / "fortran").resolve()
    preprocess_dir = (current_dir / "preprocessed_fortran/").resolve()
    output_dir = (current_dir / "matlab/interfaces/prima_matlab/").resolve()
    src_files = extract_source_files(src_dir)
    src_files.append(
        (
            current_dir / "fortran/examples/newuoa/newuoa_example_1.f90"
        ).resolve()
    )
    preprocessed_files = preprocess_sources(src_dir, src_files, preprocess_dir)
    if output_dir.exists():
        shutil.rmtree(output_dir)
    translate_sources(preprocess_dir, preprocessed_files, output_dir)
    if preprocess_dir.exists():
        shutil.rmtree(preprocess_dir)
