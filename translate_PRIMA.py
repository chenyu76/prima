#!/bin/python
import os
import re
import shutil
import subprocess
import time
from pathlib import Path
from functools import wraps
import datetime


def wrap_as_comment(text: str) -> str:
    lines = text.split("\n")
    max_col = max([len(line) for line in lines])
    lines = (
        ["%-" + "-" * max_col + "-%"]
        + ["% " + line + " " * (max_col - len(line)) + " %" for line in lines]
        + ["%-" + "-" * max_col + "-%"]
    )
    return "\n".join(lines)


def format_date_and_commit(text: str) -> str:
    today = datetime.date.today().strftime("%Y-%m-%d")
    commit = subprocess.run(
        ["git", "rev-parse", "--short=8", "main"],
        capture_output=True,
        text=True,
        check=True,
    ).stdout.strip()
    return text.format(date=today, commit=commit)


FILE_PREAMBLE = wrap_as_comment(format_date_and_commit("""
This file is translated from Fortran by 4ft2pm on {date}.
4ft2pm is a Fortran-to-Matlab translator by CHEN Yu and ZHANG Zaikun.
The fortran version is from PRIMA (https://libprima.net) with git commit {commit}.
""".strip()))

FOLDER_README = format_date_and_commit("""
This is the MATLAB of ??????
Translated from Fortran by 4ft2pm on {date}.
4ft2pm is a Fortran-to-Matlab translator by CHEN Yu and ZHANG Zaikun.
The fortran version is from PRIMA (https://libprima.net) with git commit {commit}.
""".strip())


def timer(func):
    @wraps(func)
    def wrapper(*args, **kwargs):
        print(func.__name__, end="", flush=True)
        start = time.perf_counter()
        result = func(*args, **kwargs)
        end = time.perf_counter()
        print(f": {end - start:.4f}s")
        return result

    return wrapper


@timer
def preprocess_sources(src_dir, src_files, output_dir):
    preprocessed_files = []
    for file_abs_path in src_files:
        file_rel_path = Path(os.path.relpath(file_abs_path, start=src_dir))
        source_path = src_dir / file_rel_path
        target_path = output_dir / file_rel_path
        preprocessed_files.append(target_path)

        # Create target subdirectories automatically
        target_path.parent.mkdir(parents=True, exist_ok=True)

        # .F90 -> Preprocess with gfortran
        if source_path.suffix == ".F90":
            subprocess.run(
                [
                    "gfortran",
                    "-E",
                    "-cpp",
                    "-P",
                    # "-DPRIMA_DEBUGGING=1",
                    str(source_path),
                    f"-I{src_dir}",
                    "-o",
                    str(target_path),
                ],
                check=True,
            )
            content = target_path.read_text(encoding="utf-8", errors="ignore")
            # .strip() removes all leading/trailing whitespace and newlines
            # We append a single "\n" at the end as per standard POSIX
            # file conventions
            target_path.write_text(content.strip() + "\n", encoding="utf-8")
        else:
            shutil.copy2(source_path, target_path)

    return preprocessed_files


@timer
def extract_source_files(src_dir):
    """Extracts source file paths from CMake."""
    cmake_file = src_dir / "CMakeLists.txt"

    # Parse CMakeLists.txt for source files
    content = cmake_file.read_text(encoding="utf-8", errors="ignore")
    pattern = re.compile(r"add_library\s*\((.*?)\)", re.DOTALL | re.IGNORECASE)
    matches = pattern.findall(content)

    # Extract file tokens and filter out keywords/variables
    raw_tokens = matches[0].split()
    keywords = {"SHARED", "STATIC", "INTERFACE", "MODULE", "OBJECT", "ALIAS"}
    return [
        src_dir / Path(t.strip("'\""))
        for t in raw_tokens[1:]
        if t.upper() not in keywords and not t.startswith("$")
    ]


@timer
def translate_sources(
    src_dir, preprocessed_files, output_dir, pkg_name, translator_exec="4ft2pm"
):
    return subprocess.run(
        [
            translator_exec,
            "-r",
            str(src_dir),
            "-o",
            str(output_dir),
            "--create-setup-m",
            # "--try-bit-consistency",
            # "--no-simplify",
            "--prima",
            "--as-package",
            pkg_name,
            # "--preamble",
            # FILE_PREAMBLE,
            "--max-column-width",
            "100",
        ]
    )


if __name__ == "__main__":
    current_dir = Path(__file__).parent
    src_dir = (current_dir / "fortran").resolve()
    preprocess_dir = (current_dir / "preprocessed_fortran/").resolve()
    output_dir = (current_dir / "matlab/interfaces/").resolve()
    pkg_name = "prima_mat"

    for p in [output_dir / "+fortran", output_dir / f"+{pkg_name}"]:
        if p.exists():
            shutil.rmtree(p)

    src_files = extract_source_files(src_dir)
    preprocessed_files = preprocess_sources(src_dir, src_files, preprocess_dir)
    translate_sources(preprocess_dir, preprocessed_files, output_dir, pkg_name)

    if preprocess_dir.exists():
        shutil.rmtree(preprocess_dir)
