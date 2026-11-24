#!/bin/python
import os
import re
import shutil
import subprocess
import time
import tempfile
from pathlib import Path
from functools import wraps
import datetime

from postprocess_PRIMA import postprocess


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
This file is translated from Fortran by 42mp on {date}.
42mp is a Fortran-to-Matlab translator by CHEN Yu and ZHANG Zaikun.
The fortran version is from PRIMA (https://libprima.net) with git commit {commit}.
""".strip()))

FILE_PREAMBLE = wrap_as_comment("Placeholder")

FOLDER_README = format_date_and_commit("""
This is the MATLAB of ??????
Translated from Fortran by 42mp on {date}.
42mp is a Fortran-to-Matlab translator by CHEN Yu and ZHANG Zaikun.
The fortran version is from PRIMA (https://libprima.net) with git commit {commit}.
""".strip())


def timer(func):
    @wraps(func)
    def wrapper(*args, **kwargs):
        print(func.__name__.replace("_", " "), end="", flush=True)
        start = time.perf_counter()
        result = func(*args, **kwargs)
        end = time.perf_counter()
        print(f": {end - start:.4f}s")
        return result

    return wrapper


class FortranCommentProtector:
    """
    Keep Fortran comments unchanged while CPP processes the code.

    A line in const.F90 includes `PRIMA_DEBUGGING` as a part of comments.
    This class prevent it being replaced by gfortran.
    """

    def __init__(self):
        self.comments = {}
        self.counter = 0

    def protect(self, content):
        """
        Hide free-form Fortran comments, respecting quoted strings.
        """
        lines = []
        quote = None
        for line in content.splitlines(keepends=True):
            # Do not treat the ! in #if !defined(...) as a comment.
            if line.lstrip().startswith("#"):
                lines.append(line)
                continue
            # Whole-line comments may also appear inside a continued string.
            comment_only = line.lstrip().startswith("!")
            for index, char in enumerate(line):
                if char == "!" and (comment_only or quote is None):
                    comment = line[index:].rstrip("\r\n")
                    token = f"!__PRIMA_COMMENT_{self.counter}__"
                    self.comments[token] = comment
                    self.counter += 1
                    line = line[:index] + token + line[index + len(comment) :]
                    break
                if comment_only:
                    continue
                if quote:
                    if char == quote:
                        quote = None
                elif char in ("'", '"'):
                    quote = char
            lines.append(line)
        return "".join(lines)

    def restore(self, content):
        """Replace placeholders with the original comments."""
        return re.sub(
            r"!__PRIMA_COMMENT_\d+__",
            lambda match: self.comments.get(match.group(), match.group()),
            content,
        )


@timer
def preprocess(src_dir, src_files, output_dir, debugging=0):
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
            comment_protector = FortranCommentProtector()
            protected = comment_protector.protect(
                source_path.read_text(encoding="utf-8", errors="ignore")
            )
            with tempfile.TemporaryDirectory(
                prefix="prima-preprocess-"
            ) as tmp:
                protected_path = Path(tmp) / source_path.name
                protected_path.write_text(protected, encoding="utf-8")
                subprocess.run(
                    [
                        "gfortran",
                        "-E",
                        "-cpp",
                        "-P",
                        f"-DPRIMA_DEBUGGING={debugging}",
                        str(protected_path),
                        f"-I{source_path.parent}",
                        f"-I{src_dir}",
                        "-o",
                        str(target_path),
                    ],
                    check=True,
                )
            content = target_path.read_text(encoding="utf-8", errors="ignore")
            content = comment_protector.restore(content)
            # .strip() removes all leading/trailing whitespace and newlines
            # We append a single "\n" at the end as per standard POSIX
            # file conventions
            target_path.write_text(content.strip() + "\n", encoding="utf-8")
        else:
            shutil.copy2(source_path, target_path)

    return preprocessed_files


@timer
def extract_source(src_dir):
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
def translate(options):
    return subprocess.run(
        [
            "42mp-cabal",
        ]
        + options
    )


if __name__ == "__main__":
    current_dir = Path(__file__).parent
    package_name = "prima_mat"
    src_dir = (current_dir / "fortran").resolve()
    preprocess_dir = (current_dir / "preprocessed_fortran/").resolve()
    output_dir = (current_dir / "matlab/interfaces/").resolve()
    p = output_dir / "+fortran"
    if p.exists():
        shutil.rmtree(p)
    (output_dir / "setup_mex.m").unlink(missing_ok=True)

    p = output_dir / f"+{package_name}"
    if p.exists():
        shutil.rmtree(p)

    src_files = extract_source(src_dir)
    preprocess(src_dir, src_files, preprocess_dir, debugging=1)

    prima_mode = False
    translate(
        [
            "-r",
            str(preprocess_dir),
            "-o",
            str(output_dir),
            "--create-setup-m",
            "--prima" if prima_mode else "--try-bit-consistency",
            # "--no-simplify",
            "--as-package",
            package_name,
            "--preamble",
            FILE_PREAMBLE,
            "--max-column-width",
            "100",
        ]
    )

    if preprocess_dir.exists():
        shutil.rmtree(preprocess_dir)

    if prima_mode:
        postprocess(package_name)
