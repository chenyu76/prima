#!/usr/bin/env python3
"""Postprocess the translated MATLAB package."""

from pathlib import Path
import re

DEBUG_STATE_FUNCTION = """function previous = debug_state(varargin)
    % DEBUG_STATE Get or set debug mode, returning the previous mode.
    persistent debugging
    % Persistent variables start as []; test emptiness, not existence.
    if isempty(debugging)
        debugging = false;
    end
    previous = debugging;
    if nargin > 0
        debugging = varargin{1};
    end
end
"""


def postprocess(package_name):
    generate_debug_state_function(package_name)
    replace_debugging_references(package_name)
    replace_maxpow10_references(package_name)


def replace_maxpow10_references(package_name):
    transform_package_files(
        package_name,
        lambda text: replace_variable_in_blocks(
            text,
            "consts_obj.MAXPOW10",
            "MAXPOW10 = floor(log10(realmax)) - 1;",
        ),
    )


def replace_debugging_references(package_name):
    transform_package_files(
        package_name,
        lambda text: replace_variable(
            text,
            "const_obj.DEBUGGING",
            f"{package_name}.common.debug_state()",
        ),
    )


def transform_package_files(package_name, transform):
    """Apply a text transform to the package's MATLAB files."""
    package_dir = (
        Path(__file__).parent / "matlab" / "interfaces" / f"+{package_name}"
    )
    for path in package_dir.rglob("*.m"):
        with path.open(encoding="utf-8", newline="") as source:
            text = source.read()
        updated = transform(text)
        if updated != text:
            with path.open("w", encoding="utf-8", newline="") as target:
                target.write(updated)


def replace_variable(text, variable, replacement):
    """Replace all complete variable references, ignoring case."""
    pattern = re.compile(r"(?<![\w.])" + re.escape(variable) + r"(?!\w)", re.I)
    return pattern.sub(lambda match: replacement, text)


def replace_variable_in_blocks(text, variable, statement):
    """Replace a variable with a block-local assignment's LHS.

    Insert the assignment before the first complete statement using the variable,
    once per indent block.
    Nested blocks have independent assignments; blank lines and comments do not
    change blocks. MATLAB ``...`` continuation lines belong to their first line.
    """
    pattern = re.compile(r"(?<![\w.])" + re.escape(variable) + r"(?!\w)", re.I)
    replacement = statement.split("=", 1)[0].strip()
    continuation = re.compile(r"\.\.\.\s*(?:%.*)?$")
    newline = "\r\n" if "\r\n" in text else "\n"
    lines = text.splitlines(keepends=True)
    blocks = []
    result = []
    index = 0
    while index < len(lines):
        line = lines[index]
        if not line.strip() or line.lstrip().startswith("%"):
            result.append(pattern.sub(replacement, line))
            index += 1
            continue

        indent = line[: len(line) - len(line.lstrip(" \t"))]
        width = len(indent.expandtabs())
        while blocks and blocks[-1][0] > width:
            blocks.pop()
        if not blocks or blocks[-1][0] < width:
            blocks.append((width, False))

        end = index
        while end < len(lines):
            current = lines[end]
            end += 1
            if not current.strip() or current.lstrip().startswith("%"):
                continue
            if not continuation.search(current.rstrip("\r\n")):
                break
        original = "".join(lines[index:end])
        if original.strip() == statement.strip():
            blocks[-1] = (width, True)
        if pattern.search(original) and not blocks[-1][1]:
            result.append(indent + statement.strip() + newline)
            blocks[-1] = (width, True)
        result.append(pattern.sub(replacement, original))
        index = end
    return "".join(result)


def generate_debug_state_function(package_name):
    debug_state_path = (
        Path(__file__).parent
        / "matlab"
        / "interfaces"
        / f"+{package_name}"
        / "+common"
        / "debug_state.m"
    )
    debug_state_path.parent.mkdir(parents=True, exist_ok=True)
    debug_state_path.write_text(DEBUG_STATE_FUNCTION, encoding="utf-8")


if __name__ == "__main__":
    postprocess("prima_mat")
