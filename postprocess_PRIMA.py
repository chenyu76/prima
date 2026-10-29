#!/usr/bin/env python3
"""Add runtime debugging methods to the translated package."""

from pathlib import Path

DEBUG_METHOD = """
    methods (Static)
        function value = debug_state(varargin)
            persistent debugging
            if isempty(debugging)
                debugging = false;
            end
            if nargin > 0
                debugging = varargin{1};
            end
            value = debugging;
        end
    end
"""


def postprocess(package_name):
    consts_mod_path = (
        Path(__file__).parent
        / "matlab"
        / "interfaces"
        / f"+{package_name}"
        / "+common"
        / "consts_mod.m"
    )
    text = consts_mod_path.read_text(encoding="utf-8")
    text = text.replace(
        "obj.DEBUGGING = true;",
        f"obj.DEBUGGING = {package_name}.common.consts_mod.debug_state();",
    )
    text = text.replace(
        "    methods\n",
        DEBUG_METHOD.replace("{package_name}", package_name) + "    methods\n",
        1,
    )
    consts_mod_path.write_text(text, encoding="utf-8")


if __name__ == "__main__":
    postprocess("prima_mat")
