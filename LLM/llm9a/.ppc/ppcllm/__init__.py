import os.path
from typing import List, Optional
from ppcgrader.compiler import Compiler
import ppcgrader.config
from os import path
from urllib.request import urlretrieve


class Config(ppcgrader.config.Config):
    def __init__(self):
        super().__init__(binary='llm', cfg_file=__file__, openmp=True)

    def common_flags(self, compiler: Compiler) -> Compiler:
        compiler = super().common_flags(compiler)
        compiler = compiler.add_source(
            path.join(path.dirname(__file__), 'ref.cc'))
        return compiler

    def demo_flags(self, compiler: Compiler) -> Compiler:
        return self.common_flags(compiler)

    def demo_command(self, args: List[str]) -> List[str]:
        if len(args) == 0 or len(args) == 1:
            # download model and tokenizer
            if not os.path.exists("stories15M.bin"):
                urlretrieve(
                    "https://huggingface.co/karpathy/tinyllamas/resolve/main/stories15M.bin",
                    "stories15M.bin")
            if not os.path.exists("tokenizer.bin"):
                urlretrieve(
                    "https://github.com/karpathy/llama2.c/raw/master/tokenizer.bin",
                    "tokenizer.bin")
            text = "Once upon a time, there was a test story" if len(
                args) == 0 else args[0]
            args = ["stories15M.bin", "tokenizer.bin", text]
        return [os.path.join('./', self.demo_binary)] + args

    def parse_output(self, output):
        input_data = {"seed": None, "length": None, "config": {}}
        output_data = {}
        output_errors = {
            "max_error": None,
            "location_pos": None,
            "location_tok": None,
        }
        statistics = {}

        for line in output.splitlines():
            splitted = line.split('\t')
            if splitted[0] == 'result':
                errors = {
                    'fail': True,
                    'pass': False,
                    'done': False
                }[splitted[1]]
            elif splitted[0] == 'time':
                time = float(splitted[1])
            elif splitted[0] == 'perf_wall_clock_ns':
                time = int(splitted[1]) / 1e9
                statistics[splitted[0]] = int(splitted[1])
            elif splitted[0].startswith('perf_'):
                statistics[splitted[0]] = int(splitted[1])
            elif splitted[0] in ['max_error', 'threshold']:
                output_errors[splitted[0]] = float(splitted[1])
            elif splitted[0] in ['location_tok', 'location_pos']:
                output_errors[splitted[0]] = int(splitted[1])
            elif splitted[0] in ['seed', 'length']:
                input_data[splitted[0]] = int(splitted[1])
            elif splitted[0].startswith("config."):
                input_data["config"][splitted[0][7:]] = int(splitted[1])

        return time, errors, input_data, output_data, output_errors, statistics

    def explain_terminal(self, output, color=False) -> Optional[str]:
        from .info import explain_terminal
        return explain_terminal(output, color)
