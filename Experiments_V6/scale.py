"""Scale Newick branch lengths and save the result without changing tree labels."""
import argparse
from decimal import Decimal, InvalidOperation
from pathlib import Path
import re

# Skip quoted labels and comments so colons inside them are not branch lengths.
_TOKENS = re.compile(
    r"'(?:[^']|'')*'|\[[^\]]*\]|"
    r'(?P<prefix>:\s*)(?P<length>[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?)'
)


def scale_branch_lengths(newick_tree, scale_factor):
    """Multiply explicit lengths, including scientific notation; preserve other text."""
    try:
        factor = Decimal(str(scale_factor))
    except InvalidOperation as exc:
        raise ValueError('Scale must be a finite, nonnegative number') from exc
    if not factor.is_finite() or factor < 0:
        raise ValueError('Scale must be a finite, nonnegative number')
    if not newick_tree.strip() or not newick_tree.rstrip().endswith(';'):
        raise ValueError('Newick input must end with a semicolon')

    def replace(match):
        if match.group('length') is None:
            return match.group(0)
        value = Decimal(match.group('length')) * factor
        return match.group('prefix') + str(value)

    return _TOKENS.sub(replace, newick_tree)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument('--newick', help='Quoted Newick string')
    source.add_argument('--input', type=Path, help='File containing Newick tree(s)')
    parser.add_argument('--scale', required=True, help='Nonnegative branch-length multiplier')
    parser.add_argument('--output', type=Path,
                        help='Default: <input-stem>_scale.tre, or files_scale.tre for a string')
    args = parser.parse_args()
    output = args.output or (args.input.with_name(args.input.stem + f'_{args.scale}.tre')
                             if args.input else Path('files_scale.tre'))
    if args.input and args.input.resolve() == output.resolve():
        parser.error('Output must differ from the input file')
    try:
        newick = args.input.read_text() if args.input else args.newick
        scaled = scale_branch_lengths(newick, args.scale)
        output.write_text(scaled.rstrip() + '\n')
    except (OSError, ValueError) as exc:
        parser.error(str(exc))
    print(f'Saved scaled tree(s) to {output}')


if __name__ == '__main__':
    main()
