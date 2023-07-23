#!/usr/bin/env python
r"""Python script to convert our LaTeX bibitem entries to rST.

Input argument is our Tutorial.tex file only (default stdin), output to stdout.

Looks for the \bibitem entries.
"""

import sys


def fix_lines(lines):
    r"""Loop over lines applying simple state tracking.

    e.g. Looks for this:

        \bibitem{waterman1987}
        Michael S. Waterman, Mark Eggert: ``A new algorithm for best subsequence alignments with application to tRNA-rRNA comparisons'', \textit{Journal of Molecular Biology} {\bf 197} (4): 723--728 (1987). \url{https://doi.org/10.1016/0022-2836(87)90478-5}

    Yields this:

        .. [waterman1987] Michael S. Waterman, Mark Eggert: ``A new algorithm for best subsequence alignments with application to tRNA-rRNA comparisons'', *Journal of Molecular Biology* **197** (4): 723--728 (1987). https://doi.org/10.1016/0022-2836(87)90478-5
    """  # noqa: RST214
    item = False
    for line in lines:
        if line.strip().startswith(r"\bibitem{"):
            if item:
                yield item + "\n"
            item = line.strip()[9:].split("}", 1)[0]
            if item.startswith("dehoon"):
                item = "deHoon" + item[6:]
            else:
                item = item.title()
            item = f".. [{item}]"
        elif line.strip().startswith(r"\end{thebibliography}"):
            if item:
                yield item + "\n"
            item = False
        elif item:
            line = line.replace("``", "").replace("''", "")
            line = line.replace(r"\textit{", "*")
            line = line.replace(r"{\it ", "*")
            if r"{\bf " in line:
                line = line.replace(r"} {\bf ", "* **")
                line = line.replace("} ", "** ")
                line = line.replace("}:", "**: ")
                line = line.replace("},", "**, ")
            else:
                line = line.replace("}", "*")
            line = line.replace(r"\url{", "").replace("}", "")
            line = line.replace("  ", " ")
            item += " " + line.strip()
        else:
            # sys.stderr.write("DEBUG: Ignoring this: " + line)
            pass


for f in sys.argv[1:]:
    if f == "-":
        sys.stderr.write("Parsing LaTeX on stdin\n")
        for line in fix_lines(sys.stdin):
            sys.stdout.write(line)
    else:
        sys.stderr.write("Fixing %s\n" % f)
        with open(f) as handle:
            for line in fix_lines(handle):
                sys.stdout.write(line)
