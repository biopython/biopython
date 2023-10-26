# Copyright 2010 by Eric Talevich. All rights reserved.
# Copyright 2012 by Wibowo Arindrarto. All rights reserved.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
"""Common utility functions for various Bio submodules."""


import os
from typing import Any, TypeVar, Callable, Optional, cast
from typing import Protocol

# workaround type checking method attributes from https://github.com/python/mypy/issues/2087#issuecomment-587741762

F = TypeVar("F", bound=Callable[..., object])


class _FunctionWithPrevious(Protocol[F]):
    previous: Optional[int]
    __call__: F


def function_with_previous(func: F) -> _FunctionWithPrevious[F]:
    """Decorate a function as having an attribute named 'previous'."""
    function_with_previous = cast(_FunctionWithPrevious[F], func)
    # Make sure the cast isn't a lie.
    function_with_previous.previous = None
    return function_with_previous


def find_test_dir(start_dir: Optional[str] = None) -> str:
    """Find the absolute path of Biopython's Tests directory.

    Arguments:
    start_dir -- Initial directory to begin lookup (default to current dir)

    If the directory is not found up the filesystem's root directory, an
    exception will be raised.

    """
    if not start_dir:
        # no callbacks in function signatures!
        # defaults to the current directory
        # (using __file__ would give the installed Biopython)
        start_dir = "."

    target = os.path.abspath(start_dir)
    while True:
        if os.path.isdir(os.path.join(target, "Bio")) and os.path.isdir(
            os.path.join(target, "Tests")
        ):
            # Good, we're in the Biopython root now
            return os.path.abspath(os.path.join(target, "Tests"))
        # Recurse up the tree
        # TODO - Test this on Windows
        new, tmp = os.path.split(target)
        if target == new:
            # Reached root
            break
        target = new
    raise ValueError(
        f"Not within Biopython source tree: {os.path.abspath(start_dir)!r}"
    )


def run_doctest(target_dir: Optional[str] = None, *args: Any, **kwargs: Any) -> None:
    """Run doctest for the importing module."""
    import doctest

    # default doctest options
    default_kwargs = {"optionflags": doctest.ELLIPSIS}
    kwargs.update(default_kwargs)

    cur_dir = os.path.abspath(os.curdir)

    print("Running doctests...")
    try:
        os.chdir(find_test_dir(target_dir))
        doctest.testmod(*args, **kwargs)
    finally:
        # and revert back to initial directory
        os.chdir(cur_dir)
    print("Done")


if __name__ == "__main__":
    run_doctest()
