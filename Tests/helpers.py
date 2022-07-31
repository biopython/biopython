"""Tests helpers."""
import functools
import os
import tempfile


def temporary_directory():
    """Generate temporary directory using templib."""

    def decorator(fn):
        @functools.wraps(fn)
        def wrapper(*args, **kwargs):
            with tempfile.TemporaryDirectory(
                dir=os.getcwd()
            ) as temporary_directory_path:
                kwargs.update({"temporary_directory": temporary_directory_path})
                return fn(*args, **(kwargs or {}))

        return wrapper

    return decorator
