
try:
    import pycodestyle
except ImportError:
    import pep8 as pycodestyle
import re

__version__ = '0.1.1'

BLIND_EXCEPT_REGEX = re.compile(r'(except:)')  # noqa


def check_blind_except(physical_line):
    if pycodestyle.noqa(physical_line):
        return
    match = BLIND_EXCEPT_REGEX.search(physical_line)
    if match:
        return match.start(), 'B901 blind except: statement'

check_blind_except.name = 'flake8-blind-except'
check_blind_except.version = __version__
