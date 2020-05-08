"""Wordlists loaded from package data.

We can treat them as part of the code for the imperative mood check, and
therefore we load them at import time, rather than on-demand.

"""
import re
import pkgutil
import snowballstemmer
from typing import Iterator, Dict, Set


#: Regular expression for stripping comments from the wordlists
COMMENT_RE = re.compile(r'\s*#.*')

#: Stemmer function for stemming words in English
stem = snowballstemmer.stemmer('english').stemWord


def load_wordlist(name: str) -> Iterator[str]:
    """Iterate over lines of a wordlist data file.

    `name` should be the name of a package data file within the data/
    directory.

    Whitespace and #-prefixed comments are stripped from each line.

    """
    data = pkgutil.get_data('pydocstyle', 'data/' + name)
    if data is not None:
        text = data.decode('utf8')
        for line in text.splitlines():
            line = COMMENT_RE.sub('', line).strip()
            if line:
                yield line


#: A dict mapping stemmed verbs to the imperative form
def make_imperative_verbs_dict(wordlist: Iterator[str]) -> Dict[str, Set[str]]:
    imperative_verbs = {}  # type: Dict[str, Set[str]]
    for word in wordlist:
        imperative_verbs.setdefault(stem(word), set()).add(word)
    return imperative_verbs


IMPERATIVE_VERBS = make_imperative_verbs_dict(load_wordlist('imperatives.txt'))

#: Words that are forbidden to appear as the first word in a docstring
IMPERATIVE_BLACKLIST = set(load_wordlist('imperatives_blacklist.txt'))
