# Copyright 2000-2001, Dalke Scientific Software, LLC
# Distributed under the Biopython License Agreement (see the LICENSE file).

"""Optimize an expression tree

  - remove Group nodes with no name
  - merge successive Str, single character positive Any nodes and positive
        Literals
"""

import Expression

def skip_empty_group(exp):
    while isinstance(exp, Expression.Group) and exp.name is None:
        exp = exp.expression
    return exp

# The only thing this does is remove Groups with no names
def optimize_unnamed_groups_recursive(exp):
    if isinstance(exp, Expression.MaxRepeat) or \
       isinstance(exp, Expression.Group):
        exp.expression = skip_empty_group(exp.expression)
    elif isinstance(exp, Expression.ExpressionList):
        subexps = []
        for subexp in exp.expressions:
            subexp = skip_empty_group(subexp)
            subexps.append(subexp)
        exp.expressions = subexps

    if hasattr(exp, "expression"):
        optimize_unnamed_groups_recursive(exp.expression)
    elif hasattr(exp, "expressions"):
        for exp in exp.expressions:
            optimize_unnamed_groups_recursive(exp)

def optimize_unnamed_groups(exp):
    """return an equivalent expression tree but without unnamed groups

    WARNING: has side-effect
    """

    optimize_unnamed_groups_recursive(exp)
    return skip_empty_group(exp)

def is_mergeable(exp):
    return isinstance(exp, Expression.Str) or \
           (isinstance(exp, Expression.Literal) and exp.invert == 0) or \
           (isinstance(exp, Expression.Any) and len(exp.chars) == 1 and \
            exp.invert == 0)
def get_merge_text(exp):
    if isinstance(exp, Expression.Str):
        return exp.string
    elif isinstance(exp, Expression.Literal):
        assert exp.invert == 0
        return exp.char
    elif isinstance(exp, Expression.Any):
        assert len(exp.chars) == 1 and exp.invert == 0
        return exp.chars
    raise AssertionError, "unknown exp: %s" % repr(exp)

def merge_strings(exp):
    """merge successive strings and string-like terms into a single Str

    WARNING: has side-effects
    """
    # Merge all of the children first - might promote a
    # Str("A") + (Str("B") + Str("C")) into an Str("A") + Str("BC")
    # before we do the merge into Str("ABC")
    
    if hasattr(exp, "expression"):
        exp.expression = merge_strings(exp.expression)
    elif hasattr(exp, "expressions"):
        subexps = []
        for subexp in exp.expressions:
            subexps.append(merge_strings(subexp))
        exp.expressions = subexps

    if isinstance(exp, Expression.Seq):
        all_strings = 1
        subexps = []
        for subexp in exp.expressions:
            if not subexps:
                subexps.append(subexp)
            elif not is_mergeable(subexps[-1]):
                all_strings = 0
                subexps.append(subexp)
            elif not is_mergeable(subexp):
                all_strings = 0
                subexps.append(subexp)
            else:
                # Previous and current are mergeable
                subexps[-1] = Expression.Str(get_merge_text(subexps[-1]) + \
                                             get_merge_text(subexp))
        if all_strings and subexps:
            assert len(subexps) == 1
            return subexps[0]
    
    return exp
    
def optimize(exp):
    """expression tree -> optimized expression tree

    Apply various optimizations to the expression tree.
    """
    
    exp = exp.copy()
    exp = optimize_unnamed_groups(exp)
    exp = merge_strings(exp)
    return exp
