import ast
import sys

if sys.version_info >= (3, 8):
    from importlib.metadata import version
else:
    from importlib_metadata import version


class ComprehensionChecker:
    """
    Flake8 plugin to help you write better list/set/dict comprehensions.
    """

    name = "flake8-comprehensions"
    version = version("flake8-comprehensions")

    def __init__(self, tree, *args, **kwargs):
        self.tree = tree

    messages = {
        "C400": "C400 Unnecessary generator - rewrite as a list comprehension.",
        "C401": "C401 Unnecessary generator - rewrite as a set comprehension.",
        "C402": "C402 Unnecessary generator - rewrite as a dict comprehension.",
        "C403": "C403 Unnecessary list comprehension - rewrite as a set comprehension.",
        "C404": (
            "C404 Unnecessary list comprehension - rewrite as a dict comprehension."
        ),
        "C405": "C405 Unnecessary {type} literal - ",
        "C406": "C406 Unnecessary {type} literal - ",
        "C407": "C407 Unnecessary list comprehension - '{func}' can take a generator.",
        "C408": "C408 Unnecessary {type} call - rewrite as a literal.",
        "C409": "C409 Unnecessary {type} passed to tuple() - ",
        "C410": "C410 Unnecessary {type} passed to list() - ",
        "C411": "C411 Unnecessary list call - remove the outer call to list().",
        "C412": "C412 Unnecessary {type} comprehension - 'in' can take a generator.",
        "C413": "C413 Unnecessary {outer} call around {inner}(){remediation}.",
        "C414": "C414 Unnecessary {inner} call within {outer}().",
        "C415": "C415 Unnecessary subscript reversal of iterable within {func}().",
        "C416": "C416 Unnecessary {type} comprehension - rewrite using {type}().",
    }

    def run(self):
        for node in ast.walk(self.tree):
            if isinstance(node, ast.Call) and isinstance(node.func, ast.Name):
                num_positional_args = len(node.args)

                if (
                    num_positional_args == 1
                    and isinstance(node.args[0], ast.GeneratorExp)
                    and node.func.id in ("list", "set")
                ):
                    msg_key = {"list": "C400", "set": "C401"}[node.func.id]
                    yield (
                        node.lineno,
                        node.col_offset,
                        self.messages[msg_key],
                        type(self),
                    )

                elif (
                    num_positional_args == 1
                    and isinstance(node.args[0], (ast.GeneratorExp, ast.ListComp))
                    and isinstance(node.args[0].elt, ast.Tuple)
                    and len(node.args[0].elt.elts) == 2
                    and node.func.id == "dict"
                ):
                    if isinstance(node.args[0], ast.GeneratorExp):
                        msg = "C402"
                    else:
                        msg = "C404"
                    yield (
                        node.lineno,
                        node.col_offset,
                        self.messages[msg],
                        type(self),
                    )

                elif (
                    num_positional_args == 1
                    and isinstance(node.args[0], ast.ListComp)
                    and node.func.id in ("list", "set")
                ):
                    msg_key = {"list": "C411", "set": "C403"}[node.func.id]
                    yield (
                        node.lineno,
                        node.col_offset,
                        self.messages[msg_key],
                        type(self),
                    )

                elif num_positional_args == 1 and (
                    isinstance(node.args[0], ast.Tuple)
                    and node.func.id == "tuple"
                    or isinstance(node.args[0], ast.List)
                    and node.func.id == "list"
                ):
                    suffix = "remove the outer call to {func}()."
                    msg_key = {"tuple": "C409", "list": "C410"}[node.func.id]
                    msg = self.messages[msg_key] + suffix
                    yield (
                        node.lineno,
                        node.col_offset,
                        msg.format(
                            type=type(node.args[0]).__name__.lower(), func=node.func.id
                        ),
                        type(self),
                    )

                elif (
                    num_positional_args == 1
                    and isinstance(node.args[0], (ast.Tuple, ast.List))
                    and node.func.id in ("tuple", "list", "set", "dict")
                ):
                    suffix = "rewrite as a {func} literal."
                    msg_key = {
                        "tuple": "C409",
                        "list": "C410",
                        "set": "C405",
                        "dict": "C406",
                    }[node.func.id]
                    msg = self.messages[msg_key] + suffix
                    yield (
                        node.lineno,
                        node.col_offset,
                        msg.format(
                            type=type(node.args[0]).__name__.lower(), func=node.func.id
                        ),
                        type(self),
                    )

                elif (
                    num_positional_args == 1
                    # These take 1 positional argument + some keyword arguments
                    and (
                        node.func.id
                        in ("all", "any", "frozenset", "tuple", "max", "min", "sorted",)
                    )
                    and isinstance(node.args[0], ast.ListComp)
                ):

                    yield (
                        node.lineno,
                        node.col_offset,
                        self.messages["C407"].format(func=node.func.id),
                        type(self),
                    )

                elif (
                    num_positional_args in (1, 2)
                    # These can take a second positional argument
                    and (node.func.id in ("enumerate", "sum",))
                    and isinstance(node.args[0], ast.ListComp)
                ):

                    yield (
                        node.lineno,
                        node.col_offset,
                        self.messages["C407"].format(func=node.func.id),
                        type(self),
                    )

                elif (
                    num_positional_args == 2
                    and node.func.id == "filter"
                    and isinstance(node.args[1], ast.ListComp)
                ):

                    yield (
                        node.lineno,
                        node.col_offset,
                        self.messages["C407"].format(func=node.func.id),
                        type(self),
                    )

                elif (
                    num_positional_args >= 2
                    and node.func.id == "map"
                    and any(isinstance(a, ast.ListComp) for a in node.args[1:])
                ):
                    yield (
                        node.lineno,
                        node.col_offset,
                        self.messages["C407"].format(func=node.func.id),
                        type(self),
                    )

                elif (
                    num_positional_args == 0
                    and not has_star_args(node)
                    and not has_keyword_args(node)
                    and node.func.id in ("tuple", "list", "dict")
                ):
                    yield (
                        node.lineno,
                        node.col_offset,
                        self.messages["C408"].format(type=node.func.id),
                        type(self),
                    )

                elif (
                    node.func.id in {"list", "reversed"}
                    and num_positional_args > 0
                    and isinstance(node.args[0], ast.Call)
                    and isinstance(node.args[0].func, ast.Name)
                    and node.args[0].func.id == "sorted"
                ):
                    remediation = ""
                    if node.func.id == "reversed":
                        reverse_flag_value = False
                        for keyword in node.args[0].keywords:
                            if keyword.arg != "reverse":
                                continue
                            if isinstance(keyword.value, ast.NameConstant):
                                reverse_flag_value = keyword.value.value
                            elif isinstance(keyword.value, ast.Num):
                                reverse_flag_value = bool(keyword.value.n)
                            else:
                                # Complex value
                                reverse_flag_value = None

                        if reverse_flag_value is None:
                            remediation = " - toggle reverse argument to sorted()"
                        else:
                            remediation = " - use sorted(..., reverse={!r})".format(
                                not reverse_flag_value
                            )

                    msg = self.messages["C413"].format(
                        inner=node.args[0].func.id,
                        outer=node.func.id,
                        remediation=remediation,
                    )
                    yield (
                        node.lineno,
                        node.col_offset,
                        msg,
                        type(self),
                    )

                elif (
                    num_positional_args > 0
                    and isinstance(node.args[0], ast.Call)
                    and isinstance(node.args[0].func, ast.Name)
                    and (
                        (
                            node.func.id in {"set", "sorted"}
                            and node.args[0].func.id
                            in {"list", "reversed", "sorted", "tuple"}
                        )
                        or (
                            node.func.id in {"list", "tuple"}
                            and node.args[0].func.id in {"list", "tuple"}
                        )
                        or (node.func.id == "set" and node.args[0].func.id == "set")
                    )
                ):
                    yield (
                        node.lineno,
                        node.col_offset,
                        self.messages["C414"].format(
                            inner=node.args[0].func.id, outer=node.func.id
                        ),
                        type(self),
                    )

                elif (
                    node.func.id in {"reversed", "set", "sorted"}
                    and num_positional_args > 0
                    and isinstance(node.args[0], ast.Subscript)
                    and isinstance(node.args[0].slice, ast.Slice)
                    and node.args[0].slice.lower is None
                    and node.args[0].slice.upper is None
                    and isinstance(node.args[0].slice.step, ast.UnaryOp)
                    and isinstance(node.args[0].slice.step.op, ast.USub)
                    and isinstance(node.args[0].slice.step.operand, ast.Num)
                    and node.args[0].slice.step.operand.n == 1
                ):
                    yield (
                        node.lineno,
                        node.col_offset,
                        self.messages["C415"].format(func=node.func.id),
                        type(self),
                    )

            elif isinstance(node, ast.Compare):
                if (
                    len(node.ops) == 1
                    and isinstance(node.ops[0], ast.In)
                    and len(node.comparators) == 1
                    and isinstance(
                        node.comparators[0], (ast.DictComp, ast.ListComp, ast.SetComp)
                    )
                ):
                    yield (
                        node.lineno,
                        node.col_offset,
                        self.messages["C412"].format(
                            type=comp_type[node.comparators[0].__class__]
                        ),
                        type(self),
                    )

            elif isinstance(node, (ast.ListComp, ast.SetComp)):
                if (
                    len(node.generators) == 1
                    and not node.generators[0].ifs
                    and not is_async_generator(node.generators[0])
                    and (
                        isinstance(node.elt, ast.Name)
                        and isinstance(node.generators[0].target, ast.Name)
                        and node.elt.id == node.generators[0].target.id
                    )
                ):

                    yield (
                        node.lineno,
                        node.col_offset,
                        self.messages["C416"].format(type=comp_type[node.__class__]),
                        type(self),
                    )


def has_star_args(call_node):
    return any(isinstance(a, ast.Starred) for a in call_node.args)


def has_keyword_args(call_node):
    return any(k.arg is None for k in call_node.keywords)


comp_type = {
    ast.DictComp: "dict",
    ast.ListComp: "list",
    ast.SetComp: "set",
}


if sys.version_info >= (3, 6):

    def is_async_generator(node):
        return node.is_async


else:

    def is_async_generator(node):
        return False
