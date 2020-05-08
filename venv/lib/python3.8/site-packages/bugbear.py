import ast
import builtins
from collections import namedtuple
from contextlib import suppress
from functools import lru_cache, partial
import itertools
from keyword import iskeyword
import logging
import re

import attr
import pycodestyle


__version__ = "20.1.4"

LOG = logging.getLogger("flake8.bugbear")


@attr.s(hash=False)
class BugBearChecker:
    name = "flake8-bugbear"
    version = __version__

    tree = attr.ib(default=None)
    filename = attr.ib(default="(none)")
    lines = attr.ib(default=None)
    max_line_length = attr.ib(default=79)
    visitor = attr.ib(init=False, default=attr.Factory(lambda: BugBearVisitor))
    options = attr.ib(default=None)

    def run(self):
        if not self.tree or not self.lines:
            self.load_file()
        visitor = self.visitor(filename=self.filename, lines=self.lines)
        visitor.visit(self.tree)
        for e in itertools.chain(visitor.errors, self.gen_line_based_checks()):
            if pycodestyle.noqa(self.lines[e.lineno - 1]):
                continue

            if self.should_warn(e.message[:4]):
                yield self.adapt_error(e)

    def gen_line_based_checks(self):
        """gen_line_based_checks() -> (error, error, error, ...)

        The following simple checks are based on the raw lines, not the AST.
        """
        for lineno, line in enumerate(self.lines, start=1):
            length = len(line) - 1
            if length > 1.1 * self.max_line_length:
                yield B950(lineno, length, vars=(length, self.max_line_length))

    @classmethod
    def adapt_error(cls, e):
        """Adapts the extended error namedtuple to be compatible with Flake8."""
        return e._replace(message=e.message.format(*e.vars))[:4]

    def load_file(self):
        """Loads the file in a way that auto-detects source encoding and deals
        with broken terminal encodings for stdin.

        Stolen from flake8_import_order because it's good.
        """

        if self.filename in ("stdin", "-", None):
            self.filename = "stdin"
            self.lines = pycodestyle.stdin_get_value().splitlines(True)
        else:
            self.lines = pycodestyle.readlines(self.filename)

        if not self.tree:
            self.tree = ast.parse("".join(self.lines))

    @staticmethod
    def add_options(optmanager):
        """Informs flake8 to ignore B9xx by default."""
        optmanager.extend_default_ignore(disabled_by_default)

    @lru_cache()
    def should_warn(self, code):
        """Returns `True` if Bugbear should emit a particular warning.

        flake8 overrides default ignores when the user specifies
        `ignore = ` in configuration.  This is problematic because it means
        specifying anything in `ignore = ` implicitly enables all optional
        warnings.  This function is a workaround for this behavior.

        As documented in the README, the user is expected to explicitly select
        the warnings.
        """
        if code[:2] != "B9":
            # Normal warnings are safe for emission.
            return True

        if self.options is None:
            LOG.info(
                "Options not provided to Bugbear, optional warning %s selected.", code
            )
            return True

        for i in range(2, len(code) + 1):
            if code[:i] in self.options.select:
                return True

        LOG.info(
            "Optional warning %s not present in selected warnings: %r. Not "
            "firing it at all.",
            code,
            self.options.select,
        )
        return False


def _is_identifier(arg):
    # Return True if arg is a valid identifier, per
    # https://docs.python.org/2/reference/lexical_analysis.html#identifiers

    if not isinstance(arg, ast.Str):
        return False

    return re.match(r"^[A-Za-z_][A-Za-z0-9_]*$", arg.s) is not None


def _to_name_str(node):
    # Turn Name and Attribute nodes to strings, e.g "ValueError" or
    # "pkg.mod.error", handling any depth of attribute accesses.
    if isinstance(node, ast.Name):
        return node.id
    assert isinstance(node, ast.Attribute)
    return _to_name_str(node.value) + "." + node.attr


@attr.s
class BugBearVisitor(ast.NodeVisitor):
    filename = attr.ib()
    lines = attr.ib()
    node_stack = attr.ib(default=attr.Factory(list))
    node_window = attr.ib(default=attr.Factory(list))
    errors = attr.ib(default=attr.Factory(list))
    futures = attr.ib(default=attr.Factory(set))

    NODE_WINDOW_SIZE = 4

    if False:
        # Useful for tracing what the hell is going on.

        def __getattr__(self, name):
            print(name)
            return self.__getattribute__(name)

    def visit(self, node):
        self.node_stack.append(node)
        self.node_window.append(node)
        self.node_window = self.node_window[-self.NODE_WINDOW_SIZE :]
        super().visit(node)
        self.node_stack.pop()

    def visit_ExceptHandler(self, node):
        if node.type is None:
            self.errors.append(
                B001(node.lineno, node.col_offset, vars=("bare `except:`",))
            )
        elif isinstance(node.type, ast.Tuple):
            names = [_to_name_str(e) for e in node.type.elts]
            as_ = " as " + node.name if node.name is not None else ""
            if len(names) == 0:
                vs = ("`except (){}:`".format(as_),)
                self.errors.append(B001(node.lineno, node.col_offset, vars=vs))
            elif len(names) == 1:
                self.errors.append(B013(node.lineno, node.col_offset, vars=names))
            else:
                # See if any of the given exception names could be removed, e.g. from:
                #    (MyError, MyError)  # duplicate names
                #    (MyError, BaseException)  # everything derives from the Base
                #    (Exception, TypeError)  # builtins where one subclasses another
                # but note that other cases are impractical to hande from the AST.
                # We expect this is mostly useful for users who do not have the
                # builtin exception hierarchy memorised, and include a 'shadowed'
                # subtype without realising that it's redundant.
                good = sorted(set(names), key=names.index)
                if "BaseException" in good:
                    good = ["BaseException"]
                for name, other in itertools.permutations(tuple(good), 2):
                    if issubclass(
                        getattr(builtins, name, type), getattr(builtins, other, ())
                    ):
                        if name in good:
                            good.remove(name)
                if good != names:
                    desc = good[0] if len(good) == 1 else "({})".format(", ".join(good))
                    self.errors.append(
                        B014(
                            node.lineno,
                            node.col_offset,
                            vars=(", ".join(names), as_, desc),
                        )
                    )
        self.generic_visit(node)

    def visit_UAdd(self, node):
        trailing_nodes = list(map(type, self.node_window[-4:]))
        if trailing_nodes == [ast.UnaryOp, ast.UAdd, ast.UnaryOp, ast.UAdd]:
            originator = self.node_window[-4]
            self.errors.append(B002(originator.lineno, originator.col_offset))
        self.generic_visit(node)

    def visit_Call(self, node):
        if isinstance(node.func, ast.Attribute):
            for bug in (B301, B302, B305):
                if node.func.attr in bug.methods:
                    call_path = ".".join(self.compose_call_path(node.func.value))
                    if call_path not in bug.valid_paths:
                        self.errors.append(bug(node.lineno, node.col_offset))
                    break
            else:
                self.check_for_b005(node)
        else:
            with suppress(AttributeError, IndexError):
                if (
                    node.func.id in ("getattr", "hasattr")
                    and node.args[1].s == "__call__"  # noqa: W503
                ):
                    self.errors.append(B004(node.lineno, node.col_offset))
                if (
                    node.func.id == "getattr"
                    and len(node.args) == 2  # noqa: W503
                    and _is_identifier(node.args[1])  # noqa: W503
                    and not iskeyword(node.args[1].s)  # noqa: W503
                ):
                    self.errors.append(B009(node.lineno, node.col_offset))
                elif (
                    node.func.id == "setattr"
                    and len(node.args) == 3  # noqa: W503
                    and _is_identifier(node.args[1])  # noqa: W503
                    and not iskeyword(node.args[1].s)  # noqa: W503
                ):
                    self.errors.append(B010(node.lineno, node.col_offset))

        self.generic_visit(node)

    def visit_Attribute(self, node):
        call_path = list(self.compose_call_path(node))
        if ".".join(call_path) == "sys.maxint":
            self.errors.append(B304(node.lineno, node.col_offset))
        elif len(call_path) == 2 and call_path[1] == "message":
            name = call_path[0]
            for elem in reversed(self.node_stack[:-1]):
                if isinstance(elem, ast.ExceptHandler) and elem.name == name:
                    self.errors.append(B306(node.lineno, node.col_offset))
                    break

    def visit_Assign(self, node):
        if isinstance(self.node_stack[-2], ast.ClassDef):
            # note: by hasattr below we're ignoring starred arguments, slices
            # and tuples for simplicity.
            assign_targets = {t.id for t in node.targets if hasattr(t, "id")}
            if "__metaclass__" in assign_targets:
                self.errors.append(B303(node.lineno, node.col_offset))
        elif len(node.targets) == 1:
            t = node.targets[0]
            if isinstance(t, ast.Attribute) and isinstance(t.value, ast.Name):
                if (t.value.id, t.attr) == ("os", "environ"):
                    self.errors.append(B003(node.lineno, node.col_offset))
        self.generic_visit(node)

    def visit_For(self, node):
        self.check_for_b007(node)
        self.generic_visit(node)

    def visit_Assert(self, node):
        self.check_for_b011(node)
        self.generic_visit(node)

    def visit_AsyncFunctionDef(self, node):
        self.check_for_b902(node)
        self.check_for_b006(node)
        self.generic_visit(node)

    def visit_FunctionDef(self, node):
        self.check_for_b901(node)
        self.check_for_b902(node)
        self.check_for_b006(node)
        self.generic_visit(node)

    def visit_ClassDef(self, node):
        self.check_for_b903(node)
        self.generic_visit(node)

    def visit_Try(self, node):
        self.check_for_b012(node)
        self.generic_visit(node)

    def compose_call_path(self, node):
        if isinstance(node, ast.Attribute):
            yield from self.compose_call_path(node.value)
            yield node.attr
        elif isinstance(node, ast.Name):
            yield node.id

    def check_for_b005(self, node):
        if node.func.attr not in B005.methods:
            return  # method name doesn't match

        if len(node.args) != 1 or not isinstance(node.args[0], ast.Str):
            return  # used arguments don't match the builtin strip

        call_path = ".".join(self.compose_call_path(node.func.value))
        if call_path in B005.valid_paths:
            return  # path is exempt

        s = node.args[0].s
        if len(s) == 1:
            return  # stripping just one character

        if len(s) == len(set(s)):
            return  # no characters appear more than once

        self.errors.append(B005(node.lineno, node.col_offset))

    def check_for_b006(self, node):
        for default in node.args.defaults + node.args.kw_defaults:
            if isinstance(default, B006.mutable_literals):
                self.errors.append(B006(default.lineno, default.col_offset))
            elif isinstance(default, ast.Call):
                call_path = ".".join(self.compose_call_path(default.func))
                if call_path in B006.mutable_calls:
                    self.errors.append(B006(default.lineno, default.col_offset))
                elif call_path not in B008.immutable_calls:
                    self.errors.append(B008(default.lineno, default.col_offset))

    def check_for_b007(self, node):
        targets = NameFinder()
        targets.visit(node.target)
        ctrl_names = set(filter(lambda s: not s.startswith("_"), targets.names))
        body = NameFinder()
        for expr in node.body:
            body.visit(expr)
        used_names = set(body.names)
        for name in sorted(ctrl_names - used_names):
            n = targets.names[name][0]
            self.errors.append(B007(n.lineno, n.col_offset, vars=(name,)))

    def check_for_b011(self, node):
        if isinstance(node.test, ast.NameConstant) and node.test.value is False:
            self.errors.append(B011(node.lineno, node.col_offset))

    def check_for_b012(self, node):
        def _loop(node, bad_node_types):
            if isinstance(node, (ast.AsyncFunctionDef, ast.FunctionDef)):
                return

            if isinstance(node, (ast.While, ast.For)):
                bad_node_types = (ast.Return,)

            elif isinstance(node, bad_node_types):
                self.errors.append(B012(node.lineno, node.col_offset))

            for child in ast.iter_child_nodes(node):
                _loop(child, bad_node_types)

        for child in node.finalbody:
            _loop(child, (ast.Return, ast.Continue, ast.Break))

    def walk_function_body(self, node):
        def _loop(parent, node):
            if isinstance(node, (ast.AsyncFunctionDef, ast.FunctionDef)):
                return
            yield parent, node
            for child in ast.iter_child_nodes(node):
                yield from _loop(node, child)

        for child in node.body:
            yield from _loop(node, child)

    def check_for_b901(self, node):
        if node.name == "__await__":
            return

        has_yield = False
        return_node = None

        for parent, x in self.walk_function_body(node):
            # Only consider yield when it is part of an Expr statement.
            if isinstance(parent, ast.Expr) and isinstance(
                x, (ast.Yield, ast.YieldFrom)
            ):
                has_yield = True

            if isinstance(x, ast.Return) and x.value is not None:
                return_node = x

            if has_yield and return_node is not None:
                self.errors.append(B901(return_node.lineno, return_node.col_offset))
                break

    def check_for_b902(self, node):
        if not isinstance(self.node_stack[-2], ast.ClassDef):
            return

        decorators = NameFinder()
        decorators.visit(node.decorator_list)

        if "staticmethod" in decorators.names:
            # TODO: maybe warn if the first argument is surprisingly `self` or
            # `cls`?
            return

        bases = {b.id for b in self.node_stack[-2].bases if isinstance(b, ast.Name)}
        if "type" in bases:
            if (
                "classmethod" in decorators.names
                or node.name in B902.implicit_classmethods  # noqa: W503
            ):
                expected_first_args = B902.metacls
                kind = "metaclass class"
            else:
                expected_first_args = B902.cls
                kind = "metaclass instance"
        else:
            if (
                "classmethod" in decorators.names
                or node.name in B902.implicit_classmethods  # noqa: W503
            ):
                expected_first_args = B902.cls
                kind = "class"
            else:
                expected_first_args = B902.self
                kind = "instance"

        args = node.args.args
        vararg = node.args.vararg
        kwarg = node.args.kwarg
        kwonlyargs = node.args.kwonlyargs

        if args:
            actual_first_arg = args[0].arg
            lineno = args[0].lineno
            col = args[0].col_offset
        elif vararg:
            actual_first_arg = "*" + vararg.arg
            lineno = vararg.lineno
            col = vararg.col_offset
        elif kwarg:
            actual_first_arg = "**" + kwarg.arg
            lineno = kwarg.lineno
            col = kwarg.col_offset
        elif kwonlyargs:
            actual_first_arg = "*, " + kwonlyargs[0].arg
            lineno = kwonlyargs[0].lineno
            col = kwonlyargs[0].col_offset
        else:
            actual_first_arg = "(none)"
            lineno = node.lineno
            col = node.col_offset

        if actual_first_arg not in expected_first_args:
            if not actual_first_arg.startswith(("(", "*")):
                actual_first_arg = repr(actual_first_arg)
            self.errors.append(
                B902(lineno, col, vars=(actual_first_arg, kind, expected_first_args[0]))
            )

    def check_for_b903(self, node):
        body = node.body
        if (
            body
            and isinstance(body[0], ast.Expr)  # noqa: W503
            and isinstance(body[0].value, ast.Str)  # noqa: W503
        ):
            # Ignore the docstring
            body = body[1:]

        if (
            len(body) != 1
            or not isinstance(body[0], ast.FunctionDef)  # noqa: W503
            or body[0].name != "__init__"  # noqa: W503
        ):
            # only classes with *just* an __init__ method are interesting
            return

        # all the __init__ function does is a series of assignments to attributes
        for stmt in body[0].body:
            if not isinstance(stmt, ast.Assign):
                return
            targets = stmt.targets
            if len(targets) > 1 or not isinstance(targets[0], ast.Attribute):
                return
            if not isinstance(stmt.value, ast.Name):
                return

        self.errors.append(B903(node.lineno, node.col_offset))


@attr.s
class NameFinder(ast.NodeVisitor):
    """Finds a name within a tree of nodes.

    After `.visit(node)` is called, `found` is a dict with all name nodes inside,
    key is name string, value is the node (useful for location purposes).
    """

    names = attr.ib(default=attr.Factory(dict))

    def visit_Name(self, node):
        self.names.setdefault(node.id, []).append(node)

    def visit(self, node):
        """Like super-visit but supports iteration over lists."""
        if not isinstance(node, list):
            return super().visit(node)

        for elem in node:
            super().visit(elem)
        return node


error = namedtuple("error", "lineno col message type vars")
Error = partial(partial, error, type=BugBearChecker, vars=())


B001 = Error(
    message="B001 Do not use {}, it also catches unexpected "
    "events like memory errors, interrupts, system exit, and so on.  "
    "Prefer `except Exception:`.  If you're sure what you're doing, "
    "be explicit and write `except BaseException:`."
)

B002 = Error(
    message="B002 Python does not support the unary prefix increment. Writing "
    "++n is equivalent to +(+(n)), which equals n. You meant n += 1."
)

B003 = Error(
    message="B003 Assigning to `os.environ` doesn't clear the environment. "
    "Subprocesses are going to see outdated variables, in disagreement "
    "with the current process. Use `os.environ.clear()` or the `env=` "
    "argument to Popen."
)

B004 = Error(
    message="B004 Using `hasattr(x, '__call__')` to test if `x` is callable "
    "is unreliable. If `x` implements custom `__getattr__` or its "
    "`__call__` is itself not callable, you might get misleading "
    "results. Use `callable(x)` for consistent results."
)

B005 = Error(
    message="B005 Using .strip() with multi-character strings is misleading "
    "the reader. It looks like stripping a substring. Move your "
    "character set to a constant if this is deliberate. Use "
    ".replace() or regular expressions to remove string fragments."
)
B005.methods = {"lstrip", "rstrip", "strip"}
B005.valid_paths = {}

B006 = Error(
    message="B006 Do not use mutable data structures for argument defaults.  They "
    "are created during function definition time. All calls to the function "
    "reuse this one instance of that data structure, persisting changes "
    "between them."
)
B006.mutable_literals = (ast.Dict, ast.List, ast.Set)
B006.mutable_calls = {
    "Counter",
    "OrderedDict",
    "collections.Counter",
    "collections.OrderedDict",
    "collections.defaultdict",
    "collections.deque",
    "defaultdict",
    "deque",
    "dict",
    "list",
    "set",
}
B007 = Error(
    message="B007 Loop control variable {!r} not used within the loop body. "
    "If this is intended, start the name with an underscore."
)
B008 = Error(
    message="B008 Do not perform function calls in argument defaults.  The call is "
    "performed only once at function definition time. All calls to your "
    "function will reuse the result of that definition-time function call.  If "
    "this is intended, assign the function call to a module-level variable and "
    "use that variable as a default value."
)
B008.immutable_calls = {"tuple", "frozenset"}
B009 = Error(
    message="B009 Do not call getattr with a constant attribute value, "
    "it is not any safer than normal property access."
)
B010 = Error(
    message="B010 Do not call setattr with a constant attribute value, "
    "it is not any safer than normal property access."
)
B011 = Error(
    message="B011 Do not call assert False since python -O removes these calls. "
    "Instead callers should raise AssertionError()."
)
B012 = Error(
    message="B012 return/continue/break inside finally blocks cause exceptions "
    "to be silenced. Exceptions should be silenced in except blocks. Control "
    "statements can be moved outside the finally block."
)
B013 = Error(
    message="B013 A length-one tuple literal is redundant.  "
    "Write `except {0}:` instead of `except ({0},):`."
)
B014 = Error(
    message="B014 Redundant exception types in `except ({0}){1}:`.  "
    "Write `except {2}{1}:`, which catches exactly the same exceptions."
)

# Those could be false positives but it's more dangerous to let them slip
# through if they're not.
B301 = Error(
    message="B301 Python 3 does not include `.iter*` methods on dictionaries. "
    "Remove the `iter` prefix from the method name. For Python 2 "
    "compatibility, prefer the Python 3 equivalent unless you expect "
    "the size of the container to be large or unbounded. Then use "
    "`six.iter*` or `future.utils.iter*`."
)
B301.methods = {"iterkeys", "itervalues", "iteritems", "iterlists"}
B301.valid_paths = {"six", "future.utils", "builtins"}

B302 = Error(
    message="B302 Python 3 does not include `.view*` methods on dictionaries. "
    "Remove the `view` prefix from the method name. For Python 2 "
    "compatibility, prefer the Python 3 equivalent unless you expect "
    "the size of the container to be large or unbounded. Then use "
    "`six.view*` or `future.utils.view*`."
)
B302.methods = {"viewkeys", "viewvalues", "viewitems", "viewlists"}
B302.valid_paths = {"six", "future.utils", "builtins"}

B303 = Error(
    message="B303 `__metaclass__` does nothing on Python 3. Use "
    "`class MyClass(BaseClass, metaclass=...)`. For Python 2 "
    "compatibility, use `six.add_metaclass`."
)

B304 = Error(message="B304 `sys.maxint` is not a thing on Python 3. Use `sys.maxsize`.")

B305 = Error(
    message="B305 `.next()` is not a thing on Python 3. Use the `next()` "
    "builtin. For Python 2 compatibility, use `six.next()`."
)
B305.methods = {"next"}
B305.valid_paths = {"six", "future.utils", "builtins"}

B306 = Error(
    message="B306 `BaseException.message` has been deprecated as of Python "
    "2.6 and is removed in Python 3. Use `str(e)` to access the "
    "user-readable message. Use `e.args` to access arguments passed "
    "to the exception."
)

# Warnings disabled by default.
B901 = Error(
    message="B901 Using `yield` together with `return x`. Use native "
    "`async def` coroutines or put a `# noqa` comment on this "
    "line if this was intentional."
)
B902 = Error(
    message="B902 Invalid first argument {} used for {} method. Use the "
    "canonical first argument name in methods, i.e. {}."
)
B902.implicit_classmethods = {"__new__", "__init_subclass__", "__class_getitem__"}
B902.self = ["self"]  # it's a list because the first is preferred
B902.cls = ["cls", "klass"]  # ditto.
B902.metacls = ["metacls", "metaclass", "typ", "mcs"]  # ditto.

B903 = Error(
    message="B903 Data class should either be immutable or use __slots__ to "
    "save memory. Use collections.namedtuple to generate an immutable "
    "class, or enumerate the attributes in a __slot__ declaration in "
    "the class to leave attributes mutable."
)

B950 = Error(message="B950 line too long ({} > {} characters)")

disabled_by_default = ["B901", "B902", "B903", "B950"]
