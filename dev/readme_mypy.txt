If you have functions with many arguments with complex types, you can use type aliases to try and help shorten them.

For example, instead of doing:
def broadcast_message(
        message: str,
        servers: List[Tuple[Tuple[str, int], Dict[str, str]]]) -> None:
    ...
...you could do:
ConnectionOptions = Dict[str, str]
Address = Tuple[str, int]
Server = Tuple[Address, ConnectionOptions]

def broadcast_message(message: str, servers: List[Server]) -> None:
    ...
If you'd like to type dictionaries/give types to individual fields, there's currently experimental support for TypedDicts -- you can track the issue here: https://github.com/python/mypy/issues/985.

It's currently a WIP, but I think we can probably expect to see TypedDicts added to the typing module sometime during the next few months (?). After that point, there'll probably be a little bit of lag before other 3rd party tools (like IntelliJ) support TypedDict, but it'll get there eventually.

If you'd like to see large examples of type annotations in use, I would check out Zulip and mypy's source code. Zulip also wrote a blog post talking about their experiences migrating to using mypy. I suppose typeshed might also be a good example of how to write stub files.

TypedDicts:       https://github.com/python/mypy/issues/985
Zulip:            https://github.com/zulip/zulip
Zulip post:       http://blog.zulip.org/2016/10/13/static-types-in-python-oh-mypy/
MyPy:             https://github.com/python/mypy
MyPy Cheat Sheet: http://mypy.readthedocs.io/en/latest/cheat_sheet.html
PEP 484:          https://www.python.org/dev/peps/pep-0484/
TypeShed:         https://github.com/python/typeshed
