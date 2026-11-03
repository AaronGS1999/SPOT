import re


def natural_key(text: str) -> tuple[object, ...]:
    return tuple(int(t) if t.isdigit() else t for t in re.split(r"(\d+)", text))
