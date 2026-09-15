#!/usr/bin/env python3
"""Print one scalar from a PPG12-style YAML config: cfgval.py <config> <dotted.key>

Uses PyYAML when available, otherwise a minimal fallback that understands the
flat `section:\n  key: value` layout of the config files (enough for output.*
and photonjet.* scalars). Used by the photonjet shell scripts."""
import re
import sys


def load(path):
    try:
        import yaml  # type: ignore
        with open(path) as f:
            return yaml.safe_load(f)
    except ImportError:
        pass
    doc, section = {}, None
    with open(path) as f:
        for raw in f:
            line = raw.split("#", 1)[0].rstrip()
            if not line.strip():
                continue
            m = re.match(r"^(\w+):\s*$", line)
            if m:
                section = m.group(1)
                doc.setdefault(section, {})
                continue
            m = re.match(r"^\s{2}(\w+):\s*(.+)$", line)
            if m and section is not None:
                val = m.group(2).strip()
                if val.startswith(("'", '"')) and val.endswith(("'", '"')):
                    val = val[1:-1]
                doc[section][m.group(1)] = val
    return doc


def main():
    if len(sys.argv) != 3:
        sys.exit("usage: cfgval.py <config.yaml> <dotted.key>")
    doc = load(sys.argv[1])
    cur = doc
    for k in sys.argv[2].split("."):
        if not isinstance(cur, dict) or k not in cur:
            sys.exit(f"key not found: {sys.argv[2]}")
        cur = cur[k]
    print(cur if cur is not None else "")


if __name__ == "__main__":
    main()
