from __future__ import annotations

import json

from theory.rubin_chain import reference_run


def main() -> None:
    result = reference_run()
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
