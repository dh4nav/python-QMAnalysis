import sys
import subprocess


def main():
    sys.exit(subprocess.call(
        ["sphinx-build", "-b", "html", "docs", "docs/_build/html"]))


if __name__ == "__main__":
    main()
