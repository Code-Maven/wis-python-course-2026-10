import random
import re
from pathlib import Path

STUDENTS_FILE = Path(__file__).parent / "STUDENTS.md"
ROW_PATTERN = re.compile(r"^\|\s*\[([^\]]+)\]\(")


def read_students(path=STUDENTS_FILE):
    names = []
    for line in path.read_text(encoding="utf-8").splitlines():
        match = ROW_PATTERN.match(line)
        if match:
            names.append(match.group(1))
    return names


def main():
    names = read_students()
    if not names:
        raise SystemExit("No students found")
    print(random.choice(names))


if __name__ == "__main__":
    main()
