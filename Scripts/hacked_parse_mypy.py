import os
from pathlib import Path
import re
import subprocess
from datetime import datetime

script_path = Path(__file__).absolute()
root_path = script_path.parent.parent
os.chdir(root_path)

# Get the current timestamp
timestamp = datetime.now().strftime("%d%m%y%H%M%S")

# Get the latest Git hash
git_hash = subprocess.check_output(["git", "rev-parse", "HEAD"], text=True).strip()[:7]


# Run the mypy command and capture its output
process = subprocess.Popen(
    ["mypy", "--disallow-untyped-defs", "Bio"],
    stdout=subprocess.PIPE,
    stderr=subprocess.PIPE,
    text=True
)

# Process the log
merged_lines = []
merged_lines2 = []
last_line = None

for line in process.stdout:
    line = line.strip()
    if last_line:
        match = re.match(r'(.*): note: (.*)', last_line)
        if match and re.match(r'.*:\d+:\d+:', line):
            method = match.group(2)
            out_line = f"{line.split(': error:')[0]}: {method}"
            if re.findall(r'In (member|function) "_[^\W_]', last_line):
                merged_lines2.append(out_line)
            else:
                merged_lines.append(out_line)
    last_line = line


for line in process.stderr:
    line = line.strip()
    print(line)


# Output the merged lines
output_dir = Path(".")
output_dir.joinpath(f'todo_{timestamp}_{git_hash}.txt').write_text('\n'.join(merged_lines))
output_dir.joinpath(f'todo_extra_{timestamp}_{git_hash}.txt').write_text('\n'.join(merged_lines2))
